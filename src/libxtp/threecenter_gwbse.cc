/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multiarray.h>
#include <votca/xtp/threecenter.h>

using std::flush;

namespace votca {
namespace xtp {

void TCMatrix_gwbse::Initialize(int basissize, int mmin, int mmax, int nmin,
                                int nmax) {

  // here as storage indices starting from zero
  _nmin = nmin;
  _nmax = nmax;
  _ntotal = nmax - nmin + 1;
  _mmin = mmin;
  _mmax = mmax;
  _mtotal = mmax - mmin + 1;
  _basissize = basissize;

  // vector has mtotal elements
  _matrix = std::vector<Eigen::MatrixXd>(
      _mtotal, Eigen::MatrixXd::Zero(_ntotal, _basissize));
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using either CUDA or Openmp
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {

  // Try to run the operation in a GPU, otherwise is the default Openmp
  // implementation
#if defined(USE_GPU)
  EigenCuda gpu_handle;
  try {
    _matrix = gpu_handle.right_matrix_tensor_mult(_matrix, matrix);
  } catch (const std::runtime_error& error) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp()
        << " GPU Multiplyrightwithauxmatrix failed due to: " << error.what()
        << " Using default OpenMP implementation!" << flush;
    this->MultiplyRightWithAuxMatrixOpenMP(matrix);
  }
#else
  this->MultiplyRightWithAuxMatrixOpenMP(matrix);
#endif
  return;
}

/*
 * Fill the 3-center object by looping over shells of GW basis set and
 * calling FillBlock, which calculates all 3-center overlap integrals
 * associated to a particular shell, convoluted with the DFT orbital
 * coefficients
 */
void TCMatrix_gwbse::Fill(const AOBasis& gwbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals) {
  // needed for Rebuild())
  _auxbasis = &gwbasis;
  _dftbasis = &dftbasis;
  _dft_orbitals = &dft_orbitals;

  // loop over all shells in the GW basis and get _Mmn for that shell
#pragma omp parallel for schedule(guided)  // private(_block)
  for (int is = 0; is < gwbasis.getNumofShells(); is++) {
    const AOShell& shell = gwbasis.getShell(is);
    std::vector<Eigen::MatrixXd> block;
    for (int i = 0; i < _mtotal; i++) {
      block.push_back(Eigen::MatrixXd::Zero(_ntotal, shell.getNumFunc()));
    }
    // Fill block for this shell (3-center overlap with _dft_basis +
    // multiplication with _dft_orbitals )
    FillBlock(block, shell, dftbasis, dft_orbitals);

    // put into correct position
    for (int m_level = 0; m_level < _mtotal; m_level++) {
      _matrix[m_level].block(0, shell.getStartIndex(), _ntotal,
                             shell.getNumFunc()) = block[m_level];
    }  // m-th DFT orbital
  }    // shells of GW basis set

  AOOverlap auxoverlap;
  auxoverlap.Fill(gwbasis);
  AOCoulomb auxcoulomb;
  auxcoulomb.Fill(gwbasis);
  Eigen::MatrixXd inv_sqrt = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, 5e-7);
  _removedfunctions = auxcoulomb.Removedfunctions();
  MultiplyRightWithAuxMatrix(inv_sqrt);
  return;
}

/*
 * Determines the 3-center integrals for a given shell in the GW basis
 * by calculating the 3-center overlap integral of the functions in the
 * GW shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock),
 * followed by a convolution of those with the DFT orbital coefficients
 */

void TCMatrix_gwbse::FillBlock(std::vector<Eigen::MatrixXd>& block,
                               const AOShell& auxshell, const AOBasis& dftbasis,
                               const Eigen::MatrixXd& dft_orbitals) {
  tensor3d::extent_gen extents;
  std::vector<Eigen::MatrixXd> symmstorage;
  for (int i = 0; i < auxshell.getNumFunc(); ++i) {
    symmstorage.push_back(
        Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));
  }
  const Eigen::MatrixXd dftm =
      dft_orbitals.block(0, _mmin, dft_orbitals.rows(), _mtotal);
  const Eigen::MatrixXd dftn =
      dft_orbitals.block(0, _nmin, dft_orbitals.rows(), _ntotal);
  // alpha-loop over the "left" DFT basis function
  for (int row = 0; row < dftbasis.getNumofShells(); row++) {

    const AOShell& shell_row = dftbasis.getShell(row);
    const int row_start = shell_row.getStartIndex();
    // ThreecMatrix is symmetric, restrict explicit calculation to triangular
    // matrix
    for (int col = 0; col <= row; col++) {
      const AOShell& shell_col = dftbasis.getShell(col);
      const int col_start = shell_col.getStartIndex();

      tensor3d threec_block(extents[range(0, auxshell.getNumFunc())][range(
          0, shell_row.getNumFunc())][range(0, shell_col.getNumFunc())]);
      std::fill_n(threec_block.data(), threec_block.num_elements(), 0.0);

      bool nonzero =
          FillThreeCenterRepBlock(threec_block, auxshell, shell_row, shell_col);
      if (nonzero) {
        for (int aux_c = 0; aux_c < auxshell.getNumFunc(); aux_c++) {
          for (int row_c = 0; row_c < shell_row.getNumFunc(); row_c++) {
            for (int col_c = 0; col_c < shell_col.getNumFunc(); col_c++) {
              // symmetry
              if ((col_start + col_c) > (row_start + row_c)) {
                break;
              }
              symmstorage[aux_c](row_start + row_c, col_start + col_c) =
                  threec_block[aux_c][row_c][col_c];
            }  // ROW copy
          }    // COL copy
        }      // AUX copy
      }
    }  // gamma-loop
  }    // alpha-loop
  for (int k = 0; k < auxshell.getNumFunc(); ++k) {
    const Eigen::MatrixXd& matrix = symmstorage[k];
    Eigen::MatrixXd threec_inMo =
        dftn.transpose() * matrix.selfadjointView<Eigen::Lower>() * dftm;
    for (int i = 0; i < threec_inMo.cols(); ++i) {
      block[i].col(k) = threec_inMo.col(i);
    }
  }
  return;
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using OPENMP parallelization
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrixOpenMP(
    const Eigen::MatrixXd& matrix) {
#pragma omp parallel for
  for (int i_occ = 0; i_occ < _mtotal; i_occ++) {
    Eigen::MatrixXd temp = _matrix[i_occ] * matrix;
    _matrix[i_occ] = temp;
  }
  return;
}

}  // namespace xtp
}  // namespace votca
