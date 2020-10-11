/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/threecenter.h"
using std::flush;

namespace votca {
namespace xtp {

void TCMatrix_gwbse::Initialize(Index basissize, Index mmin, Index mmax,
                                Index nmin, Index nmax) {

  // here as storage indices starting from zero
  _nmin = nmin;
  _nmax = nmax;
  _ntotal = nmax - nmin + 1;
  _mmin = mmin;
  _mmax = mmax;
  _mtotal = mmax - mmin + 1;
  _auxbasissize = basissize;

  // vector has mtotal elements
  _matrix = std::vector<Eigen::MatrixXd>(
      _mtotal, Eigen::MatrixXd::Zero(_ntotal, _auxbasissize));
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using either CUDA or Openmp
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {
  OpenMP_CUDA gemm;
  gemm.setOperators(_matrix, matrix);
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < msize(); i++) {
    gemm.MultiplyRight(_matrix[i]);
  }
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

  Fill3cMO(gwbasis, dftbasis, dft_orbitals);

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
 * GW shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock)
 */
std::vector<Eigen::MatrixXd> TCMatrix_gwbse::ComputeAO3cBlock(
    const AOShell& auxshell, const AOBasis& dftbasis) const {
  std::vector<Eigen::MatrixXd> ao3c = std::vector<Eigen::MatrixXd>(
      auxshell.getNumFunc(),
      Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));
  // alpha-loop over the "left" DFT basis function
  for (Index row = 0; row < dftbasis.getNumofShells(); row++) {

    const AOShell& shell_row = dftbasis.getShell(row);
    const Index row_start = shell_row.getStartIndex();
    // ThreecMatrix is symmetric, restrict explicit calculation to triangular
    // matrix
    for (Index col = 0; col <= row; col++) {
      const AOShell& shell_col = dftbasis.getShell(col);
      const Index col_start = shell_col.getStartIndex();

      Eigen::Tensor<double, 3> threec_block(auxshell.getNumFunc(),
                                            shell_row.getNumFunc(),
                                            shell_col.getNumFunc());
      threec_block.setZero();

      bool nonzero =
          FillThreeCenterRepBlock(threec_block, auxshell, shell_row, shell_col);
      if (nonzero) {
        for (Index aux_c = 0; aux_c < auxshell.getNumFunc(); aux_c++) {
          for (Index row_c = 0; row_c < shell_row.getNumFunc(); row_c++) {
            for (Index col_c = 0; col_c < shell_col.getNumFunc(); col_c++) {
              // symmetry
              if ((col_start + col_c) > (row_start + row_c)) {
                break;
              }
              ao3c[aux_c](row_start + row_c, col_start + col_c) =
                  threec_block(aux_c, row_c, col_c);
            }  // ROW copy
          }    // COL copy
        }      // AUX copy
      }
    }  // gamma-loop
  }    // alpha-loop

  for (Eigen::MatrixXd& mat : ao3c) {
    mat.triangularView<Eigen::Upper>() =
        mat.triangularView<Eigen::Lower>().transpose();
  }
  return ao3c;
}

void TCMatrix_gwbse::Fill3cMO(const AOBasis& gwbasis, const AOBasis& dftbasis,
                              const Eigen::MatrixXd& dft_orbitals) {

  const Eigen::MatrixXd dftm =
      dft_orbitals.block(0, _mmin, dft_orbitals.rows(), _mtotal);
  const Eigen::MatrixXd dftn =
      dft_orbitals.block(0, _nmin, dft_orbitals.rows(), _ntotal).transpose();

  OpenMP_CUDA transform;
  transform.setOperators(dftn, dftm);
#pragma omp parallel for schedule(guided)
  for (Index is = 0; is < gwbasis.getNumofShells(); is++) {
    const AOShell& shell = gwbasis.getShell(is);

    std::vector<Eigen::MatrixXd> ao3c = ComputeAO3cBlock(shell, dftbasis);

    // this is basically a transpose of AO3c and at the same time the ao->mo
    // transformation
    std::vector<Eigen::MatrixXd> block = std::vector<Eigen::MatrixXd>(
        _mtotal, Eigen::MatrixXd::Zero(_ntotal, ao3c.size()));

    Index dim = static_cast<Index>(ao3c.size());
    for (Index k = 0; k < dim; ++k) {
      transform.MultiplyLeftRight(ao3c[k]);
      for (Index i = 0; i < ao3c[k].cols(); ++i) {
        block[i].col(k) = ao3c[k].col(i);
      }
    }

    // put into correct position
    for (Index m_level = 0; m_level < _mtotal; m_level++) {
      _matrix[m_level].block(0, shell.getStartIndex(), _ntotal,
                             shell.getNumFunc()) = block[m_level];
    }  // m-th DFT orbital
  }    // shells of GW basis set
}

}  // namespace xtp
}  // namespace votca
