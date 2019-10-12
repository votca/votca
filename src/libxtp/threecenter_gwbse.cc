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

  // Try to run the operation in a Nvidia GPU, otherwise is the default Openmp
  // implementation
#if defined(USE_CUDA)
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp()
      << " Using CUDA/OpenMP for tensor matrix multiplication: " << flush;
  // Preallocated memory for all the matrices
  CudaPipeline cuda_pip;
  const Eigen::MatrixXd& head = _matrix.front();
  CudaMatrix cuma_A{head.rows(), head.cols()};
  CudaMatrix cuma_B{matrix, cuda_pip.get_stream()};
  CudaMatrix cuma_C{head.rows(), matrix.cols()};

#pragma omp parallel for schedule(dynamic)
  for (int i_occ = 0; i_occ < _mtotal; i_occ++) {
    // All the GPU communication happens through a single thread that reuses all
    // memory allocated in the GPU and it's dynamically load-balanced by OpenMP.
    // The rest of the threads use the default CPU matrix multiplication
    if (OPENMP::getThreadId() == 0) {
      cuma_A.copy_to_gpu(_matrix[i_occ]);
      cuda_pip.gemm(cuma_A, cuma_B, cuma_C);
      _matrix[i_occ] = cuma_C;
    } else {
      Eigen::MatrixXd temp = _matrix[i_occ] * matrix;
      _matrix[i_occ] = temp;
    }
  }
#else
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp()
      << " Using Default OpenMP for tensor matrix multiplication: " << flush;
  MultiplyRightWithAuxMatrixOpenMP(matrix);
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

  // If cuda is enabled the dft orbitals are sent first to the cuda gpu
  // and memory in the cuda gpu is allocated for the intermediate matrices
#if defined(USE_CUDA)
  CudaPipeline cuda_pip;
  std::array<CudaMatrix, 2> cuda_matrices =
      SendDFTMatricesToGPU(dft_orbitals, cuda_pip);
  std::array<CudaMatrix, 3> cuda_inter_matrices =
      CreateIntermediateCudaMatrices(dft_orbitals.rows());
#endif

  // loop over all shells in the GW basis and get _Mmn for that shell
#pragma omp parallel for schedule(guided)  // private(_block)
  for (int is = 0; is < gwbasis.getNumofShells(); is++) {
    const AOShell& shell = gwbasis.getShell(is);

    // Fill block for this shell (3-center overlap with _dft_basis +
    // multiplication with _dft_orbitals )
    std::vector<Eigen::MatrixXd> symmstorage =
        ComputeSymmStorage(shell, dftbasis, dft_orbitals);

    // If cuda is enable all the GPU communication happens through a single
    // thread that reuses all memory allocated in the GPU and it's dynamically
    // load-balanced by OpenMP. The remaining threads will perform the
    // convolution using the default CPU method.
    std::vector<Eigen::MatrixXd> block;
#if defined(USE_CUDA)
    if (OPENMP::getThreadId() == 0) {
      block = FillBlockCUDA(symmstorage, cuda_matrices, cuda_inter_matrices,
                            cuda_pip);
    } else {
      block = FillBlock(symmstorage, dft_orbitals);
    }
#else
    // Otherwise the convolution is performed by Eigen
    block = FillBlock(symmstorage, dft_orbitals);
#endif
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
 * GW shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock)
 */
std::vector<Eigen::MatrixXd> TCMatrix_gwbse::ComputeSymmStorage(
    const AOShell& auxshell, const AOBasis& dftbasis,
    const Eigen::MatrixXd& dft_orbitals) const {
  std::vector<Eigen::MatrixXd> symmstorage = std::vector<Eigen::MatrixXd>(
      auxshell.getNumFunc(),
      Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));
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

      Eigen::Tensor<double, 3> threec_block(auxshell.getNumFunc(),
                                            shell_row.getNumFunc(),
                                            shell_col.getNumFunc());
      threec_block.setZero();

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
                  threec_block(aux_c, row_c, col_c);
            }  // ROW copy
          }    // COL copy
        }      // AUX copy
      }
    }  // gamma-loop
  }    // alpha-loop

  return symmstorage;
}

/*
 *  Convolution of the GW shell with ALL functions with the DFT orbital
 * coefficients
 */
std::vector<Eigen::MatrixXd> TCMatrix_gwbse::FillBlock(
    const std::vector<Eigen::MatrixXd>& symmstorage,
    const Eigen::MatrixXd& dft_orbitals) const {

  std::vector<Eigen::MatrixXd> block = std::vector<Eigen::MatrixXd>(
      _mtotal, Eigen::MatrixXd::Zero(_ntotal, symmstorage.size()));
  const Eigen::MatrixXd dftm =
      dft_orbitals.block(0, _mmin, dft_orbitals.rows(), _mtotal);
  const Eigen::MatrixXd dftn =
      dft_orbitals.block(0, _nmin, dft_orbitals.rows(), _ntotal);

  int dim = static_cast<int>(symmstorage.size());
  for (auto k = 0; k < dim; ++k) {
    const Eigen::MatrixXd& matrix = symmstorage[k];
    Eigen::MatrixXd threec_inMo =
        dftn.transpose() * matrix.selfadjointView<Eigen::Lower>() * dftm;
    for (int i = 0; i < threec_inMo.cols(); ++i) {
      block[i].col(k) = threec_inMo.col(i);
    }
  }
  return block;
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

#if defined(USE_CUDA)
/*
 * Convolution of the GW shell with the DFT orbital coefficients using an Nvidia
 * GPU. The Cuda device behaves like a server that is receiving matrix-matrix
 * multiplications from a single stream (an Nvidia queue) and handle them
 * in an asynchronous way. It performs the following operations when recieving a
 * request:
 *  1. Check that there is enough space for the arrays
 *  2. Allocate memory for each matrix
 *  3. Copy the matrix to the allocated space
 *  4. Perform the matrix multiplication
 *  5. Return the result matrix
 * The Cuda device knows to which memory address it needs to copy back the
 * result. see: https://docs.nvidia.com/cuda/cublas/index.html#thread-safety2
 */
std::vector<Eigen::MatrixXd> TCMatrix_gwbse::FillBlockCUDA(
    const std::vector<Eigen::MatrixXd>& symmstorage,
    const std::array<CudaMatrix, 2>& cuda_matrices,
    std::array<CudaMatrix, 3>& cuda_inter_matrices,
    const CudaPipeline& cuda_pip) const {

  std::vector<Eigen::MatrixXd> block = std::vector<Eigen::MatrixXd>(
      _mtotal, Eigen::MatrixXd::Zero(_ntotal, symmstorage.size()));

  try {
    const CudaMatrix& cuma_A = cuda_matrices[0];
    const CudaMatrix& cuma_C = cuda_matrices[1];
    CudaMatrix& cuma_B = cuda_inter_matrices[0];
    CudaMatrix& cuma_X = cuda_inter_matrices[1];
    CudaMatrix& cuma_Y = cuda_inter_matrices[2];

    int dim = static_cast<int>(symmstorage.size());
    for (int k = 0; k < dim; ++k) {
      const Eigen::MatrixXd& matrix = symmstorage[k];
      cuma_B.copy_to_gpu(matrix.selfadjointView<Eigen::Lower>());
      cuda_pip.gemm(cuma_A, cuma_B, cuma_X);
      cuda_pip.gemm(cuma_X, cuma_C, cuma_Y);
      Eigen::MatrixXd threec_inMo = cuma_Y;
      for (int i = 0; i < threec_inMo.cols(); ++i) {
        block[i].col(k) = threec_inMo.col(i);
      }
    }
  } catch (const std::runtime_error& error) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " FillBlockCUDA failed due to: " << error.what()
        << flush;
    throw;
  }
  return block;
}

std::array<CudaMatrix, 2> TCMatrix_gwbse::SendDFTMatricesToGPU(
    const Eigen::MatrixXd& dft_orbitals, const CudaPipeline& cuda_pip) const {
  const Eigen::MatrixXd dftm =
      dft_orbitals.block(0, _mmin, dft_orbitals.rows(), _mtotal);
  const Eigen::MatrixXd dftn =
      dft_orbitals.block(0, _nmin, dft_orbitals.rows(), _ntotal);

  // Smart Pointers to the cuda arrays
  const cudaStream_t& stream = cuda_pip.get_stream();

  return {CudaMatrix{dftn.transpose(), stream}, CudaMatrix{dftm, stream}};
}

std::array<CudaMatrix, 3> TCMatrix_gwbse::CreateIntermediateCudaMatrices(
    long basissize) const {
  long mcols = _mtotal - _mmin;
  long ncols = _ntotal - _nmin;

  return {CudaMatrix{basissize, basissize}, CudaMatrix{ncols, basissize},
          CudaMatrix{ncols, mcols}};
}
#endif

}  // namespace xtp
}  // namespace votca
