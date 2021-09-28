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
#include "votca/xtp/threecenter.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/symmetric_matrix.h"

namespace votca {
namespace xtp {

void TCMatrix_gwbse::Initialize(Index basissize, Index mmin, Index mmax,
                                Index nmin, Index nmax) {

  // here as storage indices starting from zero
  nmin_ = nmin;
  nmax_ = nmax;
  ntotal_ = nmax - nmin + 1;
  mmin_ = mmin;
  mmax_ = mmax;
  mtotal_ = mmax - mmin + 1;
  auxbasissize_ = basissize;

  // vector has mtotal elements
  // largest object should be allocated in multithread fashion
  matrix_ = std::vector<Eigen::MatrixXd>(mtotal_);
#pragma omp parallel for schedule(dynamic, 4)
  for (Index i = 0; i < mtotal_; i++) {
    matrix_[i] = Eigen::MatrixXd::Zero(ntotal_, auxbasissize_);
  }
}

/*
 * Modify 3-center matrix elements consistent with use of symmetrized
 * Coulomb interaction using either CUDA or Openmp.
 */
void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {
  OpenMP_CUDA gemm;
  gemm.setOperators(matrix_, matrix);
#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index i = 0; i < msize(); i++) {
      gemm.MultiplyRight(matrix_[i], threadid);
    }
  }
}
/*
 * Fill the 3-center object by looping over shells of GW basis set and
 * calling FillBlock, which calculates all 3-center overlap integrals
 * associated to a particular shell, convoluted with the DFT orbital
 * coefficients
 */
void TCMatrix_gwbse::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals) {
  // needed for Rebuild())
  auxbasis_ = &auxbasis;
  dftbasis_ = &dftbasis;
  dft_orbitals_ = &dft_orbitals;

  Fill3cMO(auxbasis, dftbasis, dft_orbitals);

  AOOverlap auxoverlap;
  auxoverlap.Fill(auxbasis);
  AOCoulomb auxcoulomb;
  auxcoulomb.Fill(auxbasis);
  Eigen::MatrixXd inv_sqrt = auxcoulomb.Pseudo_InvSqrt_GWBSE(auxoverlap, 5e-7);
  removedfunctions_ = auxcoulomb.Removedfunctions();
  MultiplyRightWithAuxMatrix(inv_sqrt);

  return;
}

}  // namespace xtp
}  // namespace votca
