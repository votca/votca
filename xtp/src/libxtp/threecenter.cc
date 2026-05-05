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
#include "votca/xtp/aoreactionfield.h"
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
 * Helper: compute the pseudo-inverse square root of a combined
 * (V_bare + V_reac) matrix using the same overlap-metric regularisation as
 * AOCoulomb::Pseudo_InvSqrt_GWBSE.  Returns (Vm1 * Ssqrt)^T and sets
 * removedfunctions to the number of discarded eigenvalues.
 */
static Eigen::MatrixXd PseudoInvSqrtCombined(const Eigen::MatrixXd& V_combined,
                                              const Eigen::MatrixXd& S_aux,
                                              Index& removedfunctions) {
  constexpr double etol = 5e-7;
  removedfunctions = 0;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(S_aux);
  Eigen::VectorXd diag_overlap =
      Eigen::VectorXd::Zero(eo.eigenvalues().size());
  for (Index i = 0; i < diag_overlap.size(); ++i) {
    if (eo.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diag_overlap(i) = 1.0 / std::sqrt(eo.eigenvalues()(i));
    }
  }
  Eigen::MatrixXd Ssqrt = eo.eigenvectors() * diag_overlap.asDiagonal() *
                          eo.eigenvectors().transpose();

  Eigen::MatrixXd ortho = Ssqrt * V_combined * Ssqrt;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho);
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(es.eigenvalues().size());
  for (Index i = 0; i < diag.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diag(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }
  Eigen::MatrixXd Vm1 =
      es.eigenvectors() * diag.asDiagonal() * es.eigenvectors().transpose();
  return (Vm1 * Ssqrt).transpose();
}

/*
 * Fill the 3-center object by looping over shells of GW basis set and
 * calling FillBlock, which calculates all 3-center overlap integrals
 * associated to a particular shell, convoluted with the DFT orbital
 * coefficients.
 */
void TCMatrix_gwbse::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals) {
  // needed for Rebuild()
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

/*
 * Fill overload that additionally computes and adds the MM reaction-field
 * correction to the 2-centre Coulomb matrix before forming V^{-1/2}.
 *
 * The bare 2-centre Coulomb matrix V is augmented with the reaction-field
 * matrix V_reac before forming the pseudo-inverse square root:
 *
 *   V_tilde = V_bare + V_reac
 *   I^ml_mu <- M^ml_nu * (V_tilde)^{-1/2}_{nu,mu}
 *
 * so that the screened Coulomb interaction W is built from the renormalised
 * interaction v_tilde, intrinsically including MM polarisation screening.
 * See Duchemin, Jacquemin, Blase, J. Chem. Phys. 144, 164106 (2016) and
 * the VOTCA WIREs review (Tirimbó, Sundaram, Baumeier 2024) section 3.3.1.
 *
 * polar_segs must contain Thole-parametrised PolarSites.  Their induced-
 * dipole fields are used as scratch and are reset to zero on return.
 * grid_quality is passed to Vxc_Grid::GridSetup; "xfine" is recommended
 * because aux-basis functions can be more contracted than the DFT basis.
 *
 * NOTE: Rebuild() calls the plain overload and therefore does NOT re-apply
 * the reaction field.  In an evGW loop the caller should use this overload
 * for every iteration that rebuilds the 3-centre integrals.
 */
void TCMatrix_gwbse::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis,
                          const Eigen::MatrixXd& dft_orbitals,
                          const QMMolecule& qmatoms,
                          std::vector<PolarSegment>& polar_segs,
                          const std::string& grid_quality) {
  auxbasis_ = &auxbasis;
  dftbasis_ = &dftbasis;
  dft_orbitals_ = &dft_orbitals;

  Fill3cMO(auxbasis, dftbasis, dft_orbitals);

  AOOverlap auxoverlap;
  auxoverlap.Fill(auxbasis);
  AOCoulomb auxcoulomb;
  auxcoulomb.Fill(auxbasis);

  // Compute V_reac and form the combined Coulomb matrix
  AOReactionField rf;
  rf.Fill(auxbasis, qmatoms, polar_segs, grid_quality);

  Eigen::MatrixXd V_combined = auxcoulomb.Matrix() + rf.Matrix();

  // Compute pseudo-inverse square root of the combined matrix
  Index removedfunctions = 0;
  Eigen::MatrixXd inv_sqrt =
      PseudoInvSqrtCombined(V_combined, auxoverlap.Matrix(), removedfunctions);

  removedfunctions_ = removedfunctions;
  MultiplyRightWithAuxMatrix(inv_sqrt);

  return;
}

}  // namespace xtp
}  // namespace votca
