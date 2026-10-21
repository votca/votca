/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/uks_convergenceacc.h"
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/trah.h"

#include <algorithm>

#include <tuple>

namespace votca {
namespace xtp {

void UKSConvergenceAcc::Configure(const options& opt_alpha,
                                  const options& opt_beta) {
  opt_alpha_ = opt_alpha;
  opt_beta_ = opt_beta;

  nocclevels_alpha_ = opt_alpha_.numberofelectrons;
  nocclevels_beta_ = opt_beta_.numberofelectrons;

  // one shared DIIS/ADIIS history length
  diis_.setHistLength(opt_alpha_.histlength);
  trah_options_.max_iterations = opt_alpha_.davidson_max_iter;
  // adiis_.setHistLength(opt_alpha_.histlength);
}

void UKSConvergenceAcc::setLogger(Logger* log) { log_ = log; }

void UKSConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  // As ConvergenceAcc::setOverlap: removed directions are kept out of the
  // occupied and low virtual space via removed_projector_.
  S_ = &S;
  const SymmetricEigenSystem es = SymmetricEigen(S.Matrix());
  const Eigen::VectorXd& s_eig = es.values;
  Eigen::VectorXd inv_sqrt = Eigen::VectorXd::Zero(s_eig.size());
  Index removed = 0;
  for (Index i = 0; i < s_eig.size(); ++i) {
    if (s_eig(i) < etol) {
      ++removed;
    } else {
      inv_sqrt(i) = 1.0 / std::sqrt(s_eig(i));
    }
  }
  Sminusahalf = es.vectors * inv_sqrt.asDiagonal() * es.vectors.transpose();
  removed_projector_.resize(0, 0);
  if (removed > 0) {
    const Eigen::MatrixXd U = es.vectors.leftCols(removed);
    removed_projector_ = U * U.transpose();
  }
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Smallest value of AOOverlap matrix is " << s_eig(0)
      << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Removed " << removed
      << " basisfunction from inverse overlap matrix (threshold " << etol << ")"
      << std::flush;
  if (removed == 0 && s_eig(0) < 1e3 * etol) {
    XTP_LOG(Log::error, *log_)
        << TimeStamp()
        << " WARNING: the overlap matrix is nearly singular; if the SCF is "
           "unstable, raise xtpdft.overlap_tolerance above "
        << s_eig(0) << "." << std::flush;
  }
}

tools::EigenSystem UKSConvergenceAcc::SolveFockmatrix(
    const Eigen::MatrixXd& H) const {
  Eigen::MatrixXd H_ortho = Sminusahalf.transpose() * H * Sminusahalf;
  if (removed_projector_.size() > 0) {
    H_ortho += kRemovedShift * removed_projector_;
  }
  const SymmetricEigenSystem es = SymmetricEigen(H_ortho);
  if (!es.values.allFinite()) {
    throw std::runtime_error("Matrix Diagonalisation failed");
  }

  tools::EigenSystem result;
  result.eigenvalues() = es.values;
  result.eigenvectors() = Sminusahalf * es.vectors;
  return result;
}

Eigen::MatrixXd UKSConvergenceAcc::UnflattenRotation(
    const Eigen::VectorXd& v_ov, Index nao, Index nocclevels) const {
  Index nvirt = nao - nocclevels;
  Eigen::MatrixXd kappa = Eigen::MatrixXd::Zero(nao, nao);
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      double val = v_ov(i * nvirt + a);
      kappa(i, nocclevels + a) = val;
      kappa(nocclevels + a, i) = -val;
    }
  }
  return kappa;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    UKSConvergenceAcc::UnflattenCoupledRotation(const Eigen::VectorXd& v,
                                                Index nao_alpha,
                                                Index nocclevels_alpha,
                                                Index nao_beta,
                                                Index nocclevels_beta) const {
  Index n_ov_alpha = nocclevels_alpha * (nao_alpha - nocclevels_alpha);
  Eigen::VectorXd v_alpha = v.head(n_ov_alpha);
  Eigen::VectorXd v_beta = v.tail(v.size() - n_ov_alpha);
  Eigen::MatrixXd kappa_alpha =
      UnflattenRotation(v_alpha, nao_alpha, nocclevels_alpha);
  Eigen::MatrixXd kappa_beta =
      UnflattenRotation(v_beta, nao_beta, nocclevels_beta);
  return {kappa_alpha, kappa_beta};
}

Eigen::VectorXd UKSConvergenceAcc::BuildCoupledSigmaVector(
    const Eigen::VectorXd& v, const Eigen::MatrixXd& C_alpha,
    Index nocclevels_alpha, const Eigen::MatrixXd& C_beta,
    Index nocclevels_beta, const CoupledFockBuilder& coupled_fock_builder,
    double finite_diff_step) const {
  Index nao_alpha = C_alpha.rows();
  Index nvirt_alpha = nao_alpha - nocclevels_alpha;
  Index nao_beta = C_beta.rows();
  Index nvirt_beta = nao_beta - nocclevels_beta;
  Index n_ov_alpha = nocclevels_alpha * nvirt_alpha;
  Index n_ov_beta = nocclevels_beta * nvirt_beta;

  auto [kappa_alpha_trial, kappa_beta_trial] = UnflattenCoupledRotation(
      v, nao_alpha, nocclevels_alpha, nao_beta, nocclevels_beta);
  kappa_alpha_trial *= finite_diff_step;
  kappa_beta_trial *= finite_diff_step;

  // Same linearized-exponential + Lowdin-reorthonormalization approach
  // as BuildSigmaVector's own EvaluateGradientAt -- but rotates BOTH
  // channels together and builds BOTH new Fock matrices together in
  // ONE coupled_fock_builder call, rather than one channel at a time
  // with the other held fixed. This is what actually captures the
  // alpha-beta coupling: the shared Coulomb potential and the XC
  // kernel's cross-spin terms naturally see both perturbed densities
  // together inside coupled_fock_builder, rather than needing this
  // function to construct any cross-coupling term explicitly itself.
  auto EvaluateBothGradientsAt = [&](const Eigen::MatrixXd& kappa_alpha,
                                     const Eigen::MatrixXd& kappa_beta)
      -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd> {
    Eigen::MatrixXd C_alpha_rot =
        C_alpha *
        (Eigen::MatrixXd::Identity(nao_alpha, nao_alpha) + kappa_alpha);
    Eigen::MatrixXd nonortho_alpha =
        C_alpha_rot.transpose() * S_->Matrix() * C_alpha_rot;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_alpha(nonortho_alpha);
    C_alpha_rot = C_alpha_rot * es_alpha.operatorInverseSqrt();

    Eigen::MatrixXd C_beta_rot =
        C_beta * (Eigen::MatrixXd::Identity(nao_beta, nao_beta) + kappa_beta);
    Eigen::MatrixXd nonortho_beta =
        C_beta_rot.transpose() * S_->Matrix() * C_beta_rot;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_beta(nonortho_beta);
    C_beta_rot = C_beta_rot * es_beta.operatorInverseSqrt();

    Eigen::MatrixXd C_alpha_occ_rot = C_alpha_rot.leftCols(nocclevels_alpha);
    Eigen::MatrixXd D_alpha_rot = C_alpha_occ_rot * C_alpha_occ_rot.transpose();
    Eigen::MatrixXd C_beta_occ_rot = C_beta_rot.leftCols(nocclevels_beta);
    Eigen::MatrixXd D_beta_rot = C_beta_occ_rot * C_beta_occ_rot.transpose();

    SpinFock H_rot = coupled_fock_builder(D_alpha_rot, D_beta_rot);
    Eigen::MatrixXd F_MO_alpha_rot =
        C_alpha_rot.transpose() * H_rot.alpha * C_alpha_rot;
    Eigen::MatrixXd F_MO_beta_rot =
        C_beta_rot.transpose() * H_rot.beta * C_beta_rot;
    return {F_MO_alpha_rot, F_MO_beta_rot};
  };

  auto [F_MO_alpha_plus, F_MO_beta_plus] =
      EvaluateBothGradientsAt(kappa_alpha_trial, kappa_beta_trial);
  auto [F_MO_alpha_minus, F_MO_beta_minus] =
      EvaluateBothGradientsAt(-kappa_alpha_trial, -kappa_beta_trial);

  // The energy gradient is dE/dkappa_ia = -2 F_ia (see
  // EnergyGradientScale), so the Hessian-vector product is -2 times the
  // central difference of F_ia along v.
  const double scale = EnergyGradientScale / (2.0 * finite_diff_step);
  Eigen::VectorXd sigma(n_ov_alpha + n_ov_beta);
  for (Index i = 0; i < nocclevels_alpha; ++i) {
    for (Index a = 0; a < nvirt_alpha; ++a) {
      double g_plus_ia = F_MO_alpha_plus(i, nocclevels_alpha + a);
      double g_minus_ia = F_MO_alpha_minus(i, nocclevels_alpha + a);
      sigma(i * nvirt_alpha + a) = scale * (g_plus_ia - g_minus_ia);
    }
  }
  for (Index i = 0; i < nocclevels_beta; ++i) {
    for (Index a = 0; a < nvirt_beta; ++a) {
      double g_plus_ia = F_MO_beta_plus(i, nocclevels_beta + a);
      double g_minus_ia = F_MO_beta_minus(i, nocclevels_beta + a);
      sigma(n_ov_alpha + i * nvirt_beta + a) = scale * (g_plus_ia - g_minus_ia);
    }
  }
  return sigma;
}

Eigen::MatrixXd UKSConvergenceAcc::DirectMinimizationRotation(
    const Eigen::MatrixXd& H_AO, const tools::EigenSystem& MOs,
    Index nocclevels, double& predicted_energy_change) const {
  Index nao = MOs.eigenvectors().rows();
  const Eigen::MatrixXd& C = MOs.eigenvectors();
  const Eigen::VectorXd& eps = MOs.eigenvalues();

  // Orbital gradient: the occ-virt block of the MO-basis Fock matrix,
  // which vanishes exactly at self-consistency (F_MO would be block-
  // diagonal, occ-occ and virt-virt only).
  Eigen::MatrixXd F_MO = C.transpose() * H_AO * C;

  // Antisymmetric rotation generator: nonzero only in the occ-virt (and
  // virt-occ, by antisymmetry) blocks -- occ-occ/virt-virt rotations
  // would leave the density matrix (and therefore the energy) entirely
  // unchanged, so there is nothing to gain from including them.
  Eigen::MatrixXd kappa = Eigen::MatrixXd::Zero(nao, nao);
  // g_h_ratio(i,a) below stores kappa_ia's own g_ia/h_ia BEFORE the
  // whole-matrix trust-radius scaling further down -- needed again
  // afterward to compute the quadratic model's predicted energy change
  // using the FINAL (per-element- and trust-radius-clamped) kappa,
  // since g_ia and h_ia themselves do not change under that later,
  // uniform rescaling, only kappa does.
  Eigen::MatrixXd h_matrix = Eigen::MatrixXd::Zero(nao, nao);
  // Guards near-degenerate occ-virt orbital pairs from producing a
  // wildly oversized step for THAT pair specifically -- confirmed
  // directly against an independent ORCA reference run on this exact
  // system, whose own TRAH implementation repeatedly hit BOTH
  // genuinely negative gaps (occupied/virtual character swapping
  // during iteration, since occupation here is defined by column
  // index, not current energy ordering -- down to -0.635 Ha) and
  // gaps numerically indistinguishable from zero (as small as
  // 0.000002 Ha), needing its own explicit warnings and special
  // handling for both. abs() here is essential, not cosmetic: a
  // signed gap of -0.6 would otherwise pass std::max(-0.6, kMinGap)
  // as if it were the tiny POSITIVE value kMinGap, producing an
  // enormous, wrong-direction step for exactly that pair -- the
  // Hessian approximation must be treated as positive-definite
  // (a legitimate descent direction) regardless of the true, possibly
  // negative or near-zero curvature this diagonal approximation
  // cannot itself represent.
  constexpr double kMinGap = 1e-3;
  // Per-PAIR step cap, independent of and in addition to the overall,
  // whole-matrix trust-radius bound below -- confirmed necessary
  // because a single pathological pair (as above) could otherwise
  // dominate the entire rotation DIRECTION even after the whole
  // matrix is normalized to the trust radius, silently scaling every
  // other, well-behaved pair down to near-nothing while that one bad
  // pair still controls where the step actually points.
  constexpr double kMaxKappaElement = 0.1;
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = nocclevels; a < nao; ++a) {
      double gap = std::max(std::abs(eps(a) - eps(i)), kMinGap);
      // Energy gradient g_ia = -2 F_ia (EnergyGradientScale) and the
      // diagonal of the energy Hessian, 2 (eps_a - eps_i): the leading
      // term of d^2E/dkappa_ia^2 for one spin channel.
      double g_ia = EnergyGradientScale * F_MO(i, a);
      double h_ia = 2.0 * gap;
      double kappa_ia = -g_ia / h_ia;
      kappa_ia = std::clamp(kappa_ia, -kMaxKappaElement, kMaxKappaElement);
      kappa(i, a) = kappa_ia;
      kappa(a, i) = -kappa_ia;
      h_matrix(i, a) = h_ia;
    }
  }

  double knorm = kappa.norm();
  if (knorm > trust_.Radius()) {
    kappa *= (trust_.Radius() / knorm);
  }

  // Quadratic model's predicted energy change, Sum_ia[g_ia*kappa_ia
  // + 0.5*h_ia*kappa_ia^2], for the kappa actually applied (after the
  // per-element clamp and the trust-radius scaling), as needed by
  // Iterate's Fletcher accept/reject check.
  predicted_energy_change = 0.0;
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = nocclevels; a < nao; ++a) {
      double kappa_ia = kappa(i, a);
      predicted_energy_change += EnergyGradientScale * F_MO(i, a) * kappa_ia +
                                 0.5 * h_matrix(i, a) * kappa_ia * kappa_ia;
    }
  }

  // exp(kappa) ~= I + kappa is valid specifically because the trust-
  // radius bound above keeps kappa small -- re-orthonormalize
  // afterward since I+kappa is only approximately unitary, using the
  // same Lowdin (symmetric orthogonalization) approach as
  // DFTEngine::OrthogonalizeGuess.
  Eigen::MatrixXd C_new = C * (Eigen::MatrixXd::Identity(nao, nao) + kappa);
  Eigen::MatrixXd nonortho = C_new.transpose() * S_->Matrix() * C_new;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_ortho(nonortho);
  return C_new * es_ortho.operatorInverseSqrt();
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    UKSConvergenceAcc::CoupledAugmentedHessianStep(
        const Eigen::MatrixXd& H_AO_alpha, const tools::EigenSystem& MOs_alpha,
        Index nocclevels_alpha, const Eigen::MatrixXd& H_AO_beta,
        const tools::EigenSystem& MOs_beta, Index nocclevels_beta,
        const CoupledFockBuilder& coupled_fock_builder, double trust_radius,
        double& predicted_energy_change) const {
  Index nao_alpha = MOs_alpha.eigenvectors().rows();
  Index nvirt_alpha = nao_alpha - nocclevels_alpha;
  Index n_ov_alpha = nocclevels_alpha * nvirt_alpha;
  const Eigen::MatrixXd& C_alpha = MOs_alpha.eigenvectors();
  const Eigen::VectorXd& eps_alpha = MOs_alpha.eigenvalues();

  Index nao_beta = MOs_beta.eigenvectors().rows();
  Index nvirt_beta = nao_beta - nocclevels_beta;
  Index n_ov_beta = nocclevels_beta * nvirt_beta;
  const Eigen::MatrixXd& C_beta = MOs_beta.eigenvectors();
  const Eigen::VectorXd& eps_beta = MOs_beta.eigenvalues();

  Index n_ov = n_ov_alpha + n_ov_beta;

  Eigen::MatrixXd F_MO_alpha = C_alpha.transpose() * H_AO_alpha * C_alpha;
  Eigen::MatrixXd F_MO_beta = C_beta.transpose() * H_AO_beta * C_beta;

  // Energy gradient g = -2 F_ia (EnergyGradientScale) and the diagonal
  // preconditioner 2 (eps_a - eps_i), the leading diagonal term of the
  // energy Hessian that BuildCoupledSigmaVector applies; alpha block
  // first, beta block after it, as in UnflattenCoupledRotation. The
  // preconditioner has no alpha-beta cross terms; only the Davidson
  // convergence depends on it, not the step.
  Eigen::VectorXd g(n_ov);
  Eigen::VectorXd diag_h(n_ov);
  constexpr double kMinGap = 1e-3;
  for (Index i = 0; i < nocclevels_alpha; ++i) {
    for (Index a = 0; a < nvirt_alpha; ++a) {
      g(i * nvirt_alpha + a) =
          EnergyGradientScale * F_MO_alpha(i, nocclevels_alpha + a);
      double gap = std::max(
          std::abs(eps_alpha(nocclevels_alpha + a) - eps_alpha(i)), kMinGap);
      diag_h(i * nvirt_alpha + a) = 2.0 * gap;
    }
  }
  for (Index i = 0; i < nocclevels_beta; ++i) {
    for (Index a = 0; a < nvirt_beta; ++a) {
      g(n_ov_alpha + i * nvirt_beta + a) =
          EnergyGradientScale * F_MO_beta(i, nocclevels_beta + a);
      double gap = std::max(
          std::abs(eps_beta(nocclevels_beta + a) - eps_beta(i)), kMinGap);
      diag_h(n_ov_alpha + i * nvirt_beta + a) = 2.0 * gap;
    }
  }

  // Hessian products by finite differences of the coupled Fock matrices,
  // one vector at a time (two coupled Fock builds each).
  bool analytic = bool(orbital_hessian_) && !analytic_hessian_failed_;
  trah::HessianProduct hessian = [&](const Eigen::MatrixXd& V) {
    if (analytic) {
      try {
        return orbital_hessian_(MOs_alpha, MOs_beta, H_AO_alpha, H_AO_beta, V);
      } catch (const std::runtime_error& error) {
        XTP_LOG(Log::error, *log_)
            << TimeStamp() << " Analytic UKS Hessian not available ("
            << error.what() << "); using finite differences" << std::flush;
        analytic = false;
        analytic_hessian_failed_ = true;
      }
    }
    Eigen::MatrixXd HV(V.rows(), V.cols());
    for (Index k = 0; k < V.cols(); ++k) {
      HV.col(k) =
          BuildCoupledSigmaVector(V.col(k), C_alpha, nocclevels_alpha, C_beta,
                                  nocclevels_beta, coupled_fock_builder);
    }
    return HV;
  };
  const trah::StepResult step = trah::Step(g, diag_h, hessian, trust_radius,
                                           trah_options_, &trah_subspace_);
  predicted_energy_change = step.predicted;
  last_step_linear_ = step.linear;
  last_step_norm_ = step.x.norm();

  XTP_LOG(Log::warning, *log_)
      << TimeStamp() << " TRAH step: |g|=" << g.norm()
      << ", |x|=" << step.x.norm() << ", trust radius=" << trust_radius
      << ", alpha=" << step.alpha << ", micro-iterations=" << step.iterations
      << ", Hessian products=" << step.products
      << (analytic ? " (analytic)" : " (finite differences)")
      << (step.reused ? " (subspace reused)" : "")
      << ", residual=" << step.residual
      << (step.converged ? "" : " (not converged)")
      << ", predicted dE=" << predicted_energy_change << std::flush;

  return {
      trah::RotateOrbitals(C_alpha, step.x.head(n_ov_alpha), nocclevels_alpha),
      trah::RotateOrbitals(C_beta, step.x.tail(n_ov_beta), nocclevels_beta)};
}

Eigen::MatrixXd UKSConvergenceAcc::DensityMatrixGroundState_unres(
    const Eigen::MatrixXd& MOs, Index nocclevels) const {
  if (nocclevels == 0) {
    return Eigen::MatrixXd::Zero(MOs.rows(), MOs.rows());
  }
  Eigen::MatrixXd occstates = MOs.leftCols(nocclevels);
  return occstates * occstates.transpose();
}

UKSConvergenceAcc::SpinDensity UKSConvergenceAcc::DensityMatrix(
    const tools::EigenSystem& MOs_alpha,
    const tools::EigenSystem& MOs_beta) const {
  SpinDensity result;
  result.alpha = DensityMatrixGroundState_unres(MOs_alpha.eigenvectors(),
                                                nocclevels_alpha_);
  result.beta =
      DensityMatrixGroundState_unres(MOs_beta.eigenvectors(), nocclevels_beta_);
  return result;
}

void UKSConvergenceAcc::Levelshift(Eigen::MatrixXd& H,
                                   const Eigen::MatrixXd& MOs_old,
                                   const options& opt, Index nocclevels) const {
  if (opt.levelshift < 1e-9) {
    return;
  }
  Eigen::VectorXd virt = Eigen::VectorXd::Zero(H.rows());
  for (Index i = nocclevels; i < H.rows(); ++i) {
    virt(i) = opt.levelshift;
  }

  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Using levelshift:" << opt.levelshift << " Hartree"
      << std::flush;

  Eigen::MatrixXd vir = S_->Matrix() * MOs_old * virt.asDiagonal() *
                        MOs_old.transpose() * S_->Matrix();
  H += vir;
}

Eigen::MatrixXd UKSConvergenceAcc::BuildErrorMatrix(
    const Eigen::MatrixXd& dmat, const Eigen::MatrixXd& H) const {
  const Eigen::MatrixXd& S = S_->Matrix();
  return Sminusahalf.transpose() * (H * dmat * S - S * dmat * H) * Sminusahalf;
}

double UKSConvergenceAcc::CombinedError(const Eigen::MatrixXd& err_alpha,
                                        const Eigen::MatrixXd& err_beta) const {
  return std::max(err_alpha.cwiseAbs().maxCoeff(),
                  err_beta.cwiseAbs().maxCoeff());
}

UKSConvergenceAcc::SpinDensity UKSConvergenceAcc::Iterate(
    const SpinDensity& dmat, SpinFock& H, tools::EigenSystem& MOs_alpha,
    tools::EigenSystem& MOs_beta, double totE) {

  // Check the second-order step returned by the previous call, now that
  // the energy of its density (totE) is known: Fletcher's trust-radius
  // update. A rejected step is undone and, while second-order steps
  // continue, repeated at once with the smaller radius from the stored
  // Fock matrices and Davidson subspace (see ConvergenceAcc::Iterate).
  last_rejected_ = false;
  if (direct_min_pending_) {
    direct_min_pending_ = false;
    const double actual_change = totE - direct_min_pre_energy_;
    const bool accepted =
        trust_.Update(direct_min_predicted_change_, actual_change,
                      direct_min_linear_, direct_min_step_norm_);
    XTP_LOG(Log::warning, *log_)
        << TimeStamp()
        << " Second-order step check: actual dE=" << actual_change
        << ", predicted dE=" << direct_min_predicted_change_
        << ", r=" << trust_.LastRatio()
        << (accepted ? " (accepted)" : " (rejected)")
        << ", new trust radius=" << trust_.Radius() << std::flush;
    if (!accepted) {
      if (trust_.Radius() < kMinTrustRadius) {
        direct_min_floor_hit_ = true;
        trah_active_ = false;
        XTP_LOG(Log::warning, *log_)
            << TimeStamp() << " Trust radius below " << kMinTrustRadius
            << " without an accepted step; leaving second-order steps for "
               "this SCF"
            << std::flush;
      }
      MOs_alpha.eigenvectors() = direct_min_pre_MOs_alpha_;
      MOs_beta.eigenvectors() = direct_min_pre_MOs_beta_;
      MOs_alpha.eigenvalues() = direct_min_pre_MOs_alpha_energies_;
      MOs_beta.eigenvalues() = direct_min_pre_MOs_beta_energies_;
      totE_.push_back(direct_min_pre_energy_);
      last_rejected_ = true;
      usedmixing_ = false;
      if (trah_active_) {
        const SpinFock H_pre = direct_min_pre_H_;
        return SecondOrderStep(DensityMatrix(MOs_alpha, MOs_beta), H_pre,
                               MOs_alpha, MOs_beta, direct_min_pre_energy_);
      }
      return DensityMatrix(MOs_alpha, MOs_beta);
    }
  }

  // History, trimmed at ONE index for the Fock, density and error
  // histories and for DIIS's own error history: the oldest entry, or with
  // DIIS_maxout the entry with the largest error. Previously the index was
  // a stale position (never re-based after an erase, maxerror_ never
  // reset), and DIIS::Update could trim at a different index than the
  // Fock history, pairing coefficients with the wrong matrices. totE_ is
  // not trimmed: it is the chronological energy list getDeltaE reads.
  Index drop = 0;
  if (int(mathist_alpha_.size()) == opt_alpha_.histlength) {
    if (opt_alpha_.maxout) {
      drop = Index(std::max_element(errhist_.begin(), errhist_.end()) -
                   errhist_.begin());
    }
    mathist_alpha_.erase(mathist_alpha_.begin() + drop);
    mathist_beta_.erase(mathist_beta_.begin() + drop);
    dmatHist_alpha_.erase(dmatHist_alpha_.begin() + drop);
    dmatHist_beta_.erase(dmatHist_beta_.begin() + drop);
    errhist_.erase(errhist_.begin() + drop);
  }

  totE_.push_back(totE);

  if (nocclevels_alpha_ > 0 &&
      nocclevels_alpha_ < MOs_alpha.eigenvalues().size()) {
    double gap_alpha = MOs_alpha.eigenvalues()(nocclevels_alpha_) -
                       MOs_alpha.eigenvalues()(nocclevels_alpha_ - 1);
    // consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures added
    // deliberately: a real run showed level shift catastrophically
    // breaking the Davidson solve this class's own direct-minimization
    // fallback depends on (residual climbing past 395, far beyond any
    // prior failure mode) once direct-minimization engaged (originally
    // discovered against an earlier, decoupled AugmentedHessianStep
    // implementation; the same mechanism still applies to
    // CoupledAugmentedHessianStep's own diagonal preconditioner below,
    // which is built the same way). Root cause: Levelshift() modifies H
    // in place, BEFORE CoupledAugmentedHessianStep/DirectMinimizationRotation
    // ever see it -- and that preconditioner (diag_h) is built directly
    // from eps(a)-eps(i), the raw orbital energy gap, which becomes
    // systematically inflated by the shift for virtual orbitals.
    // ADIIS/DIIS do not
    // have this problem (they operate on the Fock/density matrices
    // themselves, never deriving a separate quantity like an orbital
    // gap from the shifted eigenvalues), so level shift stays active
    // for them -- this guard only disables it once direct-minimization
    // is about to take over, since level shift's own purpose
    // (discouraging occ-virt mixing during diagonalization) does not
    // even apply to a method that never diagonalizes H at all.
    if ((diiserror_ > opt_alpha_.levelshiftend &&
         opt_alpha_.levelshift > 0.0) ||
        gap_alpha < 1e-6) {
      // "- 1" margin: consecutive_adiis_failures_ is only incremented
      // LATER in this same Iterate() call (after this check runs), so
      // without the margin this would still see the PREVIOUS
      // iteration's count on the exact iteration where the threshold
      // is reached and direct-minimization actually fires -- letting
      // level shift through on precisely the iteration it needs to be
      // blocked, one iteration too late.
      if (consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures - 1 &&
          !trah_active_) {
        Levelshift(H.alpha, MOs_alpha.eigenvectors(), opt_alpha_,
                   nocclevels_alpha_);
      }
    }
  }

  if (nocclevels_beta_ > 0 &&
      nocclevels_beta_ < MOs_beta.eigenvalues().size()) {
    double gap_beta = MOs_beta.eigenvalues()(nocclevels_beta_) -
                      MOs_beta.eigenvalues()(nocclevels_beta_ - 1);
    if ((diiserror_ > opt_beta_.levelshiftend && opt_beta_.levelshift > 0.0) ||
        gap_beta < 1e-6) {
      if (consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures - 1 &&
          !trah_active_) {
        Levelshift(H.beta, MOs_beta.eigenvectors(), opt_beta_,
                   nocclevels_beta_);
      }
    }
  }

  Eigen::MatrixXd err_alpha = BuildErrorMatrix(dmat.alpha, H.alpha);
  Eigen::MatrixXd err_beta = BuildErrorMatrix(dmat.beta, H.beta);

  diiserror_ = CombinedError(err_alpha, err_beta);

  // Trailing-average trigger bookkeeping (see this class's own header
  // comment on diiserror_history_ for the full ORCA-derived reasoning)
  // -- tracked unconditionally, every iteration, regardless of what
  // this iteration goes on to do (ADIIS/DIIS/mixing/direct-
  // minimization), since the whole point is to observe the genuine,
  // realized trajectory of diiserror_ itself.
  ++total_iteration_count_;
  diiserror_history_.push_back(diiserror_);
  if (Index(diiserror_history_.size()) > kTrailingWindowSize) {
    diiserror_history_.erase(diiserror_history_.begin());
  }

  mathist_alpha_.push_back(H.alpha);
  mathist_beta_.push_back(H.beta);
  dmatHist_alpha_.push_back(dmat.alpha);
  dmatHist_beta_.push_back(dmat.beta);

  errhist_.push_back(diiserror_);

  // crucial: one shared error matrix = alpha + beta contribution
  diis_.Update(drop, err_alpha, err_beta);

  bool diis_error = false;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " DIIs error " << diiserror_ << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  // Once second-order steps have started (second_order "auto" after the
  // trigger below, or "always" from the start) they continue until the SCF
  // has converged or the trust radius has collapsed.
  if (!trah_active_ && !direct_min_floor_hit_ &&
      opt_alpha_.second_order != "off") {
    const bool failures =
        consecutive_adiis_failures_ >= kMaxConsecutiveADIISFailures;
    const bool stalled = TrailingAverageStalled();
    if (opt_alpha_.second_order == "always" || failures || stalled) {
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Switching to second-order (TRAH) steps"
          << (failures ? " ((A)DIIS failed repeatedly)" : "")
          << (stalled ? " (DIIS error stalled)" : "") << std::flush;
      trah_active_ = true;
    }
  }
  if (trah_active_) {
    return SecondOrderStep(dmat, H, MOs_alpha, MOs_beta, totE);
  }

  Eigen::MatrixXd H_guess_alpha = H.alpha;
  Eigen::MatrixXd H_guess_beta = H.beta;

  // As in ConvergenceAcc::Iterate: below DIIS_start (e.g. after a warm
  // start) no damping, a plain first step and DIIS from two entries on.
  const bool near_convergence = diiserror_ < opt_alpha_.diis_start;
  const std::size_t min_history = near_convergence ? 1 : 2;
  if ((diiserror_ < opt_alpha_.adiis_start ||
       diiserror_ < opt_alpha_.diis_start) &&
      opt_alpha_.usediis && mathist_alpha_.size() > min_history) {

    Eigen::VectorXd coeffs;

    if (diiserror_ > opt_alpha_.diis_start ||
        totE_.back() - totE_[totE_.size() - 2] > kEnergyRiseForADIIS) {
      coeffs = adiis_.CalcCoeff(dmatHist_alpha_, dmatHist_beta_, mathist_alpha_,
                                mathist_beta_);
      diis_error = !adiis_.Info() || coeffs.size() == 0;
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using ADIIS for next UKS guess" << std::flush;
    } else {
      coeffs = diis_.CalcCoeff();
      diis_error = !diis_.Info() || coeffs.size() == 0;
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " Using DIIS for next UKS guess" << std::flush;
    }

    if (diis_error) {
      ++consecutive_adiis_failures_;
      XTP_LOG(Log::warning, *log_)
          << TimeStamp() << " (A)DIIS failed using mixing instead"
          << std::flush;
      H_guess_alpha = H.alpha;
      H_guess_beta = H.beta;
    } else {
      consecutive_adiis_failures_ = 0;
      H_guess_alpha.setZero();
      H_guess_beta.setZero();
      for (Index i = 0; i < coeffs.size(); ++i) {
        if (std::abs(coeffs(i)) < 1e-8) {
          continue;
        }
        H_guess_alpha += coeffs(i) * mathist_alpha_[i];
        H_guess_beta += coeffs(i) * mathist_beta_[i];
      }
    }
  }

  MOs_alpha = SolveFockmatrix(H_guess_alpha);
  MOs_beta = SolveFockmatrix(H_guess_beta);

  SpinDensity dmatout = DensityMatrix(MOs_alpha, MOs_beta);

  if (diiserror_ > opt_alpha_.mixingend || !opt_alpha_.usediis || diis_error ||
      (mathist_alpha_.size() <= 2 && !near_convergence)) {
    // mixingend, not adiis_start -- deliberately decoupled (see the
    // options struct's own comment in convergenceacc.h): ORCA keeps
    // DampErr fully independent of DIISStart, and recommends making
    // DampErr much SMALLER for difficult systems specifically so
    // damping stays active well past the point where DIIS itself
    // starts being tried -- reusing adiis_start here could never
    // represent that independently, since it would tie "when does
    // mixing turn off" to the same value as "when does ADIIS/DIIS
    // engage at all", which are conceptually separate questions.
    usedmixing_ = true;
    // Adaptive damping (matches ORCA's own DampFac/DampMax design,
    // confirmed directly from a real ORCA log's own resolved SCF
    // settings -- see the options struct's own comment in
    // convergenceacc.h for the full reasoning): ramp LINEARLY from
    // mixingparameter (the base, e.g. 0.7) toward mixingmax (the
    // ceiling, e.g. 0.98) as consecutive_adiis_failures_ increases
    // toward kMaxConsecutiveADIISFailures, rather than applying the
    // ceiling value for the entire run regardless of whether the SCF
    // is actually struggling. Ties the ramp to the SAME signal already
    // driving the direct-minimization trigger itself, rather than
    // introducing a separate struggle metric -- consecutive_adiis_
    // failures_ resets to 0 on any successful ADIIS/DIIS step, so the
    // ramp relaxes back toward the base value just as readily as it
    // climbed.
    double ramp_fraction =
        std::min(1.0, double(consecutive_adiis_failures_) /
                          double(kMaxConsecutiveADIISFailures));
    double mixingparameter_alpha_current =
        opt_alpha_.mixingparameter +
        ramp_fraction * (opt_alpha_.mixingmax - opt_alpha_.mixingparameter);
    double mixingparameter_beta_current =
        opt_beta_.mixingparameter +
        ramp_fraction * (opt_beta_.mixingmax - opt_beta_.mixingparameter);
    dmatout.alpha = mixingparameter_alpha_current * dmat.alpha +
                    (1.0 - mixingparameter_alpha_current) * dmatout.alpha;
    dmatout.beta = mixingparameter_beta_current * dmat.beta +
                   (1.0 - mixingparameter_beta_current) * dmatout.beta;
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " Using coupled UKS mixing with adaptive alpha="
        << mixingparameter_alpha_current
        << " (base=" << opt_alpha_.mixingparameter
        << ", ceiling=" << opt_alpha_.mixingmax
        << ", ramp fraction=" << ramp_fraction << ")" << std::flush;
  } else {
    usedmixing_ = false;
  }

  return dmatout;
}

// After kAutoStartIteration iterations: true if the DIIS error has, on
// average over the last kTrailingWindowSize iterations, shrunk by less than
// a factor kMeanRatioTolerance per iteration (ORCA's AutoTRAH criterion:
// "Auto Start mean grad. ratio tolernc. 1.125"). Evaluated every iteration,
// also when the error stays above ADIIS_start and (A)DIIS is never tried.
bool UKSConvergenceAcc::TrailingAverageStalled() const {
  if (total_iteration_count_ < kAutoStartIteration ||
      Index(diiserror_history_.size()) < kTrailingWindowSize) {
    return false;
  }
  double mean_ratio = 0.0;
  Index ratio_count = 0;
  for (Index i = 1; i < Index(diiserror_history_.size()); ++i) {
    if (diiserror_history_[i - 1] > 1e-12) {
      mean_ratio += diiserror_history_[i] / diiserror_history_[i - 1];
      ++ratio_count;
    }
  }
  return ratio_count > 0 &&
         mean_ratio / double(ratio_count) > 1.0 / kMeanRatioTolerance;
}

UKSConvergenceAcc::SpinDensity UKSConvergenceAcc::SecondOrderStep(
    const SpinDensity& dmat, const SpinFock& H, tools::EigenSystem& MOs_alpha,
    tools::EigenSystem& MOs_beta, double totE) {
  // The step rotates MOs_alpha/MOs_beta, and its gradient and the Fletcher
  // check are taken relative to H and totE, which belong to dmat. If dmat
  // is not the density of these orbitals (after a damped iteration, or a
  // guess density), return their density first; the step follows on the
  // next call, from a consistent starting point.
  SpinDensity dmat_mos = DensityMatrix(MOs_alpha, MOs_beta);
  const double mismatch =
      std::max((dmat_mos.alpha - dmat.alpha).cwiseAbs().maxCoeff(),
               (dmat_mos.beta - dmat.beta).cwiseAbs().maxCoeff());
  if (mismatch > kDensityMismatchTolerance) {
    XTP_LOG(Log::warning, *log_)
        << TimeStamp()
        << " Density differs from that of the current orbitals (max "
           "deviation "
        << mismatch << "); returning the orbital density first" << std::flush;
    usedmixing_ = false;
    return dmat_mos;
  }

  direct_min_pre_energy_ = totE;
  direct_min_pre_MOs_alpha_ = MOs_alpha.eigenvectors();
  direct_min_pre_MOs_beta_ = MOs_beta.eigenvectors();
  direct_min_pre_MOs_alpha_energies_ = MOs_alpha.eigenvalues();
  direct_min_pre_MOs_beta_energies_ = MOs_beta.eigenvalues();
  direct_min_pre_H_ = H;

  Eigen::MatrixXd C_new_alpha;
  Eigen::MatrixXd C_new_beta;
  // without the linear term (diagonal path) a rejection shrinks by 0.7
  direct_min_linear_ = 0.0;
  direct_min_step_norm_ = 0.0;
  if (coupled_fock_builder_ || orbital_hessian_) {
    std::tie(C_new_alpha, C_new_beta) = CoupledAugmentedHessianStep(
        H.alpha, MOs_alpha, nocclevels_alpha_, H.beta, MOs_beta,
        nocclevels_beta_, coupled_fock_builder_, trust_.Radius(),
        direct_min_predicted_change_);
    direct_min_linear_ = last_step_linear_;
    direct_min_step_norm_ = last_step_norm_;
  } else {
    // no Fock builder (e.g. the atomic UKS guess): diagonal Hessian, one
    // channel at a time
    double predicted_alpha = 0.0;
    double predicted_beta = 0.0;
    C_new_alpha = DirectMinimizationRotation(
        H.alpha, MOs_alpha, nocclevels_alpha_, predicted_alpha);
    C_new_beta = DirectMinimizationRotation(H.beta, MOs_beta, nocclevels_beta_,
                                            predicted_beta);
    direct_min_predicted_change_ = predicted_alpha + predicted_beta;
  }
  direct_min_pending_ = true;

  MOs_alpha.eigenvectors() = C_new_alpha;
  MOs_beta.eigenvectors() = C_new_beta;
  // diagonal of the old Fock matrix in the rotated orbitals: for the
  // preconditioner and the gap checks until the next diagonalisation
  MOs_alpha.eigenvalues() =
      (C_new_alpha.transpose() * H.alpha * C_new_alpha).diagonal();
  MOs_beta.eigenvalues() =
      (C_new_beta.transpose() * H.beta * C_new_beta).diagonal();
  usedmixing_ = false;
  return DensityMatrix(MOs_alpha, MOs_beta);
}

double UKSConvergenceAcc::getDeltaE() const {
  if (totE_.size() < 2) {
    return 100.0;
  }
  return std::abs(totE_.back() - totE_[totE_.size() - 2]);
}

bool UKSConvergenceAcc::isConverged() const {
  // not on the call that undid a rejected second-order step (see
  // ConvergenceAcc::isConverged)
  return (!last_rejected_ && getDeltaE() < opt_alpha_.Econverged &&
          diiserror_ < opt_alpha_.error_converged);
}

}  // namespace xtp
}  // namespace votca