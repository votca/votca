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
  // adiis_.setHistLength(opt_alpha_.histlength);
}

void UKSConvergenceAcc::setLogger(Logger* log) { log_ = log; }

void UKSConvergenceAcc::setOverlap(AOOverlap& S, double etol) {
  S_ = &S;
  Sminusahalf = S.Pseudo_InvSqrt(etol);
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Smallest value of AOOverlap matrix is "
      << S_->SmallestEigenValue() << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Removed " << S_->Removedfunctions()
      << " basisfunction from inverse overlap matrix" << std::flush;
}

tools::EigenSystem UKSConvergenceAcc::SolveFockmatrix(
    const Eigen::MatrixXd& H) const {
  Eigen::MatrixXd H_ortho = Sminusahalf.transpose() * H * Sminusahalf;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H_ortho);

  if (es.info() != Eigen::ComputationInfo::Success) {
    throw std::runtime_error("Matrix Diagonalisation failed. DiagInfo" +
                             std::to_string(es.info()));
  }

  tools::EigenSystem result;
  result.eigenvalues() = es.eigenvalues();
  result.eigenvectors() = Sminusahalf * es.eigenvectors();
  return result;
}

Eigen::MatrixXd UKSConvergenceAcc::UnflattenRotation(const Eigen::VectorXd& v_ov,
                                                     Index nao,
                                                     Index nocclevels) const {
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

Eigen::VectorXd UKSConvergenceAcc::BuildSigmaVector(
    const Eigen::VectorXd& v_ov, const Eigen::MatrixXd& C, Index nocclevels,
    const FockBuilder& fock_builder, const Eigen::MatrixXd& g_occ_virt,
    double finite_diff_step) const {
  Index nao = C.rows();
  Index nvirt = nao - nocclevels;
  (void)g_occ_virt;  // no longer needed for a CENTRAL difference -- kept
                    // in the signature for interface stability with
                    // AugmentedHessianOperator's own construction.

  // CENTRAL, not one-sided/forward, finite-difference Hessian-vector
  // product -- confirmed necessary, not just a nicety, by this
  // class's own symmetry check (u.(H*v) vs v.(H*u), see the
  // conversation this grew out of): a forward difference's leading
  // error term is proportional to the THIRD derivative of the energy,
  // which is not symmetric between u and v the way the true Hessian
  // is, and empirically this asymmetry was found to be large (a
  // complete, ~100% relative mismatch), not a small correction on an
  // otherwise-good answer. A central difference cancels this leading,
  // odd-order error term exactly, at the cost of one extra Fock build
  // per sigma-vector evaluation (two perturbed gradient evaluations
  // instead of one, since the UNPERTURBED gradient g_occ_virt is no
  // longer needed at all for this formula).
  //
  // 1e-3, not the earlier 1e-2 -- though confirmed directly NOT to be
  // the explanation for a real, observed "stuck" pattern (predicted/
  // actual dE frozen identically across a dozen+ consecutive trust-
  // radius shrinks): a real run at 1e-2 and another at 1e-3 both got
  // stuck at the EXACT SAME trust radius (0.000465261027974) to many
  // decimal places, which a genuine step-size-floor explanation would
  // not produce (a 10x smaller step should have moved that threshold
  // by a corresponding amount). Kept at the smaller, still-reasonable
  // 1e-3 regardless, since there is no remaining reason to prefer the
  // larger value -- but the actual "stuck" mechanism is now believed
  // to be AugmentedHessianStep's own alpha_max=1000 bisection ceiling
  // saturating, not this step size; see that function's own comment.
  //
  // Now a parameter, not a hardcoded constant: needed to run the same
  // sigma-vector evaluation at two different step sizes and compare
  // them directly, testing whether finite-difference truncation vs.
  // rounding-error noise is the dominant error source here -- see the
  // conversation this grew out of. Every existing call site continues
  // to pass 1e-3 explicitly (via the default argument), so behavior is
  // unchanged unless a caller deliberately requests a different value.
  Eigen::MatrixXd kappa_trial =
      finite_diff_step * UnflattenRotation(v_ov, nao, nocclevels);

  // Same linearized-exponential + Lowdin-reorthonormalization approach
  // as DirectMinimizationRotation's own orbital update -- valid here
  // for the identical reason: finite_diff_step keeps kappa_trial small
  // regardless of how large v_ov itself is (v_ov is a Davidson trial
  // vector, not itself trust-radius bounded at this stage).
  auto EvaluateGradientAt =
      [&](const Eigen::MatrixXd& kappa) -> Eigen::MatrixXd {
    Eigen::MatrixXd C_rot = C * (Eigen::MatrixXd::Identity(nao, nao) + kappa);
    Eigen::MatrixXd nonortho = C_rot.transpose() * S_->Matrix() * C_rot;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_ortho(nonortho);
    C_rot = C_rot * es_ortho.operatorInverseSqrt();

    Eigen::MatrixXd C_occ_rot = C_rot.leftCols(nocclevels);
    Eigen::MatrixXd D_rot = C_occ_rot * C_occ_rot.transpose();

    Eigen::MatrixXd F_AO_rot = fock_builder(D_rot);
    return C_rot.transpose() * F_AO_rot * C_rot;
  };

  Eigen::MatrixXd F_MO_plus = EvaluateGradientAt(kappa_trial);
  Eigen::MatrixXd F_MO_minus = EvaluateGradientAt(-kappa_trial);

  Eigen::VectorXd sigma(nocclevels * nvirt);
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      double g_plus_ia = F_MO_plus(i, nocclevels + a);
      double g_minus_ia = F_MO_minus(i, nocclevels + a);
      sigma(i * nvirt + a) = (g_plus_ia - g_minus_ia) / (2.0 * finite_diff_step);
    }
  }
  return sigma;
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
  auto EvaluateBothGradientsAt =
      [&](const Eigen::MatrixXd& kappa_alpha, const Eigen::MatrixXd& kappa_beta)
      -> std::pair<Eigen::MatrixXd, Eigen::MatrixXd> {
    Eigen::MatrixXd C_alpha_rot =
        C_alpha * (Eigen::MatrixXd::Identity(nao_alpha, nao_alpha) + kappa_alpha);
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

  Eigen::VectorXd sigma(n_ov_alpha + n_ov_beta);
  for (Index i = 0; i < nocclevels_alpha; ++i) {
    for (Index a = 0; a < nvirt_alpha; ++a) {
      double g_plus_ia = F_MO_alpha_plus(i, nocclevels_alpha + a);
      double g_minus_ia = F_MO_alpha_minus(i, nocclevels_alpha + a);
      sigma(i * nvirt_alpha + a) =
          (g_plus_ia - g_minus_ia) / (2.0 * finite_diff_step);
    }
  }
  for (Index i = 0; i < nocclevels_beta; ++i) {
    for (Index a = 0; a < nvirt_beta; ++a) {
      double g_plus_ia = F_MO_beta_plus(i, nocclevels_beta + a);
      double g_minus_ia = F_MO_beta_minus(i, nocclevels_beta + a);
      sigma(n_ov_alpha + i * nvirt_beta + a) =
          (g_plus_ia - g_minus_ia) / (2.0 * finite_diff_step);
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
      // Approximate, diagonal Hessian: orbital energy differences --
      // the same cheap starting approximation used by most quasi-
      // Newton SCF methods (e.g. SOSCF's own initial Hessian guess).
      double h_ia = 2.0 * gap;
      double kappa_ia = -F_MO(i, a) / h_ia;
      kappa_ia = std::clamp(kappa_ia, -kMaxKappaElement, kMaxKappaElement);
      kappa(i, a) = kappa_ia;
      kappa(a, i) = -kappa_ia;
      h_matrix(i, a) = h_ia;
    }
  }

  double knorm = kappa.norm();
  if (knorm > trust_radius_current_) {
    kappa *= (trust_radius_current_ / knorm);
  }

  // Quadratic model's own predicted energy change, Sum_ia[g_ia*kappa_ia
  // + 0.5*h_ia*kappa_ia^2], using the FINAL kappa (after both the
  // per-element clamp and the whole-matrix trust-radius scaling above)
  // -- needed by Iterate's own Fletcher-style accept/reject logic.
  // g_ia = F_MO(i,a) directly, matching the sign convention kappa_ia =
  // -F_MO(i,a)/h_ia was itself built from above.
  predicted_energy_change = 0.0;
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = nocclevels; a < nao; ++a) {
      double kappa_ia = kappa(i, a);
      predicted_energy_change +=
          F_MO(i, a) * kappa_ia + 0.5 * h_matrix(i, a) * kappa_ia * kappa_ia;
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

namespace {
// Local operator matching DavidsonSolver's own MatrixReplacement
// template interface (.rows(), .diagonal(), operator*(MatrixXd)) for
// the scaled augmented Hessian, Helmich-Paris Eq. 9:
//   [[0, alpha*g^T], [alpha*g, H]]
// H*v_ov itself is never built explicitly -- only ever applied to a
// trial vector, via UKSConvergenceAcc::BuildSigmaVector's own finite-
// difference approximation (see that function's own header comment).
struct AugmentedHessianOperator {
  const Eigen::VectorXd& g;
  const Eigen::MatrixXd& C;
  Index nocclevels;
  double alpha_scale;
  const UKSConvergenceAcc::FockBuilder& fock_builder;
  const Eigen::MatrixXd& F_MO;
  const UKSConvergenceAcc* self;
  const Eigen::VectorXd& diag_h;  // approximate diagonal Hessian,
                                  // reused as the Davidson
                                  // preconditioner (Helmich-Paris
                                  // Sec. II B) -- NOT the step itself.

  Index rows() const { return 1 + g.size(); }

  Eigen::VectorXd diagonal() const {
    Eigen::VectorXd d(1 + g.size());
    d(0) = 0.0;
    d.tail(g.size()) = diag_h;
    return d;
  }

  Eigen::MatrixXd operator*(const Eigen::MatrixXd& V) const {
    Eigen::MatrixXd AV = Eigen::MatrixXd::Zero(V.rows(), V.cols());
    for (Index col = 0; col < V.cols(); ++col) {
      double v0 = V(0, col);
      Eigen::VectorXd v_ov = V.block(1, col, g.size(), 1);
      AV(0, col) = alpha_scale * g.dot(v_ov);
      Eigen::VectorXd sigma =
          self->BuildSigmaVector(v_ov, C, nocclevels, fock_builder, F_MO);
      AV.block(1, col, g.size(), 1) = alpha_scale * g * v0 + sigma;
    }
    return AV;
  }
};
}  // namespace

namespace {
// Coupled analogue of AugmentedHessianOperator: matches the SAME
// DavidsonSolver MatrixReplacement interface, but for the augmented
// Hessian over the COMBINED (alpha+beta) rotation space -- g and
// diag_h are the concatenated, both-channel vectors, and each
// operator* call rotates BOTH channels together via
// UKSConvergenceAcc::BuildCoupledSigmaVector, capturing the genuine
// alpha-beta coupling that two independent AugmentedHessianOperator
// solves cannot.
struct CoupledAugmentedHessianOperator {
  const Eigen::VectorXd& g;
  const Eigen::MatrixXd& C_alpha;
  Index nocclevels_alpha;
  const Eigen::MatrixXd& C_beta;
  Index nocclevels_beta;
  double alpha_scale;
  const UKSConvergenceAcc::CoupledFockBuilder& coupled_fock_builder;
  const UKSConvergenceAcc* self;
  const Eigen::VectorXd& diag_h;

  Index rows() const { return 1 + g.size(); }

  Eigen::VectorXd diagonal() const {
    Eigen::VectorXd d(1 + g.size());
    d(0) = 0.0;
    d.tail(g.size()) = diag_h;
    return d;
  }

  Eigen::MatrixXd operator*(const Eigen::MatrixXd& V) const {
    Eigen::MatrixXd AV = Eigen::MatrixXd::Zero(V.rows(), V.cols());
    for (Index col = 0; col < V.cols(); ++col) {
      double v0 = V(0, col);
      Eigen::VectorXd v_ov = V.block(1, col, g.size(), 1);
      AV(0, col) = alpha_scale * g.dot(v_ov);
      Eigen::VectorXd sigma = self->BuildCoupledSigmaVector(
          v_ov, C_alpha, nocclevels_alpha, C_beta, nocclevels_beta,
          coupled_fock_builder);
      AV.block(1, col, g.size(), 1) = alpha_scale * g * v0 + sigma;
    }
    return AV;
  }
};
}  // namespace

Eigen::MatrixXd UKSConvergenceAcc::AugmentedHessianStep(
    const Eigen::MatrixXd& H_AO, const tools::EigenSystem& MOs,
    Index nocclevels, const FockBuilder& fock_builder, double trust_radius,
    double& predicted_energy_change) const {
  Index nao = MOs.eigenvectors().rows();
  Index nvirt = nao - nocclevels;
  Index n_ov = nocclevels * nvirt;
  const Eigen::MatrixXd& C = MOs.eigenvectors();
  const Eigen::VectorXd& eps = MOs.eigenvalues();

  Eigen::MatrixXd F_MO = C.transpose() * H_AO * C;

  Eigen::VectorXd g(n_ov);
  Eigen::VectorXd diag_h(n_ov);
  constexpr double kMinGap = 1e-3;
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      g(i * nvirt + a) = F_MO(i, nocclevels + a);
      double gap = std::max(std::abs(eps(nocclevels + a) - eps(i)), kMinGap);
      diag_h(i * nvirt + a) = 2.0 * gap;
    }
  }

  // Temporary diagnostic (see the conversation this grew out of): a
  // residual norm in the thousands, not decreasing, needs concrete
  // magnitudes to actually locate rather than continued guessing.
  // Probes BuildSigmaVector directly on the SAME normalized direction
  // (g/||g||) used to build the Davidson initial guess below -- a
  // well-behaved Hessian-vector product on a unit-norm input should
  // itself be a reasonably-scaled vector, not thousands in magnitude.
  {
    double gnorm_diag = g.norm();
    Eigen::VectorXd probe_direction;
    if (gnorm_diag > 1e-12) {
      probe_direction = g / gnorm_diag;
    } else {
      probe_direction = Eigen::VectorXd::Unit(n_ov, 0);
    }
    Eigen::VectorXd probe_sigma =
        BuildSigmaVector(probe_direction, C, nocclevels, fock_builder, F_MO);
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " AugmentedHessianStep diagnostic: ||g||="
        << gnorm_diag << ", diag_h range=[" << diag_h.minCoeff() << ", "
        << diag_h.maxCoeff() << "], ||sigma(g/||g||)||=" << probe_sigma.norm()
        << ", max|sigma|=" << probe_sigma.cwiseAbs().maxCoeff() << std::flush;

    // Symmetry check (see the conversation this grew out of): a
    // genuinely symmetric operator must satisfy u.(H*v) == v.(H*u)
    // exactly (up to floating-point noise) for ANY u, v -- the
    // Davidson solver's own "SYMM" matrix-type assumption depends on
    // this holding for whatever BuildSigmaVector actually computes,
    // not just on the TRUE, analytic Hessian being symmetric (a
    // finite-difference approximation to it has no guarantee of
    // inheriting that property automatically). Uses two simple,
    // distinct probe directions (the normalized gradient, and the
    // first unit basis vector) rather than random vectors, so this
    // is exactly reproducible run to run.
    Eigen::VectorXd u = probe_direction;
    Eigen::VectorXd v = Eigen::VectorXd::Unit(n_ov, std::min<Index>(1, n_ov - 1));
    Eigen::VectorXd Hv = BuildSigmaVector(v, C, nocclevels, fock_builder, F_MO);
    Eigen::VectorXd Hu = probe_sigma;
    double u_dot_Hv = u.dot(Hv);
    double v_dot_Hu = v.dot(Hu);
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " AugmentedHessianStep symmetry check: u.(H*v)="
        << u_dot_Hv << ", v.(H*u)=" << v_dot_Hu
        << ", relative difference=" << std::abs(u_dot_Hv - v_dot_Hu) /
                                        std::max(std::abs(u_dot_Hv), 1e-12)
        << std::flush;

    // Step-size (truncation vs. rounding/cancellation noise) check
    // (see the conversation this grew out of): a genuine, external
    // ORCA comparison on this exact system found ORCA's own TRAH
    // (using the EXACT analytic coupled-perturbed response, not a
    // finite-difference approximation) converges where this class's
    // own AugmentedHessianStep does not -- raising the question of
    // whether the finite-difference sigma vector itself is simply too
    // noisy for this specific, difficult system, rather than the
    // deliberate alpha/beta decoupling simplification being the
    // dominant gap. If BuildSigmaVector's own result is genuinely
    // dominated by truncation error (the expected, well-behaved
    // regime), shrinking the step should make consecutive evaluations
    // agree MORE closely, converging toward some fixed answer. If
    // instead the smaller step gives a WORSE, noisier result than the
    // larger one, that specifically implicates rounding/cancellation
    // noise (from fock_builder's own integral-screening tolerance and
    // the XC grid's own finite quadrature accuracy) as dominant --
    // direct, empirical evidence for exactly the general concern that
    // Hessian-vector products via nested finite differences are
    // especially sensitive to numerical noise, since this differences
    // a quantity that is already itself a derivative.
    Eigen::VectorXd sigma_step1 = probe_sigma;  // already computed above
                                                // at the same direction
                                                // and the same default
                                                // (1e-3) step -- reused
                                                // here rather than
                                                // redundantly recomputed,
                                                // saving two Fock builds.
    Eigen::VectorXd sigma_step2 =
        BuildSigmaVector(probe_direction, C, nocclevels, fock_builder, F_MO,
                         1e-4);
    double step_relative_diff =
        (sigma_step1 - sigma_step2).norm() /
        std::max(sigma_step1.norm(), 1e-12);
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " AugmentedHessianStep step-size check: "
           "||sigma(1e-3)||="
        << sigma_step1.norm() << ", ||sigma(1e-4)||=" << sigma_step2.norm()
        << ", relative difference=" << step_relative_diff << std::flush;
  }

  // Bisection over alpha (Helmich-Paris Sec. II B, Eq. 11) to keep the
  // resulting step within the trust radius: ||kappa(alpha)||^2/alpha^2
  // <= trust_radius^2. Bounds match the paper's own default
  // [alpha_min, alpha_max] = [1, 1000] (Table I).
  double alpha_min = 1.0;
  double alpha_max = 1000.0;
  const double kOriginalAlphaMax = alpha_max;  // alpha_max itself gets
                                               // narrowed during the
                                               // bisection loop below,
                                               // so this is needed
                                               // separately to check
                                               // whether alpha_try
                                               // actually approaches
                                               // the ORIGINAL ceiling.
  Eigen::VectorXd best_kappa_flat = Eigen::VectorXd::Zero(n_ov);
  double best_mu = 0.0;

  // The paper's own, deliberately problem-tailored starting vectors
  // (Eq. 12): b0 separates out the pure level-shift direction, b1 the
  // component parallel to the gradient -- NOT the Davidson solver's
  // own generic, diagonal-based default starting guess (which has no
  // reason to already reflect this specific, bordered-matrix
  // structure, including the fixed, exact 0 in its own top-left
  // element).
  Eigen::MatrixXd initial_guess = Eigen::MatrixXd::Zero(1 + n_ov, 2);
  initial_guess(0, 0) = 1.0;
  double gnorm = g.norm();
  if (gnorm > 1e-12) {
    initial_guess.block(1, 1, n_ov, 1) = g / gnorm;
  } else {
    initial_guess(1, 1) = 1.0;
  }

  auto SolveForAlpha = [&](double alpha_try, Eigen::VectorXd& kappa_flat_out,
                          double& mu_out) {
    AugmentedHessianOperator op{g,      C,     nocclevels, alpha_try,
                               fock_builder, F_MO, this,       diag_h};
    DavidsonSolver solver(*log_);
    solver.set_matrix_type("SYMM");
    // Loosened from "normal" (1e-4): each bisection trial's own
    // Davidson solve is just one approximate step inside the outer
    // Fletcher accept/reject loop (which already independently
    // catches a bad step and retries with a smaller trust radius
    // regardless), so it does not need to fully converge on every
    // single one of up to kMaxBisectionIters trials -- a real,
    // reported cost problem confirmed directly: ~1.4s per Davidson
    // solve, times up to 20 bisection trials per spin channel, made a
    // single AugmentedHessianStep call impractically slow.
    solver.set_tolerance("loose");
    // Increased from 16 -- a real run on this system got to within
    // ~10x of the "loose" (1e-3) tolerance (residual 0.0108) by
    // iteration 13 of 15, then ran out of budget before actually
    // crossing it. Set generously higher rather than just enough to
    // barely close that specific gap, so this does not need another
    // round of incremental bumping if a slightly harder case needs a
    // few more iterations than this one did.
    solver.set_iter_max(50);
    // Real, confirmed bug: DavidsonSolver defaults max_search_space_
    // to neigen*5 = 1*5 = 5 whenever it is left unset (its own
    // constructor initializes it to 0, and solve() itself falls back
    // to neigen*5 in that case) -- far too small a subspace for this
    // (1+n_ov)-dimensional problem (761 here), and confirmed directly
    // to be why "Search Space" kept cycling 2->3->4->5->restart
    // without ever accumulating enough information to converge (0.00%
    // converged, every single trial, across every run so far).
    solver.set_max_search_space(40);
    solver.solve(op, 1, initial_guess);
    Eigen::VectorXd eigvec = solver.eigenvectors().col(0);
    mu_out = solver.eigenvalues()(0);
    double v0 = eigvec(0);
    // Guards against a degenerate eigenvector with (near-)zero
    // leading component, for which Eq. 10's own kappa(alpha) =
    // eigvec_tail/v0 normalization is ill-defined.
    if (std::abs(v0) < 1e-8) {
      kappa_flat_out = Eigen::VectorXd::Zero(g.size());
      return;
    }
    kappa_flat_out = eigvec.tail(g.size()) / v0;
  };

  // Start from alpha_min (the gentlest, least aggressive alpha*g
  // coupling), not the midpoint of [alpha_min, alpha_max] -- the paper
  // itself notes alpha ends up AT alpha_min once things are
  // well-behaved (Sec. II B), and starting the search at 500+ (the
  // naive midpoint of [1, 1000]) makes the very first Davidson solve's
  // own off-diagonal coupling alpha*g needlessly large before there is
  // any reason yet to believe that is necessary.
  double alpha_try = alpha_min;
  // Increased from 8 -- a real run showed the bisection failing to
  // actually converge within 8 iterations for at least some trials:
  // requesting trust_radius=0.098 landed on achieved step_norm=0.1374,
  // roughly 40% off target and identical to the PREVIOUS call's own
  // achieved value for trust_radius=0.14 -- i.e. hitting the iteration
  // cap without satisfying its own 1% convergence criterion at all,
  // rather than genuinely tracking the shrinking trust radius. The
  // interval-halving math (2^n shrink factor) says 8 should be enough
  // for a well-behaved, monotonic step_norm(alpha); the real data says
  // otherwise for at least part of this system's own trajectory, so
  // this trades back some of the earlier cost reduction for a bisection
  // that actually reaches its own target.
  constexpr int kMaxBisectionIters = 20;
  for (int bisection_iter = 0; bisection_iter < kMaxBisectionIters;
      ++bisection_iter) {
    Eigen::VectorXd kappa_flat;
    double mu;
    SolveForAlpha(alpha_try, kappa_flat, mu);
    double step_norm = kappa_flat.norm() / alpha_try;
    best_kappa_flat = kappa_flat;
    best_mu = mu;
    if (std::abs(step_norm - trust_radius) < 0.01 * trust_radius) {
      break;
    }
    if (step_norm > trust_radius) {
      // Step too long -- Sec. II B confirms larger alpha shrinks it.
      alpha_min = alpha_try;
    } else {
      alpha_max = alpha_try;
    }
    alpha_try = 0.5 * (alpha_min + alpha_max);
  }

  // Temporary diagnostic (see the conversation this grew out of): a
  // real run showed predicted/actual dE frozen identically across a
  // dozen+ consecutive trust-radius shrinks, at the exact same
  // threshold regardless of BuildSigmaVector's own finite-difference
  // step size -- ruling that step size out as the cause and pointing
  // instead at this bisection's own alpha_max=1000 ceiling possibly
  // saturating (once alpha itself cannot grow any further, the
  // resulting step cannot shrink any further either, no matter how
  // small trust_radius itself becomes). Prints the actual final
  // alpha_try used and how close it sits to alpha_max, to confirm or
  // rule this out directly rather than continue guessing.
  XTP_LOG(Log::warning, *log_)
      << TimeStamp() << " AugmentedHessianStep bisection diagnostic: "
         "final alpha_try="
      << alpha_try << ", original alpha_max ceiling=" << kOriginalAlphaMax
      << ", achieved step_norm=" << (best_kappa_flat.norm() / alpha_try)
      << ", requested trust_radius=" << trust_radius << std::flush;

  Eigen::MatrixXd kappa = UnflattenRotation(best_kappa_flat, nao, nocclevels);

  // Predicted energy change from the SAME quadratic model the
  // augmented-Hessian eigenvalue problem itself is built from --
  // Q(kappa)-E0 = g^T*kappa + 0.5*kappa^T*H*kappa. mu itself already
  // equals g^T*kappa + kappa^T*H*kappa (the level-shifted stationarity
  // condition, Eq. 8, dotted with kappa) at the exact solution, so
  // 0.5*(g^T*kappa + mu*||kappa||^2) is the equivalent, cheaper-to-
  // evaluate form -- avoids needing a further BuildSigmaVector call
  // just to get this number.
  predicted_energy_change =
      0.5 * (g.dot(best_kappa_flat) + best_mu * best_kappa_flat.squaredNorm());

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

  // Combined gradient and diagonal-Hessian preconditioner: alpha's own
  // block first, beta's immediately after -- matching
  // UnflattenCoupledRotation/BuildCoupledSigmaVector's own layout
  // convention. The diagonal preconditioner itself is still built
  // per-channel from each channel's OWN orbital-energy gaps (same
  // formula as AugmentedHessianStep's own diag_h) -- only the actual
  // Hessian-VECTOR product (via BuildCoupledSigmaVector, inside
  // CoupledAugmentedHessianOperator) captures the real cross-channel
  // coupling; the preconditioner is only ever an approximate guide for
  // the Davidson iteration, not the step itself, so this
  // simplification (no explicit alpha-beta cross term in the
  // preconditioner) does not undermine what this whole undertaking is
  // actually meant to fix.
  Eigen::VectorXd g(n_ov);
  Eigen::VectorXd diag_h(n_ov);
  constexpr double kMinGap = 1e-3;
  for (Index i = 0; i < nocclevels_alpha; ++i) {
    for (Index a = 0; a < nvirt_alpha; ++a) {
      g(i * nvirt_alpha + a) = F_MO_alpha(i, nocclevels_alpha + a);
      double gap = std::max(
          std::abs(eps_alpha(nocclevels_alpha + a) - eps_alpha(i)), kMinGap);
      diag_h(i * nvirt_alpha + a) = 2.0 * gap;
    }
  }
  for (Index i = 0; i < nocclevels_beta; ++i) {
    for (Index a = 0; a < nvirt_beta; ++a) {
      g(n_ov_alpha + i * nvirt_beta + a) =
          F_MO_beta(i, nocclevels_beta + a);
      double gap = std::max(
          std::abs(eps_beta(nocclevels_beta + a) - eps_beta(i)), kMinGap);
      diag_h(n_ov_alpha + i * nvirt_beta + a) = 2.0 * gap;
    }
  }

  double alpha_min = 1.0;
  double alpha_max = 1000.0;
  Eigen::VectorXd best_kappa_flat = Eigen::VectorXd::Zero(n_ov);
  double best_mu = 0.0;

  Eigen::MatrixXd initial_guess = Eigen::MatrixXd::Zero(1 + n_ov, 2);
  initial_guess(0, 0) = 1.0;
  double gnorm = g.norm();
  if (gnorm > 1e-12) {
    initial_guess.block(1, 1, n_ov, 1) = g / gnorm;
  } else {
    initial_guess(1, 1) = 1.0;
  }

  auto SolveForAlpha = [&](double alpha_try, Eigen::VectorXd& kappa_flat_out,
                          double& mu_out) {
    CoupledAugmentedHessianOperator op{
        g, C_alpha, nocclevels_alpha, C_beta, nocclevels_beta,
        alpha_try, coupled_fock_builder, this, diag_h};
    DavidsonSolver solver(*log_);
    solver.set_matrix_type("SYMM");
    solver.set_tolerance("loose");
    solver.set_iter_max(50);
    solver.set_max_search_space(40);
    solver.solve(op, 1, initial_guess);
    Eigen::VectorXd eigvec = solver.eigenvectors().col(0);
    mu_out = solver.eigenvalues()(0);
    double v0 = eigvec(0);
    if (std::abs(v0) < 1e-8) {
      kappa_flat_out = Eigen::VectorXd::Zero(g.size());
      return;
    }
    kappa_flat_out = eigvec.tail(g.size()) / v0;
  };

  double alpha_try = alpha_min;
  constexpr int kMaxBisectionIters = 20;
  for (int bisection_iter = 0; bisection_iter < kMaxBisectionIters;
      ++bisection_iter) {
    Eigen::VectorXd kappa_flat;
    double mu;
    SolveForAlpha(alpha_try, kappa_flat, mu);
    double step_norm = kappa_flat.norm() / alpha_try;
    best_kappa_flat = kappa_flat;
    best_mu = mu;
    if (std::abs(step_norm - trust_radius) < 0.01 * trust_radius) {
      break;
    }
    if (step_norm > trust_radius) {
      alpha_min = alpha_try;
    } else {
      alpha_max = alpha_try;
    }
    alpha_try = 0.5 * (alpha_min + alpha_max);
  }

  XTP_LOG(Log::warning, *log_)
      << TimeStamp() << " CoupledAugmentedHessianStep bisection diagnostic: "
         "final alpha_try="
      << alpha_try << ", achieved step_norm="
      << (best_kappa_flat.norm() / alpha_try)
      << ", requested trust_radius=" << trust_radius << std::flush;

  auto [kappa_alpha, kappa_beta] = UnflattenCoupledRotation(
      best_kappa_flat, nao_alpha, nocclevels_alpha, nao_beta, nocclevels_beta);

  // ONE, combined predicted energy change for the whole, coupled step
  // -- same formula as AugmentedHessianStep's own (Q(kappa)-E0 =
  // 0.5*(g^T*kappa + mu*||kappa||^2)), but now naturally a single
  // number for both channels together, rather than needing to be
  // summed from two separate calls the way Iterate's own decoupled
  // AugmentedHessianStep path currently does.
  predicted_energy_change =
      0.5 * (g.dot(best_kappa_flat) + best_mu * best_kappa_flat.squaredNorm());

  Eigen::MatrixXd C_alpha_new =
      C_alpha * (Eigen::MatrixXd::Identity(nao_alpha, nao_alpha) + kappa_alpha);
  Eigen::MatrixXd nonortho_alpha =
      C_alpha_new.transpose() * S_->Matrix() * C_alpha_new;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_alpha(nonortho_alpha);
  C_alpha_new = C_alpha_new * es_alpha.operatorInverseSqrt();

  Eigen::MatrixXd C_beta_new =
      C_beta * (Eigen::MatrixXd::Identity(nao_beta, nao_beta) + kappa_beta);
  Eigen::MatrixXd nonortho_beta =
      C_beta_new.transpose() * S_->Matrix() * C_beta_new;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_beta(nonortho_beta);
  C_beta_new = C_beta_new * es_beta.operatorInverseSqrt();

  return {C_alpha_new, C_beta_new};
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

  // Fletcher's trust-radius update (Helmich-Paris, J. Chem. Phys. 154,
  // 164104 (2021), Sec. II D -- confirmed directly by reading that
  // paper, not reconstructed from memory): verify whatever
  // DirectMinimizationRotation step was taken last call, now that its
  // actual effect on the energy (totE, just passed in -- computed by
  // the caller from a real Fock build on that step's own density) is
  // finally available. This CANNOT be checked within the same
  // Iterate() call that took the step, since that call has no way to
  // know what energy its own returned density will produce until the
  // caller has built a new Fock matrix from it and come back around.
  if (direct_min_pending_) {
    double actual_change = totE - direct_min_pre_energy_;
    double r = (std::abs(direct_min_predicted_change_) > 1e-14)
                  ? actual_change / direct_min_predicted_change_
                  : -1.0;  // treat a degenerate (~zero) predicted change
                          // as an outright reject, same as r<0 below --
                          // the model gave no useful information about
                          // this step at all.
    XTP_LOG(Log::warning, *log_)
        << TimeStamp() << " Direct-minimization step check: actual dE="
        << actual_change << ", predicted dE=" << direct_min_predicted_change_
        << ", r=" << r << ", trust radius=" << trust_radius_current_
        << std::flush;
    if (r < 0.0) {
      // Reject: the quadratic model was not applicable within the
      // given trust region (either the energy rose while predicted to
      // fall, or vice versa). Revert to the pre-step MOs/energy and
      // shrink the trust radius -- the NEXT call's own
      // consecutive_adiis_failures_ check will naturally retry
      // DirectMinimizationRotation from this reverted point with the
      // smaller radius, since nothing here has changed the underlying
      // (A)DIIS behavior that triggered it in the first place.
      trust_radius_current_ *= 0.7;
      // Floor tied to BuildSigmaVector's own finite-difference step
      // (kFiniteDiffStep = 1e-3, defined there): a real run showed
      // trust_radius shrinking past 1e-6 while the step, predicted
      // change, and actual change stayed EXACTLY identical every
      // time -- the bisection's own alpha_min=1 floor means the
      // gentlest achievable step cannot shrink further once alpha_min
      // itself is the binding constraint, so continuing to request an
      // even smaller trust radius changes nothing and the reject loop
      // can never resolve on its own (confirmed directly: it only
      // ended when the outer SCF's own 100-iteration budget ran out).
      // More fundamentally, a trust radius below the sigma vector's
      // own probing resolution is asking for precision the underlying
      // finite-difference model was never built to provide -- its
      // accuracy does not improve as the requested step shrinks, only
      // the requested step size does. Once hit, give up on direct
      // minimization for this SCF call rather than loop toward
      // ever-smaller radii that cannot change the outcome.
      if (trust_radius_current_ < kMinTrustRadius) {
        direct_min_floor_hit_ = true;
        XTP_LOG(Log::warning, *log_)
            << TimeStamp() << " Direct-minimization trust radius fell "
               "below its own finite-difference resolution floor ("
            << kMinTrustRadius << ") without an accepted step -- "
               "falling back to mixing instead of continuing to shrink."
            << std::flush;
      }
      MOs_alpha.eigenvectors() = direct_min_pre_MOs_alpha_;
      MOs_beta.eigenvectors() = direct_min_pre_MOs_beta_;
      MOs_alpha.eigenvalues() = direct_min_pre_MOs_alpha_energies_;
      MOs_beta.eigenvalues() = direct_min_pre_MOs_beta_energies_;
      totE_.push_back(direct_min_pre_energy_);
      direct_min_pending_ = false;
      usedmixing_ = false;
      return DensityMatrix(MOs_alpha, MOs_beta);
    } else if (r <= 0.25) {
      // Accepted, but the step was too long -- shrink for next time.
      trust_radius_current_ *= 0.7;
    } else if (r > 0.75) {
      // Accepted, and the model was a good fit -- grow for next time.
      trust_radius_current_ *= 1.2;
    }
    // 0.25 < r <= 0.75: accepted, trust radius left unchanged.
    direct_min_pending_ = false;
  }

  if (int(mathist_alpha_.size()) == opt_alpha_.histlength) {
    totE_.erase(totE_.begin() + maxerrorindex_);
    mathist_alpha_.erase(mathist_alpha_.begin() + maxerrorindex_);
    mathist_beta_.erase(mathist_beta_.begin() + maxerrorindex_);
    dmatHist_alpha_.erase(dmatHist_alpha_.begin() + maxerrorindex_);
    dmatHist_beta_.erase(dmatHist_beta_.begin() + maxerrorindex_);
  }

  totE_.push_back(totE);

  if (nocclevels_alpha_ > 0 &&
      nocclevels_alpha_ < MOs_alpha.eigenvalues().size()) {
    double gap_alpha = MOs_alpha.eigenvalues()(nocclevels_alpha_) -
                       MOs_alpha.eigenvalues()(nocclevels_alpha_ - 1);
    // consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures added
    // deliberately: a real run showed level shift catastrophically
    // breaking AugmentedHessianStep's own Davidson solve (residual
    // climbing past 395, far beyond any prior failure mode) once
    // direct-minimization engaged. Root cause: Levelshift() modifies H
    // in place, BEFORE AugmentedHessianStep/DirectMinimizationRotation
    // ever see it -- and AugmentedHessianStep's own diagonal-Hessian
    // preconditioner (diag_h) is built directly from eps(a)-eps(i),
    // the raw orbital energy gap, which becomes systematically
    // inflated by the shift for virtual orbitals. ADIIS/DIIS do not
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
      if (consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures - 1) {
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
      if (consecutive_adiis_failures_ < kMaxConsecutiveADIISFailures - 1) {
        Levelshift(H.beta, MOs_beta.eigenvectors(), opt_beta_,
                  nocclevels_beta_);
      }
    }
  }

  Eigen::MatrixXd err_alpha = BuildErrorMatrix(dmat.alpha, H.alpha);
  Eigen::MatrixXd err_beta = BuildErrorMatrix(dmat.beta, H.beta);

  diiserror_ = CombinedError(err_alpha, err_beta);

  mathist_alpha_.push_back(H.alpha);
  mathist_beta_.push_back(H.beta);
  dmatHist_alpha_.push_back(dmat.alpha);
  dmatHist_beta_.push_back(dmat.beta);

  if (opt_alpha_.maxout) {
    if (diiserror_ > maxerror_) {
      maxerror_ = diiserror_;
      maxerrorindex_ = mathist_alpha_.size() - 1;
    }
  } else {
    maxerrorindex_ = 0;
  }

  // crucial: one shared error matrix = alpha + beta contribution
  diis_.Update(maxerrorindex_, err_alpha, err_beta);

  bool diis_error = false;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " DIIs error " << diiserror_ << std::flush;
  XTP_LOG(Log::error, *log_)
      << TimeStamp() << " Delta Etot " << getDeltaE() << std::flush;

  Eigen::MatrixXd H_guess_alpha = H.alpha;
  Eigen::MatrixXd H_guess_beta = H.beta;

  if ((diiserror_ < opt_alpha_.adiis_start ||
       diiserror_ < opt_alpha_.diis_start) &&
      opt_alpha_.usediis && mathist_alpha_.size() > 2) {

    Eigen::VectorXd coeffs;

    if (diiserror_ > opt_alpha_.diis_start ||
        totE_.back() > 0.9 * totE_[totE_.size() - 2]) {
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
      if (consecutive_adiis_failures_ >= kMaxConsecutiveADIISFailures &&
          !direct_min_floor_hit_) {
        // Mirrors ORCA's own auto-TRAH trigger: after (A)DIIS has
        // visibly, repeatedly failed rather than just being slow, fall
        // back to a direct-minimization step instead of plain mixing --
        // see this class's own header comment on DirectMinimizationRotation
        // for the full reasoning and the ORCA log this was validated
        // against directly.
        XTP_LOG(Log::warning, *log_)
            << TimeStamp() << " (A)DIIS failed " << consecutive_adiis_failures_
            << " times in a row, switching to direct-minimization step"
            << std::flush;
        // Save the pre-step state so this step's actual effect can be
        // verified (and, if necessary, reverted) once its own energy
        // becomes available on the NEXT Iterate() call -- see the
        // Fletcher accept/reject check at the top of this function.
        direct_min_pre_energy_ = totE;
        direct_min_pre_MOs_alpha_ = MOs_alpha.eigenvectors();
        direct_min_pre_MOs_beta_ = MOs_beta.eigenvectors();
        direct_min_pre_MOs_alpha_energies_ = MOs_alpha.eigenvalues();
        direct_min_pre_MOs_beta_energies_ = MOs_beta.eigenvalues();

        double predicted_change_alpha = 0.0;
        double predicted_change_beta = 0.0;
        Eigen::MatrixXd C_new_alpha;
        Eigen::MatrixXd C_new_beta;
        // AugmentedHessianStep needs a Fock-builder callback (injected
        // by DFTEngine via setFockBuilderAlpha/Beta) to evaluate its
        // own finite-difference sigma vectors -- fall back to the
        // simpler, diagonal-Hessian-only DirectMinimizationRotation for
        // any caller that has not wired this up, rather than failing
        // outright.
        if (fock_builder_alpha_ && fock_builder_beta_) {
          C_new_alpha = AugmentedHessianStep(
              H.alpha, MOs_alpha, nocclevels_alpha_, fock_builder_alpha_,
              trust_radius_current_, predicted_change_alpha);
          C_new_beta = AugmentedHessianStep(
              H.beta, MOs_beta, nocclevels_beta_, fock_builder_beta_,
              trust_radius_current_, predicted_change_beta);
        } else {
          C_new_alpha = DirectMinimizationRotation(
              H.alpha, MOs_alpha, nocclevels_alpha_, predicted_change_alpha);
          C_new_beta = DirectMinimizationRotation(
              H.beta, MOs_beta, nocclevels_beta_, predicted_change_beta);
        }
        direct_min_predicted_change_ =
            predicted_change_alpha + predicted_change_beta;
        direct_min_pending_ = true;

        MOs_alpha.eigenvectors() = C_new_alpha;
        MOs_beta.eigenvectors() = C_new_beta;
        // Approximate orbital energies from the (only approximately
        // diagonal, post-rotation) MO-basis Fock matrix -- consistent
        // with the rotated orbitals themselves, and cheap to obtain;
        // used only for level-shift gap checks and diagnostics until
        // the next full diagonalization naturally supersedes them.
        MOs_alpha.eigenvalues() =
            (C_new_alpha.transpose() * H.alpha * C_new_alpha).diagonal();
        MOs_beta.eigenvalues() =
            (C_new_beta.transpose() * H.beta * C_new_beta).diagonal();

        SpinDensity dmatout_direct = DensityMatrix(MOs_alpha, MOs_beta);
        usedmixing_ = false;
        return dmatout_direct;
      }
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

  if (diiserror_ > opt_alpha_.mixingend || !opt_alpha_.usediis ||
      diis_error || mathist_alpha_.size() <= 2) {
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
    dmatout.alpha = opt_alpha_.mixingparameter * dmat.alpha +
                    (1.0 - opt_alpha_.mixingparameter) * dmatout.alpha;
    dmatout.beta = opt_beta_.mixingparameter * dmat.beta +
                   (1.0 - opt_beta_.mixingparameter) * dmatout.beta;
    XTP_LOG(Log::warning, *log_)
        << TimeStamp()
        << " Using coupled UKS mixing with alpha=" << opt_alpha_.mixingparameter
        << std::flush;
  } else {
    usedmixing_ = false;
  }

  return dmatout;
}

double UKSConvergenceAcc::getDeltaE() const {
  if (totE_.size() < 2) {
    return 100.0;
  }
  return std::abs(totE_.back() - totE_[totE_.size() - 2]);
}

bool UKSConvergenceAcc::isConverged() const {
  return (getDeltaE() < opt_alpha_.Econverged &&
          diiserror_ < opt_alpha_.error_converged);
}

}  // namespace xtp
}  // namespace votca