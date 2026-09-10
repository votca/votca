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

#ifndef VOTCA_XTP_EWALDBACKGROUND_H
#define VOTCA_XTP_EWALDBACKGROUND_H

// Standard includes
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>

// Third party includes
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>

// Local VOTCA includes
#include "votca/tools/constants.h"
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/ewaldblockjacobipreconditioner.h"
#include "votca/xtp/ewaldperiodicdipoleoperator.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/ewaldshapecorrection.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/segmentmapper.h"

/**
 * \brief ewdbgpol-equivalent calculator, built entirely on the new
 *        Ewald machinery (EwaldRegistry, EwaldRealSpaceSum,
 *        EwaldReciprocalSpaceSum, EwaldShapeCorrection,
 *        EwaldPeriodicDipoleOperator) developed alongside this class,
 *        rather than the legacy xtp/ewald/ code -- registered under a
 *        different Identify() ("ewaldbackground") so the two can coexist
 *        and be compared directly.
 *
 * Converges the self-consistent, fully polarized (neutral-state)
 * background of the whole periodic system, and writes it to an HDF5
 * checkpoint file via EwaldRegistry::WriteToCpt.
 *
 * User-facing length-related options (coulombmethod.alpha,
 * coulombmethod.k_max, realspace.r_min) are specified in nm/nm^-1, not
 * this class's own internal bohr/bohr^-1 units -- converted once in
 * ParseOptions, immediately on read. This is deliberate: everything
 * downstream of ParseOptions (member variables, logging, the actual
 * EwaldRealSpaceSum/EwaldReciprocalSpaceSum construction) stays in the
 * same bohr-based units PolarSite/EwaldRegistry already use natively, so
 * only the options-parsing boundary itself needs to know about the
 * user-facing unit choice.
 *
 * Segments may have any number of sites each -- this now matches
 * EwaldPeriodicDipoleOperator's own generalized (offset-table) vector
 * layout directly, built the same way here for the same reason: an
 * earlier version of both this class and the operator required exactly
 * one site per segment, and that limitation was lifted from the operator
 * first, validated there in isolation (see
 * test_ewaldperiodicdipoleoperatormultisite.cc), before being carried
 * through here -- so that a bug in this generalization, if there were
 * one, would be easy to isolate from the operator's own already-tested
 * multi-site behavior.
 *
 * The outer "loop until 2nd-order fields converged" structure the legacy
 * PolarBackground::Polarize() itself uses (see the project notes on that
 * class) is not reproduced here at all: EwaldPeriodicDipoleOperator's own
 * ConjugateGradient solve already *is* that convergence loop, in one
 * call, rather than a separate SOR iteration this class would need to
 * drive itself.
 */

namespace votca {
namespace xtp {

class EwaldBackground final : public QMCalculator {
 public:
  std::string Identify() const { return "ewaldbackground"; }
  bool WriteToStateFile() const { return false; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  std::string mapping_file_;
  std::string checkpoint_file_;

  double alpha_ = 0.5;      // bohr^-1 (internal). User-facing XML value
                            // for coulombmethod.alpha is nm^-1, converted
                            // in ParseOptions; only used as-is if the
                            // user explicitly sets it, otherwise
                            // overwritten in Evaluate() with a value
                            // derived from the box (see there, already
                            // in bohr^-1 at that point).
  bool alpha_explicit_ = false;
  double k_max_ = 3.0;      // bohr^-1 (internal); same nm^-1-in-XML
                            // pattern as alpha_ -- derived in Evaluate()
                            // from the (possibly also derived) alpha_,
                            // unless explicitly set.
  bool k_max_explicit_ = false;
  double thole_a_ = 0.39;
  double r_min_ = 18.897259886;  // bohr (internal) = 1.0 nm, the actual
                            // default applied in ParseOptions. User-
                            // facing XML value for realspace.r_min is
                            // nm, converted in ParseOptions; this member
                            // initializer only matters if ParseOptions
                            // somehow never runs.
  double field_tol_ = 1e-8;
  EwaldShape shape_ = EwaldShape::Cube;

  Index max_iter_ = 200;
  double pcg_tolerance_ = 1e-8;
  bool induce_ = true;  // matches legacy's own polarmethod.induce; when
                        // false, the PCG solve is skipped entirely and
                        // every site's induced dipole stays at zero
                        // (already the value it's set to before the
                        // permanent-field computation below) -- for
                        // comparing against a legacy run with induce=0,
                        // without needing a singular (zero-polarizability)
                        // PCG matrix on this side, which
                        // PolarSite::setpolarization's own
                        // eigenvalues.cwiseInverse() would produce.
  // Debug-only knob, no legacy counterpart (does not correspond to any
  // real physical option -- legacy always applies its own
  // FU12_ShapeField_At_By unconditionally). Exists to isolate whether
  // EwaldPeriodicDipoleOperator's own shape-correction term (see that
  // class's own class documentation for why it was added) is the cause
  // of a real, large PCG iteration-count regression observed once both
  // that term and Thole-damped intramolecular coupling were added in the
  // same session -- a genuine open question at the time this was added,
  // not a permanent feature. Defaults to true (matches legacy); setting
  // it false for a direct A/B comparison, everything else unchanged, is
  // the intended and only use.
  bool apply_shape_correction_to_induced_ = true;
  // Same purpose and pattern as apply_shape_correction_to_induced_
  // immediately above, but for EwaldPeriodicDipoleOperator's own
  // self_field_matrix_ term instead -- see that member's own
  // documentation.
  bool apply_self_field_correction_to_induced_ = true;
  // Same purpose and pattern as the two above, but for
  // AddIntraSegmentCoupling's own Thole-damping term instead -- see
  // EwaldPeriodicDipoleOperator's own apply_thole_damping_intramolecular_
  // documentation.
  bool apply_thole_damping_intramolecular_ = true;
  // Same purpose and pattern as the three above, but for
  // EwaldPeriodicDipoleOperator's own real_sum_.AddFieldAt call
  // (real-space intermolecular coupling) -- added after those three
  // turned out to be small contributions on both this codebase's own
  // side and legacy's own, leaving a real, substantial coupling-term
  // mismatch found this session unlocalized any further than
  // "somewhere in {thole, real-space intermolecular, reciprocal-space}
  // combined" with only those three toggles.
  bool apply_realspace_intermolecular_coupling_ = true;
  // Same purpose and pattern as apply_realspace_intermolecular_
  // coupling_ immediately above, but for EwaldPeriodicDipoleOperator's
  // own recip_sum_.AddFieldAtMany call (reciprocal-space coupling)
  // instead.
  bool apply_reciprocal_coupling_ = true;
  // Debug/experimental, no legacy counterpart -- see
  // EwaldBlockJacobiPreconditioner's own class documentation for what
  // this is for and why it exists, and its own usage-pattern note (the
  // real reason this needed its own if/else branch below rather than a
  // single runtime-switchable cg instance: the preconditioner type is a
  // compile-time Eigen::ConjugateGradient template parameter). Not yet
  // validated against a real system this size -- only its own mechanics
  // (compiles, runs, solves correctly) have been confirmed, with a
  // synthetic test operator, not this real one.
  bool use_block_jacobi_preconditioner_ = false;
  // Debug/experimental. Switches the whole solve from PCG (either
  // preconditioner) to plain weighted-Jacobi iteration (JOR) -- see
  // SolveWithJOR's own documentation for the algorithm and its direct
  // correspondence to legacy PolarBackground's own SOR. Exists
  // specifically because a real run this session found PCG's own
  // operator to be genuinely indefinite (p.A.p < 0, a direct algebraic
  // certificate, not an inference -- confirmed via
  // SolveWithIndefinitenessCheck), which CG's own convergence theory
  // cannot handle regardless of preconditioner; JOR makes no positive-
  // definiteness assumption anywhere in its own derivation.
  bool use_jor_ = false;
  // Relaxation factor for use_jor_'s own iteration. Matches legacy
  // PolarBackground's own default for this specific calculator
  // (_polar_wSOR_N = 0.35, confirmed directly against legacy's own
  // source -- NOT the more commonly cited generic 0.25 default declared
  // on APolarSite::Induce's own signature, which this calculator
  // overrides).
  double jor_omega_ = 0.35;
  // Debug/experimental. See SolveWithJOR's own declaration for what
  // this does and why it exists (literally replicating legacy's own
  // two-phase unrelaxed-then-relaxed structure, rather than relying on
  // an unverified assumption to relate this codebase's own JOR
  // recursion to legacy's).
  bool match_legacy_first_step_ = false;
  // Debug/experimental, no legacy counterpart, temporary -- exists to
  // support a direct, deterministic comparison against legacy at two
  // points in the solve where both codes compute well-defined
  // quantities: the first JOR update, mu_1 = omega * P * F_perm (no
  // induced-induced coupling yet on either side -- confirmed to match
  // legacy's own equivalent, once compared correctly, to ~1e-6), and
  // the second, where induced-induced coupling first enters (a real
  // suspicion this session, following that first match, that THIS is
  // where a genuine discrepancy lives). When set, SolveWithJOR dumps
  // mu_1 (segment id, site index, position, and dipole components --
  // both this codebase's own atomic units and legacy's own e*nm, to
  // allow a DIRECT diff against a legacy dump without a manual unit
  // conversion step first) to jor_iteration1_dump.csv after the first
  // update, continues one more iteration, dumps mu_1/mu_2/the coupling-
  // carrying residual to jor_iteration2_dump.csv after the second, then
  // stops (throws) -- deliberately not a silent early return, so this
  // can never be left on by accident and mistaken for a real solve.
  bool debug_dump_first_iteration_and_stop_ = false;
};

inline void EwaldBackground::ParseOptions(const tools::Property& opt) {
  mapping_file_ = opt.get("multipoles").as<std::string>();
  checkpoint_file_ = opt.ifExistsReturnElseReturnDefault<std::string>(
      "output.checkpoint", "ewaldbackground.hdf5");

  alpha_explicit_ = opt.exists("coulombmethod.alpha");
  if (alpha_explicit_) {
    // User-facing value is nm^-1; alpha_ is stored internally in
    // bohr^-1. alpha has units of 1/length, so converting the numeric
    // value goes the same direction as EwaldBackground's own box-derived
    // default below: bohr^-1 = nm^-1 * bohr2nm (NOT nm2bohr -- that
    // would be backwards for an inverse-length quantity; see the field
    // conversion note in ptop_dump.cc for the same kind of
    // easy-to-invert factor, worked through carefully there).
    alpha_ = opt.get("coulombmethod.alpha").as<double>() * tools::conv::bohr2nm;
  }
  // Default k_max derived from alpha_, not a fixed literal: the
  // reciprocal-space Gaussian weight exp(-k^2/4*alpha^2) means a k_max
  // that isn't scaled to alpha can generate a k-vector count many orders
  // of magnitude larger than what's actually needed for convergence, for
  // no numerical benefit -- k=6*alpha already gives exp(-9) ~ 1.2e-4, a
  // reasonable default cutoff. This was a real, costly mistake in an
  // earlier version of this class (a fixed k_max=15.0 default, entirely
  // disconnected from alpha_, caused a real run on a 1000-segment system
  // to hang for an extremely long time computing near-zero-weight
  // k-vectors) -- worth being explicit about here so it isn't
  // reintroduced by, e.g., reverting this line without reading why.
  // alpha_ itself, in turn, is resolved in Evaluate() rather than here
  // when not explicit (see there) since a good default genuinely depends
  // on the box, which isn't available yet at this point -- so k_max_'s
  // own default is also finalized there, once the (possibly derived)
  // alpha_ is known.
  k_max_explicit_ = opt.exists("coulombmethod.k_max");
  if (k_max_explicit_) {
    // Same nm^-1 (user-facing) -> bohr^-1 (internal) conversion as
    // alpha_ above.
    k_max_ = opt.get("coulombmethod.k_max").as<double>() * tools::conv::bohr2nm;
  }
  std::string shape_str = opt.ifExistsReturnElseReturnDefault<std::string>(
      "coulombmethod.shape", "cube");
  if (shape_str == "cube" || shape_str == "sphere") {
    shape_ = EwaldShape::Cube;
  } else if (shape_str == "slab") {
    shape_ = EwaldShape::Slab;
  } else {
    throw std::runtime_error(
        "EwaldBackground: coulombmethod.shape must be 'cube', 'sphere', or "
        "'slab'");
  }

  thole_a_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.thole_a", thole_a_);
  // realspace.r_min is user-facing nm; r_min_ is stored internally in
  // bohr. r_min is a plain length (not an inverse length like alpha_/
  // k_max_ above), so this conversion goes the more familiar direction:
  // bohr = nm * nm2bohr. 1.0 nm is the default here (~18.9 bohr,
  // replacing the old bare-bohr literal default of 20.0).
  double r_min_nm = opt.ifExistsReturnElseReturnDefault<double>(
      "realspace.r_min", 1.0);
  r_min_ = r_min_nm * tools::conv::nm2bohr;
  field_tol_ = opt.ifExistsReturnElseReturnDefault<double>(
      "realspace.field_tol", field_tol_);

  max_iter_ = opt.ifExistsReturnElseReturnDefault<Index>(
      "polarmethod.max_iter", max_iter_);
  pcg_tolerance_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.tolerance", pcg_tolerance_);
  induce_ = opt.ifExistsReturnElseReturnDefault<bool>("polarmethod.induce",
                                                       induce_);
  // Debug-only, undocumented on purpose (no legacy counterpart to give
  // it a natural home in the schema) -- see this member's own
  // declaration for what it's for.
  apply_shape_correction_to_induced_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_apply_shape_correction",
      apply_shape_correction_to_induced_);
  // Debug-only, undocumented on purpose, same pattern as the shape one
  // immediately above.
  apply_self_field_correction_to_induced_ =
      opt.ifExistsReturnElseReturnDefault<bool>(
          "polarmethod.debug_apply_self_field_correction",
          apply_self_field_correction_to_induced_);
  // Debug-only, undocumented on purpose, same pattern as the two above.
  apply_thole_damping_intramolecular_ =
      opt.ifExistsReturnElseReturnDefault<bool>(
          "polarmethod.debug_apply_thole_damping_intramolecular",
          apply_thole_damping_intramolecular_);
  // Debug-only, undocumented on purpose, same pattern as the three above.
  apply_realspace_intermolecular_coupling_ =
      opt.ifExistsReturnElseReturnDefault<bool>(
          "polarmethod.debug_apply_realspace_intermolecular_coupling",
          apply_realspace_intermolecular_coupling_);
  // Debug-only, undocumented on purpose, same pattern as the four above.
  apply_reciprocal_coupling_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_apply_reciprocal_coupling",
      apply_reciprocal_coupling_);
  // Debug/experimental, undocumented on purpose -- see
  // use_block_jacobi_preconditioner_'s own declaration.
  use_block_jacobi_preconditioner_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_use_block_jacobi_preconditioner",
      use_block_jacobi_preconditioner_);
  // Debug/experimental, see use_jor_'s own declaration.
  use_jor_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_use_jor", use_jor_);
  jor_omega_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.debug_jor_omega", jor_omega_);
  match_legacy_first_step_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_match_legacy_first_step", match_legacy_first_step_);
  debug_dump_first_iteration_and_stop_ =
      opt.ifExistsReturnElseReturnDefault<bool>(
          "polarmethod.debug_dump_first_iteration_and_stop",
          debug_dump_first_iteration_and_stop_);
}

// Runs the standard preconditioned CG algorithm directly (rather than via
// Eigen::ConjugateGradient) so every iteration's own curvature,
// p.dot(A*p), can be inspected. For a genuinely SPD operator this is
// mathematically guaranteed positive for every nonzero p, at every
// iteration -- CG's own convergence theory depends on it. If it is ever
// <= 0, that is a direct proof (not an inference from residual behavior)
// that the operator is NOT positive-definite, and this stops immediately
// with a clear diagnostic rather than continue computing meaningless
// iterations. This exists specifically to distinguish, for a real
// divergence found this session (residual growing past its own starting
// value with a plain PCG solve, on a system legacy's own SOR handles
// fine at the same physics), between two genuinely different
// explanations that could not be told apart from residual behavior
// alone: the true global operator being indefinite (this check would
// fire), versus a merely poorly-preconditioned but still SPD system
// (this check would never fire, even on a slowly-converging or stalled
// run). Templated on the preconditioner type for the same reason
// EwaldBackground::Evaluate's own PCG solve needs two separate branches
// below: the preconditioner type is a compile-time choice here too, and
// this avoids writing the same loop out twice.
struct PcgIndefinitenessResult {
  Eigen::VectorXd x;
  Index iterations = 0;
  double residual = 0.0;
  bool converged = false;
  // -1 if the operator was never found indefinite; otherwise the
  // 1-based iteration number at which p.dot(A*p) <= 0 was first found.
  // This is a real, unambiguous certificate whenever it fires (CG's own
  // convergence theory strictly requires p.dot(A*p) > 0 for every
  // nonzero p if the operator is genuinely SPD), but it can miss real
  // indefiniteness for a "lucky" (or, from this test's own point of
  // view, unlucky) b that happens not to excite the negative-eigenvalue
  // direction at all -- confirmed directly with a deliberately
  // constructed 3x3 counterexample (a genuinely indefinite matrix,
  // eigenvalues -1/3/5, where b=[1,1,1] is EXACTLY orthogonal to the
  // eigenvector for -1, letting CG converge to the exact right answer
  // in 2 iterations without curvature ever going non-positive). See
  // lanczos_min_eigenvalue below for a real (if still b-dependent, just
  // less narrowly so) improvement on this.
  Index indefinite_at_iteration = -1;
  double indefinite_curvature = 0.0;
  // The smallest eigenvalue of the (small, m x m, m = iterations run)
  // Lanczos tridiagonal matrix built from every alpha_k/beta_k this CG
  // run itself computed -- a standard result (the CG-Lanczos
  // correspondence; verified directly here against known eigenvalues
  // for both an SPD test matrix, exact match to 1e-14, and a genuinely
  // indefinite one with a GENERIC right-hand side, exact match there
  // too) that these Ritz values approximate the extreme eigenvalues of
  // the (preconditioned) operator using the WHOLE accumulated Krylov
  // subspace, not one iteration's own snapshot the way
  // indefinite_at_iteration is. Still not a b-independent guarantee (the
  // same 3x3 counterexample above, with its own adversarially exact
  // orthogonality, defeats this too -- confirmed directly) but
  // meaningfully more robust for any b that is not exactly orthogonal to
  // the relevant eigenvector, which a genuine physical field vector,
  // not an adversarially constructed one, essentially never is exactly.
  // NaN if fewer than 2 iterations ran (need at least 2 alphas for a
  // meaningful 2x2 tridiagonal matrix).
  double lanczos_min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
};

template <typename Preconditioner>
PcgIndefinitenessResult SolveWithIndefinitenessCheck(
    const EwaldPeriodicDipoleOperator& op, Preconditioner& precond,
    const Eigen::VectorXd& b, Index max_iter, double tol, Logger& log,
    std::chrono::steady_clock::time_point t_start) {
  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         t_start)
        .count();
  };

  PcgIndefinitenessResult result;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  Eigen::VectorXd residual_vec = b - op * x;
  const double rhs_norm2 = b.squaredNorm();
  const double threshold =
      std::max(tol * tol * rhs_norm2, std::numeric_limits<double>::min());
  double residual_norm2 = residual_vec.squaredNorm();
  double tol_error = std::sqrt(residual_norm2 / rhs_norm2);
  bool converged = (residual_norm2 < threshold);
  // Every alpha_k/beta_k this run computes, in order -- see
  // lanczos_min_eigenvalue's own documentation for what these build.
  std::vector<double> alphas;
  std::vector<double> betas;

  if (!converged) {
    Eigen::VectorXd p = precond.solve(residual_vec);
    Eigen::VectorXd z(b.size()), tmp(b.size());
    double abs_new = residual_vec.dot(p);
    Index i = 0;
    while (i < max_iter) {
      tmp.noalias() = op * p;
      const double curvature = p.dot(tmp);
      XTP_LOG(Log::info, log)
          << TimeStamp() << "   PCG iter " << (i + 1)
          << ": curvature p.A.p=" << curvature << " (" << elapsed_s() << "s)"
          << std::flush;
      if (curvature <= 0.0) {
        result.indefinite_at_iteration = i + 1;
        result.indefinite_curvature = curvature;
        XTP_LOG(Log::info, log)
            << TimeStamp() << "   PCG iter " << (i + 1)
            << ": p.A.p = " << curvature
            << " <= 0 -- the operator is NOT positive-definite (this is a "
               "direct algebraic certificate, not an inference from "
               "residual behavior). Stopping here rather than continue "
               "computing iterations that CG's own convergence theory no "
               "longer covers."
            << std::flush;
        break;
      }
      const double alpha = abs_new / curvature;
      alphas.push_back(alpha);
      x += alpha * p;
      residual_vec -= alpha * tmp;

      residual_norm2 = residual_vec.squaredNorm();
      tol_error = std::sqrt(residual_norm2 / rhs_norm2);
      if (residual_norm2 < threshold) {
        converged = true;
        ++i;
        break;
      }

      z = precond.solve(residual_vec);
      const double abs_old = abs_new;
      abs_new = residual_vec.dot(z);
      const double beta = abs_new / abs_old;
      betas.push_back(beta);
      p = z + beta * p;
      ++i;
    }
    result.iterations = i;
  } else {
    result.iterations = 0;
  }

  // Lanczos tridiagonal matrix from alphas/betas -- see
  // lanczos_min_eigenvalue's own documentation for the derivation and
  // its own direct numerical verification (both SPD and indefinite test
  // cases, exact match to the true eigenvalues in both). Needs at least
  // 2 alphas (i.e. 1 beta) for a meaningful 2x2 matrix; a 1x1 "matrix"
  // is just 1/alphas[0], not a real eigenvalue ESTIMATE of anything.
  if (alphas.size() >= 2) {
    const std::size_t m = alphas.size();
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m, m);
    T(0, 0) = 1.0 / alphas[0];
    for (std::size_t i = 1; i < m; ++i) {
      T(i, i) = 1.0 / alphas[i] + betas[i - 1] / alphas[i - 1];
      const double off = std::sqrt(betas[i - 1]) / alphas[i - 1];
      T(i, i - 1) = off;
      T(i - 1, i) = off;
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    result.lanczos_min_eigenvalue = es.eigenvalues()(0);
  }

  result.x = x;
  result.residual = tol_error;
  result.converged = converged;
  return result;
}

// Per-site polarizability blocks (P, not P^-1), used for the plain
// Jacobi-Over-Relaxation (JOR) iteration below. Confirmed directly
// against legacy PolarBackground's own source (not re-derived from
// theory alone) that this, not EwaldBlockJacobiPreconditioner's own
// larger (intramolecular-Thole-coupled) block, is the exact match to
// legacy's own per-site update: APolarSite::Induce computes
//   mu_new = (1-wSOR)*mu_old + wSOR*(-P)*(F_perm+F_induced)
// which is a textbook weighted-Jacobi step once D = P^-1 is recognized
// as this system's own diagonal block (EwaldPeriodicDipoleOperator's
// own diagonal is built from getPInv() exactly) -- D^-1 = P follows
// directly, with no re-derivation needed of the specific field-based
// expression legacy itself uses internally.
//
// Also confirmed directly (by tracing legacy's own PolarBackground::
// Evaluate) that legacy's own "SOR" is Jacobi-style, not sequential
// Gauss-Seidel-based SOR despite the shared name: legacy computes EVERY
// site's own field in full (step III.B, "(Re-)generate induction
// fields") before updating ANY site's own dipole (step III.C, "Induce
// again") -- every field used in a given iteration is built entirely
// from the PREVIOUS iteration's own induced dipoles.
class EwaldSitePolarizabilityBlocks {
 public:
  EwaldSitePolarizabilityBlocks(const EwaldRegistry& registry,
                                const std::vector<Index>& ids) {
    for (Index id : ids) {
      const PolarSegment& segment =
          registry.Get(id, EwaldChargeState::Neutral);
      for (const PolarSite& site : segment) {
        // getPInv() is this codebase's own stored quantity (see
        // EwaldPeriodicDipoleOperator's own diagonal); P itself is not
        // separately exposed anywhere, so recovering it needs a real
        // 3x3 inverse here -- cheap and exact for a non-singular
        // polarizability tensor.
        blocks_.push_back(site.getPInv().inverse());
      }
    }
  }

  Index size() const { return 3 * Index(blocks_.size()); }

  // Exposes the per-site P block directly -- added specifically to
  // support debug_dump_first_iteration_and_stop_'s own extended dump
  // (F_perm and P alongside mu_1), so a genuine mu_1 = omega*P*F_perm
  // discrepancy against legacy can be isolated to F_perm, P, or their
  // combination, rather than left as an unexplained aggregate mismatch.
  const Eigen::Matrix3d& GetBlock(Index n) const { return blocks_[n]; }

  // Returns P * v, applied per 3-DOF site block -- the D^-1 * v the
  // JOR iteration below needs (see this class's own documentation for
  // why D^-1 = P here specifically). Site ordering matches
  // EwaldPeriodicDipoleOperator's own targets_ layout exactly (same
  // registry, same ids, same per-segment site iteration order), so no
  // separate offset bookkeeping is needed here: site n's own 3 DOF are
  // always at [3n, 3n+3) for both this class and that operator.
  Eigen::VectorXd Apply(const Eigen::VectorXd& v) const {
    Eigen::VectorXd result(size());
    for (std::size_t n = 0; n < blocks_.size(); ++n) {
      result.segment<3>(3 * Index(n)) =
          blocks_[n] * v.segment<3>(3 * Index(n));
    }
    return result;
  }

 private:
  std::vector<Eigen::Matrix3d> blocks_;
};

struct JorResult {
  Eigen::VectorXd x;
  Index iterations = 0;
  double residual = 0.0;
  double max_dU = 0.0;
  double avg_dU = 0.0;
  bool converged = false;
};

// Mirrors legacy APolarSite::HistdU() exactly (confirmed directly
// against its own source, apolarsite.cc): per-site relative change in
// induced dipole moment between successive iterates, with a small-
// magnitude fallback so a near-zero dipole (e.g. this site's own first
// iteration, starting from x=0) doesn't divide by ~zero. Legacy's own
// "small" constant is documented there as "1e-10 e*nm" -- this
// codebase uses atomic units (bohr) throughout, not nm (see e.g.
// EwaldRealSpaceSum's own class documentation), so using legacy's raw
// 1e-10 unconverted would apply the wrong ABSOLUTE threshold (1e-10 nm
// and 1e-10 bohr are very different physical dipole magnitudes) --
// tools::conv::nm2bohr is this codebase's own existing conversion
// constant (already used elsewhere in this file for r_min), reused
// here rather than a separately hardcoded factor. epstol itself (see
// SolveWithJOR below) needs no such conversion: it's a ratio of two
// same-unit quantities, genuinely dimensionless regardless of which
// unit system U0/U1 happen to be expressed in.
double SiteRelativeDipoleChange(const Eigen::Vector3d& U0,
                                const Eigen::Vector3d& U1) {
  const double small = 1e-10 * tools::conv::nm2bohr;
  const double abs_U0 = U0.norm();
  const double abs_U1 = U1.norm();
  const double abs_dU = (U1 - U0).norm();
  const double abs_U = (abs_U1 > abs_U0 && small > abs_U0) ? abs_U1 : abs_U0;
  if (small > abs_U) {
    return abs_U;
  }
  return abs_dU / abs_U;
}

// Plain weighted-Jacobi iteration for A*x=b (A = EwaldPeriodicDipoleOperator),
// matching legacy PolarBackground's own SOR exactly in structure (see
// EwaldSitePolarizabilityBlocks's own documentation for the direct
// confirmation against legacy's own source). Iteration, in the general
// weighted-Jacobi form for A = D + R:
//   x_new = x_old + omega * D^-1 * (b - A*x_old)
// Unlike PCG, this makes no positive-definiteness assumption anywhere
// in its own derivation -- its convergence instead depends on the
// spectral radius of the iteration matrix (I - omega*D^-1*A) staying
// below 1, a genuinely different (and empirically, for this kind of
// system, more forgiving) condition. Verified directly, before writing
// any of this, against a small SPD test system (exact match to the
// direct solve) and against a genuinely indefinite one (the same 3x3
// counterexample used to validate this session's own indefiniteness
// checks) -- JOR converges on the indefinite case where CG structurally
// cannot, matching the entire reason for using it here.
//
// Convergence criterion matches legacy PolarBackground's own exactly
// (confirmed directly against its own source, step III.D) -- NOT the
// global relative residual ||b-Ax||/||b|| this function's own first
// version used (borrowed, incorrectly, from the PCG code path's own
// stopping rule): legacy checks a per-site relative change in induced
// dipole moment between successive iterates (SiteRelativeDipoleChange
// above) against epstol=1e-3 (legacy's own hardcoded value, not this
// codebase's own pcg_tolerance_) for EVERY site, with a looser
// avgdU < epstol*0.1 fallback that can force convergence even if a few
// individual sites remain just above epstol. The global residual is
// still computed and logged every iteration (informational only -- it
// is a genuinely different quantity, not proportional to max_dU/avg_dU
// in general, and no longer decides when this function itself stops).
//
// match_legacy_first_step, when true, uses omega=1.0 (fully unrelaxed)
// for iteration 0 only, then the configured omega for every iteration
// after -- literally replicating legacy's own two-phase structure
// (InduceDirect(), fully unrelaxed, THEN the main loop's own wSOR-
// relaxed Induce() calls) rather than applying the same omega-relaxed
// step from x=0 that legacy never does. Added after this session's own
// attempt to algebraically relate the two codebases' differing
// recursions (legacy: mu2 = mu1 - wSOR*P*FU(mu1), coefficient 1 on
// mu1, since legacy's own mu1 is already fully unrelaxed; this
// function's own un-adjusted recursion: x2 = (2-omega)*x1 +
// omega*P*C*x1, a different coefficient, purely because THIS x1 was
// only omega-relaxed to begin with) turned out to rest on an
// unverified assumption (that C, this class's own leftover-from-D
// coupling operator, equals legacy's own FU field) that was never
// independently confirmed -- this flag sidesteps that assumption
// entirely by making the two codebases' own x1/mu1 and x2/mu2 directly,
// literally comparable (same relaxation history on both sides), so a
// comparison needs only the same simple bohr<->nm dipole conversion
// already validated for x1 vs mu1, not any field-specific conversion
// or unverified operator decomposition.
JorResult SolveWithJOR(const EwaldPeriodicDipoleOperator& op,
                       const EwaldSitePolarizabilityBlocks& site_p,
                       const Eigen::VectorXd& b, Index max_iter,
                       double omega, Logger& log,
                       std::chrono::steady_clock::time_point t_start,
                       const EwaldRegistry& registry,
                       const std::vector<Index>& ids,
                       bool debug_dump_first_iteration_and_stop,
                       bool match_legacy_first_step) {
  auto elapsed_s = [&]() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         t_start)
        .count();
  };

  const double kEpsTol = 1e-3;  // legacy's own hardcoded value

  JorResult result;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  const double rhs_norm = b.norm();
  const Index n_sites = b.size() / 3;

  Index i = 0;
  for (; i < max_iter; ++i) {
    Eigen::VectorXd residual_vec = b - op * x;
    const double residual_norm = residual_vec.norm();
    const Eigen::VectorXd x_old = x;
    const double this_step_omega =
        (i == 0 && match_legacy_first_step) ? 1.0 : omega;
    x = x_old + this_step_omega * site_p.Apply(residual_vec);

    if (i == 0 && debug_dump_first_iteration_and_stop) {
      // See debug_dump_first_iteration_and_stop_'s own declaration for
      // why this exists and what it's for. x at this point is exactly
      // this_step_omega * P * F_perm (x_old is still zero here, so
      // site_p.Apply(residual_vec) IS the whole update -- no induced-
      // induced coupling has entered yet, on either this codebase's own
      // side or legacy's, since neither has performed a second field
      // evaluation). residual_vec itself, at this specific point (still
      // i==0, x_old=0), equals b exactly (residual_vec = b - op*0 = b),
      // and b is itself exactly F_perm: EwaldBackground::Evaluate's own
      // permanent-field loop writes b.segment<3>(...) = site.V()
      // directly, then immediately resets site.V() to zero -- so this
      // is the one point where dumping residual_vec captures F_perm
      // without needing to re-fetch it from anywhere. P itself (the
      // diagonal of the same 3x3 block site_p.Apply uses) is exposed
      // via GetBlock, added specifically for this dump.
      std::ofstream dump("jor_iteration1_dump.csv");
      dump << "segment_id,site_index,pos_x_bohr,pos_y_bohr,pos_z_bohr,"
              "mu_x_bohr,mu_y_bohr,mu_z_bohr,mu_x_enm,mu_y_enm,mu_z_enm,"
              "Fperm_x_bohr,Fperm_y_bohr,Fperm_z_bohr,"
              "Pxx_bohr3,Pyy_bohr3,Pzz_bohr3,"
              "Pxy_bohr3,Pxz_bohr3,Pyz_bohr3\n";
      Index n = 0;
      for (Index id : ids) {
        const PolarSegment& segment =
            registry.Get(id, EwaldChargeState::Neutral);
        Index site_index = 0;
        for (const PolarSite& site : segment) {
          const Eigen::Vector3d pos = site.getPos();
          const Eigen::Vector3d mu = x.segment<3>(3 * n);
          // bohr2nm is this codebase's own existing conversion constant
          // (tools::conv), the same one already used elsewhere in this
          // file for r_min -- reused here rather than a separately
          // hardcoded factor, for a direct diff against a legacy dump
          // without a manual conversion step first.
          const Eigen::Vector3d mu_enm = mu * tools::conv::bohr2nm;
          const Eigen::Vector3d f_perm = residual_vec.segment<3>(3 * n);
          const Eigen::Matrix3d& P = site_p.GetBlock(n);
          // Off-diagonal terms added after a real run this session
          // showed the diagonal-only reconstruction genuinely failing
          // for legacy -- this codebase's own P itself is NOT
          // guaranteed diagonal either (PolarSite applies a real
          // rotation, pinv_ = R^T * pinv_ * R, when mapping a site's
          // own polarizability into the simulation frame), so a
          // diagonal-only self-consistency check here would only be
          // reassuring, not conclusive, for the same reason.
          dump << id << "," << site_index << "," << pos.x() << ","
               << pos.y() << "," << pos.z() << "," << mu.x() << ","
               << mu.y() << "," << mu.z() << "," << mu_enm.x() << ","
               << mu_enm.y() << "," << mu_enm.z() << "," << f_perm.x()
               << "," << f_perm.y() << "," << f_perm.z() << ","
               << P(0, 0) << "," << P(1, 1) << "," << P(2, 2) << ","
               << P(0, 1) << "," << P(0, 2) << "," << P(1, 2) << "\n";
          ++site_index;
          ++n;
        }
      }
      dump.close();
      // Debug/experimental -- see DumpStagedCoupling's own declaration.
      // x here IS x1 (this loop's own i==0 update, captured before any
      // later iteration overwrites it) -- matching legacy's own FUa-FUe,
      // all computed from mu_1, not from a later, already-relaxed
      // iterate.
      op.DumpStagedCoupling(x, "new_staged_coupling_dump.csv");
      // Debug/experimental -- see DumpCachedNeighborList's own
      // declaration. Called here, right after DumpStagedCoupling's own
      // stage B (which calls AddFieldAt for this same target), to
      // capture exactly what that specific call visited/cached --
      // added after a real discrepancy was found, on the real system,
      // between DumpStagedCoupling and DumpPerPairIntermolecularField,
      // that no local reproduction attempt (small test systems,
      // various call orderings) could reproduce.
      {
        std::ofstream cache_header("new_cachedneighborlist_dump.csv");
        cache_header << "target_segment_id,source_segment_id,"
                        "translation_idx,t_x_bohr,t_y_bohr,t_z_bohr,r_bohr\n";
      }
      op.DumpCachedNeighborList(0, "new_cachedneighborlist_dump.csv");
      // Debug/experimental -- see DumpPerPairIntermolecularField's own
      // declaration. x's own dipoles are already set on every site by
      // DumpStagedCoupling's own internal setup, right above -- this
      // call relies on that already having happened, per its own
      // documented contract of not setting any dipole state itself.
      // Restricted to segment 0 to match legacy's own equivalent dump
      // (also gated to a single target, given the far finer per-pair
      // granularity here compared to the earlier per-segment
      // neighbor-list dump, which could afford to cover every target).
      {
        std::ofstream pp_header("new_perpairfield_dump.csv");
        pp_header << "target_segment_id,source_segment_id,source_site_index,"
                     "t_x_bohr,t_y_bohr,t_z_bohr,field_x_bohr,field_y_bohr,"
                     "field_z_bohr\n";
      }
      op.DumpPerPairIntermolecularField(0, "new_perpairfield_dump.csv");
      // Debug/experimental -- see DumpPerPairIntermolecularStaticField's
      // own declaration. Isolates ApplyStaticField's own contribution
      // (F_perm's own generation, and half of what AddFieldAt/
      // DumpStagedCoupling's own stage B actually compute) -- added
      // after the induced-only comparison above came back clean but
      // the CUMULATIVE FUb comparison against legacy still showed a
      // large mismatch, meaning the static half (never independently
      // re-checked on the current, real-space-minimum-image-fixed
      // code) needed the same direct treatment.
      {
        std::ofstream pps_header("new_perpairstaticfield_dump.csv");
        pps_header << "target_segment_id,source_segment_id,source_site_index,"
                      "t_x_bohr,t_y_bohr,t_z_bohr,field_x_bohr,field_y_bohr,"
                      "field_z_bohr\n";
      }
      op.DumpPerPairIntermolecularStaticField(
          0, "new_perpairstaticfield_dump.csv");
    }

    if (i == 1 && debug_dump_first_iteration_and_stop) {
      // Second dump, following a suspicion this session that the
      // induced-induced coupling term itself (not the F_perm/P terms
      // the first dump already validated to ~1e-6, matching legacy
      // exactly once that comparison was done correctly) is where a
      // real discrepancy lives. x_old here is x_1 (the already-
      // validated mu_1 = 0.35*P*F_perm); x is x_2, the result of one
      // full round of induced-induced coupling entering through
      // residual_vec = b - op*x_1, which is dumped directly (not
      // reconstructed from an assumed internal sign convention for
      // op's own diagonal/off-diagonal split) precisely so no new
      // reconstruction bug can be introduced here the way the first
      // comparison's own unit-mismatch bug was. Since the coupling
      // Ewald sums are linear in the source (induced) dipoles, and
      // x_1 = 0.35 * (legacy's own mu_1 from InduceDirect) was already
      // confirmed to ~1e-6, the coupling contribution each code's own
      // x_1 produces should relate by the same 0.35 factor if the
      // coupling term itself is correctly implemented on both sides --
      // this is the direct, checkable prediction this second dump
      // exists to test.
      std::ofstream dump("jor_iteration2_dump.csv");
      dump << "segment_id,site_index,"
              "mu1_x_bohr,mu1_y_bohr,mu1_z_bohr,"
              "mu2_x_bohr,mu2_y_bohr,mu2_z_bohr,"
              "residual2_x_bohr,residual2_y_bohr,residual2_z_bohr\n";
      Index n = 0;
      for (Index id : ids) {
        const PolarSegment& segment =
            registry.Get(id, EwaldChargeState::Neutral);
        Index site_index = 0;
        for (const PolarSite& site : segment) {
          (void)site;
          const Eigen::Vector3d mu1 = x_old.segment<3>(3 * n);
          const Eigen::Vector3d mu2 = x.segment<3>(3 * n);
          const Eigen::Vector3d res2 = residual_vec.segment<3>(3 * n);
          dump << id << "," << site_index << "," << mu1.x() << ","
               << mu1.y() << "," << mu1.z() << "," << mu2.x() << ","
               << mu2.y() << "," << mu2.z() << "," << res2.x() << ","
               << res2.y() << "," << res2.z() << "\n";
          ++site_index;
          ++n;
        }
      }
      dump.close();
      throw std::runtime_error(
          "EwaldBackground: stopped after dumping iteration-2 induced "
          "dipoles to jor_iteration2_dump.csv (debug_dump_first_"
          "iteration_and_stop was set) -- this is not a real failure, "
          "just the requested early stop.");
    }

    double max_dU = -1.0;
    double avg_dU = 0.0;
    for (Index n = 0; n < n_sites; ++n) {
      const double dU = SiteRelativeDipoleChange(x_old.segment<3>(3 * n),
                                                 x.segment<3>(3 * n));
      avg_dU += dU;
      if (dU > max_dU) {
        max_dU = dU;
      }
    }
    avg_dU /= double(n_sites);

    XTP_LOG(Log::info, log)
        << TimeStamp() << "   JOR iter " << (i + 1) << ": max_dU=" << max_dU
        << " avg_dU=" << avg_dU << " (residual=" << residual_norm / rhs_norm
        << ", informational only -- see this function's own documentation "
           "for why max_dU/avg_dU, not this, is the actual stopping "
           "criterion) ("
        << elapsed_s() << "s)" << std::flush;

    result.residual = residual_norm / rhs_norm;
    result.max_dU = max_dU;
    result.avg_dU = avg_dU;

    bool converged = (max_dU <= kEpsTol);
    if (avg_dU < kEpsTol * 0.1) {
      converged = true;
    }
    if (converged) {
      result.converged = true;
      ++i;
      break;
    }
  }
  result.iterations = i;
  result.x = x;
  return result;
}

inline bool EwaldBackground::Evaluate(Topology& top) {
  Logger log;
  log.setReportLevel(Log::info);
  log.setCommonPreface("\nEWD");
  // maverick=true means "single job running alone" in VOTCA's own
  // terminology (job-parallel calculators set this false so multiple
  // threads' log output can be buffered and interleaved cleanly) --
  // this calculator is a single xtp_run job, not job-parallel, and
  // critically, false here would silently buffer every message below
  // until a final std::cout<<log flush, which would defeat the entire
  // point of adding real-time progress output at all.
  log.setMultithreading(true);

  auto t_start = std::chrono::steady_clock::now();
  auto elapsed_s = [&](std::chrono::steady_clock::time_point since) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         since)
        .count();
  };

  XTP_LOG(Log::info, log)
      << TimeStamp() << " Starting Ewald background calculation"
      << std::flush;

  PolarMapper polmap(log);
  polmap.LoadMappingFile(mapping_file_);

  EwaldRegistry registry;
  std::vector<Index> ids;
  // offsets[n] is the global vector index where segment ids[n]'s own
  // block starts (in units of a single scalar, 3 per site); mirrors
  // EwaldPeriodicDipoleOperator's own offsets_ table exactly, since b and
  // x here must use the same layout that operator expects.
  std::vector<Index> offsets;
  ids.reserve(top.Segments().size());
  offsets.reserve(top.Segments().size() + 1);
  offsets.push_back(0);
  for (const Segment& seg : top.Segments()) {
    PolarSegment mol =
        polmap.map(seg, SegId(seg.getId(), std::string("n")));
    registry.Register(seg.getId(), EwaldChargeState::Neutral, mol);
    ids.push_back(seg.getId());
    offsets.push_back(offsets.back() + 3 * mol.size());
  }
  const Index total_size = offsets.back();

  XTP_LOG(Log::info, log)
      << TimeStamp() << " Mapped " << ids.size() << " segments, "
      << (total_size / 3) << " polarizable sites total ("
      << elapsed_s(t_start) << "s)" << std::flush;

  const Eigen::Matrix3d& box = top.getBox();
  const double volume = box.col(0).dot(box.col(1).cross(box.col(2)));

  // Standard Ewald rule of thumb: alpha ~ 3/L_min, where L_min is the
  // box's shortest lattice vector length -- chosen so erfc(alpha*r) has
  // decayed to a small value (erfc(3.0) ~ 2.2e-5) by roughly the point a
  // traditional minimum-image cutoff would sit at. EwaldRealSpaceSum's
  // own adaptive shell walk isn't actually bound by that minimum-image
  // constraint (it will walk past L_min/2 if field_tol demands it), so
  // this is a starting point rather than a hard requirement -- but it's
  // still the right scale to default to, since a fixed alpha_ literal
  // disconnected from box size is exactly the same class of mistake the
  // k_max_ default already was (see above): too large an alpha for a
  // large box pushes cost into reciprocal space (where it's cubic in
  // k_max, hence in alpha) for no numerical benefit.
  if (!alpha_explicit_) {
    const double L_min = std::min(
        {box.col(0).norm(), box.col(1).norm(), box.col(2).norm()});
    alpha_ = 3.0 / L_min;
  }
  if (!k_max_explicit_) {
    k_max_ = 6.0 * alpha_;
  }

  EwaldRealSpaceSum real_sum(box, registry, alpha_, thole_a_, r_min_,
                            field_tol_);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_, k_max_);
  EwaldShapeCorrection shape(volume, registry, shape_);

  XTP_LOG(Log::info, log)
      << TimeStamp() << " Real/reciprocal-space sums constructed ("
      << elapsed_s(t_start) << "s)" << std::flush;
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Parameters: alpha=" << alpha_ << " bohr^-1 ("
      << alpha_ * tools::conv::nm2bohr << " nm^-1), k_max=" << k_max_
      << " bohr^-1 (" << k_max_ * tools::conv::nm2bohr << " nm^-1, "
      << recip_sum.NumKVectors() << " k-vectors), thole_a=" << thole_a_
      << ", r_min=" << r_min_ << " bohr (" << r_min_ * tools::conv::bohr2nm
      << " nm), field_tol=" << field_tol_ << ", shape="
      << (shape_ == EwaldShape::Cube ? "cube" : "slab")
      << ", box volume=" << volume << " bohr^3, max_iter=" << max_iter_
      << ", pcg_tolerance=" << pcg_tolerance_
      << ", induce=" << (induce_ ? "true" : "false")
      << ", debug_apply_shape_correction="
      << (apply_shape_correction_to_induced_ ? "true" : "false")
      << ", debug_apply_self_field_correction="
      << (apply_self_field_correction_to_induced_ ? "true" : "false")
      << ", debug_apply_thole_damping_intramolecular="
      << (apply_thole_damping_intramolecular_ ? "true" : "false")
      << ", debug_apply_realspace_intermolecular_coupling="
      << (apply_realspace_intermolecular_coupling_ ? "true" : "false")
      << ", debug_apply_reciprocal_coupling="
      << (apply_reciprocal_coupling_ ? "true" : "false")
      << ", debug_use_block_jacobi_preconditioner="
      << (use_block_jacobi_preconditioner_ ? "true" : "false")
      << ", debug_use_jor=" << (use_jor_ ? "true" : "false")
      << ", debug_jor_omega=" << jor_omega_
      << ", debug_match_legacy_first_step="
      << (match_legacy_first_step_ ? "true" : "false")
      << ", debug_dump_first_iteration_and_stop="
      << (debug_dump_first_iteration_and_stop_ ? "true" : "false")
      << std::flush;

  // Permanent field only (every induced dipole is still zero at this
  // point), computed once -- this becomes b for the PCG solve below,
  // matching CalcInducedDipolesViaPCG's own b-construction pattern, but
  // with the opposite sign: b = +V here, not -V (see
  // EwaldPeriodicDipoleOperator's own class documentation for why).
  auto t_field = std::chrono::steady_clock::now();
  std::vector<std::pair<Index, PolarSite*>> targets;
  targets.reserve(std::size_t(total_size / 3));
  for (std::size_t n = 0; n < ids.size(); ++n) {
    PolarSegment& segment = registry.Get(ids[n], EwaldChargeState::Neutral);
    for (Index s = 0; s < segment.size(); ++s) {
      PolarSite& site = segment[s];
      site.setInduced_Dipole(Eigen::Vector3d::Zero());
      site.Reset();
      targets.push_back({ids[n], &site});
    }
  }
  for (const auto& entry : targets) {
    real_sum.AddFieldAt<Estatic::V>(entry.first, *entry.second,
                                    EwaldChargeState::Neutral);
  }
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Real-space permanent field done ("
      << elapsed_s(t_field) << "s)" << std::flush;

  auto t_recip = std::chrono::steady_clock::now();
  // EwaldReciprocalSpaceSum no longer takes a segment id (it never
  // excludes anything -- see its own class documentation), so only the
  // bare site pointers are needed here.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  recip_sum.AddFieldAtMany<Estatic::V>(
      recip_targets, EwaldChargeState::Neutral,
      [&](std::size_t done, std::size_t total) {
        XTP_LOG(Log::info, log)
            << TimeStamp() << "   k-space progress: " << done << "/"
            << total << " sites (" << elapsed_s(t_recip) << "s)"
            << std::flush;
      });
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Reciprocal-space permanent field done ("
      << elapsed_s(t_recip) << "s)" << std::flush;

  for (const auto& entry : targets) {
    shape.AddFieldAt<Estatic::V>(*entry.second, EwaldChargeState::Neutral);
  }
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Permanent field total (" << elapsed_s(t_field)
      << "s)" << std::flush;

  // Intramolecular static-static compensation: EwaldReciprocalSpaceSum's
  // own structure factor never excludes anything (see that class's own
  // documentation), so it unconditionally includes every same-segment
  // static-static pair's own contribution -- exactly the erf(alpha*r)/r
  // "other half" of the Ewald split (erf+erfc=1). Legacy explicitly
  // removes this same leak (see EwdInteractor::FP12_ERF_At_By's own
  // "Note the (-): This is a compensation term" comment) so that the net
  // static-static intramolecular contribution comes out to genuinely
  // zero, not the un-cancelled reciprocal-space leak alone -- see
  // EwaldRealSpaceInteractor::ApplyIntramolecularStaticCorrection's own
  // documentation for the fuller account of why (a real mechanism in
  // legacy's own code, missed on an earlier pass through it, not a new
  // design decision here). Every ordered pair within a segment
  // (i receiving j's own correction, and j receiving i's own) is
  // visited, matching legacy's own double loop.
  EwaldRealSpaceInteractor erf_interactor(alpha_, thole_a_);
  for (Index n = 0; n < Index(ids.size()); ++n) {
    PolarSegment& segment = registry.Get(ids[n], EwaldChargeState::Neutral);
    Index n_sites = segment.size();
    if (n_sites < 2) {
      continue;
    }
    for (Index i = 0; i < n_sites; ++i) {
      for (Index j = 0; j < n_sites; ++j) {
        if (i == j) {
          continue;
        }
        erf_interactor
            .ApplyIntramolecularStaticCorrection<PolarSite, Estatic::V>(
                segment[j], segment[i]);
      }
    }
  }
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Intramolecular static compensation applied ("
      << elapsed_s(t_field) << "s)" << std::flush;

  Eigen::VectorXd b(total_size);
  {
    std::size_t t = 0;
    for (std::size_t n = 0; n < ids.size(); ++n) {
      const PolarSegment& segment =
          registry.Get(ids[n], EwaldChargeState::Neutral);
      Index base = offsets[n];
      for (Index s = 0; s < segment.size(); ++s) {
        b.segment<3>(base + 3 * s) = targets[t].second->V();
        targets[t].second->Reset();
        ++t;
      }
    }
  }

  if (induce_) {
    XTP_LOG(Log::info, log)
        << TimeStamp() << " Starting " << (use_jor_ ? "JOR" : "PCG")
        << " solve (max " << max_iter_ << " iterations, tolerance "
        << pcg_tolerance_ << ")" << std::flush;
    auto t_pcg = std::chrono::steady_clock::now();

    EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape, ids,
                                   alpha_, thole_a_,
                                   apply_shape_correction_to_induced_,
                                   apply_self_field_correction_to_induced_,
                                   apply_thole_damping_intramolecular_,
                                   apply_realspace_intermolecular_coupling_,
                                   apply_reciprocal_coupling_);
    // Debug/experimental -- see DumpIntraPairThole's own declaration.
    // Called here, right after op's own construction and before the
    // solve itself, since this data depends only on geometry and
    // polarizability, not on anything the solve computes.
    if (debug_dump_first_iteration_and_stop_) {
      op.DumpIntraPairThole("new_intrapair_dump.csv");
      // Debug/experimental -- see EwaldRealSpaceSum::DumpNeighborList's
      // own declaration. Dumps every ids segment as its own target (a
      // first version of this restricted to a single segment risked
      // hiding a genuine PBC-scheme discrepancy that only shows up for
      // segments in a different geometric relationship to the box), so
      // target_segment_id is written as its own column -- matching
      // legacy's own equivalent dump's own choice, made for the same
      // reason.
      {
        std::ofstream nb_header("new_neighborlist_dump.csv");
        nb_header << "target_segment_id,source_segment_id,r_bohr,t_x_bohr,"
                     "t_y_bohr,t_z_bohr\n";
      }
      for (Index id : ids) {
        const PolarSegment& segment =
            registry.Get(id, EwaldChargeState::Neutral);
        real_sum.DumpNeighborListAppend(id, segment[0],
                                       EwaldChargeState::Neutral,
                                       "new_neighborlist_dump.csv");
      }
    }
    Eigen::VectorXd x;
    Index iterations;
    double residual;
    bool converged;
    Index indefinite_at_iteration = -1;
    double indefinite_curvature = 0.0;
    double lanczos_min_eigenvalue = std::numeric_limits<double>::quiet_NaN();

    // use_jor_ is an outer, first choice: JOR and PCG are two entirely
    // different algorithms (see SolveWithJOR's own documentation for
    // why JOR is the one actually chosen once the operator was directly
    // confirmed indefinite), not two variants of the same solve -- JOR
    // has no preconditioner-type branching of its own (EwaldSitePolarizabilityBlocks
    // is the one and only D^-1 it uses, matching legacy exactly), so it
    // never enters the PCG-specific branches below at all.
    if (use_jor_) {
      EwaldSitePolarizabilityBlocks site_p(registry, ids);
      auto result =
          SolveWithJOR(op, site_p, b, max_iter_, jor_omega_, log, t_pcg,
                      registry, ids, debug_dump_first_iteration_and_stop_,
                      match_legacy_first_step_);
      x = result.x;
      iterations = result.iterations;
      residual = result.residual;
      converged = result.converged;

      XTP_LOG(Log::info, log)
          << TimeStamp() << " JOR finished after " << iterations
          << " iterations, max_dU=" << result.max_dU
          << " avg_dU=" << result.avg_dU << " (residual=" << residual
          << ", informational) (" << elapsed_s(t_pcg) << "s)" << std::flush;
      if (!converged) {
        throw std::runtime_error(
            "EwaldBackground: JOR did not converge within max_iter");
      }
    } else {

    // The preconditioner type is a compile-time choice here too (see
    // SolveWithIndefinitenessCheck's own documentation for why that
    // function is templated on it) -- these two branches genuinely
    // instantiate two different specializations, they cannot share one
    // call. EwaldBlockJacobiPreconditioner builds itself fully in its
    // own constructor (see its own class documentation); Eigen's own
    // DiagonalPreconditioner does not -- it needs an explicit compute(op)
    // call first, unlike the constructor-based pattern
    // EwaldBlockJacobiPreconditioner itself uses. Getting this backwards
    // (assuming both work the same way) would be a real, silent
    // correctness bug -- solve() on an uncompute()'d DiagonalPreconditioner
    // does not throw, it just returns nonsense -- so this is deliberately
    // NOT written as a single shared code path that "just" swaps the
    // preconditioner type.
    if (use_block_jacobi_preconditioner_) {
      EwaldBlockJacobiPreconditioner precond(registry, ids, alpha_, thole_a_);
      auto result = SolveWithIndefinitenessCheck(op, precond, b, max_iter_,
                                                 pcg_tolerance_, log, t_pcg);
      x = result.x;
      iterations = result.iterations;
      residual = result.residual;
      converged = result.converged;
      indefinite_at_iteration = result.indefinite_at_iteration;
      indefinite_curvature = result.indefinite_curvature;
      lanczos_min_eigenvalue = result.lanczos_min_eigenvalue;
    } else {
      Eigen::DiagonalPreconditioner<double> precond;
      precond.compute(op);
      auto result = SolveWithIndefinitenessCheck(op, precond, b, max_iter_,
                                                 pcg_tolerance_, log, t_pcg);
      x = result.x;
      iterations = result.iterations;
      residual = result.residual;
      converged = result.converged;
      indefinite_at_iteration = result.indefinite_at_iteration;
      indefinite_curvature = result.indefinite_curvature;
      lanczos_min_eigenvalue = result.lanczos_min_eigenvalue;
    }

    XTP_LOG(Log::info, log)
        << TimeStamp() << " PCG finished after " << iterations
        << " iterations, residual " << residual << " (" << elapsed_s(t_pcg)
        << "s)" << std::flush;
    if (!std::isnan(lanczos_min_eigenvalue)) {
      XTP_LOG(Log::info, log)
          << TimeStamp() << " Lanczos min eigenvalue estimate (from this "
             "run's own alpha/beta coefficients, using the whole "
             "accumulated Krylov subspace -- see "
             "PcgIndefinitenessResult::lanczos_min_eigenvalue's own "
             "documentation for what this can and cannot guarantee): "
          << lanczos_min_eigenvalue << std::flush;
      if (lanczos_min_eigenvalue < 0.0) {
        XTP_LOG(Log::info, log)
            << TimeStamp()
            << " This is NEGATIVE -- strong evidence (not an absolute "
               "guarantee; see this run's own documentation) that the "
               "operator is not positive-definite."
            << std::flush;
      }
    }

    if (indefinite_at_iteration >= 0) {
      throw std::runtime_error(
          "EwaldBackground: PCG's own operator was found to be NOT "
          "positive-definite at iteration " +
          std::to_string(indefinite_at_iteration) + " (p.A.p = " +
          std::to_string(indefinite_curvature) +
          " <= 0) -- this is a direct algebraic certificate, not an "
          "inference from residual behavior.");
    }

    if (!converged) {
      throw std::runtime_error(
          "EwaldBackground: PCG did not converge (max_iter reached, "
          "operator was never found indefinite along the way)");
    }
    }

    for (std::size_t n = 0; n < ids.size(); ++n) {
      PolarSegment& segment =
          registry.Get(ids[n], EwaldChargeState::Neutral);
      Index base = offsets[n];
      for (Index s = 0; s < segment.size(); ++s) {
        segment[s].setInduced_Dipole(x.segment<3>(base + 3 * s));
      }
    }
  } else {
    // Matches legacy's own polarmethod.induce=0 behaviour: every site's
    // induced dipole stays at zero (already set before the
    // permanent-field computation above) -- see induce_'s own class
    // documentation for why this is a separate code path rather than a
    // zero-polarizability run through the same PCG solve.
    XTP_LOG(Log::info, log)
        << TimeStamp()
        << " polarmethod.induce is false: skipping the PCG solve, "
           "induced dipoles left at zero"
        << std::flush;
  }

  auto t_cpt = std::chrono::steady_clock::now();
  // Write the permanent field (b) back into every site's V() before the
  // checkpoint write, unconditionally regardless of induce_ -- neither
  // code path above leaves V() holding it otherwise: the PCG path
  // overwrites V() repeatedly during its own iterations (each
  // EwaldPeriodicDipoleOperator::multiply() call internally
  // Reset()s and rewrites every site's V() via RawMultiply), leaving
  // whatever the *last* iteration's trial field happened to be, not the
  // permanent field; the induce_=false path leaves V() at the zero
  // Reset() left it just after b was built above. Restoring b here is
  // what makes the permanent field actually recoverable from the
  // checkpoint at all -- without this, hdf5_dump has no way to report
  // anything but 0.0 for it (see that tool's own former limitation note,
  // now resolved by this).
  {
    for (std::size_t n = 0; n < ids.size(); ++n) {
      PolarSegment& segment = registry.Get(ids[n], EwaldChargeState::Neutral);
      Index base = offsets[n];
      for (Index s = 0; s < segment.size(); ++s) {
        segment[s].V() = b.segment<3>(base + 3 * s);
      }
    }
  }

  CheckpointFile cpf(checkpoint_file_, CheckpointAccessLevel::CREATE);
  CheckpointWriter w = cpf.getWriter();
  registry.WriteToCpt(w);

  XTP_LOG(Log::info, log)
      << TimeStamp() << " Checkpoint written to " << checkpoint_file_
      << " (" << elapsed_s(t_cpt) << "s)" << std::flush;
  XTP_LOG(Log::info, log)
      << TimeStamp() << " Ewald background calculation done, total "
      << elapsed_s(t_start) << "s" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDBACKGROUND_H
