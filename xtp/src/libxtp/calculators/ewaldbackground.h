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
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/tools/constants.h"
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/ewaldblockjacobipreconditioner.h"
#include "votca/xtp/ewaldperiodicdipoleoperator.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/ewaldparameters.h"
#include "votca/xtp/ewaldshapecorrection.h"
#include "votca/xtp/ewaldsolvers.h"
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
  // Dimensionless: the real-space distance cutoff is
  // screening_factor_/alpha. See EwaldRealSpaceSum's own constructor for
  // what it controls and what raising or lowering it costs. Unitless, so
  // unlike r_min/alpha/k_max it needs no nm<->bohr conversion.
  double screening_factor_ = 6.0;
  EwaldShape shape_ = EwaldShape::Cube;

  Index max_iter_ = 200;
  double pcg_tolerance_ = 1e-8;
  bool induce_ = true;  // matches legacy's own polarmethod.induce; when
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
  screening_factor_ = opt.ifExistsReturnElseReturnDefault<double>(
      "realspace.screening_factor", screening_factor_);
  if (screening_factor_ <= 0.0) {
    throw std::runtime_error(
        "EwaldBackground: realspace.screening_factor must be positive");
  }

  max_iter_ = opt.ifExistsReturnElseReturnDefault<Index>(
      "polarmethod.max_iter", max_iter_);
  pcg_tolerance_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.tolerance", pcg_tolerance_);
  induce_ = opt.ifExistsReturnElseReturnDefault<bool>("polarmethod.induce",
                                                       induce_);
  // Debug-only, undocumented on purpose (no legacy counterpart to give
  // it a natural home in the schema) -- see this member's own
  // declaration for what it's for.
  // Debug-only, undocumented on purpose, same pattern as the shape one
  // immediately above.
  // Debug-only, undocumented on purpose, same pattern as the two above.
  // Debug-only, undocumented on purpose, same pattern as the three above.
  // Debug-only, undocumented on purpose, same pattern as the four above.
  // Debug/experimental, undocumented on purpose -- see
  // use_block_jacobi_preconditioner_'s own declaration.
  use_block_jacobi_preconditioner_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.use_block_jacobi_preconditioner",
      use_block_jacobi_preconditioner_);
  // Debug/experimental, see use_jor_'s own declaration.
  use_jor_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.use_jor", use_jor_);
  jor_omega_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.jor_omega", jor_omega_);
  match_legacy_first_step_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_match_legacy_first_step", match_legacy_first_step_);
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
                            field_tol_, /*shell_width=*/0.945,
                            /*n_max=*/15, screening_factor_);
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
      << " nm), field_tol=" << field_tol_
      << ", screening_factor=" << screening_factor_ << " (real-space cutoff "
      << screening_factor_ / alpha_ << " bohr)"
      << ", shape="
      << (shape_ == EwaldShape::Cube ? "cube" : "slab")
      << ", box volume=" << volume << " bohr^3, max_iter=" << max_iter_
      << ", pcg_tolerance=" << pcg_tolerance_
      << ", induce=" << (induce_ ? "true" : "false")
      << ", use_block_jacobi_preconditioner="
      << (use_block_jacobi_preconditioner_ ? "true" : "false")
      << ", use_jor=" << (use_jor_ ? "true" : "false")
      << ", jor_omega=" << jor_omega_
      << ", debug_match_legacy_first_step="
      << (match_legacy_first_step_ ? "true" : "false")
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
            << total << " k-vectors (" << elapsed_s(t_recip) << "s)"
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
  // EwaldRealSpaceInteractor::ApplyErfStaticFieldCorrection's own
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
            .ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(
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
                                   alpha_, thole_a_);
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
    {
      // Phase breakdown of the matvec -- see
      // EwaldPeriodicDipoleOperator::RawMultiplyTimings' own
      // declaration. n_calls exceeds the reported iteration count by
      // one, since baseline_ = RawMultiply(0) is built in the operator's
      // own constructor. The first call also builds the real-space
      // neighbour cache, so it is much more expensive than the rest and
      // inflates the real_space per-call average; read the per-call
      // figures as an upper bound on steady-state cost.
      const auto& tm = op.Timings();
      const double n = double(std::max<Index>(tm.n_calls, 1));
      XTP_LOG(Log::info, log)
          << TimeStamp() << " RawMultiply phase breakdown over "
          << tm.n_calls << " calls (total " << tm.total() << "s):"
          << std::flush;
      auto line = [&](const char* nm, double t) {
        XTP_LOG(Log::info, log)
            << TimeStamp() << (boost::format("   %1$-12s %2$8.3fs total  "
                                             "%3$7.3fs/call  %4$5.1f%%") %
                               nm % t % (t / n) %
                               (tm.total() > 0 ? 100.0 * t / tm.total() : 0.0))
                                  .str()
            << std::flush;
      };
      line("setup", tm.setup);
      line("real_space", tm.real_space);
      line("reciprocal", tm.reciprocal);
      line("shape", tm.shape);
      line("assemble", tm.assemble);
      line("intra", tm.intra);

      // Neighbour-list size -- the real cost driver behind the
      // real_space line above. See EwaldRealSpaceSum::NeighborStats.
      const auto ns = real_sum.GetNeighborStats();
      XTP_LOG(Log::info, log)
          << TimeStamp()
          << (boost::format(
                  "   neighbours: %1$.1f entries/target over %2$d targets "
                  "(%3$d kept, %4$d culled beyond %5$.1f bohr = %6$.1f%%)") %
              ns.entries_per_target() % ns.targets % ns.entries % ns.culled %
              real_sum.RealSpaceCutoff() % (100.0 * ns.culled_fraction()))
                 .str()
          << std::flush;
    }
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

  // The convergence parameters travel with the converged state. A job
  // that later embeds a foreground in this background must use the same
  // alpha, k_max and shape -- alpha in particular decides how the
  // interaction is split between the real- and reciprocal-space sums, so
  // a different value is not the same physics, and nothing about the
  // mismatch would be visible at run time. See EwaldParameters.
  {
    EwaldParameters params;
    params.alpha = alpha_;
    params.k_max = k_max_;
    params.r_min = r_min_;
    params.field_tol = field_tol_;
    params.thole_a = thole_a_;
    params.screening_factor = screening_factor_;
    params.shape = shape_;
    params.box = box;
    CheckpointWriter wp = w.openChild("ewald_parameters");
    params.WriteToCpt(wp);
  }

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
