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
#include "votca/xtp/ewaldparameters.h"
#include "votca/xtp/ewaldperiodicdipoleoperator.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/ewaldshapecorrection.h"
#include "votca/xtp/ewaldsolvers.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/segmentmapper.h"

/**
 * \brief ewdbgpol-equivalent calculator built on EwaldRegistry,
 *        EwaldRealSpaceSum, EwaldReciprocalSpaceSum, EwaldShapeCorrection
 *        and EwaldPeriodicDipoleOperator instead of the legacy xtp/ewald/
 *        code. Registered as "ewaldbackground" so that both can coexist.
 *
 * Converges the self-consistent, fully polarized (neutral-state)
 * background of the whole periodic system, and writes it to an HDF5
 * checkpoint file via EwaldRegistry::WriteToCpt.
 *
 * Length-related options (coulombmethod.alpha, coulombmethod.k_max,
 * realspace.r_min) are given in nm or nm^-1 and converted to bohr or
 * bohr^-1 in ParseOptions; everything downstream uses atomic units.
 *
 * Segments may have any number of sites; b and x use the variable-width
 * layout of EwaldPeriodicDipoleOperator.
 *
 * The outer "loop until 2nd-order fields converged" of legacy
 * PolarBackground::Polarize() is not reproduced: the iterative solve
 * (PCG or JOR) is that loop.
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

  double alpha_ = 0.5;  // bohr^-1. XML value (coulombmethod.alpha) is
                        // in nm^-1. Unless set explicitly, derived from
                        // the box in Evaluate().
  bool alpha_explicit_ = false;
  double k_max_ = 3.0;  // bohr^-1, XML value in nm^-1. Unless set
                        // explicitly, derived from alpha_ in Evaluate().
  bool k_max_explicit_ = false;
  double thole_a_ = 0.39;
  double r_min_ = 18.897259886;  // bohr (= 1.0 nm, the default applied
                                 // in ParseOptions). XML value
                                 // (realspace.r_min) is in nm.
  double field_tol_ = 1e-8;
  // Dimensionless: the real-space distance cutoff is
  // screening_factor_/alpha (see the EwaldRealSpaceSum constructor).
  double screening_factor_ = 6.0;
  EwaldShape shape_ = EwaldShape::Cube;

  Index max_iter_ = 200;
  double pcg_tolerance_ = 1e-8;
  bool induce_ = true;  // legacy polarmethod.induce; if false, induced
                        // dipoles stay zero and no solve is run.
  // Debug/experimental, no legacy counterpart: use
  // EwaldBlockJacobiPreconditioner instead of the diagonal one.
  bool use_block_jacobi_preconditioner_ = false;
  // Debug/experimental. Use weighted-Jacobi iteration (SolveWithJOR,
  // equivalent to legacy PolarBackground's SOR) instead of PCG. JOR does
  // not require a positive-definite operator, whereas CG fails if
  // SolveWithIndefinitenessCheck finds p.A.p <= 0.
  bool use_jor_ = false;
  // JOR relaxation factor. Legacy PolarBackground default
  // (_polar_wSOR_N = 0.35), which overrides the 0.25 default of
  // APolarSite::Induce.
  double jor_omega_ = 0.35;
  // Debug/experimental. Unrelaxed first JOR step, replicating legacy's
  // iteration history (see SolveWithJOR).
  bool match_legacy_first_step_ = false;
};

inline void EwaldBackground::ParseOptions(const tools::Property& opt) {
  mapping_file_ = opt.get("multipoles").as<std::string>();
  checkpoint_file_ = opt.ifExistsReturnElseReturnDefault<std::string>(
      "output.checkpoint", "ewaldbackground.hdf5");

  alpha_explicit_ = opt.exists("coulombmethod.alpha");
  if (alpha_explicit_) {
    // Inverse length: value in bohr^-1 = value in nm^-1 * bohr2nm.
    alpha_ = opt.get("coulombmethod.alpha").as<double>() * tools::conv::bohr2nm;
  }
  // Unless given explicitly, alpha_ and k_max_ are derived in Evaluate(),
  // where the box is known; k_max_ is always a multiple of alpha_, since
  // the reciprocal-space weight is exp(-k^2/4 alpha^2). (A fixed k_max
  // unrelated to alpha produced orders of magnitude too many k-vectors.)
  k_max_explicit_ = opt.exists("coulombmethod.k_max");
  if (k_max_explicit_) {
    // nm^-1 -> bohr^-1, as for alpha_.
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

  thole_a_ = opt.ifExistsReturnElseReturnDefault<double>("polarmethod.thole_a",
                                                         thole_a_);
  // Length: bohr = nm * nm2bohr. Default 1.0 nm (~18.9 bohr).
  double r_min_nm =
      opt.ifExistsReturnElseReturnDefault<double>("realspace.r_min", 1.0);
  r_min_ = r_min_nm * tools::conv::nm2bohr;
  field_tol_ = opt.ifExistsReturnElseReturnDefault<double>(
      "realspace.field_tol", field_tol_);
  screening_factor_ = opt.ifExistsReturnElseReturnDefault<double>(
      "realspace.screening_factor", screening_factor_);
  if (screening_factor_ <= 0.0) {
    throw std::runtime_error(
        "EwaldBackground: realspace.screening_factor must be positive");
  }

  max_iter_ = opt.ifExistsReturnElseReturnDefault<Index>("polarmethod.max_iter",
                                                         max_iter_);
  pcg_tolerance_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.tolerance", pcg_tolerance_);
  induce_ =
      opt.ifExistsReturnElseReturnDefault<bool>("polarmethod.induce", induce_);
  // Debug/experimental options, intentionally not in the XML schema
  // (see the member declarations).
  use_block_jacobi_preconditioner_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.use_block_jacobi_preconditioner",
      use_block_jacobi_preconditioner_);
  use_jor_ = opt.ifExistsReturnElseReturnDefault<bool>("polarmethod.use_jor",
                                                       use_jor_);
  jor_omega_ = opt.ifExistsReturnElseReturnDefault<double>(
      "polarmethod.jor_omega", jor_omega_);
  match_legacy_first_step_ = opt.ifExistsReturnElseReturnDefault<bool>(
      "polarmethod.debug_match_legacy_first_step", match_legacy_first_step_);
}

inline bool EwaldBackground::Evaluate(Topology& top) {
  Logger log;
  log.setReportLevel(Log::info);
  log.setCommonPreface("\nEWD");
  // Single (non-job-parallel) calculator: print log messages immediately
  // instead of buffering them until the end, so progress is visible.
  log.setMultithreading(true);

  auto t_start = std::chrono::steady_clock::now();
  auto elapsed_s = [&](std::chrono::steady_clock::time_point since) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                         since)
        .count();
  };

  XTP_LOG(Log::info, log) << TimeStamp()
                          << " Starting Ewald background calculation"
                          << std::flush;

  PolarMapper polmap(log);
  polmap.LoadMappingFile(mapping_file_);

  EwaldRegistry registry;
  std::vector<Index> ids;
  // offsets[n]: scalar index where segment ids[n] starts (3 per site),
  // the same layout as EwaldPeriodicDipoleOperator::offsets_.
  std::vector<Index> offsets;
  ids.reserve(top.Segments().size());
  offsets.reserve(top.Segments().size() + 1);
  offsets.push_back(0);
  for (const Segment& seg : top.Segments()) {
    PolarSegment mol = polmap.map(seg, SegId(seg.getId(), std::string("n")));
    registry.Register(seg.getId(), EwaldChargeState::Neutral, mol);
    ids.push_back(seg.getId());
    offsets.push_back(offsets.back() + 3 * mol.size());
  }
  const Index total_size = offsets.back();

  XTP_LOG(Log::info, log) << TimeStamp() << " Mapped " << ids.size()
                          << " segments, " << (total_size / 3)
                          << " polarizable sites total (" << elapsed_s(t_start)
                          << "s)" << std::flush;

  const Eigen::Matrix3d& box = top.getBox();
  const double volume = box.col(0).dot(box.col(1).cross(box.col(2)));

  // alpha and k_max are derived together by balancing the cost of the two
  // sums; the result is independent of alpha (tested). Both cutoffs are
  // fixed multiples of alpha, so the accuracy does not depend on it:
  //
  //   real space  EwaldRealSpaceSum truncates at screening_factor/alpha,
  //               leaving a tail of order erfc(screening_factor).
  //   reciprocal  the Gaussian weight is exp(-k^2/4 alpha^2), so a cutoff
  //               at s_k*alpha leaves exp(-s_k^2/4).
  //
  // so with N sites in a cell of volume V the work at ONE target is
  //
  //   N_real = (N/V)(4pi/3)(s_r/alpha)^3      falling as alpha^-3
  //   N_k    = (4pi/3)(s_k alpha)^3 V/(8pi^3)   rising as alpha^3
  //
  // and d/dalpha [A alpha^-3 + B alpha^3] = 0 gives alpha^6 = A/B, i.e.
  //
  //   alpha = sqrt(2 pi) (N^(1/6) / V^(1/3)) sqrt(s_r/s_k)
  //
  // the standard Ewald result, whose N^(1/6) is what makes the method
  // O(N^(3/2)) rather than O(N^2).
  //
  // The former defaults alpha = 3/L_min, k_max = 6*alpha fixed the
  // reciprocal sum at ~98 k-vectors for any cubic cell while the
  // real-space cutoff (2L) grew with the system (thiophene: 301,691 terms
  // per target at 1000 segments, 1,508,063 at 5000). They also truncated
  // the two halves at very different accuracies (erfc(6) = 2e-17 vs
  // exp(-9) = 1.2e-4).
  //
  // screening_factor s_r is therefore the single accuracy parameter:
  // eps = erfc(s_r) and s_k = 2*sqrt(-ln eps) give the same reciprocal
  // tail (s_k = 12.39 for s_r = 6). With the balanced alpha this is 3x
  // (216 methane) to 62x (400,000 sites) faster than the old defaults.
  //
  // Setting <alpha> or <k_max> explicitly still overrides either half.
  if (!alpha_explicit_ || !k_max_explicit_) {
    const double s_r = screening_factor_;
    // erfc underflows to 0 for s_r beyond ~27; clamp to keep s_k finite.
    const double eps = std::max(std::erfc(s_r), 1e-300);
    const double s_k = 2.0 * std::sqrt(-std::log(eps));

    if (!alpha_explicit_) {
      const double n_sites = double(total_size / 3);
      alpha_ = std::sqrt(2.0 * tools::conv::Pi) * std::pow(n_sites, 1.0 / 6.0) /
               std::cbrt(volume) * std::sqrt(s_r / s_k);

      // Not clamped against r_min_: r_min_ only bounds the shell search
      // by translation magnitude, while the cutoff culls by pair
      // separation, so accuracy is erfc(s_r) either way. An r_min_
      // beyond the cutoff only costs time; this is reported below.
    }
    if (!k_max_explicit_) {
      k_max_ = s_k * alpha_;
    }
  }

  EwaldRealSpaceSum real_sum(box, registry, alpha_, thole_a_, r_min_,
                             field_tol_, /*shell_width=*/0.945,
                             /*n_max=*/15, screening_factor_);
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_, k_max_);
  EwaldShapeCorrection shape(volume, registry, shape_);

  XTP_LOG(Log::info, log) << TimeStamp()
                          << " Real/reciprocal-space sums constructed ("
                          << elapsed_s(t_start) << "s)" << std::flush;
  XTP_LOG(Log::info, log) << TimeStamp() << " Ewald split: alpha=" << alpha_
                          << " bohr^-1 (" << alpha_ * tools::conv::nm2bohr
                          << " nm^-1), k_max=" << k_max_ << " bohr^-1 ("
                          << k_max_ * tools::conv::nm2bohr << " nm^-1, "
                          << recip_sum.NumKVectors()
                          << " k-vectors), real-space cutoff "
                          << screening_factor_ / alpha_
                          << " bohr (screening_factor=" << screening_factor_
                          << ")" << std::flush;
  // With the derived alpha the two halves should be of comparable size;
  // a large imbalance indicates a hand-set alpha or an unsuitable cell.
  {
    const double n_sites = double(total_size / 3);
    const double n_real_est = (n_sites / volume) *
                              (4.0 * tools::conv::Pi / 3.0) *
                              std::pow(screening_factor_ / alpha_, 3);
    XTP_LOG(Log::info, log)
        << TimeStamp()
        << " Ewald balance (terms per target, estimated): " << "real "
        << n_real_est << " vs reciprocal " << recip_sum.NumKVectors()
        << std::flush;

    // r_min_ beyond the cutoff forces the search through shells whose
    // pairs are all culled: no effect on accuracy, only on cost.
    const double cutoff = screening_factor_ / alpha_;
    if (r_min_ > cutoff) {
      XTP_LOG(Log::info, log)
          << TimeStamp() << " NOTE: r_min (" << r_min_
          << " bohr) lies beyond the real-space cutoff (" << cutoff
          << " bohr). The shell search will walk shells that contribute "
             "nothing before it is allowed to stop. Accuracy is unaffected "
             "-- it is erfc(screening_factor) either way -- but lowering "
             "r_min to at most the cutoff would save that work."
          << std::flush;
    }
  }
  XTP_LOG(Log::info, log) << TimeStamp() << " Cell: volume=" << volume
                          << " bohr^3, shape="
                          << (shape_ == EwaldShape::Cube ? "cube" : "slab")
                          << ", r_min=" << r_min_ << " bohr ("
                          << r_min_ * tools::conv::bohr2nm
                          << " nm), field_tol=" << field_tol_
                          << ", thole_a=" << thole_a_ << std::flush;
  XTP_LOG(Log::info, log) << TimeStamp() << " Induction: induce="
                          << (induce_ ? "true" : "false")
                          << ", max_iter=" << max_iter_
                          << ", tolerance=" << pcg_tolerance_
                          << ", solver=" << (use_jor_ ? "JOR" : "PCG")
                          << (use_jor_ ? ""
                                       : (use_block_jacobi_preconditioner_
                                              ? " (block-Jacobi)"
                                              : " (unpreconditioned)"))
                          << std::flush;
  if (use_jor_) {
    XTP_LOG(Log::info, log)
        << TimeStamp() << " JOR omega=" << jor_omega_ << std::flush;
  }
  if (match_legacy_first_step_) {
    XTP_LOG(Log::info, log)
        << TimeStamp() << " debug_match_legacy_first_step is ON" << std::flush;
  }

  // Permanent field (all induced dipoles zero), computed once; it is the
  // right-hand side b = +V of the solve (opposite sign to
  // CalcInducedDipolesViaPCG, see EwaldPeriodicDipoleOperator).
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
  // Builds every neighbour list in parallel and applies it; the operator
  // constructed below reuses these lists.
  real_sum.AddFieldAtMany<Estatic::V>(targets, EwaldChargeState::Neutral);
  XTP_LOG(Log::info, log) << TimeStamp() << " Real-space permanent field done ("
                          << elapsed_s(t_field) << "s)" << std::flush;

  auto t_recip = std::chrono::steady_clock::now();
  // The reciprocal-space sum excludes nothing, so it needs no segment ids.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  recip_sum.AddFieldAtMany<Estatic::V>(
      recip_targets, EwaldChargeState::Neutral,
      [&](std::size_t done, std::size_t total) {
        // Debug level: the k-vector count is already logged above.
        XTP_LOG(Log::debug, log)
            << TimeStamp() << "   k-space progress: " << done << "/" << total
            << " k-vectors (" << elapsed_s(t_recip) << "s)" << std::flush;
      });
  XTP_LOG(Log::info, log) << TimeStamp()
                          << " Reciprocal-space permanent field done ("
                          << elapsed_s(t_recip) << "s)" << std::flush;

  // The shape field -(4 pi / 3V) M (cube) is the same at every site, but
  // EwaldShapeCorrection::AddFieldAt sums M over the registry on every
  // call (O(N^2) per site, ~110 s for 70k sites). Evaluate once on a
  // probe and add, as EwaldPeriodicDipoleOperator does.
  auto t_shape = std::chrono::steady_clock::now();
  {
    PolarSite shape_probe(-1, "X", Eigen::Vector3d::Zero());
    shape.AddFieldAt<Estatic::V>(shape_probe, EwaldChargeState::Neutral);
    const Eigen::Vector3d shape_field = shape_probe.V();
    for (const auto& entry : targets) {
      entry.second->V() += shape_field;
    }
  }
  XTP_LOG(Log::info, log) << TimeStamp() << " Shape field done ("
                          << elapsed_s(t_shape) << "s)" << std::flush;
  // Total time since the start of the real-space pass.
  XTP_LOG(Log::info, log) << TimeStamp() << " Permanent field total ("
                          << elapsed_s(t_field) << "s)" << std::flush;
  auto t_intra = std::chrono::steady_clock::now();

  // Intramolecular static-static compensation: the reciprocal-space sum
  // excludes nothing, so it contains the erf(alpha*r)/r part of every
  // same-segment permanent pair. As legacy EwdInteractor::FP12_ERF_At_By
  // does, this is subtracted so that the net intramolecular permanent
  // interaction is zero (see
  // EwaldRealSpaceInteractor::ApplyErfStaticFieldCorrection). All ordered
  // pairs are visited, as in legacy's double loop.
  //
  // This includes i == j: the erf self-field of a permanent dipole is
  // -(4/3)alpha^3/sqrt(pi) * mu (zero for a charge) and would otherwise
  // bias the induced dipoles by an alpha^3-dependent amount. It is the
  // permanent-field counterpart of SelfFieldMatrix() in
  // EwaldPeriodicDipoleOperator; the r < kCoincidenceTol branch of
  // ApplyErfStaticFieldCorrection handles this limit.
  EwaldRealSpaceInteractor erf_interactor(alpha_, thole_a_);
  for (Index n = 0; n < Index(ids.size()); ++n) {
    PolarSegment& segment = registry.Get(ids[n], EwaldChargeState::Neutral);
    Index n_sites = segment.size();
    for (Index i = 0; i < n_sites; ++i) {
      for (Index j = 0; j < n_sites; ++j) {
        erf_interactor.ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(
            segment[j], segment[i]);
      }
    }
  }
  XTP_LOG(Log::info, log) << TimeStamp()
                          << " Intramolecular static compensation applied ("
                          << elapsed_s(t_intra) << "s)" << std::flush;

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

    // JOR uses EwaldSitePolarizabilityBlocks as D^-1 (as legacy does) and
    // has no preconditioner choice.
    if (use_jor_) {
      EwaldSitePolarizabilityBlocks site_p(registry, ids);
      auto result = SolveWithJOR(op, site_p, b, max_iter_, jor_omega_, log,
                                 t_pcg, match_legacy_first_step_);
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

      // The preconditioner type is a template parameter, hence two
      // branches. EwaldBlockJacobiPreconditioner is fully built by its
      // constructor, whereas DiagonalPreconditioner requires compute(op)
      // before use.
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
        // Matvec phase breakdown (RawMultiplyTimings). n_calls includes
        // the baseline_ call in the operator constructor.
        const auto& tm = op.Timings();
        const double n = double(std::max<Index>(tm.n_calls, 1));
        XTP_LOG(Log::info, log)
            << TimeStamp() << " RawMultiply phase breakdown over " << tm.n_calls
            << " calls (total " << tm.total() << "s):" << std::flush;
        auto line = [&](const char* nm, double t) {
          XTP_LOG(Log::info, log)
              << TimeStamp()
              << (boost::format("   %1$-12s %2$8.3fs total  "
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

        // Neighbour-list size, which determines the real_space cost
        // (EwaldRealSpaceSum::NeighborStats).
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
        XTP_LOG(Log::info, log) << TimeStamp() << " Lanczos min eigenvalue: "
                                << lanczos_min_eigenvalue << std::flush;
        if (lanczos_min_eigenvalue < 0.0) {
          XTP_LOG(Log::info, log) << TimeStamp()
                                  << " NEGATIVE -- the operator is likely not "
                                     "positive-definite."
                                  << std::flush;
        }
      }

      if (indefinite_at_iteration >= 0) {
        throw std::runtime_error(
            "EwaldBackground: PCG operator is NOT positive-definite "
            "(found at iteration " +
            std::to_string(indefinite_at_iteration) +
            ", p.A.p = " + std::to_string(indefinite_curvature) + " <= 0).");
      }

      if (!converged) {
        throw std::runtime_error(
            "EwaldBackground: PCG did not converge (max_iter reached, "
            "operator was never found indefinite along the way)");
      }
    }

    for (std::size_t n = 0; n < ids.size(); ++n) {
      PolarSegment& segment = registry.Get(ids[n], EwaldChargeState::Neutral);
      Index base = offsets[n];
      for (Index s = 0; s < segment.size(); ++s) {
        segment[s].setInduced_Dipole(x.segment<3>(base + 3 * s));
      }
    }
  } else {
    // As legacy polarmethod.induce = 0: induced dipoles stay zero (set
    // before the permanent-field computation).
    XTP_LOG(Log::info, log)
        << TimeStamp()
        << " polarmethod.induce is false: skipping the PCG solve, "
           "induced dipoles left at zero"
        << std::flush;
  }

  auto t_cpt = std::chrono::steady_clock::now();
  // Store the permanent field b in V() before writing the checkpoint, so
  // that it can be read back (e.g. by hdf5_dump). Otherwise V() holds the
  // last trial field of the solve, or zero if induce_ is false.
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

  // Store the convergence parameters with the converged state: a job
  // embedding a foreground in this background must use the same alpha,
  // k_max and shape, since alpha determines how the interaction is split
  // between the real- and reciprocal-space sums. See EwaldParameters.
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

  XTP_LOG(Log::info, log) << TimeStamp() << " Checkpoint written to "
                          << checkpoint_file_ << " (" << elapsed_s(t_cpt)
                          << "s)" << std::flush;
  XTP_LOG(Log::info, log) << TimeStamp()
                          << " Ewald background calculation done, total "
                          << elapsed_s(t_start) << "s" << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDBACKGROUND_H
