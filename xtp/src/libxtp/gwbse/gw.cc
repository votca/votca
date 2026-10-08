/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

// Standard includes
#include <chrono>
#include <fstream>
#include <iostream>

// Local VOTCA includes
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/anderson_mixing.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/gw.h"
#include "votca/xtp/newton_rapson.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigmafactory.h"

namespace votca {
namespace xtp {

void GW::configure(const options& opt) {
  opt_ = opt;

  // Normalize legacy and new grid-search settings once at configuration time.
  // From this point on, the solver operates only on the decoupled canonical
  // controls, independent of how the user specified them in XML.
  qp_solver::NormalizeGridSearchOptions(opt_);
  qptotal_ = opt_.qpmax - opt_.qpmin + 1;
  rpa_.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa_.setOutOfWindowShift(opt_.out_of_window_shift);
  sigma_ = SigmaFactory().Create(opt_.sigma_integration, Mmn_, rpa_);
  sigma_->setTimings(&timings_);
  rpa_.setTimings(&timings_);
  Sigma_base::options sigma_opt;
  sigma_opt.homo = opt_.homo;
  sigma_opt.qpmax = opt_.qpmax;
  sigma_opt.qpmin = opt_.qpmin;
  sigma_opt.rpamin = opt_.rpamin;
  sigma_opt.rpamax = opt_.rpamax;
  sigma_opt.eta = opt_.eta;
  sigma_opt.alpha = opt_.alpha;
  sigma_opt.quadrature_scheme = opt_.quadrature_scheme;
  sigma_opt.order = opt_.order;
  sigma_->configure(sigma_opt);
  Sigma_x_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Sigma_c_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Sigma_reac_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
}

void GW::setReactionField(const Eigen::MatrixXd& R) {
  reaction_field_ = R;
  dressing_ = (R.size() > 0) ? EnvironmentScreening::DressingMatrix(R)
                             : Eigen::MatrixXd();
}

void GW::DressForEnvironment() {
  if (dressing_.size() > 0 && !Mmn_.Dressed()) {
    Mmn_.DressAuxIndex(dressing_);
    XTP_LOG(Log::info, log_)
        << TimeStamp()
        << " Dressed 3c integrals with (1+R)^1/2 for environment screening"
        << std::flush;
  }
}

void GW::PrintReactionFieldShifts() const {
  const Index homo_rel = opt_.homo - opt_.qpmin;
  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << " Environment reaction field, static COH+SEX, diagonal [eV]:"
      << std::flush;
  for (Index i = 0; i < qptotal_; i++) {
    std::string level = "  Level";
    if (i == homo_rel) {
      level = "  HOMO ";
    } else if (i == homo_rel + 1) {
      level = "  LUMO ";
    }
    XTP_LOG(Log::info, log_)
        << level
        << (boost::format(" = %1$4d Sigma_reac = %2$+1.6f") % (i + opt_.qpmin) %
            (Sigma_reac_(i, i) * tools::conv::hrt2ev))
               .str()
        << std::flush;
  }
  if (homo_rel >= 0 && homo_rel + 1 < qptotal_) {
    XTP_LOG(Log::error, log_)
        << (boost::format("  HOMO %1$+1.4f eV, LUMO %2$+1.4f eV, gap %3$+1.4f "
                          "eV") %
            (Sigma_reac_(homo_rel, homo_rel) * tools::conv::hrt2ev) %
            (Sigma_reac_(homo_rel + 1, homo_rel + 1) * tools::conv::hrt2ev) %
            ((Sigma_reac_(homo_rel + 1, homo_rel + 1) -
              Sigma_reac_(homo_rel, homo_rel)) *
             tools::conv::hrt2ev))
               .str()
        << std::flush;
  }
}

double GW::CalcHomoLumoShift(Eigen::VectorXd frequencies) const {
  double DFTgap = dft_energies_(opt_.homo + 1) - dft_energies_(opt_.homo);
  double QPgap = frequencies(opt_.homo + 1 - opt_.qpmin) -
                 frequencies(opt_.homo - opt_.qpmin);
  return QPgap - DFTgap;
}

Eigen::MatrixXd GW::getHQP() const {
  return Sigma_x_ + Sigma_c_ + Sigma_reac_ - vxc_ +
         Eigen::MatrixXd(
             dft_energies_.segment(opt_.qpmin, qptotal_).asDiagonal());
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> GW::DiagonalizeQPHamiltonian()
    const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQP());
  PrintQP_Energies(qpdiag.eigenvalues());
  return qpdiag;
}

void GW::PrintGWA_Energies() const {
  Eigen::VectorXd gwa_energies = getGWAResults();
  double shift = CalcHomoLumoShift(gwa_energies);

  XTP_LOG(Log::error, log_)
      << (boost::format(
              "  ====== Perturbative quasiparticle energies (Hartree) ====== "))
             .str()
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("   DeltaHLGap = %1$+1.6f Hartree") % shift).str()
      << std::flush;

  for (Index i = 0; i < qptotal_; i++) {
    std::string level = "  Level";
    if ((i + opt_.qpmin) == opt_.homo) {
      level = "  HOMO ";
    } else if ((i + opt_.qpmin) == opt_.homo + 1) {
      level = "  LUMO ";
    }

    // With an environment, its static part of the correlation gets a
    // column of its own, S-R; S-C stays the dynamic correlation.
    const std::string env =
        (reaction_field_.size() > 0)
            ? (boost::format(" S-R = %1$+1.4f") % Sigma_reac_(i, i)).str()
            : std::string();
    XTP_LOG(Log::error, log_)
        << level
        << (boost::format(" = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = "
                          "%4$+1.4f S-C = %5$+1.4f") %
            (i + opt_.qpmin) % dft_energies_(i + opt_.qpmin) % vxc_(i, i) %
            Sigma_x_(i, i) % Sigma_c_(i, i))
               .str()
        << env << (boost::format(" GWA = %1$+1.4f") % gwa_energies(i)).str()
        << std::flush;
  }
  return;
}

void GW::PrintQSGW_Composition(double threshold) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << "  ====== QSGW QP state composition (dominant KS contributions) ======"
      << std::flush;

  for (Index n = 0; n < qptotal_; n++) {
    const Index abs_n = n + opt_.qpmin;
    std::string level = "  Level";
    if (abs_n == opt_.homo)
      level = "  HOMO ";
    else if (abs_n == opt_.homo + 1)
      level = "  LUMO ";

    // Collect contributions above threshold, sorted by weight descending
    Eigen::VectorXd weights = qsgw_rotation_.col(n).cwiseAbs2();
    std::vector<std::pair<double, Index>> contribs;
    for (Index m = 0; m < qptotal_; m++) {
      if (weights(m) >= threshold) {
        contribs.push_back({weights(m), m + opt_.qpmin});
      }
    }
    std::sort(contribs.begin(), contribs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    std::string line;
    for (const auto& [w, ks] : contribs) {
      std::string ks_label;
      if (ks == opt_.homo)
        ks_label = "(HOMO)";
      else if (ks == opt_.homo + 1)
        ks_label = "(LUMO)";

      line += (boost::format("  %1$5.1f%% KS_%2$d%3$s") % (100.0 * w) % ks %
               ks_label)
                  .str();
    }
    XTP_LOG(Log::error, log_)
        << level << (boost::format(" = %1$4d:") % abs_n).str() << line
        << std::flush;
  }
}

void GW::PrintQSGW_Energies(const std::string& seed_label,
                            const Eigen::VectorXd& seed_energies,
                            const Eigen::VectorXd& qsgw_energies) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " QSGW quasiparticle energies" << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format(
              "  ====== QSGW quasiparticle energies (Hartree) ====== "))
             .str()
      << std::flush;

  for (Index i = 0; i < qptotal_; i++) {
    std::string level = "  Level";
    if ((i + opt_.qpmin) == opt_.homo) {
      level = "  HOMO ";
    } else if ((i + opt_.qpmin) == opt_.homo + 1) {
      level = "  LUMO ";
    }
    XTP_LOG(Log::error, log_)
        << level
        << (boost::format(" = %1$4d %2$s = %3$+1.6f  QSGW = %4$+1.6f") %
            (i + opt_.qpmin) % seed_label % seed_energies(i) % qsgw_energies(i))
               .str()
        << std::flush;
  }
}

void GW::PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const {
  Eigen::VectorXd gwa_energies = getGWAResults();
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Full quasiparticle Hamiltonian  " << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format(
              "  ====== Diagonalized quasiparticle energies (Hartree) "
              "====== "))
             .str()
      << std::flush;
  for (Index i = 0; i < qptotal_; i++) {
    std::string level = "  Level";
    if ((i + opt_.qpmin) == opt_.homo) {
      level = "  HOMO ";
    } else if ((i + opt_.qpmin) == opt_.homo + 1) {
      level = "  LUMO ";
    }
    XTP_LOG(Log::error, log_)
        << level
        << (boost::format(" = %1$4d PQP = %2$+1.6f DQP = %3$+1.6f ") %
            (i + opt_.qpmin) % gwa_energies(i) % qp_diag_energies(i))
               .str()
        << std::flush;
  }
  return;
}

Eigen::VectorXd GW::ScissorShift_DFTlevel(
    const Eigen::VectorXd& dft_energies) const {
  Eigen::VectorXd shifted_energies = dft_energies;
  shifted_energies.segment(opt_.homo + 1, dft_energies.size() - opt_.homo - 1)
      .array() += opt_.shift;
  return shifted_energies;
}

void GW::CalculateGWPerturbation() {
  timings_.Reset();
  qsgw_final_energies_.resize(0);
  {
    auto t = timings_.Measure("exchange Sigma_x");
    Sigma_x_ = (1 - opt_.ScaHFX) * sigma_->CalcExchangeMatrix();
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Calculated Hartree exchange contribution"
      << std::flush;
  // The environment's static correlation, W - v restricted to v_reac:
  // part of Sigma_c by its physics, held apart because it is frequency
  // independent. Not scaled by ScaHFX: the DFT starting point has no
  // reaction field to double count.
  if (reaction_field_.size() > 0) {
    auto t = timings_.Measure("reaction field Sigma_reac");
    Sigma_reac_ = sigma_->CalcReactionFieldMatrix(reaction_field_);
    PrintReactionFieldShifts();
  }
  // From here on the RPA and Sigma_c see u = v + v_reac.
  {
    auto t = timings_.Measure("dress Mmn for the environment");
    DressForEnvironment();
  }
  // dftenergies has size aobasissize
  // rpaenergies/Mmn have size rpatotal
  // gwaenergies/frequencies have size qptotal
  // homo index is relative to dft_energies
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Scissor shifting DFT energies by: " << opt_.shift
      << " Hrt" << std::flush;
  Eigen::VectorXd dft_shifted_energies = ScissorShift_DFTlevel(dft_energies_);
  rpa_.setRPAInputEnergies(
      dft_shifted_energies.segment(opt_.rpamin, opt_.rpamax - opt_.rpamin + 1));
  Eigen::VectorXd frequencies =
      dft_shifted_energies.segment(opt_.qpmin, qptotal_);

  Anderson mixing_;
  mixing_.Configure(opt_.gw_mixing_order, opt_.gw_mixing_alpha);
  // The QP energies of the last QP solve and the RPA energies it saw (G and
  // W); after the solve, the RPA energies and, with mixing, frequencies are
  // already the input of the next iteration.
  Eigen::VectorXd qp_solution = frequencies;
  Eigen::VectorXd rpa_energies_solved = rpa_.getRPAInputEnergies();

  for (Index i_gw = 0; i_gw < opt_.gw_sc_max_iterations; ++i_gw) {
    gw_sc_iteration_ = i_gw;
    if (opt_.reset_3c > 0 && i_gw % opt_.reset_3c == 0 && i_gw != 0) {
      auto t = timings_.Measure("rebuild Mmn");
      Mmn_.Rebuild();
      XTP_LOG(Log::info, log_)
          << TimeStamp() << " Rebuilding 3c integrals" << std::flush;
      DressForEnvironment();  // Rebuild leaves them bare
    }
    {
      auto t = timings_.Measure("screening: other");
      sigma_->PrepareScreening();
    }
    if (const std::string summary = sigma_->ScreeningSummary();
        !summary.empty()) {
      const Log::Level level =
          (i_gw == 0) ? Log::Level::error : Log::Level::info;
      XTP_LOG(level, log_) << TimeStamp() << " " << summary << std::flush;
    }
    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Calculated screening via RPA" << std::flush;
    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Solving QP equations " << std::flush;
    if (opt_.gw_mixing_order > 0 && i_gw > 0) {
      mixing_.UpdateInput(frequencies);
    }

    {
      auto t = timings_.Measure("QP equations (Sigma_c diagonal)");
      frequencies = SolveQP(frequencies);
    }
    qp_solution = frequencies;
    rpa_energies_solved = rpa_.getRPAInputEnergies();

    if (opt_.gw_sc_max_iterations > 1) {
      Eigen::VectorXd rpa_energies_old = rpa_.getRPAInputEnergies();
      if (opt_.gw_mixing_order > 0 && i_gw > 0) {
        if (opt_.gw_mixing_order == 1) {
          XTP_LOG(Log::debug, log_)
              << "GWSC using linear mixing with alpha: " << opt_.gw_mixing_alpha
              << std::flush;
        } else {
          XTP_LOG(Log::debug, log_)
              << "GWSC using Anderson mixing with history "
              << opt_.gw_mixing_order << ", alpha: " << opt_.gw_mixing_alpha
              << std::flush;
        }
        mixing_.UpdateOutput(frequencies);
        Eigen::VectorXd mixed_frequencies = mixing_.MixHistory();
        rpa_.UpdateRPAInputEnergies(dft_energies_, mixed_frequencies,
                                    opt_.qpmin);
        frequencies = mixed_frequencies;

      } else {
        XTP_LOG(Log::debug, log_) << "GWSC using plain update " << std::flush;
        rpa_.UpdateRPAInputEnergies(dft_energies_, frequencies, opt_.qpmin);
      }

      for (int i = 0; i < frequencies.size(); i++) {
        XTP_LOG(Log::debug, log_)
            << "... GWSC iter " << i_gw << " state " << i << " "
            << std::setprecision(9) << frequencies(i) << std::flush;
      }

      XTP_LOG(Log::info, log_)
          << TimeStamp() << " GW_Iteration:" << i_gw
          << " Shift[Hrt]:" << CalcHomoLumoShift(frequencies) << std::flush;
      if (Converged(rpa_.getRPAInputEnergies(), rpa_energies_old,
                    opt_.gw_sc_limit)) {
        XTP_LOG(Log::info, log_) << TimeStamp() << " Converged after "
                                 << i_gw + 1 << " GW iterations." << std::flush;
        break;
      } else if (i_gw == opt_.gw_sc_max_iterations - 1) {
        XTP_LOG(Log::error, log_)
            << TimeStamp()
            << " WARNING! GW-self-consistency cycle not converged after "
            << opt_.gw_sc_max_iterations << " iterations." << std::flush;
        XTP_LOG(Log::error, log_)
            << TimeStamp() << "      Run continues. Inspect results carefully!"
            << std::flush;
        break;
      }
    }
  }
  {
    // Sigma_c of the last QP solve, at its solution and with the energies in
    // G it saw, so that the PQP energies solve their QP equation. The RPA
    // energies updated after the solve (input of a next iteration) stay.
    PrintCompetingRoots(Log::error);
    auto t = timings_.Measure("QP equations (Sigma_c diagonal)");
    const Eigen::VectorXd rpa_energies_next = rpa_.getRPAInputEnergies();
    rpa_.setRPAInputEnergies(rpa_energies_solved);
    Sigma_c_.diagonal() = sigma_->CalcCorrelationDiag(qp_solution);
    rpa_.setRPAInputEnergies(rpa_energies_next);
  }
  PrintGWA_Energies();
  timings_.Report(log_, Log::error, "GW timing summary");
}

Eigen::VectorXd GW::getGWAResults() const {
  // When QSGW with virtual trimming has run, return the pre-computed merged
  // energy vector (QSGW for the trimmed window, seed for excluded virtuals).
  // This avoids recomputing from sigma matrices which were zeroed post-loop.
  if (qsgw_final_energies_.size() > 0) {
    return qsgw_final_energies_;
  }
  return Sigma_x_.diagonal() + Sigma_c_.diagonal() + Sigma_reac_.diagonal() -
         vxc_.diagonal() + dft_energies_.segment(opt_.qpmin, qptotal_);
}

void GW::PrintCompetingRoots(Log::Level level) const {
  Index n = 0;
  for (const auto& c : competing_roots_) {
    n += c.competing ? 1 : 0;
  }
  if (n == 0) {
    return;
  }
  XTP_LOG(level, log_) << TimeStamp() << " " << n
                       << " QP level(s) with two roots of comparable weight "
                          "(second Z at least half the largest):"
                       << std::flush;
  for (Index i = 0; i < Index(competing_roots_.size()); ++i) {
    const auto& c = competing_roots_[std::size_t(i)];
    if (c.competing) {
      XTP_LOG(level, log_)
          << (boost::format("   Level = %1$4d chosen %2$+1.6f (Z = %3$.3f), "
                            "other %4$+1.6f (Z = %5$.3f)") %
              (i + opt_.qpmin) % c.omega % c.Z % c.omega_alt % c.Z_alt)
                 .str()
          << std::flush;
    }
  }
}

Eigen::VectorXd GW::SolveQP(const Eigen::VectorXd& frequencies) const {
  sigma_->ResetDiagEvalCounter();
  competing_roots_.assign(std::size_t(qptotal_), qp_solver::CompetingRoots{});
  Eigen::VectorXd env = Eigen::VectorXd::Zero(qptotal_);

  const Eigen::VectorXd intercepts =
      dft_energies_.segment(opt_.qpmin, qptotal_) + Sigma_x_.diagonal() +
      Sigma_reac_.diagonal() - vxc_.diagonal();

  Eigen::VectorXd frequencies_new = frequencies;
  Eigen::Array<bool, Eigen::Dynamic, 1> converged =
      Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(qptotal_);

  QPStats total_stats;

  // HOMO-LUMO midpoint of the current QP guesses (DFT energies where a
  // level is outside the QP range), for SolveQP_Grid's restricted window.
  {
    const Index h = opt_.homo - opt_.qpmin;
    const double e_homo =
        (h >= 0 && h < qptotal_) ? frequencies[h] : dft_energies_(opt_.homo);
    const double e_lumo = (h + 1 >= 0 && h + 1 < qptotal_)
                              ? frequencies[h + 1]
                              : dft_energies_(opt_.homo + 1);
    qp_midgap_ = 0.5 * (e_homo + e_lumo);
  }

#ifdef _OPENMP
  Index use_threads =
      OPENMP::getMaxThreads() > qptotal_ ? qptotal_ : OPENMP::getMaxThreads();
#else
  Index use_threads = 1;
#endif

#pragma omp parallel for schedule(dynamic) num_threads(use_threads)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {

    double initial_f = frequencies[gw_level];
    double intercept = intercepts[gw_level];
    boost::optional<double> newf;
    QPStats local_stats;

    if (opt_.qp_solver == "fixedpoint") {
      newf = SolveQP_FixedPoint(intercept, initial_f, gw_level, &local_stats);
    }
    if (newf) {
      frequencies_new[gw_level] = newf.value();
      converged[gw_level] = true;
    } else {
      newf = SolveQP_Grid(intercept, initial_f, gw_level, &local_stats);
      if (newf) {
        frequencies_new[gw_level] = newf.value();
        converged[gw_level] = true;
      } else {
        newf =
            SolveQP_Linearisation(intercept, initial_f, gw_level, &local_stats);
        if (newf) {
          frequencies_new[gw_level] = newf.value();
        }
      }
    }

#pragma omp critical
    { total_stats.Add(local_stats); }
  }

  if (!converged.all()) {
    std::vector<Index> states;
    for (Index s = 0; s < converged.size(); s++) {
      if (!converged[s]) {
        states.push_back(s);
      }
    }
    IndexParser rp;
    XTP_LOG(Log::error, log_) << TimeStamp() << " Not converged PQP states are:"
                              << rp.CreateIndexString(states) << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Increase the grid search interval" << std::flush;
  }

  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Sigma diagonal evaluations in SolveQP: "
                           << sigma_->GetDiagEvalCounter() << std::flush;

  XTP_LOG(Log::info, log_) << TimeStamp() << " QP diagnostics: " << "scan="
                           << total_stats.sigma_scan_calls
                           << " refine=" << total_stats.sigma_refine_calls
                           << " deriv_sigma="
                           << total_stats.sigma_derivative_calls
                           << " other=" << total_stats.sigma_other_calls
                           << " total_sigma=" << total_stats.TotalSigmaCalls()
                           << " unique_omega="
                           << total_stats.sigma_unique_frequencies
                           << " repeat_sigma=" << total_stats.sigma_repeat_calls
                           << " deriv_calls=" << total_stats.deriv_calls
                           << std::flush;

  return frequencies_new;
}

boost::optional<double> GW::SolveQP_Linearisation(double intercept0,
                                                  double frequency0,
                                                  Index gw_level,
                                                  QPStats* stats) const {
  boost::optional<double> newf = boost::none;

  QPFunc fqp(gw_level, *sigma_.get(), intercept0);

  double sigma = fqp.sigma(frequency0, EvalStage::Other);
  double dsigma_domega = fqp.deriv(frequency0);
  double Z = 1.0 - dsigma_domega;
  if (std::abs(Z) > 1e-9) {
    newf = frequency0 + (intercept0 - frequency0 + sigma) / Z;
  }

  if (stats != nullptr) {
    *stats = fqp.GetStats();
  }
  return newf;
}

boost::optional<double> GW::SolveQP_Grid_Windowed_Adaptive(
    double intercept0, double frequency0, Index gw_level, double left_limit,
    double right_limit, bool allow_rejected_return, QPStats* stats) const {

  QPFunc fqp(gw_level, *sigma_.get(), intercept0);

  qp_solver::SolverOptions solver_opt;

  // Pass only canonical search controls into the shared grid-search utility.
  // RKS and UKS intentionally use the same search semantics; only the sigma
  // evaluator differs between the two implementations.
  solver_opt.g_sc_limit = opt_.g_sc_limit;
  solver_opt.qp_bisection_max_iter = opt_.g_sc_max_iterations;
  solver_opt.qp_full_window_half_width = opt_.qp_full_window_half_width;
  solver_opt.qp_dense_spacing = opt_.qp_dense_spacing;
  solver_opt.qp_adaptive_shell_width = opt_.qp_adaptive_shell_width;
  solver_opt.qp_adaptive_shell_count = opt_.qp_adaptive_shell_count;
  solver_opt.prefer_nearest_root =
      opt_.qp_root_continuity && gw_sc_iteration_ > 0;

  QPWindowDiagnostics wdiag;
  std::vector<QPRootCandidate> accepted_roots;
  std::vector<QPRootCandidate> rejected_roots;

  bool use_brent = false;
  if (opt_.qp_root_finder == "brent") {
    use_brent = true;
  }

  auto result = qp_solver::SolveQP_Grid_Windowed(
      fqp, frequency0, left_limit, right_limit, gw_sc_iteration_, solver_opt,
      &wdiag, &accepted_roots, &rejected_roots, use_brent);

  if (Log::current_level > Log::error) {
#pragma omp critical
    {
      if (!accepted_roots.empty() || !rejected_roots.empty()) {
        XTP_LOG(Log::info, log_)
            << " Adaptive scan qplevel:" << gw_level << " center="
            << std::max(left_limit, std::min(right_limit, frequency0))
            << " shells=" << wdiag.shells_explored
            << " first_interval_shell=" << wdiag.first_interval_shell
            << " first_accepted_shell=" << wdiag.first_accepted_shell
            << " intervals=" << wdiag.intervals_found << std::flush;
      }
    }
  }

  if (stats != nullptr) {
    *stats = fqp.GetStats();
  }

  if (!accepted_roots.empty()) {
    if (result) {
      RecordCompetingRoots(gw_level, accepted_roots, result.value());
    }
    return result;
  }

  if (!rejected_roots.empty() && !allow_rejected_return) {
    if (Log::current_level > Log::error) {
#pragma omp critical
      {
        XTP_LOG(Log::info, log_)
            << " Adaptive scan qplevel:" << gw_level
            << " produced only rejected roots in [" << left_limit << ", "
            << right_limit << "], forcing wider retry" << std::flush;
      }
    }
    return boost::none;
  }

  return result;
}

boost::optional<double> GW::SolveQP_Grid_Windowed_Dense(
    double intercept0, double frequency0, Index gw_level, double left_limit,
    double right_limit, bool allow_rejected_return, QPStats* stats) const {

  QPFunc fqp(gw_level, *sigma_.get(), intercept0);

  qp_solver::SolverOptions solver_opt;

  // The dense path performs a uniform sign-change scan over the requested
  // interval. Its spacing is qp_dense_spacing, not the adaptive shell width.
  solver_opt.g_sc_limit = opt_.g_sc_limit;
  solver_opt.qp_bisection_max_iter = opt_.g_sc_max_iterations;
  solver_opt.qp_full_window_half_width = opt_.qp_full_window_half_width;
  solver_opt.qp_dense_spacing = opt_.qp_dense_spacing;
  solver_opt.qp_adaptive_shell_width = opt_.qp_adaptive_shell_width;
  solver_opt.qp_adaptive_shell_count = opt_.qp_adaptive_shell_count;
  solver_opt.prefer_nearest_root =
      opt_.qp_root_continuity && gw_sc_iteration_ > 0;

  const bool use_brent = (opt_.qp_root_finder == "brent");

  std::vector<QPRootCandidate> accepted_roots;
  std::vector<QPRootCandidate> rejected_roots;
  std::vector<std::pair<double, double>> roots;  // omega : Z

  if (left_limit < right_limit) {
    double freq_prev = left_limit;
    double targ_prev = fqp.value(freq_prev, EvalStage::Scan);

    // Dense scanning is the robustness path: it uniformly samples the entire
    // interval to avoid missing sign changes that a coarser adaptive shell
    // search might step over.
    const Index n_steps = std::max<Index>(
        2, static_cast<Index>(
               std::ceil((right_limit - left_limit) / opt_.qp_dense_spacing)) +
               1);

    auto node = [&](Index i_node) {
      return (i_node == n_steps - 1)
                 ? right_limit
                 : std::min(right_limit,
                            left_limit + static_cast<double>(i_node) *
                                             opt_.qp_dense_spacing);
    };
    {
      // all scan nodes in one batch
      std::vector<double> nodes;
      for (Index i_node = 1; i_node < n_steps; ++i_node) {
        nodes.push_back(node(i_node));
      }
      qp_solver::PrefetchIfAvailable(fqp, nodes, 0);
    }

    for (Index i_node = 1; i_node < n_steps; ++i_node) {
      const double freq = node(i_node);

      const double targ = fqp.value(freq, EvalStage::Scan);

      if (targ_prev * targ < 0.0) {
        auto cand =
            qp_solver::RefineQPInterval(freq_prev, targ_prev, freq, targ, fqp,
                                        frequency0, solver_opt, use_brent);
        if (cand) {
          roots.emplace_back(cand->omega, cand->Z);
          if (cand->accepted) {
            accepted_roots.push_back(*cand);
          } else {
            rejected_roots.push_back(*cand);
          }
        }
      }

      freq_prev = freq;
      targ_prev = targ;
    }
  }

  if (Log::current_level > Log::error) {
#pragma omp critical
    {
      if (accepted_roots.empty() && rejected_roots.empty()) {
        XTP_LOG(Log::info, log_)
            << " Dense scan qplevel:" << gw_level << " found no roots in ["
            << left_limit << ", " << right_limit << "]" << std::flush;
      } else {
        XTP_LOG(Log::info, log_)
            << " Dense scan qplevel:" << gw_level << " roots (omega:Z)\n\t\t";
        for (const auto& root : roots) {
          XTP_LOG(Log::info, log_) << std::setprecision(5) << root.first << ":"
                                   << root.second << " ";
        }
        XTP_LOG(Log::info, log_) << std::flush;
      }
    }
  }

  if (stats != nullptr) {
    *stats = fqp.GetStats();
  }

  if (!accepted_roots.empty()) {
    const double best = qp_solver::SelectRoot(accepted_roots, solver_opt).omega;
    RecordCompetingRoots(gw_level, accepted_roots, best);
    return best;
  }

  if (!rejected_roots.empty()) {
    if (!allow_rejected_return) {
      if (Log::current_level > Log::error) {
#pragma omp critical
        {
          XTP_LOG(Log::info, log_)
              << " Dense scan qplevel:" << gw_level
              << " produced only rejected roots in [" << left_limit << ", "
              << right_limit << "], forcing wider retry" << std::flush;
        }
      }
      return boost::none;
    }

    auto least_bad = std::max_element(
        rejected_roots.begin(), rejected_roots.end(),
        [](const QPRootCandidate& a, const QPRootCandidate& b) {
          return qp_solver::ScoreRoot(a) < qp_solver::ScoreRoot(b);
        });
    return least_bad->omega;
  }

  return boost::none;
}

boost::optional<double> GW::SolveQP_Grid_Windowed(
    double intercept0, double frequency0, Index gw_level, double left_limit,
    double right_limit, bool allow_rejected_return, QPStats* stats) const {

  if (opt_.qp_grid_search_mode == "adaptive") {
    return SolveQP_Grid_Windowed_Adaptive(intercept0, frequency0, gw_level,
                                          left_limit, right_limit,
                                          allow_rejected_return, stats);
  }

  if (opt_.qp_grid_search_mode == "dense") {
    return SolveQP_Grid_Windowed_Dense(intercept0, frequency0, gw_level,
                                       left_limit, right_limit,
                                       allow_rejected_return, stats);
  }

  if (opt_.qp_grid_search_mode == "adaptive_with_dense_fallback") {
    QPStats total_stats;
    auto adaptive = SolveQP_Grid_Windowed_Adaptive(
        intercept0, frequency0, gw_level, left_limit, right_limit,
        allow_rejected_return, &total_stats);
    if (adaptive) {
      if (stats != nullptr) {
        *stats = total_stats;
      }
      return adaptive;
    }

    QPStats dense_stats;
    auto dense = SolveQP_Grid_Windowed_Dense(
        intercept0, frequency0, gw_level, left_limit, right_limit,
        allow_rejected_return, &dense_stats);
    total_stats.Add(dense_stats);

    if (Log::current_level > Log::error) {
#pragma omp critical
      {
        XTP_LOG(Log::info, log_)
            << " Adaptive QP scan failed for qplevel:" << gw_level
            << ", retrying dense grid scan in [" << left_limit << ", "
            << right_limit << "]" << std::flush;
      }
    }

    if (stats != nullptr) {
      *stats = total_stats;
    }
    return dense;
  }

  throw std::runtime_error("Unknown gw.qp_grid_search_mode '" +
                           opt_.qp_grid_search_mode + "'");
}

boost::optional<double> GW::SolveQP_Grid(double intercept0, double frequency0,
                                         Index gw_level, QPStats* stats) const {
  // The full QP search window is now controlled explicitly and no longer
  // inferred from the dense-grid spacing.
  const double range = opt_.qp_full_window_half_width;

  const double full_left_limit = frequency0 - range;
  const double full_right_limit = frequency0 + range;

  double restricted_left_limit = full_left_limit;
  double restricted_right_limit = full_right_limit;

  bool use_restricted_window = false;

  if (opt_.qp_restrict_search) {
    const Index mo_level = gw_level + opt_.qpmin;
    const bool is_occupied = (mo_level <= opt_.homo);

    // The fixed bounds (-qp_zero_margin for occupied levels,
    // qp_virtual_min_energy for virtual ones) assume the level itself lies
    // inside them. An environment can move it outside -- an embedded LUMO
    // below -0.1 Ha -- and then every restricted scan fails and falls back
    // to the full dense window. In that case only, the bound moves to the
    // HOMO-LUMO midpoint, which separates occupied from virtual levels
    // wherever they are.
    if (is_occupied) {
      const double bound = (frequency0 < -opt_.qp_zero_margin)
                               ? -opt_.qp_zero_margin
                               : qp_midgap_;
      restricted_right_limit = std::min(full_right_limit, bound);
    } else {
      const double bound = (frequency0 > opt_.qp_virtual_min_energy)
                               ? opt_.qp_virtual_min_energy
                               : qp_midgap_;
      restricted_left_limit = std::max(full_left_limit, bound);
    }

    const double tol = 1e-12;
    use_restricted_window =
        (std::abs(restricted_left_limit - full_left_limit) > tol) ||
        (std::abs(restricted_right_limit - full_right_limit) > tol);
  }

  if (use_restricted_window && restricted_left_limit < restricted_right_limit) {
    auto restricted =
        SolveQP_Grid_Windowed(intercept0, frequency0, gw_level,
                              restricted_left_limit, restricted_right_limit,
                              false,  // accepted roots only
                              stats);
    if (restricted) {
      return restricted;
    }

    if (Log::current_level > Log::error) {
#pragma omp critical
      {
        XTP_LOG(Log::info, log_)
            << " Restricted QP search failed for qplevel:" << gw_level
            << " in window [" << restricted_left_limit << ", "
            << restricted_right_limit << "], forcing full dense window ["
            << full_left_limit << ", " << full_right_limit << "]" << std::flush;
      }
    }

    // Important: do NOT go back to adaptive here.
    // The restricted window already failed to produce an accepted root.
    // Force the old robust full-window dense search.
    return SolveQP_Grid_Windowed_Dense(intercept0, frequency0, gw_level,
                                       full_left_limit, full_right_limit, true,
                                       stats);
  }

  return SolveQP_Grid_Windowed(intercept0, frequency0, gw_level,
                               full_left_limit, full_right_limit, true, stats);
}

boost::optional<double> GW::SolveQP_FixedPoint(double intercept0,
                                               double frequency0,
                                               Index gw_level,
                                               QPStats* stats) const {
  boost::optional<double> newf = boost::none;
  QPFunc f(gw_level, *sigma_.get(), intercept0);
  NewtonRapson<QPFunc> newton = NewtonRapson<QPFunc>(
      opt_.g_sc_max_iterations, opt_.g_sc_limit, opt_.qp_solver_alpha);
  double freq_new = newton.FindRoot(f, frequency0);
  if (newton.getInfo() == NewtonRapson<QPFunc>::success) {
    newf = freq_new;
  }
  if (stats != nullptr) {
    *stats = f.GetStats();
  }
  return newf;
}

bool GW::Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                   double epsilon) const {
  Index state = 0;
  bool energies_converged = true;
  double diff_max = (e1 - e2).cwiseAbs().maxCoeff(&state);
  if (diff_max > epsilon) {
    energies_converged = false;
  }
  XTP_LOG(Log::info, log_) << TimeStamp() << " E_diff max=" << diff_max
                           << " StateNo:" << state << std::flush;
  return energies_converged;
}

void GW::CalculateHQP() {
  if (!sigma_->HasOffDiagonal()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " WARNING: no off-diagonal Sigma_c for sigma "
        << "integration " << opt_.sigma_integration
        << "; the DQP Hamiltonian has no off-diagonal correlation"
        << std::flush;
  }
  const auto start = std::chrono::steady_clock::now();
  Eigen::VectorXd diag_backup = Sigma_c_.diagonal();
  Sigma_c_ = sigma_->CalcCorrelationOffDiag(getGWAResults());
  Sigma_c_.diagonal() = diag_backup;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Sigma_c off-diagonal ("
      << qptotal_ * (qptotal_ - 1) / 2 << " pairs) took "
      << std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
             .count()
      << " s" << std::flush;
}

// =============================================================================
// GW::CalculateQSGW
//
// Quasiparticle self-consistent GW: iterates the rotation of the
// three-centre integrals until the QP energies are converged.
//
// Index conventions (same as evGW throughout):
//   gw_level  : 0-based index within the QP window  [0, qptotal_)
//   MO level  : absolute MO index = gw_level + opt_.qpmin
//
// In the DFT-MO basis of the window
//   H_QSGW = (diag(e_DFT) - V_xc) + U tilde_Sigma U^T,
// with U the current QP orbitals (columns, phi_i = sum_k U_ki psi_k) and
//   tilde_Sigma_ij = Sigma_x,ij + 1/2 [Sigma_c,ij(E_i) + Sigma_c,ij(E_j)]
// in the QP basis, E_i the QP energy of phi_i. tilde_Sigma is evaluated
// with the integrals rotated to the QP orbitals (both orbital indices), so
// external indices, internal sums, energies and screening all refer to the
// same orbitals. CalcCorrelationOffDiag returns the symmetrised off-diagonal
// part, CalcCorrelationDiag the diagonal.
// =============================================================================
void GW::CalculateQSGW() {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Starting QSGW self-consistency loop  " << std::flush;

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " QSGW start: qptotal_=" << qptotal_
      << " vxc rows=" << vxc_.rows() << " cols=" << vxc_.cols()
      << " Sigma_x rows=" << Sigma_x_.rows() << std::flush;

  // Initialise: start from evGW/G0W0 seed energies.
  const Eigen::VectorXd e_qp_full = getGWAResults();  // full [qpmin,qpmax]
  qsgw_seed_energies_ = e_qp_full;

  // ── Virtual-level threshold: auto-trim qpmax ─────────────────────────────
  // Scan virtual levels from LUMO upward. Exclude the first virtual whose
  // perturbative correction |e_QP - e_DFT| exceeds qsgw_max_virt_correction.
  // Occupied states are never excluded — large core corrections are physical.
  Index qsgw_qpmax = opt_.qpmax;
  {
    const Index lumo_local = opt_.homo - opt_.qpmin + 1;
    bool trimmed = false;
    for (Index n = lumo_local; n < qptotal_; n++) {
      const double corr =
          std::abs(e_qp_full(n) - dft_energies_(opt_.qpmin + n));
      if (corr > opt_.qsgw_max_virt_correction) {
        qsgw_qpmax = opt_.qpmin + n - 1;
        trimmed = true;
        XTP_LOG(Log::error, log_)
            << TimeStamp() << "  QSGW virtual threshold: level "
            << (opt_.qpmin + n)
            << " correction = " << corr * votca::tools::conv::hrt2ev
            << " eV exceeds "
            << opt_.qsgw_max_virt_correction * votca::tools::conv::hrt2ev
            << " eV. Trimming QSGW window to [" << opt_.qpmin << ","
            << qsgw_qpmax << "]." << std::flush;
        break;
      }
    }
    if (!trimmed) {
      XTP_LOG(Log::error, log_)
          << TimeStamp() << "  QSGW virtual threshold: all corrections within "
          << opt_.qsgw_max_virt_correction * votca::tools::conv::hrt2ev
          << " eV. Using full window [" << opt_.qpmin << "," << qsgw_qpmax
          << "]." << std::flush;
    }
  }

  const Index qsgw_qptotal = qsgw_qpmax - opt_.qpmin + 1;
  const bool window_trimmed = (qsgw_qpmax < opt_.qpmax);

  // If trimmed, reconfigure sigma for the reduced window
  if (window_trimmed) {
    Sigma_base::options sigma_opt;
    sigma_opt.homo = opt_.homo;
    sigma_opt.qpmin = opt_.qpmin;
    sigma_opt.qpmax = qsgw_qpmax;
    sigma_opt.rpamin = opt_.rpamin;
    sigma_opt.rpamax = opt_.rpamax;
    sigma_opt.eta = opt_.eta;
    sigma_opt.quadrature_scheme = opt_.quadrature_scheme;
    sigma_opt.order = opt_.order;
    sigma_opt.alpha = opt_.alpha;
    sigma_->configure(sigma_opt);
    Sigma_x_ = Eigen::MatrixXd::Zero(qsgw_qptotal, qsgw_qptotal);
    Sigma_c_ = Eigen::MatrixXd::Zero(qsgw_qptotal, qsgw_qptotal);
  }

  Eigen::VectorXd e_qp = e_qp_full.head(qsgw_qptotal);
  qsgw_rotation_ = Eigen::MatrixXd::Identity(qsgw_qptotal, qsgw_qptotal);

  // QSGW requires a local (GGA/LDA) DFT starting point.
  // With a hybrid functional, the HF exchange embedded in the DFT eigenvalues
  // introduces a starting-point dependence that is not removed by the QSGW
  // self-consistency loop. See Faleev, van Schilfgaarde & Kotani PRL 2004
  // and Betzinger, Friedrich & Blugel PRB 2010 for discussion.
  if (opt_.ScaHFX > 0.0) {
    throw std::runtime_error(
        "GW::CalculateQSGW: QSGW is not compatible with hybrid DFT starting "
        "points (ScaHFX = " +
        std::to_string(opt_.ScaHFX) +
        "). Use a pure GGA or LDA functional as the DFT starting point.");
  }

  if (opt_.sigma_integration == "cda") {
    throw std::runtime_error(
        "GW::CalculateQSGW: QSGW is not supported with the CDA sigma "
        "integration method. The Gaussian quadrature grid used by CDA becomes "
        "ill-conditioned as the QP wavefunctions rotate away from the DFT-MO "
        "basis, causing numerical divergence. Use sigma_integration=ppm or "
        "sigma_integration=exact instead.");
  }

  // H0 = diag(e_DFT) - V_xc in the DFT-MO basis of the window (the full
  // V_xc matrix: H_KS is diagonal there, H_KS - V_xc is not)
  Eigen::MatrixXd H0 = -vxc_.topLeftCorner(qsgw_qptotal, qsgw_qptotal);
  H0.diagonal() += dft_energies_.segment(opt_.qpmin, qsgw_qptotal);

  // Store a copy of the original (unrotated) Mmn_ QP-window slices.
  // At each iteration we restore Mmn_ to the DFT-MO state and apply the full
  // rotation qsgw_rotation_ in one shot.
  // This is necessary because dU (eigenvectors of H_QSGW, always in DFT-MO
  // basis) gives the total rotation from DFT-MOs directly, not an incremental
  // rotation relative to the previous QP basis.
  // Store ALL m-slices in the full RPA window, not just the QP window.
  // PrepareScreening (PPM, CDA) calls MultiplyRightWithAuxMatrix which
  // overwrites all slices including those outside the QP window (e.g. the
  // core level when qpmin > rpamin). If we only restore QP-window slices,
  // out-of-window slices accumulate PPM-basis corruption across iterations,
  // causing the persistent oscillation seen when qpmin > rpamin.
  // Dressed for the correlation part, as in CalculateGWPerturbation --
  // the caller may have refilled the integrals since. Exchange and the
  // reaction-field self-energy are computed with kernels that undo it.
  DressForEnvironment();
  const Index mtotal = Mmn_.msize();
  std::vector<Eigen::MatrixXd> Mmn_orig(mtotal);
  for (Index m = 0; m < mtotal; m++) {
    Mmn_orig[m] = Mmn_[m];
  }
  // The auxiliary frame the snapshot is in (PrepareScreening of the seed
  // run may have rotated it), restored with it below. The bare Coulomb
  // terms are blind to it; an environment reaction field is not.
  const Eigen::MatrixXd aux_frame_orig = Mmn_.AuxFrame();

  // Fixed-point problem in tilde_Sigma (DFT-MO basis of the window):
  //   S_in -> H = H0 + S_in -> (E, U) -> S_out = F(S_in),
  // solved with Anderson/Pulay mixing of the pairs (S_in, S_out); energies
  // and orbitals always come from the same H. Converged when the residual
  // max|S_out - S_in| and the change of the QP energies are below
  // qsgw_sc_limit; the result is the eigensystem of H0 + S_out.
  // Reuses gw_mixing_order and gw_mixing_alpha.
  Anderson qsgw_mixer;
  qsgw_mixer.Configure(opt_.gw_mixing_order, opt_.gw_mixing_alpha);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << "  QSGW mixer: order=" << opt_.gw_mixing_order
      << " alpha=" << opt_.gw_mixing_alpha << std::flush;

  // Start consistently: G and W with the seed QP energies in iteration 0
  // (after G0W0 the RPA energies are still the DFT ones).
  rpa_.UpdateRPAInputEnergies(dft_energies_, e_qp, opt_.qpmin);

  const Index nflat = qsgw_qptotal * qsgw_qptotal;
  auto diagonalise = [&](const Eigen::VectorXd& S, Index iter) {
    const Eigen::MatrixXd H = H0 + Eigen::Map<const Eigen::MatrixXd>(
                                       S.data(), qsgw_qptotal, qsgw_qptotal);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info() != Eigen::ComputationInfo::Success) {
      throw std::runtime_error(
          "GW::CalculateQSGW: diagonalisation of H_QSGW failed at iter " +
          std::to_string(iter));
    }
    return es;
  };

  Eigen::VectorXd S_in;  // the tilde_Sigma the current (e_qp, U) came from
  bool converged = false;
  for (Index iter = 0; iter < opt_.qsgw_max_iterations; iter++) {

    // ── Step 1: integrals in the current QP orbitals ────────────────────────
    for (Index m = 0; m < mtotal; m++) {
      Mmn_[m] = Mmn_orig[m];
    }
    Mmn_.RestoreAuxFrame(aux_frame_orig);
    if (iter > 0) {
      Mmn_.RotateOrbitals(qsgw_rotation_, opt_.qpmin, qsgw_qpmax);
    }

    // Environment reaction field in the current QP basis: it depends on
    // the wavefunctions, not on the energies, so it is recomputed with the
    // rotation, like the exchange.
    Eigen::MatrixXd sigma_reac;
    if (reaction_field_.size() > 0) {
      sigma_reac = sigma_->CalcReactionFieldMatrix(reaction_field_);
    }

    // Recompute screening W and Sigma in current QP basis.
    sigma_->PrepareScreening();

    // Exchange: symmetric, frequency-independent.
    Sigma_x_ = sigma_->CalcExchangeMatrix();

    // ── Step 2: S_out, the symmetrised static self-energy in the QP basis
    // (e_qp(i) is the energy of QP orbital i), then in the DFT-MO basis
    Eigen::MatrixXd tilde_Sigma =
        Sigma_x_ + sigma_->CalcCorrelationOffDiag(e_qp);
    if (sigma_reac.size() > 0) {
      tilde_Sigma += sigma_reac;  // static, symmetric: no symmetrization
    }
    tilde_Sigma.diagonal() += sigma_->CalcCorrelationDiag(e_qp);
    tilde_Sigma = qsgw_rotation_ * tilde_Sigma * qsgw_rotation_.transpose();
    const Eigen::VectorXd S_out =
        Eigen::Map<const Eigen::VectorXd>(tilde_Sigma.data(), nflat);

    // ── Step 3: residual and next input ──────────────────────────────────────
    Eigen::VectorXd S_next;
    double residual = std::numeric_limits<double>::infinity();
    if (iter == 0) {
      // the seed state is not of the form H0 + S_in: start from its output
      S_next = S_out;
    } else {
      residual = (S_out - S_in).cwiseAbs().maxCoeff();
      qsgw_mixer.UpdateInput(S_in);
      qsgw_mixer.UpdateOutput(S_out);
      S_next = qsgw_mixer.MixHistory();
    }

    // ── Step 4: QP energies and orbitals of the next input ───────────────────
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
        diagonalise(S_next, iter);
    const double diff_max = (es.eigenvalues() - e_qp).cwiseAbs().maxCoeff();
    XTP_LOG(Log::error, log_)
        << TimeStamp() << "  QSGW iter " << iter
        << "  max|dE_QP| = " << diff_max * votca::tools::conv::hrt2ev
        << " eV  max|Sigma_out - Sigma_in| = "
        << residual * votca::tools::conv::hrt2ev << " eV" << std::flush;

    if (residual < opt_.qsgw_sc_limit && diff_max < opt_.qsgw_sc_limit) {
      XTP_LOG(Log::error, log_) << TimeStamp() << "  QSGW converged in "
                                << iter + 1 << " iterations." << std::flush;
      converged = true;
      // the result: the eigensystem of H0 + S_out, the output of the last
      // evaluation (it differs from the input by less than the residual)
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_out =
          diagonalise(S_out, iter);
      e_qp = es_out.eigenvalues();
      qsgw_rotation_ = es_out.eigenvectors();
      rpa_.UpdateRPAInputEnergies(dft_energies_, e_qp, opt_.qpmin);
      break;
    }

    S_in = S_next;
    e_qp = es.eigenvalues();
    qsgw_rotation_ = es.eigenvectors();
    rpa_.UpdateRPAInputEnergies(dft_energies_, e_qp, opt_.qpmin);
  }
  if (!converged) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << "  WARNING: QSGW did not converge in "
        << opt_.qsgw_max_iterations << " iterations. Inspect results carefully."
        << std::flush;
  }

  // The QSGW energies are the eigenvalues of H_QSGW, the QP orbitals its
  // eigenvectors; getGWAResults() returns the former. Sigma_x_ holds the
  // last iteration's exchange in the QP basis of that iteration.
  qsgw_final_energies_ = e_qp;

  // If the virtual window was trimmed, merge converged QSGW energies for
  // [qpmin, qsgw_qpmax] with perturbative seed energies for (qsgw_qpmax,
  // qpmax]. Pad rotation matrix with identity for excluded virtuals so that
  // getQSGWRotation() returns a qptotal_ x qptotal_ matrix for BSE hookup.
  if (window_trimmed) {
    const Index n_excluded = qptotal_ - qsgw_qptotal;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << "  QSGW: merging converged energies [" << opt_.qpmin
        << "," << qsgw_qpmax << "] with perturbative seed energies for "
        << n_excluded << " excluded virtual(s) [" << qsgw_qpmax + 1 << ","
        << opt_.qpmax << "]" << std::flush;

    // Pad rotation with identity for excluded virtual levels
    Eigen::MatrixXd U_full = Eigen::MatrixXd::Identity(qptotal_, qptotal_);
    U_full.topLeftCorner(qsgw_qptotal, qsgw_qptotal) = qsgw_rotation_;
    qsgw_rotation_ = U_full;

    // Build merged energy vector: QSGW for trimmed window, seed for rest.
    // Store in qsgw_final_energies_ so getGWAResults() returns the correct
    // full-window result without needing to recompute from sigma matrices.
    Eigen::VectorXd e_merged = e_qp_full;
    e_merged.head(qsgw_qptotal) = e_qp;
    qsgw_final_energies_ = e_merged;

    // Update RPA with merged energies so RPAInputEnergies() is consistent.
    rpa_.UpdateRPAInputEnergies(dft_energies_, e_merged, opt_.qpmin);

    // Restore Sigma_x_ and Sigma_c_ to full qptotal_ x qptotal_ size so that
    // any downstream code doesn't encounter a size mismatch with
    // vxc_.diagonal() (which has qptotal_ elements). getGWAResults() bypasses
    // these matrices via qsgw_final_energies_, so zeroing them is safe.
    Sigma_x_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
    Sigma_c_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);

    // Restore the full-size sigma configuration.
    Sigma_base::options sigma_opt_full;
    sigma_opt_full.homo = opt_.homo;
    sigma_opt_full.qpmin = opt_.qpmin;
    sigma_opt_full.qpmax = opt_.qpmax;
    sigma_opt_full.rpamin = opt_.rpamin;
    sigma_opt_full.rpamax = opt_.rpamax;
    sigma_opt_full.eta = opt_.eta;
    sigma_opt_full.quadrature_scheme = opt_.quadrature_scheme;
    sigma_opt_full.order = opt_.order;
    sigma_opt_full.alpha = opt_.alpha;
    sigma_->configure(sigma_opt_full);
  }

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " QSGW loop complete." << std::flush;
}

void GW::PlotSigma(std::string filename, Index steps, double spacing,
                   std::string states) const {

  Eigen::VectorXd frequencies =
      rpa_.getRPAInputEnergies().segment(opt_.qpmin - opt_.rpamin, qptotal_);

  std::vector<Index> state_inds;
  IndexParser rp;
  std::vector<Index> parsed_states = rp.CreateIndexVector(states);
  for (Index gw_level : parsed_states) {
    if (gw_level >= opt_.qpmin && gw_level <= opt_.qpmax) {
      state_inds.push_back(gw_level);
    }
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " PQP(omega) written to '" << filename
      << "' for states " << rp.CreateIndexString(state_inds) << std::flush;

  const Index num_states = state_inds.size();

  const Eigen::VectorXd intercept =
      dft_energies_.segment(opt_.qpmin, qptotal_) + Sigma_x_.diagonal() +
      Sigma_reac_.diagonal() - vxc_.diagonal();
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(steps, 2 * num_states);
#pragma omp parallel for schedule(dynamic)
  for (Index grid_point = 0; grid_point < steps; grid_point++) {
    const double offset =
        ((double)grid_point - ((double)(steps - 1) / 2.0)) * spacing;
    for (Index i = 0; i < num_states; i++) {
      const Index gw_level = state_inds[i];
      const double omega = frequencies(gw_level) + offset;
      double sigma = sigma_->CalcCorrelationDiagElement(gw_level, omega);
      mat(grid_point, 2 * i) = omega;
      mat(grid_point, 2 * i + 1) = sigma + intercept[gw_level];
    }
  }

  std::ofstream out;
  out.open(filename);
  for (Index i = 0; i < num_states; i++) {
    const Index gw_level = state_inds[i];
    out << boost::format("#%1$somega_%2$d\tE_QP(omega)_%2$d") %
               (i == 0 ? "" : "\t") % gw_level;
  }
  out << std::endl;
  boost::format numFormat("%+1.6f");
  Eigen::IOFormat matFormat(Eigen::StreamPrecision, 0, "\t", "\n");
  out << numFormat % mat.format(matFormat) << std::endl;
  out.close();
}

}  // namespace xtp
}  // namespace votca
