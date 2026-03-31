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
#include <fstream>
#include <iostream>

// Local VOTCA includes
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/anderson_mixing.h"
#include "votca/xtp/gw.h"
#include "votca/xtp/newton_rapson.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigmafactory.h"

namespace votca {
namespace xtp {

void GW::configure(const options& opt) {
  opt_ = opt;
  qptotal_ = opt_.qpmax - opt_.qpmin + 1;
  rpa_.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  sigma_ = SigmaFactory().Create(opt_.sigma_integration, Mmn_, rpa_);
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

}

double GW::CalcHomoLumoShift(Eigen::VectorXd frequencies) const {
  double DFTgap = dft_energies_(opt_.homo + 1) - dft_energies_(opt_.homo);
  double QPgap = frequencies(opt_.homo + 1 - opt_.qpmin) -
                 frequencies(opt_.homo - opt_.qpmin);
  return QPgap - DFTgap;
}

Eigen::MatrixXd GW::getHQP() const {
  return Sigma_x_ + Sigma_c_ - vxc_ +
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

    XTP_LOG(Log::error, log_)
        << level
        << (boost::format(" = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = "
                          "%4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") %
            (i + opt_.qpmin) % dft_energies_(i + opt_.qpmin) % vxc_(i, i) %
            Sigma_x_(i, i) % Sigma_c_(i, i) % gwa_energies(i))
               .str()
        << std::flush;
  }
  return;
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

  Sigma_x_ = (1 - opt_.ScaHFX) * sigma_->CalcExchangeMatrix();
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Calculated Hartree exchange contribution"
      << std::flush;
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

  for (Index i_gw = 0; i_gw < opt_.gw_sc_max_iterations; ++i_gw) {
    gw_sc_iteration_ = i_gw;
    if (i_gw % opt_.reset_3c == 0 && i_gw != 0) {
      Mmn_.Rebuild();
      XTP_LOG(Log::info, log_)
          << TimeStamp() << " Rebuilding 3c integrals" << std::flush;
    }
    sigma_->PrepareScreening();
    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Calculated screening via RPA" << std::flush;
    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Solving QP equations " << std::flush;
    if (opt_.gw_mixing_order > 0 && i_gw > 0) {
      mixing_.UpdateInput(frequencies);
    }

    frequencies = SolveQP(frequencies);

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
  Sigma_c_.diagonal() = sigma_->CalcCorrelationDiag(frequencies);
  PrintGWA_Energies();
}

Eigen::VectorXd GW::getGWAResults() const {
  return Sigma_x_.diagonal() + Sigma_c_.diagonal() - vxc_.diagonal() +
         dft_energies_.segment(opt_.qpmin, qptotal_);
}

Eigen::VectorXd GW::SolveQP(const Eigen::VectorXd& frequencies) const {
  sigma_->ResetDiagEvalCounter();
  Eigen::VectorXd env = Eigen::VectorXd::Zero(qptotal_);

  const Eigen::VectorXd intercepts =
      dft_energies_.segment(opt_.qpmin, qptotal_) + Sigma_x_.diagonal() -
      vxc_.diagonal();

  Eigen::VectorXd frequencies_new = frequencies;
  Eigen::Array<bool, Eigen::Dynamic, 1> converged =
      Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(qptotal_);

  QPStats total_stats;

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
    {
      total_stats.Add(local_stats);
    }
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

  XTP_LOG(Log::info, log_) << TimeStamp() << " QP diagnostics: "
                           << "scan=" << total_stats.sigma_scan_calls
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
  solver_opt.g_sc_limit = opt_.g_sc_limit;
  solver_opt.qp_bisection_max_iter = opt_.g_sc_max_iterations;
  solver_opt.qp_grid_spacing = opt_.qp_grid_spacing;
  solver_opt.qp_grid_steps = opt_.qp_grid_steps;

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
  solver_opt.g_sc_limit = opt_.g_sc_limit;
  solver_opt.qp_bisection_max_iter = opt_.g_sc_max_iterations;
  solver_opt.qp_grid_spacing = opt_.qp_grid_spacing;
  solver_opt.qp_grid_steps = opt_.qp_grid_steps;

  const bool use_brent = (opt_.qp_root_finder == "brent");

  std::vector<QPRootCandidate> accepted_roots;
  std::vector<QPRootCandidate> rejected_roots;
  std::vector<std::pair<double, double>> roots;  // omega : Z

  if (left_limit < right_limit) {
    double freq_prev = left_limit;
    double targ_prev = fqp.value(freq_prev, EvalStage::Scan);

    const Index n_steps =
        std::max<Index>(2, static_cast<Index>(
                               std::ceil((right_limit - left_limit) /
                                         opt_.qp_grid_spacing)) +
                               1);

    for (Index i_node = 1; i_node < n_steps; ++i_node) {
      const double freq =
          (i_node == n_steps - 1)
              ? right_limit
              : std::min(right_limit,
                         left_limit + static_cast<double>(i_node) *
                                          opt_.qp_grid_spacing);

      const double targ = fqp.value(freq, EvalStage::Scan);

      if (targ_prev * targ < 0.0) {
        auto cand = qp_solver::RefineQPInterval(freq_prev, targ_prev, freq, targ,
                                                fqp, frequency0, solver_opt,
                                                use_brent);
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
        XTP_LOG(Log::info, log_) << " Dense scan qplevel:" << gw_level
                                 << " roots (omega:Z)\n\t\t";
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
    auto best =
        std::max_element(accepted_roots.begin(), accepted_roots.end(),
                         [](const QPRootCandidate& a,
                            const QPRootCandidate& b) {
                           return qp_solver::ScoreRoot(a) <
                                  qp_solver::ScoreRoot(b);
                         });
    return best->omega;
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

    auto least_bad =
        std::max_element(rejected_roots.begin(), rejected_roots.end(),
                         [](const QPRootCandidate& a,
                            const QPRootCandidate& b) {
                           return qp_solver::ScoreRoot(a) <
                                  qp_solver::ScoreRoot(b);
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
    auto adaptive =
        SolveQP_Grid_Windowed_Adaptive(intercept0, frequency0, gw_level,
                                       left_limit, right_limit,
                                       allow_rejected_return, &total_stats);
    if (adaptive) {
      if (stats != nullptr) {
        *stats = total_stats;
      }
      return adaptive;
    }

    QPStats dense_stats;
    auto dense = SolveQP_Grid_Windowed_Dense(intercept0, frequency0, gw_level,
                                             left_limit, right_limit,
                                             allow_rejected_return,
                                             &dense_stats);
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
  const double range =
      opt_.qp_grid_spacing * double(opt_.qp_grid_steps - 1) / 2.0;

  const double full_left_limit = frequency0 - range;
  const double full_right_limit = frequency0 + range;

  double restricted_left_limit = full_left_limit;
  double restricted_right_limit = full_right_limit;

  bool use_restricted_window = false;

  if (opt_.qp_restrict_search) {
    const Index mo_level = gw_level + opt_.qpmin;
    const bool is_occupied = (mo_level <= opt_.homo);

    if (is_occupied) {
      restricted_right_limit = std::min(full_right_limit, -opt_.qp_zero_margin);
    } else {
      restricted_left_limit =
          std::max(full_left_limit, opt_.qp_virtual_min_energy);
    }

    const double tol = 1e-12;
    use_restricted_window =
        (std::abs(restricted_left_limit - full_left_limit) > tol) ||
        (std::abs(restricted_right_limit - full_right_limit) > tol);
  }

  if (use_restricted_window && restricted_left_limit < restricted_right_limit) {
    auto restricted = SolveQP_Grid_Windowed(intercept0, frequency0, gw_level,
                                            restricted_left_limit,
                                            restricted_right_limit,
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
            << restricted_right_limit
            << "], forcing full dense window ["
            << full_left_limit << ", " << full_right_limit << "]"
            << std::flush;
      }
    }

    // Important: do NOT go back to adaptive here.
    // The restricted window already failed to produce an accepted root.
    // Force the old robust full-window dense search.
    return SolveQP_Grid_Windowed_Dense(intercept0, frequency0, gw_level,
                                       full_left_limit, full_right_limit,
                                       true, stats);
  }

  return SolveQP_Grid_Windowed(intercept0, frequency0, gw_level,
                               full_left_limit, full_right_limit,
                               true, stats);
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
  Eigen::VectorXd diag_backup = Sigma_c_.diagonal();
  Sigma_c_ = sigma_->CalcCorrelationOffDiag(getGWAResults());
  Sigma_c_.diagonal() = diag_backup;
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
      dft_energies_.segment(opt_.qpmin, qptotal_) + Sigma_x_.diagonal() -
      vxc_.diagonal();
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
