/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 */

#include <fstream>
#include <iostream>
#include <limits>

#include "votca/xtp/IndexParser.h"
#include "votca/xtp/anderson_mixing.h"
#include "votca/xtp/gw_uks.h"
#include "votca/xtp/newton_rapson.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/sigmafactory_uks.h"

#include "self_energy_evaluators/sigma_ppm_uks.h"

namespace votca {
namespace xtp {

GW_UKS::GW_UKS(Logger& log, TCMatrix_gwbse_spin& Mmn,
               const Eigen::MatrixXd& vxc_alpha,
               const Eigen::MatrixXd& vxc_beta,
               const Eigen::VectorXd& dft_energies_alpha,
               const Eigen::VectorXd& dft_energies_beta)
    : log_(log),
      Mmn_(Mmn),
      vxc_alpha_(vxc_alpha),
      vxc_beta_(vxc_beta),
      dft_energies_alpha_(dft_energies_alpha),
      dft_energies_beta_(dft_energies_beta),
      rpa_(log, Mmn) {}

void GW_UKS::configure(const options& opt) {
  opt_ = opt;
  qptotal_ = opt_.qpmax - opt_.qpmin + 1;
  rpa_.configure(opt_.homo_alpha, opt_.homo_beta, opt_.rpamin, opt_.rpamax);

  sigma_alpha_ = SigmaFactory_UKS().Create(opt_.sigma_integration, Mmn_, rpa_,
                                           TCMatrix::SpinChannel::Alpha);
  sigma_beta_ = SigmaFactory_UKS().Create(opt_.sigma_integration, Mmn_, rpa_,
                                          TCMatrix::SpinChannel::Beta);

  if (!sigma_alpha_ || !sigma_beta_) {
    throw std::runtime_error(
        "GW_UKS currently supports only implemented unrestricted sigma "
        "evaluators. At present this is 'ppm', 'cda', and 'exact'.");
  }

  Sigma_base_UKS::options sigma_opt;
  sigma_opt.homo = opt_.homo_alpha;
  sigma_opt.qpmax = opt_.qpmax;
  sigma_opt.qpmin = opt_.qpmin;
  sigma_opt.rpamin = opt_.rpamin;
  sigma_opt.rpamax = opt_.rpamax;
  sigma_opt.eta = opt_.eta;
  sigma_opt.alpha = opt_.alpha;
  sigma_opt.quadrature_scheme = opt_.quadrature_scheme;
  sigma_opt.order = opt_.order;
  sigma_alpha_->configure(sigma_opt);
  sigma_opt.homo = opt_.homo_beta;
  sigma_beta_->configure(sigma_opt);

  Sigma_x_alpha_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Sigma_x_beta_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Sigma_c_alpha_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  Sigma_c_beta_ = Eigen::MatrixXd::Zero(qptotal_, qptotal_);

  auto print_spin_ranges = [&](Spin spin, Index homo) {
    const Index lumo = homo + 1;

    const Index n_occ_rpa =
        (opt_.rpamin <= homo) ? (homo - opt_.rpamin + 1) : 0;
    const Index n_virt_rpa =
        (lumo <= opt_.rpamax) ? (opt_.rpamax - lumo + 1) : 0;

    const Index qp_occ_min = opt_.qpmin;
    const Index qp_occ_max = std::min(opt_.qpmax, homo);
    const Index n_occ_qp =
        (qp_occ_min <= qp_occ_max) ? (qp_occ_max - qp_occ_min + 1) : 0;

    const Index qp_virt_min = std::max(opt_.qpmin, lumo);
    const Index qp_virt_max = opt_.qpmax;
    const Index n_virt_qp =
        (qp_virt_min <= qp_virt_max) ? (qp_virt_max - qp_virt_min + 1) : 0;

    XTP_LOG(Log::error, log_)
        << TimeStamp() << " UKS " << SpinName(spin)
        << " effective ranges: HOMO=" << homo << " LUMO=" << lumo
        << " | RPA occ[" << opt_.rpamin << ":" << homo << "]"
        << " virt[" << lumo << ":" << opt_.rpamax << "]"
        << " (#occ=" << n_occ_rpa << ", #virt=" << n_virt_rpa << ")"
        << " | GW occ[" << qp_occ_min << ":" << qp_occ_max << "]"
        << " virt[" << qp_virt_min << ":" << qp_virt_max << "]"
        << " (#occ=" << n_occ_qp << ", #virt=" << n_virt_qp << ")"
        << std::flush;
  };

  print_spin_ranges(Spin::Alpha, opt_.homo_alpha);
  print_spin_ranges(Spin::Beta, opt_.homo_beta);

  if (opt_.sigma_integration == "ppm") {
    auto* ppm_a = dynamic_cast<Sigma_PPM_UKS*>(sigma_alpha_.get());
    auto* ppm_b = dynamic_cast<Sigma_PPM_UKS*>(sigma_beta_.get());

    if (ppm_a == nullptr || ppm_b == nullptr) {
      throw std::runtime_error(
          "GW_UKS: expected Sigma_PPM_UKS evaluators for "
          "sigma_integration=ppm");
    }

    ppm_a->SetSharedPPM(ppm_);
    ppm_b->SetSharedPPM(ppm_);
  }
}

const Eigen::VectorXd& GW_UKS::DftEnergies(Spin spin) const {
  return (spin == Spin::Alpha) ? dft_energies_alpha_ : dft_energies_beta_;
}
const Eigen::MatrixXd& GW_UKS::Vxc(Spin spin) const {
  return (spin == Spin::Alpha) ? vxc_alpha_ : vxc_beta_;
}
Eigen::MatrixXd& GW_UKS::SigmaX(Spin spin) {
  return (spin == Spin::Alpha) ? Sigma_x_alpha_ : Sigma_x_beta_;
}
Eigen::MatrixXd& GW_UKS::SigmaC(Spin spin) {
  return (spin == Spin::Alpha) ? Sigma_c_alpha_ : Sigma_c_beta_;
}
const Eigen::MatrixXd& GW_UKS::SigmaX(Spin spin) const {
  return (spin == Spin::Alpha) ? Sigma_x_alpha_ : Sigma_x_beta_;
}
const Eigen::MatrixXd& GW_UKS::SigmaC(Spin spin) const {
  return (spin == Spin::Alpha) ? Sigma_c_alpha_ : Sigma_c_beta_;
}
Sigma_base_UKS& GW_UKS::SigmaEvaluator(Spin spin) {
  return (spin == Spin::Alpha) ? *sigma_alpha_ : *sigma_beta_;
}
const Sigma_base_UKS& GW_UKS::SigmaEvaluator(Spin spin) const {
  return (spin == Spin::Alpha) ? *sigma_alpha_ : *sigma_beta_;
}
Index GW_UKS::Homo(Spin spin) const {
  return (spin == Spin::Alpha) ? opt_.homo_alpha : opt_.homo_beta;
}
const char* GW_UKS::SpinName(Spin spin) const {
  return (spin == Spin::Alpha) ? "alpha" : "beta";
}

std::string GW_UKS::LevelLabel(Spin spin, Index level) const {
  if (level == Homo(spin)) {
    return "  HOMO ";
  } else if (level == Homo(spin) + 1) {
    return "  LUMO ";
  } else {
    return "  Level";
  }
}

const char* GW_UKS::OccupationTag(Spin spin, Index level) const {
  return (level <= Homo(spin)) ? "occ" : "virt";
}

Eigen::VectorXd GW_UKS::ScissorShift_DFTlevel(
    const Eigen::VectorXd& dft_energies, Index homo) const {
  Eigen::VectorXd shifted_energies = dft_energies;
  shifted_energies.segment(homo + 1, dft_energies.size() - homo - 1).array() +=
      opt_.shift;
  return shifted_energies;
}

double GW_UKS::CalcSpinHomoLumoShift(const Eigen::VectorXd& frequencies,
                                     Spin spin) const {
  const Index homo = Homo(spin);
  const Eigen::VectorXd& dft = DftEnergies(spin);
  const double DFTgap = dft(homo + 1) - dft(homo);
  const double QPgap =
      frequencies(homo + 1 - opt_.qpmin) - frequencies(homo - opt_.qpmin);
  return QPgap - DFTgap;
}

void GW_UKS::CalculateGWPerturbation() {
  Sigma_x_alpha_ = (1 - opt_.ScaHFX) * sigma_alpha_->CalcExchangeMatrix();
  Sigma_x_beta_ = (1 - opt_.ScaHFX) * sigma_beta_->CalcExchangeMatrix();
  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << " Calculated spin-resolved Hartree exchange contribution"
      << std::flush;

  Eigen::VectorXd dft_shifted_alpha =
      ScissorShift_DFTlevel(dft_energies_alpha_, opt_.homo_alpha);
  Eigen::VectorXd dft_shifted_beta =
      ScissorShift_DFTlevel(dft_energies_beta_, opt_.homo_beta);

  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << " Scissor shifting alpha/beta DFT energies by: " << opt_.shift
      << " Hrt" << std::flush;

  rpa_.setRPAInputEnergies(
      dft_shifted_alpha.segment(opt_.rpamin, opt_.rpamax - opt_.rpamin + 1),
      dft_shifted_beta.segment(opt_.rpamin, opt_.rpamax - opt_.rpamin + 1));

  Eigen::VectorXd frequencies_alpha =
      dft_shifted_alpha.segment(opt_.qpmin, qptotal_);
  Eigen::VectorXd frequencies_beta =
      dft_shifted_beta.segment(opt_.qpmin, qptotal_);

  Anderson mixing_alpha;
  Anderson mixing_beta;
  mixing_alpha.Configure(opt_.gw_mixing_order, opt_.gw_mixing_alpha);
  mixing_beta.Configure(opt_.gw_mixing_order, opt_.gw_mixing_alpha);

  for (Index i_gw = 0; i_gw < opt_.gw_sc_max_iterations; ++i_gw) {
    gw_sc_iteration_ = i_gw;
    if (i_gw % opt_.reset_3c == 0 && i_gw != 0) {
      Mmn_.alpha.Rebuild();
      Mmn_.beta.Rebuild();
      XTP_LOG(Log::info, log_)
          << TimeStamp() << " Rebuilding alpha/beta 3c integrals" << std::flush;
    }

    if (opt_.sigma_integration == "ppm") {
      ppm_.PPM_construct_parameters(rpa_);
      sigma_alpha_->PrepareScreening();
      sigma_beta_->PrepareScreening();
    } else {
      sigma_alpha_->PrepareScreening();
      sigma_beta_->PrepareScreening();
    }

    XTP_LOG(Log::info, log_)
        << TimeStamp() << " Calculated unrestricted screening via RPA"
        << std::flush;

    if (opt_.gw_mixing_order > 0 && i_gw > 0) {
      mixing_alpha.UpdateInput(frequencies_alpha);
      mixing_beta.UpdateInput(frequencies_beta);
    }

    frequencies_alpha = SolveQP(Spin::Alpha, frequencies_alpha);
    frequencies_beta = SolveQP(Spin::Beta, frequencies_beta);

    if (opt_.gw_sc_max_iterations > 1) {
      Eigen::VectorXd rpa_alpha_old = rpa_.getRPAInputEnergiesAlpha();
      Eigen::VectorXd rpa_beta_old = rpa_.getRPAInputEnergiesBeta();

      if (opt_.gw_mixing_order > 0 && i_gw > 0) {
        mixing_alpha.UpdateOutput(frequencies_alpha);
        mixing_beta.UpdateOutput(frequencies_beta);
        frequencies_alpha = mixing_alpha.MixHistory();
        frequencies_beta = mixing_beta.MixHistory();
      }

      rpa_.UpdateRPAInputEnergies(dft_energies_alpha_, dft_energies_beta_,
                                  frequencies_alpha, frequencies_beta,
                                  opt_.qpmin);

      XTP_LOG(Log::info, log_)
          << TimeStamp() << " GW_Iteration:" << i_gw << " Shift_alpha[Hrt]:"
          << CalcSpinHomoLumoShift(frequencies_alpha, Spin::Alpha)
          << " Shift_beta[Hrt]:"
          << CalcSpinHomoLumoShift(frequencies_beta, Spin::Beta) << std::flush;

      const bool converged_alpha = Converged(rpa_.getRPAInputEnergiesAlpha(),
                                             rpa_alpha_old, opt_.gw_sc_limit);
      const bool converged_beta = Converged(rpa_.getRPAInputEnergiesBeta(),
                                            rpa_beta_old, opt_.gw_sc_limit);
      if (converged_alpha && converged_beta) {
        XTP_LOG(Log::info, log_)
            << TimeStamp() << " Converged after " << i_gw + 1
            << " unrestricted GW iterations." << std::flush;
        break;
      } else if (i_gw == opt_.gw_sc_max_iterations - 1) {
        XTP_LOG(Log::error, log_)
            << TimeStamp()
            << " WARNING! UKS GW-self-consistency cycle not converged after "
            << opt_.gw_sc_max_iterations << " iterations." << std::flush;
        break;
      }
    }
  }

  Sigma_c_alpha_.diagonal() =
      sigma_alpha_->CalcCorrelationDiag(frequencies_alpha);
  Sigma_c_beta_.diagonal() = sigma_beta_->CalcCorrelationDiag(frequencies_beta);
  PrintGWA_Energies(Spin::Alpha);
  PrintGWA_Energies(Spin::Beta);
}

Eigen::VectorXd GW_UKS::getGWAResultsAlpha() const {
  return Sigma_x_alpha_.diagonal() + Sigma_c_alpha_.diagonal() -
         vxc_alpha_.diagonal() +
         dft_energies_alpha_.segment(opt_.qpmin, qptotal_);
}
Eigen::VectorXd GW_UKS::getGWAResultsBeta() const {
  return Sigma_x_beta_.diagonal() + Sigma_c_beta_.diagonal() -
         vxc_beta_.diagonal() +
         dft_energies_beta_.segment(opt_.qpmin, qptotal_);
}
const Eigen::VectorXd& GW_UKS::RPAInputEnergiesAlpha() const {
  return rpa_.getRPAInputEnergiesAlpha();
}
const Eigen::VectorXd& GW_UKS::RPAInputEnergiesBeta() const {
  return rpa_.getRPAInputEnergiesBeta();
}

Eigen::VectorXd GW_UKS::SolveQP(Spin spin,
                                const Eigen::VectorXd& frequencies) const {
  const Eigen::VectorXd intercepts =
      DftEnergies(spin).segment(opt_.qpmin, qptotal_) +
      SigmaX(spin).diagonal() - Vxc(spin).diagonal();

  Eigen::VectorXd frequencies_new = frequencies;

  Index use_threads = qptotal_;
#ifdef _OPENMP
  use_threads =
      OPENMP::getMaxThreads() > qptotal_ ? qptotal_ : OPENMP::getMaxThreads();
#endif

#pragma omp parallel for schedule(dynamic) num_threads(use_threads)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    const double initial_f = frequencies[gw_level];
    const double intercept = intercepts[gw_level];

    boost::optional<double> newf;
    QPStats fixed_stats;
    QPStats grid_stats;
    QPStats lin_stats;

    if (opt_.qp_solver == "fixedpoint") {
      newf = SolveQP_FixedPoint(spin, intercept, initial_f, gw_level,
                                &fixed_stats);
    }

    if (newf) {
      frequencies_new[gw_level] = newf.value();
    } else {
      newf = SolveQP_Grid(spin, intercept, initial_f, gw_level, &grid_stats);
      if (newf) {
        frequencies_new[gw_level] = newf.value();
      } else {
        newf = SolveQP_Linearisation(spin, intercept, initial_f, gw_level,
                                     &lin_stats);
        if (newf) {
          frequencies_new[gw_level] = newf.value();
        }
      }
    }

    if (Log::current_level > Log::error) {
#pragma omp critical
      {
        QPStats total_stats;
        total_stats.Add(fixed_stats);
        total_stats.Add(grid_stats);
        total_stats.Add(lin_stats);

        const Index mo_level = opt_.qpmin + gw_level;

        XTP_LOG(Log::info, log_)
            << " QP stats " << LevelLabel(spin, mo_level) << " mo=" << mo_level
            << " spin=" << SpinName(spin)
            << " sigma_calls=" << total_stats.TotalSigmaCalls()
            << " scan=" << total_stats.sigma_scan_calls
            << " refine=" << total_stats.sigma_refine_calls
            << " other=" << total_stats.sigma_other_calls
            << " deriv_calls=" << total_stats.deriv_calls
            << " unique_omega=" << total_stats.sigma_unique_frequencies
            << " repeat_omega=" << total_stats.sigma_repeat_calls << std::flush;
      }
    }
  }

  return frequencies_new;
}

boost::optional<double> GW_UKS::SolveQP_Linearisation(Spin spin,
                                                      double intercept0,
                                                      double frequency0,
                                                      Index gw_level,
                                                      QPStats* stats) const {
  boost::optional<double> newf = boost::none;

  QPFunc fqp(gw_level, SigmaEvaluator(spin), intercept0);

  const double sigma = fqp.sigma(frequency0, EvalStage::Other);
  const double dsigma_domega = fqp.deriv(frequency0);
  const double Z = 1.0 - dsigma_domega;

  if (std::abs(Z) > 1e-9) {
    newf = frequency0 + (intercept0 - frequency0 + sigma) / Z;
  }

  if (stats != nullptr) {
    *stats = fqp.GetStats();
  }

  return newf;
}

boost::optional<double> GW_UKS::SolveQP_Grid(Spin spin, double intercept0,
                                             double frequency0, Index gw_level,
                                             QPStats* stats) const {
  const double range =
      opt_.qp_grid_spacing * double(opt_.qp_grid_steps - 1) / 2.0;

  const double full_left_limit = frequency0 - range;
  const double full_right_limit = frequency0 + range;

  double restricted_left_limit = full_left_limit;
  double restricted_right_limit = full_right_limit;

  bool use_restricted_window = false;

  if (opt_.qp_restrict_search) {
    const Index mo_level = gw_level + opt_.qpmin;
    const bool is_occupied = (mo_level <= Homo(spin));

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
    auto restricted = SolveQP_Grid_Windowed(spin, intercept0, frequency0,
                                            gw_level, restricted_left_limit,
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
            << " Restricted QP search failed for "
            << LevelLabel(spin, opt_.qpmin + gw_level) << " in window ["
            << restricted_left_limit << ", " << restricted_right_limit
            << "], forcing full dense window [" << full_left_limit << ", "
            << full_right_limit << "]" << std::flush;
      }
    }

    return SolveQP_Grid_Windowed_Dense(spin, intercept0, frequency0, gw_level,
                                       full_left_limit, full_right_limit,
                                       true, stats);
  }

  return SolveQP_Grid_Windowed(spin, intercept0, frequency0, gw_level,
                               full_left_limit, full_right_limit,
                               true, stats);
  }

boost::optional<double> GW_UKS::SolveQP_Grid_Windowed_Adaptive(
    Spin spin, double intercept0, double frequency0, Index gw_level,
    double left_limit, double right_limit, bool allow_rejected_return,
    QPStats* stats) const {

  QPFunc fqp(gw_level, SigmaEvaluator(spin), intercept0);

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
      if (accepted_roots.empty() && rejected_roots.empty()) {
        XTP_LOG(Log::info, log_)
            << " No roots found for " << LevelLabel(spin, opt_.qpmin + gw_level)
            << " (center="
            << std::max(left_limit, std::min(right_limit, frequency0))
            << ", shells=" << wdiag.shells_explored << ")" << std::flush;
      } else {
        XTP_LOG(Log::info, log_)
            << " Adaptive scan " << LevelLabel(spin, opt_.qpmin + gw_level)
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
            << " Adaptive scan " << LevelLabel(spin, opt_.qpmin + gw_level)
            << " produced only rejected roots in [" << left_limit << ", "
            << right_limit << "], forcing wider retry" << std::flush;
      }
    }
    return boost::none;
  }

  return result;
}

boost::optional<double> GW_UKS::SolveQP_Grid_Windowed_Dense(
    Spin spin, double intercept0, double frequency0, Index gw_level,
    double left_limit, double right_limit, bool allow_rejected_return,
    QPStats* stats) const {

  QPFunc fqp(gw_level, SigmaEvaluator(spin), intercept0);

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
            << " Dense scan " << LevelLabel(spin, opt_.qpmin + gw_level)
            << " found no roots in [" << left_limit << ", " << right_limit
            << "]" << std::flush;
      } else {
        XTP_LOG(Log::info, log_)
            << " Dense scan " << LevelLabel(spin, opt_.qpmin + gw_level)
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
              << " Dense scan " << LevelLabel(spin, opt_.qpmin + gw_level)
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

boost::optional<double> GW_UKS::SolveQP_Grid_Windowed(
    Spin spin, double intercept0, double frequency0, Index gw_level,
    double left_limit, double right_limit, bool allow_rejected_return,
    QPStats* stats) const {

  if (opt_.qp_grid_search_mode == "adaptive") {
    return SolveQP_Grid_Windowed_Adaptive(spin, intercept0, frequency0,
                                          gw_level, left_limit, right_limit,
                                          allow_rejected_return, stats);
  }

  if (opt_.qp_grid_search_mode == "dense") {
    return SolveQP_Grid_Windowed_Dense(spin, intercept0, frequency0, gw_level,
                                       left_limit, right_limit, 
                                       allow_rejected_return, stats);
  }

  if (opt_.qp_grid_search_mode == "adaptive_with_dense_fallback") {
    QPStats total_stats;
    auto adaptive =
        SolveQP_Grid_Windowed_Adaptive(spin, intercept0, frequency0, gw_level,
                                       left_limit, right_limit,
                                       allow_rejected_return, &total_stats);
    if (adaptive) {
      if (stats != nullptr) {
        *stats = total_stats;
      }
      return adaptive;
    }

    QPStats dense_stats;
    auto dense =
        SolveQP_Grid_Windowed_Dense(spin, intercept0, frequency0, gw_level,
                                    left_limit, right_limit, 
                                    allow_rejected_return, &dense_stats);
    total_stats.Add(dense_stats);

    if (Log::current_level > Log::error) {
#pragma omp critical
      {
        XTP_LOG(Log::info, log_)
            << " Adaptive QP scan failed for "
            << LevelLabel(spin, opt_.qpmin + gw_level)
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

boost::optional<double> GW_UKS::SolveQP_FixedPoint(Spin spin, double intercept0,
                                                   double frequency0,
                                                   Index gw_level,
                                                   QPStats* stats) const {
  boost::optional<double> newf = boost::none;

  QPFunc f(gw_level, SigmaEvaluator(spin), intercept0);
  NewtonRapson<QPFunc> newton = NewtonRapson<QPFunc>(
      opt_.g_sc_max_iterations, opt_.g_sc_limit, opt_.qp_solver_alpha);

  const double freq_new = newton.FindRoot(f, frequency0);
  if (newton.getInfo() == NewtonRapson<QPFunc>::success) {
    newf = freq_new;
  }

  if (stats != nullptr) {
    *stats = f.GetStats();
  }

  return newf;
}

bool GW_UKS::Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                       double epsilon) const {
  Index state = 0;
  const double diff_max = (e1 - e2).cwiseAbs().maxCoeff(&state);
  XTP_LOG(Log::info, log_) << TimeStamp() << " E_diff max=" << diff_max
                           << " StateNo:" << state << std::flush;
  return diff_max <= epsilon;
}

void GW_UKS::CalculateHQP() {
  Eigen::VectorXd diag_backup_alpha = Sigma_c_alpha_.diagonal();
  Eigen::VectorXd diag_backup_beta = Sigma_c_beta_.diagonal();
  Sigma_c_alpha_ = sigma_alpha_->CalcCorrelationOffDiag(getGWAResultsAlpha());
  Sigma_c_beta_ = sigma_beta_->CalcCorrelationOffDiag(getGWAResultsBeta());
  Sigma_c_alpha_.diagonal() = diag_backup_alpha;
  Sigma_c_beta_.diagonal() = diag_backup_beta;
}

Eigen::MatrixXd GW_UKS::getHQPAlpha() const {
  return Sigma_x_alpha_ + Sigma_c_alpha_ - vxc_alpha_ +
         Eigen::MatrixXd(
             dft_energies_alpha_.segment(opt_.qpmin, qptotal_).asDiagonal());
}
Eigen::MatrixXd GW_UKS::getHQPBeta() const {
  return Sigma_x_beta_ + Sigma_c_beta_ - vxc_beta_ +
         Eigen::MatrixXd(
             dft_energies_beta_.segment(opt_.qpmin, qptotal_).asDiagonal());
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
    GW_UKS::DiagonalizeQPHamiltonianAlpha() const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQPAlpha());
  PrintQP_Energies(Spin::Alpha, qpdiag.eigenvalues());
  return qpdiag;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
    GW_UKS::DiagonalizeQPHamiltonianBeta() const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQPBeta());
  PrintQP_Energies(Spin::Beta, qpdiag.eigenvalues());
  return qpdiag;
}

void GW_UKS::PrintGWA_Energies(Spin spin) const {
  const Eigen::VectorXd gwa =
      (spin == Spin::Alpha) ? getGWAResultsAlpha() : getGWAResultsBeta();

  const Index homo = Homo(spin);
  const Eigen::VectorXd& dft = DftEnergies(spin);

  XTP_LOG(Log::error, log_)
      << (boost::format("  ====== Perturbative %1% quasiparticle energies "
                        "(Hartree) ======") %
          SpinName(spin))
             .str()
      << std::flush;

  if (opt_.qpmin <= homo && homo + 1 <= opt_.qpmax) {
    const double dft_gap = dft(homo + 1) - dft(homo);
    const double qp_gap = gwa(homo + 1 - opt_.qpmin) - gwa(homo - opt_.qpmin);

    XTP_LOG(Log::error, log_)
        << (boost::format("   %1% DFT gap = %2$+1.6f Hartree   %1% GWA gap = "
                          "%3$+1.6f Hartree   Delta = %4$+1.6f Hartree") %
            SpinName(spin) % dft_gap % qp_gap % (qp_gap - dft_gap))
               .str()
        << std::flush;
  }

  for (Index i = 0; i < qptotal_; i++) {
    const Index level = i + opt_.qpmin;
    XTP_LOG(Log::error, log_)
        << LevelLabel(spin, level)
        << (boost::format(" = %1$4d (%2%) %3% DFT = %4$+1.4f VXC = %5$+1.4f "
                          "S-X = %6$+1.4f S-C = %7$+1.4f GWA = %8$+1.4f") %
            level % OccupationTag(spin, level) % SpinName(spin) %
            DftEnergies(spin)(level) % Vxc(spin)(i, i) % SigmaX(spin)(i, i) %
            SigmaC(spin)(i, i) % gwa(i))
               .str()
        << std::flush;
  }
}

void GW_UKS::PrintQP_Energies(Spin spin,
                              const Eigen::VectorXd& qp_diag_energies) const {
  const Eigen::VectorXd gwa =
      (spin == Spin::Alpha) ? getGWAResultsAlpha() : getGWAResultsBeta();

  const Index homo = Homo(spin);
  const Eigen::VectorXd& dft = DftEnergies(spin);

  XTP_LOG(Log::error, log_) << TimeStamp() << " Full " << SpinName(spin)
                            << " quasiparticle Hamiltonian" << std::flush;

  if (opt_.qpmin <= homo && homo + 1 <= opt_.qpmax) {
    const double dft_gap = dft(homo + 1) - dft(homo);
    const double pqp_gap = gwa(homo + 1 - opt_.qpmin) - gwa(homo - opt_.qpmin);
    const double dqp_gap = qp_diag_energies(homo + 1 - opt_.qpmin) -
                           qp_diag_energies(homo - opt_.qpmin);

    XTP_LOG(Log::error, log_)
        << (boost::format("   %1% DFT gap = %2$+1.6f Hartree   %1% PQP gap = "
                          "%3$+1.6f Hartree   %1% DQP gap = %4$+1.6f Hartree") %
            SpinName(spin) % dft_gap % pqp_gap % dqp_gap)
               .str()
        << std::flush;
  }

  for (Index i = 0; i < qptotal_; i++) {
    const Index level = i + opt_.qpmin;
    XTP_LOG(Log::error, log_)
        << LevelLabel(spin, level)
        << (boost::format(" = %1$4d (%2%) %3% PQP = %4$+1.6f DQP = %5$+1.6f") %
            level % OccupationTag(spin, level) % SpinName(spin) % gwa(i) %
            qp_diag_energies(i))
               .str()
        << std::flush;
  }
}

}  // namespace xtp
}  // namespace votca
