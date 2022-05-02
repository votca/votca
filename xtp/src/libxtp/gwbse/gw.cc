/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
  const Eigen::VectorXd intercepts =
      dft_energies_.segment(opt_.qpmin, qptotal_) + Sigma_x_.diagonal() -
      vxc_.diagonal();
  Eigen::VectorXd frequencies_new = frequencies;
  Eigen::Array<bool, Eigen::Dynamic, 1> converged =
      Eigen::Array<bool, Eigen::Dynamic, 1>::Zero(qptotal_);
#ifdef _OPENMP
  Index use_threads =
      OPENMP::getMaxThreads() > qptotal_ ? qptotal_ : OPENMP::getMaxThreads();
#endif
#pragma omp parallel for schedule(dynamic) num_threads(use_threads)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {

    double initial_f = frequencies[gw_level];
    double intercept = intercepts[gw_level];
    boost::optional<double> newf;
    if (opt_.qp_solver == "fixedpoint") {
      newf = SolveQP_FixedPoint(intercept, initial_f, gw_level);
    }
    if (newf) {
      frequencies_new[gw_level] = newf.value();
      converged[gw_level] = true;
    } else {
      newf = SolveQP_Grid(intercept, initial_f, gw_level);
      if (newf) {
        frequencies_new[gw_level] = newf.value();
        converged[gw_level] = true;
      } else {
        newf = SolveQP_Linearisation(intercept, initial_f, gw_level);
        if (newf) {
          frequencies_new[gw_level] = newf.value();
        }
      }
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
  return frequencies_new;
}

boost::optional<double> GW::SolveQP_Linearisation(double intercept0,
                                                  double frequency0,
                                                  Index gw_level) const {
  boost::optional<double> newf = boost::none;

  double sigma = sigma_->CalcCorrelationDiagElement(gw_level, frequency0);
  double dsigma_domega =
      sigma_->CalcCorrelationDiagElementDerivative(gw_level, frequency0);
  double Z = 1.0 - dsigma_domega;
  if (std::abs(Z) > 1e-9) {
    newf = frequency0 + (intercept0 - frequency0 + sigma) / Z;
  }
  return newf;
}

boost::optional<double> GW::SolveQP_Grid(double intercept0, double frequency0,
                                         Index gw_level) const {
  std::vector<std::pair<double, double>> roots;
  const double range =
      opt_.qp_grid_spacing * double(opt_.qp_grid_steps - 1) / 2.0;
  boost::optional<double> newf = boost::none;
  double freq_prev = frequency0 - range;
  QPFunc fqp(gw_level, *sigma_.get(), intercept0);
  double targ_prev = fqp.value(freq_prev);
  double qp_energy = 0.0;
  double gradient_max = std::numeric_limits<double>::max();
  bool pole_found = false;
  for (Index i_node = 1; i_node < opt_.qp_grid_steps; ++i_node) {
    double freq = freq_prev + opt_.qp_grid_spacing;
    double targ = fqp.value(freq);
    if (targ_prev * targ < 0.0) {  // Sign change
      double f = SolveQP_Bisection(freq_prev, targ_prev, freq, targ, fqp);
      double gradient = fqp.deriv(f);
      double qp_weight = -1.0 / gradient;
      roots.push_back(std::make_pair(f, qp_weight));
      if (std::abs(gradient) < gradient_max) {
        gradient_max = std::abs(gradient);
        qp_energy = f;
        pole_found = true;
      }
    }
    freq_prev = freq;
    targ_prev = targ;
  }
  if (Log::current_level > Log::error) {
#pragma omp critical
    {
      if (!pole_found) {
        XTP_LOG(Log::info, log_)
            << " No roots found for qplevel:" << gw_level << std::flush;
      } else {
        XTP_LOG(Log::info, log_) << " Roots found for qplevel:" << gw_level
                                 << " (qpenergy:qpweight)\n\t\t";
        for (auto& root : roots) {
          XTP_LOG(Log::info, log_) << std::setprecision(5) << root.first << ":"
                                   << root.second << " ";
        }
        XTP_LOG(Log::info, log_) << "Root chosen " << qp_energy << std::flush;
      }
    }
  }

  if (pole_found) {
    newf = qp_energy;
  }
  return newf;
}

boost::optional<double> GW::SolveQP_FixedPoint(double intercept0,
                                               double frequency0,
                                               Index gw_level) const {
  boost::optional<double> newf = boost::none;
  QPFunc f(gw_level, *sigma_.get(), intercept0);
  NewtonRapson<QPFunc> newton = NewtonRapson<QPFunc>(
      opt_.g_sc_max_iterations, opt_.g_sc_limit, opt_.qp_solver_alpha);
  double freq_new = newton.FindRoot(f, frequency0);
  if (newton.getInfo() == NewtonRapson<QPFunc>::success) {
    newf = freq_new;
  }
  return newf;
}

// https://en.wikipedia.org/wiki/Bisection_method
double GW::SolveQP_Bisection(double lowerbound, double f_lowerbound,
                             double upperbound, double f_upperbound,
                             const QPFunc& f) const {

  if (f_lowerbound * f_upperbound > 0) {
    throw std::runtime_error(
        "Bisection needs a postive and negative function value");
  }
  double zero = 0.0;
  while (true) {
    double c = 0.5 * (lowerbound + upperbound);
    if (std::abs(upperbound - lowerbound) < opt_.g_sc_limit) {
      zero = c;
      break;
    }
    double y_c = f.value(c);
    if (std::abs(y_c) < opt_.g_sc_limit) {
      zero = c;
      break;
    }
    if (y_c * f_lowerbound > 0) {
      lowerbound = c;
      f_lowerbound = y_c;
    } else {
      upperbound = c;
      f_upperbound = y_c;
    }
  }
  return zero;
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
