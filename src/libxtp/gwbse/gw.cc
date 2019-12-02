/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_exact.h"
#include "votca/xtp/sigma_ppm.h"
#include <votca/xtp/gw.h>

namespace votca {
namespace xtp {

void GW::configure(const options& opt) {
  _opt = opt;
  _qptotal = _opt.qpmax - _opt.qpmin + 1;
  _rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
  if (_opt.sigma_integration == "exact") {
    _sigma = std::make_unique<Sigma_Exact>(Sigma_Exact(_Mmn, _rpa));
  } else if (_opt.sigma_integration == "ppm") {
    _sigma = std::make_unique<Sigma_PPM>(Sigma_PPM(_Mmn, _rpa));
  }
  Sigma_base::options sigma_opt;
  sigma_opt.homo = _opt.homo;
  sigma_opt.qpmax = _opt.qpmax;
  sigma_opt.qpmin = _opt.qpmin;
  sigma_opt.rpamin = _opt.rpamin;
  sigma_opt.rpamax = _opt.rpamax;
  sigma_opt.eta = _opt.eta;
  _sigma->configure(sigma_opt);
  _Sigma_x = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  _Sigma_c = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
}

double GW::CalcHomoLumoShift(Eigen::VectorXd frequencies) const {
  double DFTgap = _dft_energies(_opt.homo + 1) - _dft_energies(_opt.homo);
  double QPgap = frequencies(_opt.homo + 1 - _opt.qpmin) -
                 frequencies(_opt.homo - _opt.qpmin);
  return QPgap - DFTgap;
}

Eigen::MatrixXd GW::getHQP() const {
  return _Sigma_x + _Sigma_c - _vxc +
         Eigen::MatrixXd(
             _dft_energies.segment(_opt.qpmin, _qptotal).asDiagonal());
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> GW::DiagonalizeQPHamiltonian()
    const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQP());
  PrintQP_Energies(qpdiag.eigenvalues());
  return qpdiag;
}

void GW::PrintGWA_Energies() const {

  double shift = CalcHomoLumoShift(_gwa_energies);

  XTP_LOG(Log::error, _log)
      << (boost::format(
              "  ====== Perturbative quasiparticle energies (Hartree) ====== "))
             .str()
      << std::flush;
  XTP_LOG(Log::error, _log)
      << (boost::format("   DeltaHLGap = %1$+1.6f Hartree") % shift).str()
      << std::flush;

  for (Index i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }

    XTP_LOG(Log::error, _log)
        << level
        << (boost::format(" = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = "
                          "%4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") %
            (i + _opt.qpmin) % _dft_energies(i + _opt.qpmin) % _vxc(i, i) %
            _Sigma_x(i, i) % _Sigma_c(i, i) % _gwa_energies(i))
               .str()
        << std::flush;
  }
  return;
}

void GW::PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const {
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Full quasiparticle Hamiltonian  " << std::flush;
  XTP_LOG(Log::error, _log)
      << (boost::format(
              "  ====== Diagonalized quasiparticle energies (Hartree) "
              "====== "))
             .str()
      << std::flush;
  for (Index i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }
    XTP_LOG(Log::error, _log)
        << level
        << (boost::format(" = %1$4d PQP = %2$+1.4f DQP = %3$+1.4f ") %
            (i + _opt.qpmin) % _gwa_energies(i) % qp_diag_energies(i))
               .str()
        << std::flush;
  }
  return;
}

Eigen::VectorXd GW::ScissorShift_DFTlevel(
    const Eigen::VectorXd& dft_energies) const {
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Scissor shifting DFT energies by: " << _opt.shift
      << " Hrt" << std::flush;
  Eigen::VectorXd shifted_energies = dft_energies;
  shifted_energies.segment(_opt.homo + 1, dft_energies.size() - _opt.homo - 1)
      .array() += _opt.shift;
  return shifted_energies;
}

void GW::CalculateGWPerturbation() {

  _Sigma_x = (1 - _opt.ScaHFX) * _sigma->CalcExchangeMatrix();
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Calculated Hartree exchange contribution"
      << std::flush;
  // dftenergies has size aobasissize
  // rpaenergies/Mmn have size rpatotal
  // gwaenergies/frequencies have size qptotal
  // homo index is relative to dft_energies
  Eigen::VectorXd dft_shifted_energies = ScissorShift_DFTlevel(_dft_energies);
  Eigen::VectorXd rpa_energies =
      dft_shifted_energies.segment(_opt.rpamin, _opt.rpamax - _opt.rpamin + 1);
  _rpa.setRPAInputEnergies(rpa_energies);
  Eigen::VectorXd frequencies =
      dft_shifted_energies.segment(_opt.qpmin, _qptotal);
  for (Index i_gw = 0; i_gw < _opt.gw_sc_max_iterations; ++i_gw) {

    if (i_gw % _opt.reset_3c == 0 && i_gw != 0) {
      _Mmn.Rebuild();
      XTP_LOG(Log::info, _log)
          << TimeStamp() << " Rebuilding 3c integrals" << std::flush;
    }
    _sigma->PrepareScreening();
    XTP_LOG(Log::info, _log)
        << TimeStamp() << " Calculated screening via RPA" << std::flush;
    if (_opt.qp_solver == "grid") {
      frequencies = SolveQP_Grid(frequencies);
      XTP_LOG(Log::error, _log)
          << TimeStamp() << " Solved QP equation on QP grid" << std::flush;
    } else if (_opt.qp_solver == "fixedpoint") {
      frequencies = SolveQP_FixedPoint(frequencies);
      XTP_LOG(Log::error, _log)
          << TimeStamp() << " Solved QP equation self-consistently"
          << std::flush;
    }
    _gwa_energies = frequencies;
    _Sigma_c.diagonal() = _sigma->CalcCorrelationDiag(_gwa_energies);
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Calculated correlation diagonal" << std::flush;
    Eigen::VectorXd rpa_energies_old = _rpa.getRPAInputEnergies();
    _rpa.UpdateRPAInputEnergies(_dft_energies, frequencies, _opt.qpmin);
    XTP_LOG(Log::info, _log)
        << TimeStamp() << " GW_Iteration:" << i_gw
        << " Shift[Hrt]:" << CalcHomoLumoShift(_gwa_energies) << std::flush;
    if (Converged(_rpa.getRPAInputEnergies(), rpa_energies_old,
                  _opt.gw_sc_limit)) {
      XTP_LOG(Log::error, _log) << TimeStamp() << " Converged after "
                                << i_gw + 1 << " GW iterations." << std::flush;
      break;
    } else if (i_gw == _opt.gw_sc_max_iterations - 1 &&
               _opt.gw_sc_max_iterations > 1) {
      XTP_LOG(Log::error, _log)
          << TimeStamp()
          << " WARNING! GW-self-consistency cycle not converged after "
          << _opt.gw_sc_max_iterations << " iterations." << std::flush;
      XTP_LOG(Log::error, _log)
          << TimeStamp() << "      Run continues. Inspect results carefully!"
          << std::flush;
      break;
    } else {
      double alpha = 0.0;
      rpa_energies = (1 - alpha) * rpa_energies + alpha * rpa_energies_old;
    }
  }

  PrintGWA_Energies();
}

Eigen::VectorXd GW::SolveQP_Grid(const Eigen::VectorXd& frequencies) const {
  const Index qptotal = _opt.qpmax - _opt.qpmin + 1;
  const double range =
      _opt.qp_grid_spacing * (double)(_opt.qp_grid_steps - 1) / 2.0;
  const Eigen::VectorXd intercept =
      _dft_energies.segment(_opt.qpmin, _qptotal) + _Sigma_x.diagonal() -
      _vxc.diagonal();
  Eigen::VectorXd frequencies_new = frequencies;
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal; ++gw_level) {
    const double frequency0 = frequencies[gw_level];
    const double intercept0 = intercept[gw_level];
    double freq_prev = frequency0 - range;
    double sigc_prev = _sigma->CalcCorrelation(gw_level, gw_level, freq_prev);
    double targ_prev = sigc_prev + intercept0 - freq_prev;
    double qp_energy = 0.0;
    double pole_weight_max = -1.0;
    for (Index i_node = 1; i_node < _opt.qp_grid_steps; ++i_node) {
      double freq = freq_prev + _opt.qp_grid_spacing;
      double sigc = _sigma->CalcCorrelation(gw_level, gw_level, freq);
      double targ = sigc + intercept0 - freq;
      if (targ_prev * targ < 0.0) {  // Sign change
        double dsigc_dfreq = (sigc - sigc_prev) / _opt.qp_grid_spacing;
        double dtarg_dfreq = (targ - targ_prev) / _opt.qp_grid_spacing;
        // Calculate fixed-point estimate of the pole (=root)
        double pole = freq_prev - targ_prev / dtarg_dfreq;
        // Calculate pole weight Z \in (0, 1)
        double pole_weight = 1.0 / (1.0 - dsigc_dfreq);
        if (pole_weight >= 1e-5 && pole_weight > pole_weight_max) {
          qp_energy = pole;
          pole_weight_max = pole_weight;
        }
      }
      freq_prev = freq;
      sigc_prev = sigc;
      targ_prev = targ;
    }
    if (pole_weight_max >= 0.0) {
      frequencies_new[gw_level] = qp_energy;
    }
  }
  return frequencies_new;
}

Eigen::VectorXd GW::SolveQP_FixedPoint(
    const Eigen::VectorXd& frequencies) const {
  const Index qptotal = _opt.qpmax - _opt.qpmin + 1;
  const Eigen::VectorXd intercept =
      _dft_energies.segment(_opt.qpmin, _qptotal) + _Sigma_x.diagonal() -
      _vxc.diagonal();
  Eigen::VectorXd frequencies_new = frequencies;
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal; ++gw_level) {
    const double intercept0 = intercept[gw_level];
    double freq = frequencies[gw_level];
    for (Index g_iter = 0; g_iter < _opt.g_sc_max_iterations; ++g_iter) {
      double freq_prev = freq;
      double sigc = _sigma->CalcCorrelation(gw_level, gw_level, freq);
      // Update frequency via fixed-point iteration
      freq = sigc + intercept0;
      if (std::abs(freq - freq_prev) < _opt.g_sc_limit) {
        break;  // Converged
      } else if (_opt.g_sc_max_iterations != 0 &&
                 g_iter == _opt.g_sc_max_iterations - 1) {
        break;  // Not converged
      }
    }
    frequencies_new[gw_level] = freq;
  }
  return frequencies_new;
}

bool GW::Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                   double epsilon) const {
  Index state = 0;
  bool energies_converged = true;
  double diff_max = (e1 - e2).cwiseAbs().maxCoeff(&state);
  if (diff_max > epsilon) {
    energies_converged = false;
  }
  XTP_LOG(Log::info, _log) << TimeStamp() << " E_diff max=" << diff_max
                           << " StateNo:" << state << std::flush;
  return energies_converged;
}

void GW::CalculateHQP() {
  _rpa.UpdateRPAInputEnergies(_dft_energies, _gwa_energies, _opt.qpmin);
  Eigen::VectorXd diag_backup = _Sigma_c.diagonal();
  _Sigma_c = _sigma->CalcCorrelationOffDiag(_gwa_energies);
  _Sigma_c.diagonal() = diag_backup;
}

}  // namespace xtp
}  // namespace votca
