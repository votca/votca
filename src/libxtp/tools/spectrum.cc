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

#include "spectrum.h"
#include <boost/math/constants/constants.hpp>
#include <votca/tools/constants.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

void Spectrum::Initialize(tools::Property& options) {

  // update options with the VOTCASHARE defaults
  UpdateWithDefaults(options, "xtp");
  std::string key = "options." + Identify();

  // orbitals file or pure DFT output
  _orbfile = options.get(key + ".input").as<std::string>();
  _output_file = options.get(key + ".output").as<std::string>();
  _n_pt = options.get(key + ".points").as<int>();
  _lower = options.get(key + ".lower").as<double>();
  _upper = options.get(key + ".upper").as<double>();
  _fwhm = options.get(key + ".fwhm").as<double>();
  _spectrum_type = options.get(key + ".type").as<std::string>();
  _minexc = options.get(key + ".minexc").as<int>();
  exc_lambda = options.get(key + ".maxexc").as<int>();
  _shiftby = options.get(key + ".shift").as<double>();

  return;
}

bool Spectrum::Evaluate() {

  _log.setReportLevel(logDEBUG);
  _log.setMultithreading(true);

  _log.setPreface(logINFO, "\n... ...");
  _log.setPreface(logERROR, "\n... ...");
  _log.setPreface(logWARNING, "\n... ...");
  _log.setPreface(logDEBUG, "\n... ...");

  XTP_LOG(logDEBUG, _log) << "Calculating absorption spectrum plot " << _orbfile
                          << std::flush;

  Orbitals orbitals;
  // load the QM data from serialized orbitals object
  XTP_LOG(logDEBUG, _log) << " Loading QM data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  // check if orbitals contains singlet energies and transition dipoles
  if (!orbitals.hasBSESinglets()) {
    throw std::runtime_error(
        "BSE singlet energies not stored in QM data file!");
  }

  if (!orbitals.hasTransitionDipoles()) {
    throw std::runtime_error(
        "BSE transition dipoles not stored in QM data file!");
  }

  const Eigen::VectorXd& BSESingletEnergies = orbitals.BSESingletEnergies();
  const std::vector<Eigen::Vector3d>& TransitionDipoles =
      orbitals.TransitionDipoles();
  std::vector<double> osc = orbitals.Oscillatorstrengths();

  int n_exc = exc_lambda - _minexc + 1;

  if (exc_lambda > int(TransitionDipoles.size())) {
    XTP_LOG(logDEBUG, _log)
        << " Transition dipoles for some excitations missing! " << std::flush;
    exc_lambda = int(TransitionDipoles.size());
  }

  XTP_LOG(logDEBUG, _log) << " Considering " << n_exc
                          << " excitation with max energy "
                          << BSESingletEnergies(exc_lambda) *
                                 tools::conv::hrt2ev
                          << " eV / min wave length "
                          << evtonm(BSESingletEnergies[exc_lambda - 1] *
                                    tools::conv::hrt2ev)
                          << " nm" << std::flush;

  /*
   *
   * For a single excitation, broaden by Lineshape function L(v-W)
   *    eps(v) = f * L(v-W)
   *
   * where
   *       v: energy
   *       f: oscillator strength in dipole-length gauge
   *       W: excitation energy
   *
   * Lineshape function depend on FWHM and can be
   *
   *      Gaussian
   *          L(v-W) = 1/(sqrt(2pi)sigma) * exp(-0.5 (v-W)^2/sigma^2
   *
   *
   *            with sigma: derived from FWHM (FWHM/2.3548)
   *
   *     Lorentzian
   *          L(v-W) = 1/pi * 0.5 FWHM/( (v-w)^2 + 0.25*FWHM^2 )
   *
   * Full spectrum is superposition of individual spectra.
   *
   *  Alternatively, one can calculate the imaginary part of the
   *  frequency-dependent dielectric function
   *
   *   IM(eps(v)) ~ 1/v^2 * W^2 * |td|^2 * L(v-W)
   *              = 1/v^2 * W   * f      * L(v-W)
   *
   *
   */

  std::ofstream ofs(_output_file, std::ofstream::out);

  if (_spectrum_type == "energy") {
    _fwhm /= tools::conv::hrt2ev;
    ofs << "# E(eV)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (int i_pt = 0; i_pt <= _n_pt; i_pt++) {

      double e =
          (_lower + i_pt * (_upper - _lower) / _n_pt) / tools::conv::hrt2ev;

      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;
      double eps_TruncLorentzian = 0.0;
      double imeps_TruncLorentzian = 0.0;

      for (int i_exc = _minexc; i_exc <= exc_lambda; i_exc++) {
        eps_Gaussian +=
            osc[i_exc] *
            Gaussian(e,
                     BSESingletEnergies(i_exc) + _shiftby / tools::conv::hrt2ev,
                     _fwhm);
        imeps_Gaussian += osc[i_exc] * BSESingletEnergies(i_exc) *
                          Gaussian(e, BSESingletEnergies(i_exc), _fwhm);
        eps_Lorentzian +=
            osc[i_exc] * Lorentzian(e, BSESingletEnergies(i_exc), _fwhm);
        imeps_Lorentzian += osc[i_exc] * BSESingletEnergies(i_exc) *
                            Lorentzian(e, BSESingletEnergies(i_exc), _fwhm);
        eps_TruncLorentzian +=
            osc[i_exc] *
            TruncatedLorentzian(e, BSESingletEnergies(i_exc), _fwhm);
        imeps_TruncLorentzian +=
            osc[i_exc] * BSESingletEnergies(i_exc) *
            TruncatedLorentzian(e, BSESingletEnergies(i_exc), _fwhm);
      }

      ofs << e * tools::conv::hrt2ev << "    " << eps_Gaussian << "   "
          << imeps_Gaussian << "   " << eps_Lorentzian << "   "
          << imeps_Lorentzian << "  " << eps_TruncLorentzian << "   "
          << imeps_TruncLorentzian << std::endl;
    }

    XTP_LOG(logDEBUG, _log)
        << " Spectrum in energy range from  " << _lower << " to " << _upper
        << " eV and with broadening of FWHM " << _fwhm * tools::conv::hrt2ev
        << " eV written to file  " << _output_file << std::flush;
  }

  if (_spectrum_type == "wavelength") {

    ofs << "# lambda(nm)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (int i_pt = 0; i_pt <= _n_pt; i_pt++) {

      double lambda = (_lower + i_pt * (_upper - _lower) / _n_pt);
      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;
      double eps_TruncLorentzian = 0.0;
      double imeps_TruncLorentzian = 0.0;

      for (int i_exc = _minexc; i_exc <= exc_lambda; i_exc++) {
        double exc_lambda =
            nmtoev(BSESingletEnergies(i_exc) * tools::conv::hrt2ev + _shiftby);
        eps_Gaussian += osc[i_exc] * Gaussian(lambda, exc_lambda, _fwhm);
        imeps_Gaussian +=
            osc[i_exc] * exc_lambda * Gaussian(lambda, exc_lambda, _fwhm);
        eps_Lorentzian += osc[i_exc] * Lorentzian(lambda, exc_lambda, _fwhm);
        imeps_Lorentzian +=
            osc[i_exc] * exc_lambda * Lorentzian(lambda, exc_lambda, _fwhm);
        eps_TruncLorentzian +=
            osc[i_exc] * TruncatedLorentzian(lambda, exc_lambda, _fwhm);
        imeps_TruncLorentzian += osc[i_exc] * BSESingletEnergies(i_exc) *
                                 TruncatedLorentzian(lambda, exc_lambda, _fwhm);
      }

      ofs << lambda << "    " << eps_Gaussian << "   " << imeps_Gaussian
          << "   " << eps_Lorentzian << "   " << imeps_Lorentzian << "   "
          << eps_TruncLorentzian << "   " << imeps_TruncLorentzian << std::endl;
    }
    XTP_LOG(logDEBUG, _log)
        << " Spectrum in wavelength range from  " << _lower << " to " << _upper
        << " nm and with broadening of FWHM " << _fwhm
        << " nm written to file  " << _output_file << std::flush;
  }

  ofs.close();
  return true;
}

double Spectrum::TruncatedLorentzian(double x, double center, double fwhm) {

  double result;
  double abs_diff = std::abs(x - center);
  if (abs_diff > 0.5 * fwhm && abs_diff < fwhm) {
    result = 1.0 / (0.25 * fwhm * fwhm) -
             1.0 / (pow(abs_diff - fwhm, 2) + 0.25 * fwhm * fwhm);
  } else if (abs_diff < 0.5 * fwhm) {
    result = 1.0 / (std::pow(x - center, 2) + 0.25 * fwhm * fwhm);
  } else {
    result = 0.0;
  }
  return 0.5 * fwhm * result / boost::math::constants::pi<double>();
}

double Spectrum::Lorentzian(double x, double center, double fwhm) {
  return 0.5 * fwhm / (std::pow(x - center, 2) + 0.25 * fwhm * fwhm) /
         boost::math::constants::pi<double>();
}

double Spectrum::Gaussian(double x, double center, double fwhm) {
  // FWHM = 2*sqrt(2 ln2) sigma = 2.3548 sigma
  double sigma = fwhm / 2.3548;
  return std::exp(-0.5 * std::pow((x - center) / sigma, 2)) / sigma /
         sqrt(2.0 * boost::math::constants::pi<double>());
}

double Spectrum::evtonm(double eV) { return 1241.0 / eV; }

double Spectrum::evtoinvcm(double eV) { return 8065.73 * eV; }

double Spectrum::nmtoinvcm(double nm) { return 1241.0 * 8065.73 / nm; }

double Spectrum::invcmtonm(double invcm) { return 1.0e7 / invcm; }

double Spectrum::nmtoev(double nm) { return 1241.0 / nm; }

}  // namespace xtp
}  // namespace votca
