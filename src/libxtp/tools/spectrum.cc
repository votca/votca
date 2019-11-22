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

  std::string key = "options." + Identify();

  // orbitals file or pure DFT output
  _orbfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".orbitals");
  _output_file = options.ifExistsReturnElseReturnDefault<std::string>(
      key + ".output", _output_file);
  _n_pt =
      options.ifExistsReturnElseReturnDefault<Index>(key + ".points", _n_pt);
  _lower = options.get(key + ".lower").as<double>();
  _upper = options.ifExistsReturnElseThrowRuntimeError<double>(key + ".upper");
  _fwhm = options.ifExistsReturnElseThrowRuntimeError<double>(key + ".fwhm");

  if (options.exists(key + ".type")) {
    std::vector<std::string> choices = {"energy", "wavelength"};
    _spectrum_type =
        options.ifExistsAndinListReturnElseThrowRuntimeError<std::string>(
            key + ".type", choices);
  }
  _minexc =
      options.ifExistsReturnElseReturnDefault<Index>(key + ".minexc", _minexc);
  _maxexc =
      options.ifExistsReturnElseReturnDefault<Index>(key + ".maxexc", _maxexc);
  _shiftby =
      options.ifExistsReturnElseReturnDefault<double>(key + ".shift", _shiftby);

  return;
}

bool Spectrum::Evaluate() {
  OPENMP::setMaxThreads(_nThreads);
  _log.setReportLevel(logDEBUG);
  _log.setMultithreading(true);

  _log.setPreface(logINFO, "\n... ...");
  _log.setPreface(logERROR, "\n... ...");
  _log.setPreface(logWARNING, "\n... ...");
  _log.setPreface(logDEBUG, "\n... ...");

  XTP_LOG_SAVE(logDEBUG, _log)
      << "Calculating absorption spectrum plot " << _orbfile << std::flush;

  Orbitals orbitals;
  // load the QM data from serialized orbitals object
  XTP_LOG_SAVE(logDEBUG, _log)
      << " Loading QM data from " << _orbfile << std::flush;
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

  const Eigen::VectorXd BSESingletEnergies =
      orbitals.BSESinglets().eigenvalues() * tools::conv::hrt2ev;
  const std::vector<Eigen::Vector3d>& TransitionDipoles =
      orbitals.TransitionDipoles();
  Eigen::VectorXd osc = orbitals.Oscillatorstrengths();

  if (_maxexc > Index(TransitionDipoles.size())) {
    _maxexc = Index(TransitionDipoles.size()) - 1;
  }

  Index n_exc = _maxexc - _minexc + 1;
  XTP_LOG_SAVE(logDEBUG, _log)
      << " Considering " << n_exc << " excitation with max energy "
      << BSESingletEnergies(_maxexc) << " eV / min wave length "
      << evtonm(BSESingletEnergies[_maxexc - 1]) << " nm" << std::flush;

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
    ofs << "# E(eV)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (Index i_pt = 0; i_pt <= _n_pt; i_pt++) {

      double e = (_lower + double(i_pt) * (_upper - _lower) / double(_n_pt));

      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;

      for (Index i_exc = _minexc; i_exc <= _maxexc; i_exc++) {
        eps_Gaussian +=
            osc[i_exc] *
            Gaussian(e, BSESingletEnergies(i_exc) + _shiftby, _fwhm);
        imeps_Gaussian += osc[i_exc] * BSESingletEnergies(i_exc) *
                          Gaussian(e, BSESingletEnergies(i_exc), _fwhm);
        eps_Lorentzian +=
            osc[i_exc] * Lorentzian(e, BSESingletEnergies(i_exc), _fwhm);
        imeps_Lorentzian += osc[i_exc] * BSESingletEnergies(i_exc) *
                            Lorentzian(e, BSESingletEnergies(i_exc), _fwhm);
      }

      ofs << e << "    " << eps_Gaussian << "   " << imeps_Gaussian << "   "
          << eps_Lorentzian << "   " << imeps_Lorentzian << std::endl;
    }

    XTP_LOG_SAVE(logDEBUG, _log)
        << " Spectrum in energy range from  " << _lower << " to " << _upper
        << " eV and with broadening of FWHM " << _fwhm
        << " eV written to file  " << _output_file << std::flush;
  }

  if (_spectrum_type == "wavelength") {

    ofs << "# lambda(nm)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (Index i_pt = 0; i_pt <= _n_pt; i_pt++) {

      double lambda =
          (_lower + double(i_pt) * (_upper - _lower) / double(_n_pt));
      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;

      for (Index i_exc = _minexc; i_exc <= _maxexc; i_exc++) {
        double exc_lambda = nmtoev(BSESingletEnergies(i_exc) + _shiftby);
        eps_Gaussian += osc[i_exc] * Gaussian(lambda, exc_lambda, _fwhm);
        imeps_Gaussian +=
            osc[i_exc] * exc_lambda * Gaussian(lambda, exc_lambda, _fwhm);
        eps_Lorentzian += osc[i_exc] * Lorentzian(lambda, exc_lambda, _fwhm);
        imeps_Lorentzian +=
            osc[i_exc] * exc_lambda * Lorentzian(lambda, exc_lambda, _fwhm);
      }

      ofs << lambda << "    " << eps_Gaussian << "   " << imeps_Gaussian
          << "   " << eps_Lorentzian << "   " << imeps_Lorentzian << std::endl;
    }
    XTP_LOG_SAVE(logDEBUG, _log)
        << " Spectrum in wavelength range from  " << _lower << " to " << _upper
        << " nm and with broadening of FWHM " << _fwhm
        << " nm written to file  " << _output_file << std::flush;
  }

  ofs.close();
  return true;
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
