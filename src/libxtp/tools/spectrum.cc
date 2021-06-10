/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Third party includes
#include <boost/math/constants/constants.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include <votca/xtp/orbitals.h>

// Local private VOTCA includes
#include "spectrum.h"

namespace votca {
namespace xtp {

void Spectrum::ParseOptions(const tools::Property& options) {

  // orbitals file or pure DFT output
  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".orbitals", job_name_ + ".orb");

  output_file_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + " spectrum_.dat");

  n_pt_ = options.get(".points").as<Index>();
  lower_ = options.get(".lower").as<double>();
  upper_ = options.get(".upper").as<double>();
  fwhm_ = options.get(".fwhm").as<double>();

  spectrum_type_ = options.get(".type").as<std::string>();
  minexc_ = options.get(".minexc").as<Index>();
  maxexc_ = options.get(".maxexc").as<Index>();
  shiftby_ = options.get(".shift").as<double>();
}

bool Spectrum::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  XTP_LOG(Log::error, log_)
      << "Calculating absorption spectrum plot " << orbfile_ << std::flush;

  Orbitals orbitals;
  // load the QM data from serialized orbitals object
  XTP_LOG(Log::error, log_)
      << " Loading QM data from " << orbfile_ << std::flush;
  orbitals.ReadFromCpt(orbfile_);

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

  if (maxexc_ > Index(TransitionDipoles.size())) {
    maxexc_ = Index(TransitionDipoles.size()) - 1;
  }

  Index n_exc = maxexc_ - minexc_ + 1;
  XTP_LOG(Log::error, log_)
      << " Considering " << n_exc << " excitation with max energy "
      << BSESingletEnergies(maxexc_) << " eV / min wave length "
      << evtonm(BSESingletEnergies[maxexc_ - 1]) << " nm" << std::flush;

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

  std::ofstream ofs(output_file_, std::ofstream::out);

  if (spectrum_type_ == "energy") {
    ofs << "# E(eV)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (Index i_pt = 0; i_pt <= n_pt_; i_pt++) {

      double e = (lower_ + double(i_pt) * (upper_ - lower_) / double(n_pt_));

      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;

      for (Index i_exc = minexc_; i_exc <= maxexc_; i_exc++) {
        eps_Gaussian +=
            osc[i_exc] *
            Gaussian(e, BSESingletEnergies(i_exc) + shiftby_, fwhm_);
        imeps_Gaussian += osc[i_exc] * BSESingletEnergies(i_exc) *
                          Gaussian(e, BSESingletEnergies(i_exc), fwhm_);
        eps_Lorentzian +=
            osc[i_exc] * Lorentzian(e, BSESingletEnergies(i_exc), fwhm_);
        imeps_Lorentzian += osc[i_exc] * BSESingletEnergies(i_exc) *
                            Lorentzian(e, BSESingletEnergies(i_exc), fwhm_);
      }

      ofs << e << "    " << eps_Gaussian << "   " << imeps_Gaussian << "   "
          << eps_Lorentzian << "   " << imeps_Lorentzian << std::endl;
    }

    XTP_LOG(Log::error, log_)
        << " Spectrum in energy range from  " << lower_ << " to " << upper_
        << " eV and with broadening of FWHM " << fwhm_
        << " eV written to file  " << output_file_ << std::flush;
  }

  if (spectrum_type_ == "wavelength") {

    ofs << "# lambda(nm)    epsGaussian    IM(eps)Gaussian   epsLorentz    "
           "Im(esp)Lorentz\n";
    for (Index i_pt = 0; i_pt <= n_pt_; i_pt++) {

      double lambda =
          (lower_ + double(i_pt) * (upper_ - lower_) / double(n_pt_));
      double eps_Gaussian = 0.0;
      double imeps_Gaussian = 0.0;
      double eps_Lorentzian = 0.0;
      double imeps_Lorentzian = 0.0;

      for (Index i_exc = minexc_; i_exc <= maxexc_; i_exc++) {
        double exc_lambda = nmtoev(BSESingletEnergies(i_exc) + shiftby_);
        eps_Gaussian += osc[i_exc] * Gaussian(lambda, exc_lambda, fwhm_);
        imeps_Gaussian +=
            osc[i_exc] * exc_lambda * Gaussian(lambda, exc_lambda, fwhm_);
        eps_Lorentzian += osc[i_exc] * Lorentzian(lambda, exc_lambda, fwhm_);
        imeps_Lorentzian +=
            osc[i_exc] * exc_lambda * Lorentzian(lambda, exc_lambda, fwhm_);
      }

      ofs << lambda << "    " << eps_Gaussian << "   " << imeps_Gaussian
          << "   " << eps_Lorentzian << "   " << imeps_Lorentzian << std::endl;
    }
    XTP_LOG(Log::error, log_)
        << " Spectrum in wavelength range from  " << lower_ << " to " << upper_
        << " nm and with broadening of FWHM " << fwhm_
        << " nm written to file  " << output_file_ << std::flush;
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
