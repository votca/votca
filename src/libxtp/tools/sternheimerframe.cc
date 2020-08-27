/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes

// Local private VOTCA includes
#include "sternheimerframe.h"

namespace votca {
namespace xtp {

void SternheimerFrame::Initialize(const tools::Property &user_options) {

  // sternheimer option
  if (options.exists(key + ".sternheimer")) {
    XTP_LOG(Log::error, *_pLog) << " Running Sternheimer" << flush;
    _gwopt.omegain = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegain", _gwopt.omegain);
    _gwopt.omegafin = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegafin", _gwopt.omegafin);
    _gwopt.step = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.step", _gwopt.step);
    _gwopt.imshift = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.imshift", _gwopt.imshift);
    _gwopt.resolution = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.resolution", _gwopt.resolution);
    _gwopt.do_precalc_fxc = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.do_precalc_fxc", _gwopt.do_precalc_fxc);
    _gwopt.calculation = options.ifExistsReturnElseReturnDefault<std::string>(
        key + ".sternheimer.calculation", _gwopt.calculation);
    _gwopt.spatialgridtype =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.spatialgridtype", _gwopt.spatialgridtype);
    _gwopt.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level", _gwopt.level);
    _gwopt.quadrature_scheme =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.quadrature_scheme", _gwopt.quadrature_scheme);
    _gwopt.quadrature_order = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.quadrature_order", _gwopt.quadrature_order);
    _gwopt.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level", _gwopt.level);
    XTP_LOG(Log::error, *_pLog)
        << " Omega initial: " << _gwopt.omegain << flush;
    XTP_LOG(Log::error, *_pLog) << " Omega final: " << _gwopt.omegafin << flush;
    XTP_LOG(Log::error, *_pLog) << " Step: " << _gwopt.step << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Imaginary shift: " << _gwopt.imshift << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Resolution: " << _gwopt.resolution << flush;
    XTP_LOG(Log::error, *_pLog)
        << " GW-Sternheimer level: " << _gwopt.level << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Calculation: " << _gwopt.calculation << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Quadrature: " << _gwopt.quadrature_scheme
        << " Order: " << _gwopt.quadrature_order << flush;
  }
};
bool SternheimerFrame::Evaluate() {

     if (_do_Sternheimer) {
    const double ev2hrt = 1 / votca::tools::conv::hrt2ev;

    Sternheimer sternheimer(_orbitals, _pLog);

    sternheimer.setUpMatrices();

    Sternheimer::options_sternheimer opt;
    opt.start_frequency_grid = _gwopt.omegain;
    opt.end_frequency_grid = _gwopt.omegafin;
    opt.number_of_frequency_grid_points = _gwopt.step;
    opt.imaginary_shift_pade_approx = _gwopt.imshift;
    opt.do_precalc_fxc = _gwopt.do_precalc_fxc;
    opt.number_output_grid_points = _gwopt.resolution;
    opt.numerical_Integration_grid_type = _gwopt.spatialgridtype;
    opt.level = _gwopt.level;
    opt.quadrature_order = _gwopt.quadrature_order;
    opt.quadrature_scheme = _gwopt.quadrature_scheme;
    // opt.level = _gwopt.level;

    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Started Sternheimer " << flush;
    if (_gwopt.calculation == "polarizability") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Started Sternheimer Polarizability" << flush;
      sternheimer.configurate(opt);
      std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();
      sternheimer.printIsotropicAverage(polar);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Finished Sternheimer Polarizability" << flush;
    }
    if (_gwopt.calculation == "gradient") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Started Sternheimer Energy Gradient" << flush;
      sternheimer.configurate(opt);
      std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();
      sternheimer.printHellmannFeynmanForces(EPC);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Finished Sternheimer Energy Gradient" << flush;
    }
    if (_gwopt.calculation == "mogradient") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Started Sternheimer MO Energy Gradient" << flush;
      sternheimer.configurate(opt);
      for (Index n = 0; n < _orbitals.MOs().eigenvalues().size(); ++n) {

        std::vector<Eigen::Vector3cd> EPC = sternheimer.MOEnergyGradient(n, n);
        sternheimer.printMOEnergyGradient(EPC, n, n);
      }
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Finished Sternheimer MO Energy Gradient" << flush;
    }
    if (_gwopt.calculation == "koopmanalpha") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp()
          << " Started Sternheimer Koopman's relaxation coefficients" << flush;
      sternheimer.configurate(opt);
      for (Index n = 0; n < _orbitals.MOs().eigenvalues().size(); ++n) {

        std::complex<double> alpha = sternheimer.KoopmanRelaxationCoeff(n, 2.0);
        sternheimer.printKoopmanRelaxationCoeff(alpha, n);
      }
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp()
          << " Finished Sternheimer Koopman's relaxation coefficients" << flush;
    }

    if (_gwopt.calculation == "koopman") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Started Sternheimer Koopman's compliant" << flush;
      sternheimer.configurate(opt);
      for (Index n = 0; n < _orbitals.MOs().eigenvalues().size(); ++n) {

        double f = 1.0;  // For occupied states
        if (n > _orbitals.getHomo()) {
          f = 0.0;  // for unoccupied state
        }
        std::complex<double> alpha = sternheimer.KoopmanRelaxationCoeff(n, 1);
        std::complex<double> correction = sternheimer.KoopmanCorrection(n, f);
        std::cout << correction + _orbitals.MOs().eigenvalues()(n) << std::endl;
        sternheimer.printKoopman(alpha, correction, n);
      }
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Finished Sternheimer Koopman's compliant"
          << flush;
    }

    if (_gwopt.calculation == "gwsternheimer") {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Started Sternheimer GW Hey Ho" << flush;
      sternheimer.configurate(opt);
      sternheimer.printGW(opt.level);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Finished Sternheimer GW" << flush;
    }

    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Finished Sternheimer" << flush;

  }


}

}  // namespace xtp
}  // namespace votca