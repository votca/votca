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

    tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);


    XTP_LOG(Log::error, *_log) << " Running Sternheimer" << flush;
    XTP_LOG(Log::error, *_log) << " Started parsing input parameters" << flush;

    _options.omegain = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegain");
    _options.omegafin = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegafin");
    _options.step = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.step");
    _options.imshift = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.imshift");
    _options.resolution = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.resolution");
    _options.do_precalc_fxc = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.do_precalc_fxc");
    _options.calculation = options.ifExistsReturnElseReturnDefault<std::string>(
        key + ".sternheimer.calculation");
    _options.spatialgridtype =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.spatialgridtype");
    _options.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level");
    _options.quadrature_scheme =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.quadrature_scheme");
    _options.quadrature_order = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.quadrature_order");
    _options.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level");

    XTP_LOG(Log::error, *_log) << " Task:" <<  _opt.calculation << flush;

    XTP_LOG(Log::error, *_log)
        << " Omega initial: " << _opt.omegain << flush;
    XTP_LOG(Log::error, *_log) << " Omega final: " << _opt.omegafin << flush;
    XTP_LOG(Log::error, *_log) << " Step: " << _opt.step << flush;
    XTP_LOG(Log::error, *_log)
        << " Imaginary shift: " << _opt.imshift << flush;
    XTP_LOG(Log::error, *_log)
        << " Resolution: " << _opt.resolution << flush;
    XTP_LOG(Log::error, *_log)
        << " GW-Sternheimer level: " << _opt.level << flush;
    XTP_LOG(Log::error, *_log)
        << " Calculation: " << _opt.calculation << flush;
    XTP_LOG(Log::error, *_log)
        << " Quadrature: " << _opt.quadrature_scheme
        << " Order: " << _opt.quadrature_order << flush;
};
bool SternheimerFrame::Evaluate() {

   OPENMP::setMaxThreads(_nThreads);

  _log.setReportLevel(Log::current_level);

  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (_do_guess) {
    XTP_LOG(Log::error, _log)
        << "Reading guess from " << _guess_file << std::flush;
    orbitals.ReadFromCpt(_guess_file);
  } else {
    XTP_LOG(Log::error, _log)
        << "Reading structure from " << _xyzfile << std::flush;
    orbitals.QMAtoms().LoadFromFile(_xyzfile);
  }


  Sternheimer sternheimer(orbitals, _log);

  sternheimer.setUpMatrices();

  XTP_LOG(Log::error, *_log)
      << TimeStamp() << " Started Sternheimer " << flush;
  if (_opt.calculation == "polarizability") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer Polarizability" << flush;
    sternheimer.configurate(opt);
    std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();
    sternheimer.printIsotropicAverage(polar);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Polarizability" << flush;
  }
  if (_opt.calculation == "gradient") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer Energy Gradient" << flush;
    sternheimer.configurate(opt);
    std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();
    sternheimer.printHellmannFeynmanForces(EPC);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Energy Gradient" << flush;
  }
  if (_opt.calculation == "mogradient") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer MO Energy Gradient" << flush;
    sternheimer.configurate(opt);
    for (Index n = 0; n < _orbitals.MOs().eigenvalues().size(); ++n) {

      std::vector<Eigen::Vector3cd> EPC = sternheimer.MOEnergyGradient(n, n);
      sternheimer.printMOEnergyGradient(EPC, n, n);
    }
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer MO Energy Gradient" << flush;
  }
  if (_opt.calculation == "koopmanalpha") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp()
        << " Started Sternheimer Koopman's relaxation coefficients" << flush;
    sternheimer.configurate(opt);
    for (Index n = 0; n < _orbitals.MOs().eigenvalues().size(); ++n) {

      std::complex<double> alpha = sternheimer.KoopmanRelaxationCoeff(n, 2.0);
      sternheimer.printKoopmanRelaxationCoeff(alpha, n);
    }
    XTP_LOG(Log::error, *_log)
        << TimeStamp()
        << " Finished Sternheimer Koopman's relaxation coefficients" << flush;
  }

  if (_opt.calculation == "koopman") {
    XTP_LOG(Log::error, *_log)
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
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Koopman's compliant" << flush;
  }

  if (_opt.calculation == "gwsternheimer") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer GW Hey Ho" << flush;
    sternheimer.configurate(opt);
    sternheimer.printGW(opt.level);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer GW" << flush;
  }

  XTP_LOG(Log::error, *_log)
      << TimeStamp() << " Finished Sternheimer" << flush;

  return true;
}

}  // namespace xtp
}  // namespace votca