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
#include "votca/xtp/gwbseengine.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/staticregion.h"


// Local private VOTCA includes
#include "sternheimerframe.h"

using std::flush;

namespace votca {
namespace xtp {

void SternheimerFrame::Initialize(const tools::Property &user_options) {

    std::string key = "sternheimer";

    tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

    _xyzfile = options.ifExistsReturnElseReturnDefault<std::string>(
      ".molecule", _job_name + ".xyz");


    XTP_LOG(Log::error, *_log) << " Running Sternheimer" << flush;
    XTP_LOG(Log::error, *_log) << " Started parsing input parameters" << flush;

    _options.start_frequency_grid = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegain",_options.start_frequency_grid);
    _options.end_frequency_grid = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.omegafin",_options.end_frequency_grid);
    _options.number_of_frequency_grid_points = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.step",_options.number_of_frequency_grid_points);
    _options.imaginary_shift_pade_approx = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.imshift",_options.imaginary_shift_pade_approx);
    _options.number_output_grid_points = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.resolution",_options.number_output_grid_points);
    _options.do_precalc_fxc = options.ifExistsReturnElseReturnDefault<double>(
        key + ".sternheimer.do_precalc_fxc",_options.do_precalc_fxc);
    _options.calculation = options.ifExistsReturnElseReturnDefault<std::string>(
        key + ".sternheimer.calculation",_options.calculation);
    _options.numerical_Integration_grid_type =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.spatialgridtype",_options.numerical_Integration_grid_type);
    _options.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level",_options.level);
    _options.quadrature_scheme =
        options.ifExistsReturnElseReturnDefault<std::string>(
            key + ".sternheimer.quadrature_scheme",_options.quadrature_scheme);
    _options.quadrature_order = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.quadrature_order",_options.quadrature_order);
    _options.level = options.ifExistsReturnElseReturnDefault<Index>(
        key + ".sternheimer.level",_options.level);

    XTP_LOG(Log::error, *_log) << " Task:" <<  _options.calculation << flush;

    XTP_LOG(Log::error, *_log)
        << " Omega initial: " << _options.start_frequency_grid << flush;
    XTP_LOG(Log::error, *_log) << " Omega final: " << _options.end_frequency_grid << flush;
    XTP_LOG(Log::error, *_log) << " Step: " << _options.number_of_frequency_grid_points << flush;
    XTP_LOG(Log::error, *_log)
        << " Imaginary shift: " << _options.imaginary_shift_pade_approx << flush;
    XTP_LOG(Log::error, *_log)
        << " Resolution: " << _options.number_output_grid_points << flush;
    XTP_LOG(Log::error, *_log)
        << " GW-Sternheimer level: " << _options.level << flush;
    XTP_LOG(Log::error, *_log)
        << " Calculation: " << _options.calculation << flush;
    XTP_LOG(Log::error, *_log)
        << " Quadrature: " << _options.quadrature_scheme
        << " Order: " << _options.quadrature_order << flush;
};
bool SternheimerFrame::Evaluate() {

   OPENMP::setMaxThreads(_nThreads);

  _log->setReportLevel(Log::current_level);

  _log->setMultithreading(true);
  _log->setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (_do_guess) {
    XTP_LOG(Log::error, *_log)
        << "Reading guess from " << _guess_file << std::flush;
    orbitals.ReadFromCpt(_guess_file);
  } else {
    XTP_LOG(Log::error, *_log)
        << "Reading structure from " << _xyzfile << std::flush;
    orbitals.QMAtoms().LoadFromFile(_xyzfile);
  }


  Sternheimer sternheimer(orbitals, _log);

  sternheimer.configurate(_options);

  sternheimer.setUpMatrices();

  XTP_LOG(Log::error, *_log)
      << TimeStamp() << " Started Sternheimer " << flush;
  if (_options.calculation == "polarizability") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer Polarizability" << flush;
    std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();
    sternheimer.printIsotropicAverage(polar);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Polarizability" << flush;
  }
  if (_options.calculation == "gradient") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer Energy Gradient" << flush;
    std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();
    sternheimer.printHellmannFeynmanForces(EPC);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Energy Gradient" << flush;
  }
  if (_options.calculation == "mogradient") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer MO Energy Gradient" << flush;
    for (Index n = 0; n < orbitals.MOs().eigenvalues().size(); ++n) {

      std::vector<Eigen::Vector3cd> EPC = sternheimer.MOEnergyGradient(n, n);
      sternheimer.printMOEnergyGradient(EPC, n, n);
    }
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer MO Energy Gradient" << flush;
  }
  if (_options.calculation == "koopmanalpha") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp()
        << " Started Sternheimer Koopman's relaxation coefficients" << flush;
    for (Index n = 0; n < orbitals.MOs().eigenvalues().size(); ++n) {

      std::complex<double> alpha = sternheimer.KoopmanRelaxationCoeff(n, 2.0);
      sternheimer.printKoopmanRelaxationCoeff(alpha, n);
    }
    XTP_LOG(Log::error, *_log)
        << TimeStamp()
        << " Finished Sternheimer Koopman's relaxation coefficients" << flush;
  }

  if (_options.calculation == "koopman") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer Koopman's compliant" << flush;
    for (Index n = 0; n < orbitals.MOs().eigenvalues().size(); ++n) {

      double f = 1.0;  // For occupied states
      if (n > orbitals.getHomo()) {
        f = 0.0;  // for unoccupied state
      }
      std::complex<double> alpha = sternheimer.KoopmanRelaxationCoeff(n, 1);
      std::complex<double> correction = sternheimer.KoopmanCorrection(n, f);
      std::cout << correction + orbitals.MOs().eigenvalues()(n) << std::endl;
      sternheimer.printKoopman(alpha, correction, n);
    }
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer Koopman's compliant" << flush;
  }

  if (_options.calculation == "gwsternheimer") {
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Started Sternheimer GW Hey Ho" << flush;
    sternheimer.printGW(_options.level);
    XTP_LOG(Log::error, *_log)
        << TimeStamp() << " Finished Sternheimer GW" << flush;
  }

  XTP_LOG(Log::error, *_log)
      << TimeStamp() << " Finished Sternheimer" << flush;

  return true;
}

}  // namespace xtp
}  // namespace votca