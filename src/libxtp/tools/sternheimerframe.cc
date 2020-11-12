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
#include "votca/xtp/padeapprox.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/staticregion.h"

// Local private VOTCA includes
#include "sternheimerframe.h"

using std::flush;

namespace votca {
namespace xtp {

void SternheimerFrame::ParseOptions(const tools::Property &user_options) {

  std::string key = "sternheimer";

  _log.setReportLevel(Log::current_level);

  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  tools::Property options =
      LoadDefaultsAndUpdateWithUserOptions("xtp", user_options);

  _orbfile = options.get(".orb_file").as<std::string>();

  XTP_LOG(Log::error, _log) << " Running Sternheimer" << flush;

  _options.start_frequency_grid = options.get(".start_frequency").as<double>();
  _options.end_frequency_grid = options.get(".end_frequency").as<double>();
  _options.number_of_frequency_grid_points = options.get(".steps").as<Index>();
  _options.imaginary_shift_pade_approx = options.get(".imshift").as<double>();
  _options.number_output_grid_points = options.get(".resolution").as<Index>();
  _options.do_precalc_fxc = options.get(".do_precalc_fxc").as<bool>();
  _options.task = options.get(".task").as<std::string>();
  _options.numerical_Integration_grid_type =
      options.get(".fxc_integration_grid").as<std::string>();
  _options.max_iterations_sc_sternheimer = options.get(".max_iter").as<Index>();
  _options.tolerance_sc_sternheimer = options.get(".tolerance").as<double>();
  _options.max_mixing_history = options.get(".max_hist").as<Index>();

  XTP_LOG(Log::error, _log) << " Task: " << _options.task << flush;

  XTP_LOG(Log::error, _log) << flush;

  if (_options.task == "polarizability") {
    XTP_LOG(Log::error, _log)
        << " Omega initial: " << _options.start_frequency_grid << flush;
    XTP_LOG(Log::error, _log)
        << " Omega final: " << _options.end_frequency_grid << flush;
    XTP_LOG(Log::error, _log)
        << " Steps: " << _options.number_of_frequency_grid_points << flush;
    XTP_LOG(Log::error, _log)
        << " Resolution: " << _options.number_output_grid_points << flush;
  }
  if (_options.task == "polarizability") {
    XTP_LOG(Log::error, _log)
        << " Imaginary shift: " << _options.imaginary_shift_pade_approx
        << flush;
  }
};
bool SternheimerFrame::Run() {

  OPENMP::setMaxThreads(_nThreads);

  // set logger

  _log.setReportLevel(Log::error);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Reading from orbitals from file: " << _orbfile
      << flush;

  // Get orbitals object
  Orbitals orbitals;

  orbitals.ReadFromCpt(_orbfile);

  XTP_LOG(Log::error, _log) << " Orbital data: " << flush;
  XTP_LOG(Log::error, _log)
      << " Basis size: " << orbitals.getBasisSetSize() << flush;
  XTP_LOG(Log::error, _log)
      << " XC Functional: " << orbitals.getXCFunctionalName() << flush;
  XTP_LOG(Log::error, _log)
      << " Basis name: " << orbitals.getDFTbasisName() << flush;

  Sternheimer sternheimer(orbitals, &_log);

  sternheimer.configurate(_options);

  sternheimer.setUpMatrices();

  std::string outfile = _options.task + ".dat";

  XTP_LOG(Log::error, _log) << " Output file: " << outfile << flush;

  std::ofstream ofs(outfile, std::ofstream::out);

  XTP_LOG(Log::error, _log) << TimeStamp() << " Started Sternheimer " << flush;

  if (_options.task == "polarizability") {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Started Sternheimer Polarizability" << flush;
    std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();
    std::vector<std::complex<double>> grid = sternheimer.BuildGrid(
        _options.start_frequency_grid, _options.end_frequency_grid,
        _options.number_output_grid_points, 0);
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Calculation complete" << flush;
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Writing output to " << outfile << flush;
    ofs << "#Freq (ev) \t polarizability_isotropic_average" << std::endl;
    for (Index i = 0; i < polar.size(); i++) {
      ofs << real(grid.at(i)) * votca::tools::conv::hrt2ev << "\t"
          << real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
                 real(polar.at(i)(0, 0)) / 3
          << std::endl;
    }

    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Finished Sternheimer Polarizability" << flush;
  }

  if (_options.task == "gradient") {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Started Sternheimer Energy Gradient" << flush;

    QMMolecule mol = orbitals.QMAtoms();

    std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();

    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Calculation complete" << flush;
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Writing output to " << outfile << flush;

    ofs << "\n"
        << "#Atom_Type "
        << "Atom_Index "
        << "Gradient x y z " << std::endl;
    for (int i = 0; i < EPC.size(); i++) {
      ofs << mol.at(i).getElement() << " " << i << " " << EPC[i][0].real()
          << " " << EPC[i][1].real() << " " << EPC[i][2].real() << std::endl;
    }
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Finished Sternheimer Energy Gradient" << flush;
  }

  if (_options.task == "mogradient") {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Started Sternheimer MO Energy Gradient" << flush;

    QMMolecule mol = orbitals.QMAtoms();

    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Calculation complete" << flush;
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Writing output to " << outfile << flush;

    ofs << "MO index "
        << "Atom_Type "
        << "Atom_Index "
        << "Gradient x y z " << std::endl;

    for (Index n = 0; n < orbitals.MOs().eigenvalues().size(); ++n) {

      std::vector<Eigen::Vector3cd> EPC = sternheimer.MOEnergyGradient(n, n);

      for (int i = 0; i < EPC.size(); i++) {
        ofs << n << " " << mol.at(i).getElement() << " " << i << " "
            << EPC[i][0].real() << " " << EPC[i][1].real() << " "
            << EPC[i][2].real() << std::endl;
      }
    }
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Finished Sternheimer MO Energy Gradient" << flush;
  }

  ofs.close();

  XTP_LOG(Log::error, _log) << TimeStamp() << " Finished Sternheimer" << flush;

  return true;
}

}  // namespace xtp
}  // namespace votca