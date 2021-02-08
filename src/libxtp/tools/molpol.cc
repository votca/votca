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

// Local VOTCA includes
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmpackage.h"
#include "votca/xtp/qmpackagefactory.h"

// Local private VOTCA includes
#include "molpol.h"

namespace votca {
namespace xtp {

void MolPol::ParseOptions(const tools::Property& options) {

  std::string mps_input = options.ifExistsReturnElseReturnDefault<std::string>(
      ".mpsinput", _job_name + ".mps");

  _input.LoadFromFile(mps_input);
  _mps_output = options.ifExistsReturnElseReturnDefault<std::string>(
      ".mpsoutput", _job_name + "_polar.mps");
  _polar_options = options.get(".options_polar");

  // polar targer or qmpackage logfile
  const std::string& mode = options.get("mode").as<std::string>();

  if (mode == "file") {
    Eigen::VectorXd target_vec =
        options.get(".target_polarisability").as<Eigen::VectorXd>();
    if (target_vec.size() != 6) {
      throw std::runtime_error(
          "ERROR <options.molpol.target> "
          " should have this format: pxx pxy pxz pyy pyz pzz");
    }
    target_vec *= std::pow(tools::conv::ang2bohr, 3);
    _polarization_target(0, 0) = target_vec(0);
    _polarization_target(1, 0) = target_vec(1);
    _polarization_target(0, 1) = target_vec(1);
    _polarization_target(2, 0) = target_vec(2);
    _polarization_target(0, 2) = target_vec(2);
    _polarization_target(1, 1) = target_vec(3);
    _polarization_target(2, 1) = target_vec(4);
    _polarization_target(1, 2) = target_vec(4);
    _polarization_target(2, 2) = target_vec(5);
  } else {
    std::string qm_package = options.get(".qmpackage").as<std::string>();
    std::string log_file = options.get(".logfile").as<std::string>();
    Logger log;
    log.setPreface(Log::info, "\n ...");
    log.setPreface(Log::error, "\n ...");
    log.setPreface(Log::warning, "\n ...");
    log.setPreface(Log::debug, "\n ...");
    log.setReportLevel(Log::current_level);
    log.setMultithreading(true);

    // Set-up QM package
    XTP_LOG(Log::error, log)
        << "Using package <" << qm_package << ">" << std::flush;
    QMPackageFactory factory;
    std::unique_ptr<QMPackage> qmpack = factory.Create(qm_package);
    qmpack->setLog(&log);
    qmpack->setRunDir(".");
    qmpack->setLogFileName(log_file);
    _polarization_target = qmpack->GetPolarizability();
  }

  Eigen::VectorXd default_weights = Eigen::VectorXd::Ones(_input.size());
  _weights = options.ifExistsReturnElseReturnDefault<Eigen::VectorXd>(
      ".weights", default_weights);

  _tolerance_pol = options.get(".tolerance").as<double>();
  _max_iter = options.get(".iterations").as<Index>();
}

Eigen::Vector3d MolPol::CalcInducedDipole(
    const PolarSegment& input, const Eigen::Vector3d& ext_field) const {
  Logger log;
  log.setMultithreading(false);
  log.setCommonPreface("\n ...");

  log.setReportLevel(Log::current_level);

  PolarRegion pol(0, log);
  pol.Initialize(_polar_options);
  pol.push_back(input);

  PolarSegment& workmol = pol[0];
  for (PolarSite& site : workmol) {
    site.V() = ext_field;
  }
  std::vector<std::unique_ptr<Region>> empty;  // pol interacts with nobody else
  pol.Evaluate(empty);
  Eigen::Vector3d induced_dipole = Eigen::Vector3d::Zero();
  for (const PolarSite& site : workmol) {
    induced_dipole += site.Induced_Dipole();
  }
  if (Log::current_level > 0) {
    std::cout << log;
  }
  return induced_dipole;
}

Eigen::Matrix3d MolPol::CalcClassicalPol(const PolarSegment& input) const {

  double eVnm_2_hrtbohr = tools::conv::ev2hrt / tools::conv::nm2bohr;
  double fieldstrength = (0.1 * eVnm_2_hrtbohr);
  Eigen::Matrix3d polarization = Eigen::Matrix3d::Zero();
  Eigen::Vector3d ext_field = fieldstrength * Eigen::Vector3d::UnitX();
  // central differences scheme
  Eigen::Vector3d dxplus = CalcInducedDipole(input, ext_field);
  Eigen::Vector3d dxminus = CalcInducedDipole(input, -ext_field);
  polarization.col(0) = dxplus - dxminus;
  ext_field = fieldstrength * Eigen::Vector3d::UnitY();
  Eigen::Vector3d dyplus = CalcInducedDipole(input, ext_field);
  Eigen::Vector3d dyminus = CalcInducedDipole(input, -ext_field);
  polarization.col(1) = dyplus - dyminus;
  ext_field = fieldstrength * Eigen::Vector3d::UnitZ();
  Eigen::Vector3d dzplus = CalcInducedDipole(input, ext_field);
  Eigen::Vector3d dzminus = CalcInducedDipole(input, -ext_field);
  polarization.col(2) = dzplus - dzminus;

  return -polarization / (2 * fieldstrength);
}

void MolPol::Printpolarization(const Eigen::Matrix3d& result) const {
  std::cout << std::endl << "First principle polarization [A^3]" << std::flush;
  double conversion = std::pow(tools::conv::bohr2ang, 3);
  std::cout << std::endl << _polarization_target * conversion << std::flush;
  std::cout << std::endl << "Scaled classical polarization [A^3]" << std::flush;
  std::cout << std::endl << result * conversion << std::flush;
  std::cout << std::endl
            << "EigenValues classical polarization [A^3]" << std::flush;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es2;
  es2.computeDirect(result, Eigen::EigenvaluesOnly);
  Eigen::Matrix3d diag = es2.eigenvalues().asDiagonal();
  std::cout << std::endl << diag * conversion << std::flush;
}

bool MolPol::Run() {

  PolarSegment polar = _input;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(_polarization_target, Eigen::EigenvaluesOnly);
  const double pol_volume_target = std::pow(es.eigenvalues().prod(), 1.0 / 3.0);
  for (Index iter = 0; iter < _max_iter; iter++) {

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es2;
    Eigen::Matrix3d pol = CalcClassicalPol(polar);
    es2.computeDirect(pol, Eigen::EigenvaluesOnly);
    const double pol_volume_iter =
        std::pow(es2.eigenvalues().prod(), 1.0 / 3.0);
    double scale = pol_volume_target / pol_volume_iter - 1;
    std::cout << "\nIteration " << iter + 1 << " of " << _max_iter
              << " target:" << pol_volume_target
              << " current:" << pol_volume_iter << std::endl;

    for (Index i = 0; i < polar.size(); i++) {
      Eigen::Matrix3d local_pol = polar[i].getpolarization();
      polar[i].setpolarization(local_pol * (1 + scale * _weights[i]));
    }

    if (std::abs(scale) < _tolerance_pol) {
      std::cout << std::endl
                << "... ... Iterative refinement : *CONVERGED*" << std::flush;
      std::cout << std::endl
                << "... ... Scaling coefficient  : " << scale << std::flush;
      polar.WriteMPS(_mps_output, "MOLPOL (OPTIMIZED)");
      Printpolarization(pol);
      break;
    } else if (iter == (_max_iter - 1)) {
      std::cout << std::endl
                << "... ... Iterative refinement : *FAILED*" << std::flush;
      std::cout << std::endl
                << "... ... ERROR Convergence not achieved. "
                << "Check your input mps-file, target polarizability <target> "
                << std::flush;
      Printpolarization(pol);
    }
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
