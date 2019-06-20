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

#include "molpol.h"
#include "votca/xtp/polarregion.h"

namespace votca {
namespace xtp {

void MolPol::Initialize(tools::Property& options) {
  std::string key = "options." + Identify();
  std::string mps_input =
      options.ifExistsReturnElseThrowRuntimeError<std::string>(key +
                                                               ".mpsinput");

  _input.LoadFromFile(mps_input);
  _mps_output = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".mpsoutput");
  _polar_options = options.get(key);

  Eigen::VectorXd target_vec =
      options.ifExistsReturnElseThrowRuntimeError<Eigen::VectorXd>(key +
                                                                   ".target");
  if (target_vec.size() != 6) {
    throw std::runtime_error(
        "ERROR <options.molpol.target> "
        " should have this format: pxx pxy pxz pyy pyz pzz");
  }
  target_vec *= std::pow(tools::conv::ang2bohr, 3);
  _polarisation_target(0, 0) = target_vec(0);
  _polarisation_target(1, 0) = target_vec(1);
  _polarisation_target(0, 1) = target_vec(1);
  _polarisation_target(2, 0) = target_vec(2);
  _polarisation_target(0, 2) = target_vec(2);
  _polarisation_target(1, 1) = target_vec(3);
  _polarisation_target(2, 1) = target_vec(4);
  _polarisation_target(1, 2) = target_vec(4);
  _polarisation_target(2, 2) = target_vec(5);

  Eigen::VectorXd default_weights = Eigen::VectorXd::Ones(_input.size());
  _weights = options.ifExistsReturnElseReturnDefault<Eigen::VectorXd>(
      key + ".weights", default_weights);

  _tolerance_pol = options.ifExistsReturnElseReturnDefault<double>(
      key + ".tolerance", _tolerance_pol);

  _max_iter = options.ifExistsReturnElseReturnDefault<int>(key + ".iterations",
                                                           _max_iter);
}

Eigen::Vector3d MolPol::Polarize(PolarRegion& pol,
                                 const Eigen::Vector3d& ext_field) const {
  pol.Reset();
  std::vector<std::unique_ptr<Region>> empty;
  PolarSegment& workmol = pol[0];
  for (PolarSite& site : workmol) {
    site.V().segment<3>(1) = ext_field;
  }
  pol.Evaluate(empty);
  Eigen::Vector3d induced_dipole = Eigen::Vector3d::Zero();
  for (const PolarSite& site : workmol) {
    induced_dipole += site.Induced_Dipole();
  }
  return induced_dipole;
}

Eigen::Matrix3d MolPol::CalcClassicalPol(const PolarSegment& input) const {
  Logger log;
  log.setMultithreading(false);
  log.setPreface(logINFO, (boost::format("\n ...")).str());
  log.setPreface(logERROR, (boost::format("\n ...")).str());
  log.setPreface(logWARNING, (boost::format("\n ...")).str());
  log.setPreface(logDEBUG, (boost::format("\n ...")).str());
  if (tools::globals::verbose) {
    log.setReportLevel(logDEBUG);
  }
  PolarRegion pol(0, log);
  pol.Initialize(_polar_options);
  pol.push_back(input);

  double eVnm_2_hrtbohr = tools::conv::ev2hrt / tools::conv::nm2bohr;
  Eigen::Matrix3d polarisation = Eigen::Matrix3d::Zero();
  Eigen::Vector3d zero = Polarize(pol, Eigen::Vector3d::Zero());
  Eigen::Vector3d ext_field = 0.1 * eVnm_2_hrtbohr * Eigen::Vector3d::UnitX();
  polarisation.col(0) = Polarize(pol, ext_field);
  ext_field = 0.1 * eVnm_2_hrtbohr * Eigen::Vector3d::UnitY();
  polarisation.col(1) = Polarize(pol, ext_field);
  ext_field = 0.1 * eVnm_2_hrtbohr * Eigen::Vector3d::UnitZ();
  polarisation.col(2) = Polarize(pol, ext_field);
  polarisation.colwise() -= zero;
  std::cout << log;
  return -polarisation / (0.1 * eVnm_2_hrtbohr);
}

void MolPol::PrintPolarisation(const Eigen::Matrix3d& result) const {
  std::cout << std::endl << "First principle polarisation [A^3]" << std::flush;
  double conversion = std::pow(tools::conv::bohr2ang, 3);
  std::cout << std::endl << _polarisation_target * conversion << std::flush;
  std::cout << std::endl << "Scaled classical polarisation [A^3]" << std::flush;
  std::cout << std::endl << result * conversion << std::flush;
  std::cout << std::endl
            << "EigenValues classical polarisation [A^3]" << std::flush;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es2;
  es2.computeDirect(result, Eigen::EigenvaluesOnly);
  Eigen::Matrix3d diag = es2.eigenvalues().asDiagonal();
  std::cout << std::endl << diag * conversion << std::flush;
}

bool MolPol::Evaluate() {

  PolarSegment polar = _input;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(_polarisation_target, Eigen::EigenvaluesOnly);
  const double pol_volume_target = std::pow(es.eigenvalues().prod(), 1.0 / 3.0);
  for (int iter = 0; iter < _max_iter; iter++) {

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es2;
    Eigen::Matrix3d pol = CalcClassicalPol(polar);
    es2.computeDirect(pol, Eigen::EigenvaluesOnly);
    const double pol_volume_iter =
        std::pow(es2.eigenvalues().prod(), 1.0 / 3.0);
    double scale = pol_volume_target / pol_volume_iter - 1;
    std::cout << "\nIteration " << iter + 1 << " of " << _max_iter
              << " target:" << pol_volume_target
              << " current:" << pol_volume_iter << std::endl;

    if (std::abs(scale) < _tolerance_pol) {
      std::cout << std::endl
                << "... ... Iterative refinement : *CONVERGED*" << std::flush;
      std::cout << std::endl
                << "... ... Scaling coefficient  : " << scale << std::flush;
      polar.WriteMPS(_mps_output, "MOLPOL (OPTIMIZED)");
      PrintPolarisation(pol);
      break;
    } else if (iter == (_max_iter - 1)) {
      std::cout << std::endl
                << "... ... Iterative refinement : *FAILED*" << std::flush;
      std::cout << std::endl
                << "... ... ERROR Convergence not achieved. "
                << "Check your input mps-file, target polarizability <target> "
                << std::flush;
      PrintPolarisation(pol);
    }

    for (int i = 0; i < polar.size(); i++) {
      PolarSite& site = polar[i];
      Eigen::Matrix3d local_pol = site.getPolarisation();
      site.setPolarisation(local_pol * std::pow(1 + scale * _weights[i], 2));
    }
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
