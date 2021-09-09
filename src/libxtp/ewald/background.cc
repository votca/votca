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

#include "background.h"
#include "votca/xtp/tictoc.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace votca {
namespace xtp {

Background::Background(Logger& log, const Eigen::Matrix3d& uc_matrix,
                       const EwaldOptions options,
                       std::vector<PolarSegment>& polar_background)
    : log_(log), unit_cell_(uc_matrix), options_(options) {
  // Convert site data to a cartesian representation
  for (const PolarSegment& pseg : polar_background) {
    EwdSegment eseg(pseg);
    ewald_background_.push_back(eseg);
  }
  // Place atoms in the simulation box
  for (EwdSegment& seg : ewald_background_) {
    for (EwdSite& site : seg) {
      site.updatePos(unit_cell_.placeCoordInBox(site.getPos()));
    }
    // Should be deleted when merged and after compare with ctp is done
    seg.calcPos();
  }
}

void Background::Polarize() {

  TicToc timer;

  Index systemSize = computeSystemSize(ewald_background_);

  XTP_LOG(Log::error, log_) << unit_cell_ << std::endl;
  XTP_LOG(Log::error, log_) << "Setup real and reciprocal space" << std::endl;
  RSpace rspace(options_, unit_cell_, ewald_background_, log_);
  KSpace kspace(options_, unit_cell_, ewald_background_, log_);

  XTP_LOG(Log::error, log_)
      << "Compute real space permanent fields" << std::endl;
  rspace.computeStaticField();
  XTP_LOG(Log::error, log_)
      << "Compute reciprocal space permanent fields" << std::endl;
  kspace.computeStaticField();
  kspace.computeShapeField();
  kspace.computeIntraMolecularCorrection();

  std::ofstream infile2;
  infile2.open("staticFieldXTP.txt");
  infile2 << "id x y z q Ex Ey Ez" << std::endl;
  for (const auto& seg : ewald_background_) {
    for (const auto& site : seg) {
      infile2 << site << std::endl;
    }
  }

  /*******************************************************
   * We turn the problem into a linear problem (A + B)x = b
   * b = the permanent field
   * x = are the induced dipoles
   * A = the "interaction" tensor that turns x into corresponding induced field
   *     at that location
   * B = The inverse polarization matrix (a super matrix containing on it's
   *     diagonal the 3x3 inverse polarization matrices)
   * *****************************************************/

  // Get static field from the sites and convert it into a big 1D vector
  // The  same for the initial guess

  XTP_LOG(Log::error, log_) << "Compute the initial guess" << std::endl;
  Eigen::VectorXd staticField = Eigen::VectorXd::Zero(systemSize);
  Eigen::VectorXd initialGuess = Eigen::VectorXd::Zero(systemSize);
  Index index = 0;
  for (auto& seg : ewald_background_) {
    for (auto& site : seg) {
      site.induceDirect();  // compute induced dipole based on static field
      Eigen::Vector3d E = site.getStaticField();
      Eigen::Vector3d induced_dipole = site.getInducedDipole();
      staticField.segment<3>(index) = E;
      initialGuess.segment<3>(index) = induced_dipole;
      index += 3;
    }
  }

  XTP_LOG(Log::error, log_) << "Done with static fields, elapsed time: "
                            << timer.elapsedTimeAsString() << std::endl;

  XTP_LOG(Log::error, log_) << "Create dipole interaction matrix" << std::endl;

  // Set up the dipole interaction matrix
  Eigen::MatrixXd inducedDipoleInteraction(systemSize, systemSize);
  inducedDipoleInteraction.fill(0);
  rspace.addInducedDipoleInteractionTo(inducedDipoleInteraction);
  kspace.addInducedDipoleInteractionTo(inducedDipoleInteraction);
  kspace.addShapeCorrectionTo(inducedDipoleInteraction);
  kspace.addSICorrectionTo(inducedDipoleInteraction);

  XTP_LOG(Log::error, log_) << "Add inverse polarization matrix" << std::endl;
  // Add  the inverse polarization
  Index diagIndex = 0;
  for (auto& seg : ewald_background_) {
    for (auto& site : seg) {
      inducedDipoleInteraction.block<3, 3>(diagIndex, diagIndex) +=
          site.getPolarizationMatrix().inverse();
      diagIndex += 3;
    }
  }

  Eigen::VectorXd inducedField = inducedDipoleInteraction * initialGuess;

  XTP_LOG(Log::error, log_) << "Done setting up matrices, elapsed time: "
                            << timer.elapsedTimeAsString() << std::endl;
  XTP_LOG(Log::error, log_)
      << "Starting the conjugate gradient solver" << std::endl;

  // Solving the linear system
  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setMaxIterations(100);
  cg.setTolerance(1e-3);
  cg.compute(inducedDipoleInteraction);
  Eigen::VectorXd x = cg.solveWithGuess(staticField, initialGuess);

  // Give some solver output to the user
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " CG: #iterations: " << cg.iterations()
      << ", estimated error: " << cg.error() << std::endl;
  if (cg.info() == Eigen::ComputationInfo::NoConvergence) {
    XTP_LOG(Log::error, log_) << "PCG iterations did not converge" << std::endl;
  }
  if (cg.info() == Eigen::ComputationInfo::NumericalIssue) {
    XTP_LOG(Log::error, log_) << "PCG had a numerical issue" << std::endl;
  }
  XTP_LOG(Log::error, log_)
      << "Done solving, elapsed time: " << timer.elapsedTimeAsString()
      << std::endl;

  // TODO: Write the results x back to the ewald sites
}




Index Background::computeSystemSize(
    std::vector<EwdSegment>& segments) const {
  Index systemSize = 0;
  for (const auto& seg : segments) {
    systemSize += 3 * seg.size();
  }
  return systemSize;
}

}  // namespace xtp
}  // namespace votca
