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
#include "votca/xtp/polarregion.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace votca {
namespace xtp {

Background::Background(Logger& log, const Eigen::Matrix3d& uc_matrix,
                       const EwaldOptions options,
                       std::vector<PolarSegment>& polar_background)
    : log_(log),
      unit_cell_(uc_matrix),
      options_(options),
      rspace(log),
      kspace(log) {
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
  rspace.Initialize(options_, unit_cell_, ewald_background_);
  kspace.Initialize(options_, unit_cell_, ewald_background_);
}

void Background::Polarize() {

  TicToc timer;

  Index systemSize = computeSystemSize(ewald_background_);

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

  // Copy paste the solution back in to the sites
  index = 0;
  for (auto& seg : ewald_background_) {
    for (auto& site : seg) {
      site.setInducedDipole(x.segment<3>(index));
      index += 3;
    }
  }
}

void Background::writeToStateFile(std::string& state_file) {
  // Write the polarization state including all segments
  CheckpointFile cpf(state_file, CheckpointAccessLevel::CREATE);
  CheckpointWriter a = cpf.getWriter();
  CheckpointWriter ww = a.openChild("polar_background");
  for (EwdSegment& seg : ewald_background_) {
    CheckpointWriter www =
        ww.openChild("background_" + std::to_string(seg.getId()));
    seg.WriteToCpt(www);
  }

  // Write the ewald options
  CheckpointWriter w2 = a.openChild("ewald_options");
  w2(options_.alpha, "alpha");
  w2(options_.k_cutoff, "k_cutoff");
  w2(options_.r_cutoff, "r_cutoff");
  w2(options_.sharpness, "sharpness");
  switch (options_.shape) {
    case Shape::cube:
      w2(std::string("cube"), "shape");
      break;
    case Shape::xyslab:
      w2(std::string("xyslab"), "shape");
      break;
    case Shape::sphere:
      w2(std::string("sphere"), "shape");
      break;
    default:
      throw std::runtime_error("Shape not recognized!");
  }
  // Write the unit cell info
  CheckpointWriter w3 = a.openChild("unit_cell");
  w3(unit_cell_.getMatrix(), "unit_cell_matrix");
}

void Background::readFromStateFile(const std::string& state_file) {
  // Read sites and multipoles etc.
  CheckpointFile cpf(state_file, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("polar_background");
  for (Index i = 0; i < r.getNumDataSets(); ++i) {
    CheckpointReader rr = r.openChild("background_" + std::to_string(i));
    ewald_background_.push_back(EwdSegment(rr));
  }
  // ewald options
  CheckpointReader r2 = cpf.getReader("ewald_options");
  r2(options_.alpha, "alpha");
  r2(options_.k_cutoff, "k_cutoff");
  r2(options_.r_cutoff, "r_cutoff");
  r2(options_.sharpness, "sharpness");
  std::string shape;
  r2(shape, "shape");
  if (shape == "cube") {
    options_.shape = Shape::cube;
  } else if (shape == "xyslab") {
    options_.shape = Shape::xyslab;
  } else if (shape == "sphere") {
    options_.shape = Shape::sphere;
  } else {
    throw std::runtime_error("Unknown shape in state file.");
  }
  // unit cell info
  CheckpointReader r3 = cpf.getReader("unit_cell");
  Eigen::Matrix3d uc;
  r3(uc, "unit_cell_matrix");
  unit_cell_.reinitialize(uc);
}

Index Background::computeSystemSize(std::vector<EwdSegment>& segments) const {
  Index systemSize = 0;
  for (const auto& seg : segments) {
    systemSize += 3 * seg.size();
  }
  return systemSize;
}


void Background::ApplyBackgroundFields(
    std::vector<std::unique_ptr<votca::xtp::Region>>& regions,
    const std::vector<std::vector<SegId>>& region_seg_ids) {
  // Sanity check
  if (region_seg_ids.size() > 2) {
    throw std::runtime_error(
        "You requested a calculation with more than two inner regions. This is "
        "not allowed, the ewald background can only be used with a polar inner "
        "region or with a qm region inside a polar region.");
  }
  // Because we implement it step wise
  if (region_seg_ids.size() > 1) {
    throw std::runtime_error(
        "The qm region inside a polar region inside a ewald background is not "
        "yet implemented.");
  }

  // Apply background fields to sites in the polarization cloud
  if (region_seg_ids.size() == 1) {  // i.e. only a polariation cloud
    PolarRegion* pCloud = dynamic_cast<PolarRegion*>(regions[0].get());
    std::vector<SegId> pCloud_original_ids = region_seg_ids[0];
    for (Index i = 0; i < pCloud->size(); i++){
      // (*pCloud)[i] will be the ith segment in pCloud
      bgFieldAtSegment((*pCloud)[i], pCloud_original_ids);
    }
  }

}

void Background::bgFieldAtSegment(PolarSegment& seg, std::vector<SegId> pCloud_indices){ ; 
}

}  // namespace xtp
}  // namespace votca
