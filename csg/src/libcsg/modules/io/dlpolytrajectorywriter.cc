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

// Standard includes
#include <iomanip>
#include <string>
#include <filesystem>

// Local private VOTCA includes
#include "dlpolytrajectorywriter.h"

namespace votca {
namespace csg {

using namespace std;

void DLPOLYTrajectoryWriter::Open(string file, bool bAppend)
// open/create a dlpoly configuration or trajectory file
// NOTE: allowed file naming - <name>.dlpc or <name>.dlph (convention:
// ".dlpc"="CONFIG_CGV", ".dlph"="HISTORY_CGV")
{
  if (bAppend) {
    throw std::runtime_error(
        "Error: appending to dlpoly files not implemented");
  }

  std::filesystem::path filepath(file.c_str());
  string out_name = "HISTORY_CGV";

  if (!filepath.has_extension()) {

    throw std::ios_base::failure("Error on creating dlpoly file '" + file +
                                 "' - extension is expected, .dlph or .dlpc");

  } else if (filepath.extension() == ".dlpc") {

    isConfig_ = true;
    out_name = "CONFIG_CGV";

  } else if (filepath.extension() == ".dlph") {

    isConfig_ = false;

  } else {
    throw std::ios_base::failure("Error on creating dlpoly file '" + file +
                                 "' - wrong extension, use .dlph or .dlpc");
  }

  if (!filepath.has_stem()) {
    if (filepath.parent_path().string().size() == 0) {
      fname_ = out_name;
    } else {
      fname_ = filepath.parent_path().string() + "/" + out_name;
    }
  } else {
    fname_ = file;
  }

  fl_.open(fname_);
  if (!fl_.is_open()) {
    throw std::ios_base::failure("Error on creating dlpoly file '" + fname_ +
                                 "'");
  }
}

void DLPOLYTrajectoryWriter::Close() { fl_.close(); }

void DLPOLYTrajectoryWriter::Write(Topology *conf) {
  static Index nstep = 1;
  const double scale = 10.0;  // nm -> A factor
  Index mavecs = 0;
  Index mpbct = 0;

  if (conf->HasForce() && conf->HasVel()) {
    mavecs = 2;
  } else if (conf->HasVel()) {
    mavecs = 1;
  }

  if (conf->getBoxType() == BoundaryCondition::typeOrthorhombic) {
    mpbct = 2;
  }
  if (conf->getBoxType() == BoundaryCondition::typeTriclinic) {
    mpbct = 3;
  }

  if (isConfig_) {
    double energy = 0.0;
    fl_ << "From VOTCA with love" << endl;
    fl_ << setw(10) << mavecs << setw(10) << mpbct << setw(10)
        << conf->BeadCount() << setw(20) << energy << endl;
    Eigen::Matrix3d m = conf->getBox();
    for (Index i = 0; i < 3; i++) {
      fl_ << fixed << setprecision(10) << setw(20) << m(i, 0) * scale
          << setw(20) << m(i, 1) * scale << setw(20) << m(i, 2) * scale << endl;
    }

  } else {
    static double dstep = 0.0;
    if (nstep == 1) {
      fl_ << "From VOTCA with love" << endl;
      fl_ << setw(10) << mavecs << setw(10) << mpbct << setw(10)
          << conf->BeadCount() << endl;
      dstep = conf->getTime() / (double)(conf->getStep());
    }

    fl_ << "timestep" << setprecision(9) << setw(10) << conf->getStep()
        << setw(10) << conf->BeadCount() << setw(10) << mavecs << setw(10)
        << mpbct;
    fl_ << setprecision(9) << setw(12) << dstep << setw(12) << conf->getTime()
        << endl;

    Eigen::Matrix3d m = conf->getBox();
    for (Index i = 0; i < 3; i++) {
      fl_ << setprecision(12) << setw(20) << m(i, 0) * scale << setw(20)
          << m(i, 1) * scale << setw(20) << m(i, 2) * scale << endl;
    }
  }

  for (Index i = 0; i < conf->BeadCount(); i++) {
    Bead *bead = conf->getBead(i);

    // AB: DL_POLY needs bead TYPE, not name!

    if (isConfig_) {
      fl_ << setw(8) << left << bead->getType() << right << setw(10) << i + 1
          << endl;
    } else {
      fl_ << setw(8) << left << bead->getType() << right << setw(10) << i + 1;
      fl_ << setprecision(6) << setw(12) << bead->getMass() << setw(12)
          << bead->getQ() << setw(12) << "   0.0" << endl;
    }

    // alternative with atom NAME & fixed floating point format (in case the
    // need arises)
    // fl_ << setw(8) << left << bead->getName() << right << setw(10) << i+1;
    // fl_ << fixed << setprecision(6) << setw(12) << bead->getMass() <<
    // setw(12)
    //<< bead->getQ() << "   0.0" << endl;

    // nm -> Angs
    fl_ << resetiosflags(std::ios::fixed) << setprecision(12) << setw(20)
        << bead->getPos().x() * scale;
    fl_ << setw(20) << bead->getPos().y() * scale << setw(20)
        << bead->getPos().z() * scale << endl;

    if (mavecs > 0) {
      if (!bead->HasVel()) {
        throw std::ios_base::failure(
            "Error: dlpoly frame is supposed to contain velocities, but bead "
            "does not have v-data");
      }

      // nm -> Angs
      fl_ << setprecision(12) << setw(20) << bead->getVel().x() * scale
          << setw(20);
      fl_ << bead->getVel().y() * scale << setw(20)
          << bead->getVel().z() * scale << endl;

      if (mavecs > 1) {
        if (!bead->HasF()) {
          throw std::ios_base::failure(
              "Error: dlpoly frame is supposed to contain forces, but bead "
              "does not have f-data");
        }

        // nm -> Angs
        fl_ << setprecision(12) << setw(20) << bead->getF().x() * scale
            << setw(20);
        fl_ << bead->getF().y() * scale << setw(20) << bead->getF().z() * scale
            << endl;
      }
    }
  }
  nstep++;
}

}  // namespace csg
}  // namespace votca
