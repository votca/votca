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

// Standard includes
#include <cmath> /* ceil */

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/grid.h"

namespace votca {
namespace xtp {

void Grid::printGridtoxyzfile(std::string filename) {
  // unit is Angstrom in xyz file
  std::ofstream points;
  points.open(filename, std::ofstream::out);
  points << gridpoints_.size() << std::endl;
  points << std::endl;
  for (const auto& point : gridpoints_) {
    points << "X " << point.x() * tools::conv::bohr2ang << " "
           << point.y() * tools::conv::bohr2ang << " "
           << point.z() * tools::conv::bohr2ang << std::endl;
  }
  points.close();
  return;
}

void Grid::setupgrid(const QMMolecule& Atomlist) {

  tools::Elements elements;
  std::pair<Eigen::Vector3d, Eigen::Vector3d> extension =
      Atomlist.CalcSpatialMinMax();
  Eigen::Array3d min = extension.first.array();
  Eigen::Array3d max = extension.second.array();
  Eigen::Array3d doublesteps = (max - min + 2 * padding_) / gridspacing_;
  Eigen::Array<votca::Index, 3, 1> steps = (doublesteps.ceil()).cast<Index>();

  // needed to symmetrize grid around molecule
  Eigen::Array3d padding =
      (doublesteps - steps.cast<double>()) * gridspacing_ * 0.5 + padding_;
  Eigen::Array3d minpos = min - padding;
  for (Index i = 0; i <= steps.x(); i++) {
    double x = minpos.x() + double(i) * gridspacing_;
    for (Index j = 0; j <= steps.y(); j++) {
      double y = minpos.y() + double(j) * gridspacing_;
      for (Index k = 0; k <= steps.z(); k++) {
        double z = minpos.z() + double(k) * gridspacing_;
        bool is_valid = false;
        Eigen::Vector3d gridpos(x, y, z);
        for (const QMAtom& atom : Atomlist) {
          const Eigen::Vector3d& atompos = atom.getPos();
          double distance2 = (gridpos - atompos).squaredNorm();
          double atomcutoff =
              elements.getVdWChelpG(atom.getElement()) * tools::conv::ang2bohr;
          if (distance2 < (atomcutoff * atomcutoff)) {
            is_valid = false;
            break;
          } else if (distance2 < (cutoff_ * cutoff_)) {
            is_valid = true;
          }
        }
        if (is_valid) {
          gridpoints_.push_back(gridpos);
        }
      }
    }
  }

  gridvalues_ = Eigen::VectorXd::Zero(gridpoints_.size());
  return;
}

}  // namespace xtp
}  // namespace votca
