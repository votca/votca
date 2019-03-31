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

#include <math.h> /* ceil */
#include <votca/tools/constants.h>
#include <votca/xtp/grid.h>

namespace votca {
namespace xtp {
using namespace tools;

void Grid::printGridtoxyzfile(std::string filename) {
  // unit is Angstrom in xyz file
  std::ofstream points;
  points.open(filename.c_str(), std::ofstream::out);
  points << _gridpoints.size() << endl;
  points << endl;
  for (const auto& point : _gridpoints) {
    points << "X " << point.getX() * conv::bohr2ang << " "
           << point.getY() * conv::bohr2ang << " "
           << point.getZ() * conv::bohr2ang << endl;
  }
  points.close();
  return;
}

void Grid::setupgrid(std::vector<QMAtom*>& Atomlist) {

  Elements elements;

  double xmin = std::numeric_limits<double>::max();
  double ymin = xmin;
  double zmin = xmin;

  double xmax = std::numeric_limits<double>::min();
  double ymax = xmax;
  double zmax = xmax;
  double xtemp, ytemp, ztemp;

  for (const QMAtom* atom : Atomlist) {
    const tools::vec& pos = atom->getPos();
    xtemp = pos.getX();
    ytemp = pos.getY();
    ztemp = pos.getZ();
    if (xtemp < xmin) xmin = xtemp;
    if (xtemp > xmax) xmax = xtemp;
    if (ytemp < ymin) ymin = ytemp;
    if (ytemp > ymax) ymax = ytemp;
    if (ztemp < zmin) zmin = ztemp;
    if (ztemp > zmax) zmax = ztemp;
  }

  vec lowerbound = vec(xmin - _padding, ymin - _padding, zmin - _padding);
  vec upperbound = vec(xmax + _padding, ymax + _padding, zmax + _padding);
  vec steps = (upperbound - lowerbound) / _gridspacing;
  int xsteps = int(ceil(steps.getX()));
  int ysteps = int(ceil(steps.getY()));
  int zsteps = int(ceil(steps.getZ()));

  // needed to symmetrize grid around molecule
  double padding_x = (steps.getX() - xsteps) * _gridspacing * 0.5 + _padding;
  double padding_y = (steps.getY() - ysteps) * _gridspacing * 0.5 + _padding;
  double padding_z = (steps.getZ() - zsteps) * _gridspacing * 0.5 + _padding;

  for (int i = 0; i <= xsteps; i++) {
    double x = xmin - padding_x + i * _gridspacing;
    for (int j = 0; j <= ysteps; j++) {
      double y = ymin - padding_y + j * _gridspacing;
      for (int k = 0; k <= zsteps; k++) {
        double z = zmin - padding_z + k * _gridspacing;
        bool is_valid = false;
        vec gridpos = vec(x, y, z);
        for (const QMAtom* atom : Atomlist) {
          vec atompos = atom->getPos();
          double distance2 = (gridpos - atompos) * (gridpos - atompos);
          double atomcutoff =
              elements.getVdWChelpG(atom->getType()) * tools::conv::ang2bohr;
          if (distance2 < (atomcutoff * atomcutoff)) {
            is_valid = false;
            break;
          } else if (distance2 < (_cutoff * _cutoff))
            is_valid = true;
        }
        if (is_valid) {
          _gridpoints.push_back(gridpos);
        }
      }
    }
  }

  _gridvalues = Eigen::VectorXd::Zero(_gridpoints.size());
  return;
}

}  // namespace xtp
}  // namespace votca
