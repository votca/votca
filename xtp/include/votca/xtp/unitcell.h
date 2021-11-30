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

#pragma once
#ifndef VOTCA_XTP_UNITCELL_H
#define VOTCA_XTP_UNITCELL_H

#include <array>
#include <iomanip>
#include <vector>
// Local VOTCA includes
#include "bgsegment.h"

namespace votca {
namespace xtp {
class UnitCell {
 public:
  UnitCell(Eigen::Matrix3d _box) { initialize(_box); }

  UnitCell() : UnitCell(Eigen::Matrix3d::Identity()) { ; }

  ~UnitCell() = default;

  double getVolume() const { return cell_volume; }

  const Eigen::Matrix3d& getMatrix() const { return cell_matrix; }

  void reinitialize(Eigen::Matrix3d& mat) { initialize(mat); }

  void initialize(Eigen::Matrix3d& mat) {
    cell_matrix = mat;
    Eigen::Vector3d L1 = cell_matrix.col(0);
    Eigen::Vector3d L2 = cell_matrix.col(1);
    Eigen::Vector3d L3 = cell_matrix.col(2);

    bool is_gromacs_triclinic_box = L1[1] == 0 && L1[2] == 0 && L2[2] == 0 &&
                                    L1[0] > 0 && L2[1] > 0 && L3[2] > 0 &&
                                    L2[0] < 0.5 * L1[0] &&
                                    L3[0] < 0.5 * L1[0] && L3[1] < 0.5 * L2[1];
    if (!is_gromacs_triclinic_box) {
      std::cout << "a: " << L1.transpose() << std::endl;
      std::cout << "b: " << L2.transpose() << std::endl;
      std::cout << "c: " << L3.transpose() << std::endl;
      throw std::runtime_error(
          "The simulation box is not a triclinic box in the GROMACS format,\n"
          "the box vectors (a,b,c) should satisfy:\n"
          "a_y = a_z = b_z = 0\n"
          "a_x > 0, b_y > 0, c_z > 0\n"
          "b_x < 0.5 a_x, c_x < 0.5 a_x, c_y < 0.5 b_y\n");
    }
    cell_matrix_inv = cell_matrix.inverse();
    cell_volume = L1.dot(L2.cross(L3));
    centerOfBox = L1 * 0.5 + L2 * 0.5 + L3 * 0.5;
  }

  Eigen::Matrix3d getInverseMatrix() const {
    return 2 * boost::math::constants::pi<double>() * cell_matrix_inv;
  }

  Eigen::Vector3d getKVector(Index nx, Index ny, Index nz) const {
    Eigen::Vector3d n((double)nx, (double)ny, (double)nz);
    return getKVector(n);
  }

  Eigen::Vector3d getKVector(Eigen::Vector3d n) const {
    return 2 * boost::math::constants::pi<double>() * cell_matrix_inv * n;
  }

  Eigen::Vector3d getLVector(Index nx, Index ny, Index nz) const {
    Eigen::Vector3d n((double)nx, (double)ny, (double)nz);
    return getLVector(n);
  }

  Eigen::Vector3d getLVector(Eigen::Vector3d n) const {
    return cell_matrix * n;
  }

  Eigen::Vector3d minImage(const Eigen::Vector3d v1,
                           const Eigen::Vector3d v2) const {
    Eigen::Vector3d r_tp = v1 - v2;
    Eigen::Vector3d r_dp =
        r_tp - cell_matrix.col(2) * std::round(r_tp.z() / cell_matrix(2, 2));
    Eigen::Vector3d r_sp =
        r_dp - cell_matrix.col(1) * std::round(r_dp.y() / cell_matrix(1, 1));
    return r_sp - cell_matrix.col(0) * std::round(r_sp.x() / cell_matrix(0, 0));
  }

  Eigen::Vector3d placeCoordInBox(const Eigen::Vector3d pos) {
    return minImage(pos, centerOfBox) + centerOfBox;
  }

  Eigen::Vector3d minImage(const BGSegment seg1, const BGSegment seg2) const {
    return minImage(seg1.getPos(), seg2.getPos());
  }

  Eigen::Vector3d minImage(const BGSite site1, const BGSite site2) const {
    return minImage(site1.getPos(), site2.getPos());
  }

  friend std::ostream& operator<<(std::ostream& out, const UnitCell cell) {
    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
                                 ", ", " ", "[", "]", "", "");
    out << std::setprecision(8) << "RSpace Cell (nm)  : "
        << (tools::conv::bohr2nm * cell.getMatrix()).format(CommaInitFmt)
        << std::endl;
    out << std::setprecision(5) << "KSpace Cell (nm-1): "
        << ((1.0 / tools::conv::bohr2nm) * cell.getInverseMatrix())
               .format(CommaInitFmt);
    return out;
  }

  /**
   * Calculates maximum cutoff for which the minimal image convention still
   * makes sense.
   */
  double maxRCutOff() const {
    Eigen::Vector3d AxB = cell_matrix.col(0).cross(cell_matrix.col(1));
    Eigen::Vector3d BxC = cell_matrix.col(1).cross(cell_matrix.col(2));
    Eigen::Vector3d CxA = cell_matrix.col(2).cross(cell_matrix.col(0));

    double Wa = std::abs(cell_matrix.col(0).dot(BxC)) / BxC.norm();
    double Wb = std::abs(cell_matrix.col(1).dot(CxA)) / CxA.norm();
    double Wc = std::abs(cell_matrix.col(2).dot(AxB)) / AxB.norm();

    return 0.5 * std::min(std::min(Wa, Wb), Wc);
  }

  std::array<Index, 3> getNrOfRealSpaceCopiesForCutOff(double cutoff) {
    std::array<Index, 3> res;
    for (Index i = 0; i < 3; ++i) {
      res[i] = static_cast<Index>(
                   std::ceil(cutoff / cell_matrix.col(i).norm() - 0.5)) +
               1;
    }
    return res;
  }

 private:
  double cell_volume;
  Eigen::Vector3d centerOfBox;
  Eigen::Matrix3d cell_matrix;
  Eigen::Matrix3d cell_matrix_inv;

};  // namespace xtp
}  // namespace xtp
}  // namespace votca

#endif