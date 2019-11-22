/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#pragma once
#ifndef VOTCA_CSG_BOUNDARYCONDITION_H
#define VOTCA_CSG_BOUNDARYCONDITION_H

#include <votca/tools/eigen.h>

namespace votca {
namespace csg {

class BoundaryCondition {

 public:
  virtual ~BoundaryCondition() = default;

  /**
   * set the simulation box
   * \param box triclinic box matrix
   */
  void setBox(const Eigen::Matrix3d &box) { _box = box; };

  /**
   * get the simulation box
   * \return triclinic box matrix
   */
  const Eigen::Matrix3d &getBox() { return _box; };

  /**
   * get the volume of the box
   * \return box volume as double
   */
  virtual double BoxVolume();

  /**
   * get shortest connection vector between r_i and r_j with respect to the
   * (periodic) box \return shortest distance vector
   */
  virtual Eigen::Vector3d BCShortestConnection(
      const Eigen::Vector3d &r_i, const Eigen::Vector3d &r_j) const = 0;

  enum eBoxtype { typeAuto = 0, typeTriclinic, typeOrthorhombic, typeOpen };
  virtual eBoxtype getBoxType() = 0;

 protected:
  Eigen::Matrix3d _box;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BOUNDARYCONDITION_H
