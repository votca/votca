/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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
#include <sstream>
#include <string>

// Local VOTCA includes
#include "votca/csg/bead.h"
#include "votca/csg/interaction.h"
#include "votca/csg/topology.h"

namespace votca {
namespace csg {

using namespace std;

void Interaction::RebuildName() {
  std::stringstream s;
  if (_mol != -1) {
    s << "molecule " << _mol;
  }
  if (!_group.empty()) {
    s << ":" << _group;
    if (_group_id != -1) {
      s << " " << _group_id;
    }
  }
  if (_index != -1) {
    s << ":index " << _index;
  }
  _name = s.str();
}

double IBond::EvaluateVar(const Topology &top) const {
  return top.getDist(_beads[0], _beads[1]).norm();
}

Eigen::Vector3d IBond::Grad(const Topology &top, Index bead) const {
  Eigen::Vector3d r = top.getDist(_beads[0], _beads[1]);
  r.normalize();
  return (bead == 0) ? -r : r;
}

double IAngle::EvaluateVar(const Topology &top) const {
  Eigen::Vector3d v1(top.getDist(_beads[1], _beads[0]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));
  return std::acos(v1.dot(v2) / sqrt(v1.squaredNorm() * v2.squaredNorm()));
}

Eigen::Vector3d IAngle::Grad(const Topology &top, Index bead) const {
  Eigen::Vector3d v1(top.getDist(_beads[1], _beads[0]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));

  double acos_prime =
      1.0 / (sqrt(1 - std::pow(v1.dot(v2), 2) /
                          (v1.squaredNorm() * v2.squaredNorm())));
  switch (bead) {
    case (0):
      return acos_prime *
             (-v2 / (v1.norm() * v2.norm()) +
              (v1.dot(v2) * v1) / (v1.squaredNorm() * v2.squaredNorm()));
    case (1):
      return acos_prime *
             ((v1 + v2) / (v1.norm() * v2.norm()) -
              (v1.dot(v2)) * (v2.squaredNorm() * v1 + v1.squaredNorm() * v2) /
                  (std::pow(v1.norm(), 3) * std::pow(v2.norm(), 3)));
    case (2):
      return acos_prime * (-v1 / (v1.norm() * v2.norm())) +
             (v1.dot(v2) * v2 / (v1.norm() * std::pow(v2.norm(), 3)));
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}

double IDihedral::EvaluateVar(const Topology &top) const {
  Eigen::Vector3d v1(top.getDist(_beads[0], _beads[1]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));
  Eigen::Vector3d v3(top.getDist(_beads[2], _beads[3]));
  Eigen::Vector3d n1 = v1.cross(v2);  // calculate the normal vector
  Eigen::Vector3d n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  return sign *
         std::acos(n1.dot(n2) / sqrt(n1.squaredNorm() * n2.squaredNorm()));
}

Eigen::Vector3d IDihedral::Grad(const Topology &top, Index bead) const {
  Eigen::Vector3d v1(top.getDist(_beads[0], _beads[1]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));
  Eigen::Vector3d v3(top.getDist(_beads[2], _beads[3]));
  Eigen::Vector3d n1, n2;
  n1 = v1.cross(v2);  // calculate the normal vector
  n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  Eigen::Vector3d returnvec;

  Eigen::Matrix3d e = Eigen::Matrix3d::Identity();

  double acos_prime =
      sign * (-1.0 / (sqrt(1 - std::pow(n1.dot(n2), 2) /
                                   (n1.squaredNorm() * n2.squaredNorm()))));
  switch (bead) {
    case (0): {
      for (Index i = 0; i < 3; i++) {
        returnvec[i] = n2.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n1.dot(v2.cross(e.col(i))) /
                           (n2.norm() * std::pow(n1.norm(), 3));
      }
      return acos_prime * returnvec;
    }
    case (1): {
      for (Index i = 0; i < 3; i++) {
        returnvec[i] =
            (n1.dot(v3.cross(e.col(i))) +
             n2.dot(e.col(i).cross(v1) + e.col(i).cross(v2))) /
                (n1.norm() * n2.norm()) -
            n1.dot(n2) * ((n1.dot(e.col(i).cross(v1) + e.col(i).cross(v2))) /
                              (n2.norm() * std::pow(n1.norm(), 3)) +
                          n2.dot(v3.cross(e.col(i))) /
                              (n1.norm() * std::pow(n2.norm(), 3)));
      }
      return acos_prime * returnvec;
    };
    case (2): {
      for (Index i = 0; i < 3; i++) {
        returnvec[i] =
            (n1.dot(e.col(i).cross(v2) + e.col(i).cross(v3)) +
             n2.dot(v1.cross(e.col(i)))) /
                (n1.norm() * n2.norm()) -
            n1.dot(n2) * (n1.dot(v1.cross(e.col(i))) /
                              (n2.norm() * std::pow(n1.norm(), 3)) +
                          (n2.dot(e.col(i).cross(v2) + e.col(i).cross(v3))) /
                              (n1.norm() * std::pow(n2.norm(), 3)));
      }
      return acos_prime * returnvec;
    };
    case (3): {  //
      for (Index i = 0; i < 3; i++) {
        returnvec[i] = n1.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n2.dot(v2.cross(e.col(i))) /
                           (n1.norm() * std::pow(n2.norm(), 3));
      }
      return acos_prime * returnvec;
    };
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}

}  // namespace csg
}  // namespace votca
