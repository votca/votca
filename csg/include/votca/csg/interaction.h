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
#ifndef VOTCA_CSG_INTERACTION_H
#define VOTCA_CSG_INTERACTION_H

// Standard includes
#include <sstream>
#include <string>

// Local VOTCA includes
#include "bead.h"
#include "topology.h"

namespace votca {
namespace csg {

/**
    \brief base class for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction {
 public:
  Interaction() = default;

  virtual ~Interaction() = default;
  virtual double EvaluateVar(const Topology &top) = 0;

  std::string getName() const { return name_; }

  void setGroup(const std::string &group) {
    group_ = group;
    RebuildName();
  }
  const std::string &getGroup() const {
    assert(group_.compare("") != 0);
    return group_;
  }

  // the group id is set by topology, when interaction is added to it
  // \todo if the group name is changed later, group id should be updated by
  // topology
  Index getGroupId() {
    assert(group_id_ != -1);
    return group_id_;
  }
  void setGroupId(Index id) { group_id_ = id; }

  void setIndex(const Index &index) {
    index_ = index;
    RebuildName();
  }
  const Index &getIndex() const {
    assert(index_ != -1);
    return index_;
  }

  void setMolecule(const Index &mol) {
    mol_ = mol;
    RebuildName();
  }
  const Index &getMolecule() const {
    assert(mol_ != -1);
    return mol_;
  }

  virtual Eigen::Vector3d Grad(const Topology &top, Index bead) = 0;
  Index BeadCount() const { return beads_.size(); }
  Index getBeadId(Index bead) const {
    assert(bead > -1 && boost::lexical_cast<size_t>(bead) < beads_.size());
    return beads_[bead];
  }

 protected:
  Index index_ = -1;
  std::string group_ = "";
  Index group_id_ = -1;
  std::string name_ = "";
  Index mol_ = -1;
  std::vector<Index> beads_;

  void RebuildName();
};

inline void Interaction::RebuildName() {
  std::stringstream s;
  if (mol_ != -1) {
    {
      s << "molecule " << mol_;
    }
  }
  if (!group_.empty()) {
    s << ":" << group_;
    if (group_id_ != -1) {
      s << " " << group_id_;
    }
  }
  if (index_ != -1) {
    {
      s << ":index " << index_;
    }
  }
  name_ = s.str();
}

/**
    \brief bond interaction
*/
class IBond : public Interaction {
 public:
  IBond(Index bead1, Index bead2) {
    beads_.resize(2);
    beads_[0] = bead1;
    beads_[1] = bead2;
  }

  IBond(std::list<Index> &beads) {
    assert(beads.size() >= 2);
    beads_.resize(2);
    for (Index i = 0; i < 2; ++i) {
      beads_[i] = beads.front();
      beads.pop_front();
    }
  }
  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) override;

 private:
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  IAngle(Index bead1, Index bead2, Index bead3) {
    beads_.resize(3);
    beads_[0] = bead1;
    beads_[1] = bead2;
    beads_[2] = bead3;
  }
  IAngle(std::list<Index> &beads) {
    assert(beads.size() >= 3);
    beads_.resize(3);
    for (Index i = 0; i < 3; ++i) {
      beads_[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) override;

 private:
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  IDihedral(Index bead1, Index bead2, Index bead3, Index bead4) {
    beads_.resize(4);
    beads_[0] = bead1;
    beads_[1] = bead2;
    beads_[2] = bead3;
    beads_[3] = bead4;
  }
  IDihedral(std::list<Index> &beads) {
    assert(beads.size() >= 4);
    beads_.resize(4);
    for (Index i = 0; i < 4; ++i) {
      beads_[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, Index bead) override;

 private:
};

inline double IBond::EvaluateVar(const Topology &top) {
  return top.getDist(beads_[0], beads_[1]).norm();
}

inline Eigen::Vector3d IBond::Grad(const Topology &top, Index bead) {
  Eigen::Vector3d r = top.getDist(beads_[0], beads_[1]);
  r.normalize();
  return (bead == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(const Topology &top) {
  Eigen::Vector3d v1(top.getDist(beads_[1], beads_[0]));
  Eigen::Vector3d v2(top.getDist(beads_[1], beads_[2]));
  return std::acos(v1.dot(v2) / sqrt(v1.squaredNorm() * v2.squaredNorm()));
}

inline Eigen::Vector3d IAngle::Grad(const Topology &top, Index bead) {
  Eigen::Vector3d v1(top.getDist(beads_[1], beads_[0]));
  Eigen::Vector3d v2(top.getDist(beads_[1], beads_[2]));

  double acos_prime =
      1.0 / (sqrt(1 - std::pow(v1.dot(v2), 2) /
                          (v1.squaredNorm() * v2.squaredNorm())));
  switch (bead) {
    case (0):
      return acos_prime *
             (-v2 / (v1.norm() * v2.norm()) +
              (v1.dot(v2) * v1) / (v1.squaredNorm() * v2.squaredNorm()));
      break;
    case (1):
      return acos_prime *
             ((v1 + v2) / (v1.norm() * v2.norm()) -
              (v1.dot(v2)) * (v2.squaredNorm() * v1 + v1.squaredNorm() * v2) /
                  (std::pow(v1.norm(), 3) * std::pow(v2.norm(), 3)));
      break;
    case (2):
      return acos_prime * (-v1 / (v1.norm() * v2.norm())) +
             (v1.dot(v2) * v2 / (v1.norm() * std::pow(v2.norm(), 3)));
      break;
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}

inline double IDihedral::EvaluateVar(const Topology &top) {
  Eigen::Vector3d v1(top.getDist(beads_[0], beads_[1]));
  Eigen::Vector3d v2(top.getDist(beads_[1], beads_[2]));
  Eigen::Vector3d v3(top.getDist(beads_[2], beads_[3]));
  Eigen::Vector3d n1 = v1.cross(v2);  // calculate the normal vector
  Eigen::Vector3d n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  return sign *
         std::acos(n1.dot(n2) / sqrt(n1.squaredNorm() * n2.squaredNorm()));
}

inline Eigen::Vector3d IDihedral::Grad(const Topology &top, Index bead) {
  Eigen::Vector3d v1(top.getDist(beads_[0], beads_[1]));
  Eigen::Vector3d v2(top.getDist(beads_[1], beads_[2]));
  Eigen::Vector3d v3(top.getDist(beads_[2], beads_[3]));
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
      break;
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
      break;
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
      break;
    };
    case (3): {  //
      for (Index i = 0; i < 3; i++) {
        returnvec[i] = n1.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n2.dot(v2.cross(e.col(i))) /
                           (n1.norm() * std::pow(n2.norm(), 3));
      }
      return acos_prime * returnvec;
      break;
    };
  }
  // should never reach this
  assert(false);
  return Eigen::Vector3d::Zero();
}
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_INTERACTION_H
