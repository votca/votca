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

#ifndef _VOTCA_CSG_INTERACTION_H
#define _VOTCA_CSG_INTERACTION_H

#include "bead.h"
#include "topology.h"
#include <sstream>
#include <string>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

/**
    \brief base calss for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction {
 public:
  Interaction() : _index(-1), _group(""), _group_id(-1), _name(""), _mol(-1){};

  virtual ~Interaction() = default;
  virtual double EvaluateVar(const Topology &top) = 0;

  std::string getName() const { return _name; }

  void setGroup(const std::string &group) {
    _group = group;
    RebuildName();
  }
  const std::string &getGroup() const {
    assert(_group.compare("") != 0);
    return _group;
  }

  // the group id is set by topology, when interaction is added to it
  // \todo if the group name is changed later, group id should be updated by
  // topology
  int getGroupId() {
    assert(_group_id != -1);
    return _group_id;
  }
  void setGroupId(int id) { _group_id = id; }

  void setIndex(const int &index) {
    _index = index;
    RebuildName();
  }
  const int &getIndex() const {
    assert(_index != -1);
    return _index;
  }

  void setMolecule(const int &mol) {
    _mol = mol;
    RebuildName();
  }
  const int &getMolecule() const {
    assert(_mol != -1);
    return _mol;
  }

  virtual Eigen::Vector3d Grad(const Topology &top, int bead) = 0;
  int BeadCount() { return _beads.size(); }
  int getBeadId(int bead) {
    assert(bead > -1 && boost::lexical_cast<size_t>(bead) < _beads.size());
    return _beads[bead];
  }

 protected:
  int _index;
  std::string _group;
  int _group_id;
  std::string _name;
  int _mol;
  std::vector<int> _beads;

  void RebuildName();
};

inline void Interaction::RebuildName() {
  std::stringstream s;
  if (_mol != -1) s << "molecule " << _mol;
  if (!_group.empty()) {
    s << ":" << _group;
    if (_group_id != -1) {
      s << " " << _group_id;
    }
  }
  if (_index != -1) s << ":index " << _index;
  _name = s.str();
}

/**
    \brief bond interaction
*/
class IBond : public Interaction {
 public:
  IBond(int bead1, int bead2) {
    _beads.resize(2);
    _beads[0] = bead1;
    _beads[1] = bead2;
  }

  IBond(std::list<int> &beads) {
    assert(beads.size() >= 2);
    _beads.resize(2);
    for (int i = 0; i < 2; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }
  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, int bead) override;

 private:
};

/**
    \brief angle interaction
*/
class IAngle : public Interaction {
 public:
  IAngle(int bead1, int bead2, int bead3) {
    _beads.resize(3);
    _beads[0] = bead1;
    _beads[1] = bead2;
    _beads[2] = bead3;
  }
  IAngle(std::list<int> &beads) {
    assert(beads.size() >= 3);
    _beads.resize(3);
    for (int i = 0; i < 3; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, int bead) override;

 private:
};

/**
    \brief dihedral interaction
*/
class IDihedral : public Interaction {
 public:
  IDihedral(int bead1, int bead2, int bead3, int bead4) {
    _beads.resize(4);
    _beads[0] = bead1;
    _beads[1] = bead2;
    _beads[2] = bead3;
    _beads[3] = bead4;
  }
  IDihedral(std::list<int> &beads) {
    assert(beads.size() >= 4);
    _beads.resize(4);
    for (int i = 0; i < 4; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top) override;
  Eigen::Vector3d Grad(const Topology &top, int bead) override;

 private:
};

inline double IBond::EvaluateVar(const Topology &top) {
  return top.getDist(_beads[0], _beads[1]).norm();
}

inline Eigen::Vector3d IBond::Grad(const Topology &top, int bead) {
  Eigen::Vector3d r = top.getDist(_beads[0], _beads[1]);
  r.normalize();
  return (bead == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(const Topology &top) {
  Eigen::Vector3d v1(top.getDist(_beads[1], _beads[0]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));
  return std::acos(v1.dot(v2) / sqrt(v1.squaredNorm() * v2.squaredNorm()));
}

inline Eigen::Vector3d IAngle::Grad(const Topology &top, int bead) {
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
  Eigen::Vector3d v1(top.getDist(_beads[0], _beads[1]));
  Eigen::Vector3d v2(top.getDist(_beads[1], _beads[2]));
  Eigen::Vector3d v3(top.getDist(_beads[2], _beads[3]));
  Eigen::Vector3d n1 = v1.cross(v2);  // calculate the normal vector
  Eigen::Vector3d n2 = v2.cross(v3);  // calculate the normal vector
  double sign = (v1.dot(n2) < 0) ? -1 : 1;
  return sign *
         std::acos(n1.dot(n2) / sqrt(n1.squaredNorm() * n2.squaredNorm()));
}

inline Eigen::Vector3d IDihedral::Grad(const Topology &top, int bead) {
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
      for (int i = 0; i < 3; i++) {
        returnvec[i] = n2.dot(v2.cross(e.col(i))) / (n1.norm() * n2.norm()) -
                       n1.dot(n2) * n1.dot(v2.cross(e.col(i))) /
                           (n2.norm() * std::pow(n1.norm(), 3));
      }
      return acos_prime * returnvec;
      break;
    }
    case (1): {
      for (int i = 0; i < 3; i++) {
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
      for (int i = 0; i < 3; i++) {
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
      for (int i = 0; i < 3; i++) {
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

#endif  // _VOTCA_CSG_INTERACTION_H
