/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

namespace votca {
namespace csg {
using namespace votca::tools;
using namespace std;

/**
    \brief base calss for all interactions

    This is the base class for all interactions.

    \todo double names/groups right, add molecules!!
*/
class Interaction {
public:
  Interaction() : _index(-1), _group(""), _group_id(-1), _name(""), _mol(-1){};

  virtual ~Interaction() {}
  virtual double EvaluateVar(const Topology &top) = 0;

  string getName() const { return _name; }

  void setGroup(const string &group) {
    _group = group;
    RebuildName();
  }
  const string &getGroup() const {
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

  virtual vec Grad(const Topology &top, int bead) = 0;
  int BeadCount() { return _beads.size(); }
  int getBeadId(int bead) {
    assert(bead > -1 && boost::lexical_cast<size_t>(bead) < _beads.size());
    return _beads[bead];
  }

protected:
  int _index;
  string _group;
  int _group_id;
  string _name;
  int _mol;
  vector<int> _beads;

  void RebuildName();
};

inline void Interaction::RebuildName() {
  stringstream s;
  if (_mol != -1)
    s << "molecule " << _mol;
  if (!_group.empty()) {
    s << ":" << _group;
    if (_group_id != -1) {
      s << " " << _group_id;
    }
  }
  if (_index != -1)
    s << ":index " << _index;
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

  IBond(list<int> &beads) {
    assert(beads.size() >= 2);
    _beads.resize(2);
    for (int i = 0; i < 2; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }
  double EvaluateVar(const Topology &top);
  vec Grad(const Topology &top, int bead);

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
  IAngle(list<int> &beads) {
    assert(beads.size() >= 3);
    _beads.resize(3);
    for (int i = 0; i < 3; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top);
  vec Grad(const Topology &top, int bead);

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
  IDihedral(list<int> &beads) {
    assert(beads.size() >= 4);
    _beads.resize(4);
    for (int i = 0; i < 4; ++i) {
      _beads[i] = beads.front();
      beads.pop_front();
    }
  }

  double EvaluateVar(const Topology &top);
  vec Grad(const Topology &top, int bead);

private:
};

inline double IBond::EvaluateVar(const Topology &top) {
  return abs(top.getDist(_beads[0], _beads[1]));
}

inline vec IBond::Grad(const Topology &top, int bead) {
  vec r = top.getDist(_beads[0], _beads[1]);
  r.normalize();
  return (bead == 0) ? -r : r;
}

inline double IAngle::EvaluateVar(const Topology &top) {
  vec v1(top.getDist(_beads[1], _beads[0]));
  vec v2(top.getDist(_beads[1], _beads[2]));
  return acos(v1 * v2 / sqrt((v1 * v1) * (v2 * v2)));
}

inline vec IAngle::Grad(const Topology &top, int bead) {
  vec v1(top.getDist(_beads[1], _beads[0]));
  vec v2(top.getDist(_beads[1], _beads[2]));

  double acos_prime =
      1.0 /
      (sqrt(1 -
            (v1 * v2) * (v1 * v2) / (abs(v1) * abs(v2) * abs(v1) * abs(v2))));
  switch (bead) {
  case (0):
    return acos_prime *
           (-v2 / (abs(v1) * abs(v2)) +
            (v1 * v2) * v1 / (abs(v2) * abs(v1) * abs(v1) * abs(v1)));
    break;
  case (1):
    return acos_prime *
           ((v1 + v2) / (abs(v1) * abs(v2)) -
            (v1 * v2) * ((v2 * v2) * v1 + (v1 * v1) * v2) /
                (abs(v1) * abs(v1) * abs(v1) * abs(v2) * abs(v2) * abs(v2)));
    break;
  case (2):
    return acos_prime *
           (-v1 / (abs(v1) * abs(v2)) +
            (v1 * v2) * v2 / (abs(v1) * abs(v2) * abs(v2) * abs(v2)));
    break;
  }
  // should never reach this
  assert(false);
  return vec(0, 0, 0);
}

inline double IDihedral::EvaluateVar(const Topology &top) {
  vec v1(top.getDist(_beads[0], _beads[1]));
  vec v2(top.getDist(_beads[1], _beads[2]));
  vec v3(top.getDist(_beads[2], _beads[3]));
  vec n1, n2;
  n1 = v1 ^ v2; // calculate the normal vector
  n2 = v2 ^ v3; // calculate the normal vector
  double sign = (v1 * n2 < 0) ? -1 : 1;
  return sign * acos(n1 * n2 / sqrt((n1 * n1) * (n2 * n2)));
}

inline vec IDihedral::Grad(const Topology &top, int bead) {
  vec v1(top.getDist(_beads[0], _beads[1]));
  vec v2(top.getDist(_beads[1], _beads[2]));
  vec v3(top.getDist(_beads[2], _beads[3]));
  vec n1, n2;
  n1 = v1 ^ v2; // calculate the normal vector
  n2 = v2 ^ v3; // calculate the normal vector
  double sign = (v1 * n2 < 0) ? -1 : 1;
  vec returnvec;                             // vector to return
  double returnvec0, returnvec1, returnvec2; // components of the return vector
  vec e0(1, 0, 0); // unit vector pointing in x-direction
  vec e1(0, 1, 0); // unit vector pointing in y-direction
  vec e2(0, 0, 1); // unit vector pointing in z-direction

  double acos_prime =
      (-1.0 / (sqrt(1 -
                    (n1 * n2) * (n1 * n2) /
                        (abs(n1) * abs(n2) * abs(n1) * abs(n2))))) *
      sign;
  switch (bead) {
  case (0): { //
    returnvec0 = acos_prime * ((n2 * (v2 ^ e0)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n1 * (v2 ^ e0))) /
                                   (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
    returnvec1 = acos_prime * ((n2 * (v2 ^ e1)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n1 * (v2 ^ e1))) /
                                   (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
    returnvec2 = acos_prime * ((n2 * (v2 ^ e2)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n1 * (v2 ^ e2))) /
                                   (abs(n1) * abs(n1) * abs(n1) * abs(n2)));
    returnvec.setX(returnvec0);
    returnvec.setY(returnvec1);
    returnvec.setZ(returnvec2);
    return returnvec;
    break;
  }
  case (1): { //
    returnvec0 =
        acos_prime *
        ((n1 * (v3 ^ e0) + n2 * ((e0 ^ v1) + (e0 ^ v2))) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * ((e0 ^ v1) + (e0 ^ v2))) /
               (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * (v3 ^ e0)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec1 =
        acos_prime *
        ((n1 * (v3 ^ e1) + n2 * ((e1 ^ v1) + (e1 ^ v2))) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * ((e1 ^ v1) + (e1 ^ v2))) /
               (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * (v3 ^ e1)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec2 =
        acos_prime *
        ((n1 * (v3 ^ e2) + n2 * ((e2 ^ v1) + (e2 ^ v2))) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * ((e2 ^ v1) + (e2 ^ v2))) /
               (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * (v3 ^ e2)) / (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec.setX(returnvec0);
    returnvec.setY(returnvec1);
    returnvec.setZ(returnvec2);
    return returnvec;
    break;
  };
  case (2): { //
    returnvec0 =
        acos_prime *
        ((n1 * ((e0 ^ v2) + (e0 ^ v3)) + n2 * (v1 ^ e0)) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * (v1 ^ e0)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * ((e0 ^ v2) + (e0 ^ v3))) /
               (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec1 =
        acos_prime *
        ((n1 * ((e1 ^ v2) + (e1 ^ v3)) + n2 * (v1 ^ e1)) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * (v1 ^ e1)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * ((e1 ^ v2) + (e1 ^ v3))) /
               (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec2 =
        acos_prime *
        ((n1 * ((e2 ^ v2) + (e2 ^ v3)) + n2 * (v1 ^ e2)) / (abs(n1) * abs(n2)) -
         ((n1 * n2) *
          ((n1 * (v1 ^ e2)) / (abs(n1) * abs(n1) * abs(n1) * abs(n2)) +
           (n2 * ((e2 ^ v2) + (e2 ^ v3))) /
               (abs(n1) * abs(n2) * abs(n2) * abs(n2)))));
    returnvec.setX(returnvec0);
    returnvec.setY(returnvec1);
    returnvec.setZ(returnvec2);
    return returnvec;
    break;
  };
  case (3): { //
    returnvec0 = acos_prime * ((n1 * (v2 ^ e0)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n2 * (v2 ^ e0))) /
                                   (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
    returnvec1 = acos_prime * ((n1 * (v2 ^ e1)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n2 * (v2 ^ e1))) /
                                   (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
    returnvec2 = acos_prime * ((n1 * (v2 ^ e2)) / (abs(n1) * abs(n2)) -
                               ((n1 * n2) * (n2 * (v2 ^ e2))) /
                                   (abs(n1) * abs(n2) * abs(n2) * abs(n2)));
    returnvec.setX(returnvec0);
    returnvec.setY(returnvec1);
    returnvec.setZ(returnvec2);
    return returnvec;
    break;
  };
  }
  // should never reach this
  assert(false);
  return vec(0, 0, 0);
}
}
}

#endif // _VOTCA_CSG_INTERACTION_H
