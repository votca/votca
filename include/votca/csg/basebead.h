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

#ifndef _VOTCA_CSG_BASEBEAD_H
#define _VOTCA_CSG_BASEBEAD_H

#include <votca/csg/topologyitem.h>

#include <votca/tools/identity.h>
#include <votca/tools/name.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class Molecule;
class BeadType;

/**
 * \brief information about a base bead
 *
 * The Base Bead class describes the core functionality of an atom or a coarse
 * grained bead. It stores information like the id, the name, the mass, the
 * charge and the residue it belongs to and the position
 *
 **/
class BaseBead : public TopologyItem, public Name, public Identity<int> {
public:
  /**
   * destructor
   */
  virtual ~BaseBead() {}

  /**
   * get the bead type
   * \return bead type object
   */
  virtual const BeadType *getType() const { return _type; }

  /**
   * set the bead type
   * \param bead type object
   */
  virtual void setType(BeadType *type) { _type = type; }

  /**
   * get the mass of the base bead
   * \return base bead mass
   */
  virtual const double &getMass() const { return _mass; }

  /**
   * get the charge of the base bead
   * \return base bead charge
   */
  virtual const double &getQ() const { return _q; }

  /**
   * set the mass of the base bead
   * \param m base bead mass
   */
  virtual void setMass(const double &m) { _mass = m; }

  /**
   * set the charge of the base bead
   * \param q base bead charge
   */
  virtual void setQ(const double &q) { _q = q; }

  /**
   * set the position of the base bead
   * \param r base bead position
   */
  virtual void setPos(const vec &r);

  /**
   * get the position of the base bead
   * \return base bead position
   */
  virtual const vec &getPos() const;

  /**
   * direct access (read/write) to the position of the base bead
   * \return reference to position
   */
  virtual vec &Pos() {
    assert(_bPos);
    return _pos;
  }

  /** does this configuration store positions? */
  bool HasPos() { return _bPos; }

  /**
   * molecule the base bead belongs to
   * \return Molecule object
   */
  virtual Molecule *getMolecule() { return _mol; }

  virtual void setMolecule(Molecule *mol);

protected:
  BaseBead()
      : _type(nullptr), _mol(nullptr), _mass(0.0), _q(0.0), _bPos(false),
        TopologyItem(nullptr){};

  BeadType *_type;
  Molecule *_mol;

  double _mass;
  double _q;
  vec _pos;

  bool _bPos;
};

inline void BaseBead::setMolecule(Molecule *mol) { _mol = mol; }

inline void BaseBead::setPos(const vec &r) {
  _bPos = true;
  _pos = r;
}

inline const vec &BaseBead::getPos() const {
  assert(_bPos);
  return _pos;
}
}
}

#endif // _VOTCA_CSG_BASEBEAD_H
