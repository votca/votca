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

#include <memory>

#include <votca/csg/topologyitem.h>
#include <votca/csg/moleculeitem.h>

#include <votca/tools/identity.h>
#include <votca/tools/name.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BeadType;

/**
 * \brief information about a base bead
 *
 * The Base Bead class describes the core functionality of an atom or a coarse
 * grained bead. It stores information like the id, the name, the mass, the
 * charge and the residue it belongs to and the position
 *
 **/
class BaseBead : public TopologyItem,
                 public MoleculeItem,
                 public virtual Name,
                 public virtual Identity<int> {
public:
  /**
   * destructor
   */
  virtual ~BaseBead() {}

  /**
   * get the bead type
   * \return const bead type pointer
   */
  virtual const std::weak_ptr<BeadType> getType() const { return _type;}

  /**
   * set the bead type
   * \param bead type object
   */
  virtual void setType(std::weak_ptr<BeadType> type) { _type = type; }

  /**
   * get the bead type
   * \return - non constant bead type pointer
   */
  virtual std::weak_ptr<BeadType> Type() const { return _type; }

  /**
   * get the name of the bead type
   * \return - string indicates the name or throw error if bead type is not 
   * accesible.
   **/
  std::string getBeadTypeName();

  /**
   * get the id of the bead type
   * \return - int indicated the id or throw an error if bead type is not
   * accessible
   **/
  int getBeadTypeId();

  /**
   * get the mass of the base bead
   * \return - base bead mass
   */
  virtual const double &getMass() const { return _mass; }

  /**
   * set the mass of the base bead
   * \param - base bead mass
   */
  virtual void setMass(const double &m) { _mass = m; }

  /**
   * set the position of the base bead
   * \param - base bead position
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
  bool HasPos() const { return _bPos; }

  /** set has position to true */
  void HasPos(bool true_or_false) { _bPos = true_or_false; }

protected:
  BaseBead()
      : TopologyItem(nullptr), MoleculeItem(nullptr), 
      _mass(0.0), _bPos(false){};


  std::weak_ptr<BeadType> _type;

  double _mass;
  vec _pos;

  bool _bPos;
};

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
