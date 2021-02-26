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

#ifndef VOTCA_CSG_MAP_H
#define VOTCA_CSG_MAP_H
#pragma once

// Standard includes
#include <memory>
#include <vector>

// VOTCA includes
#include <votca/tools/eigen.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "boundarycondition.h"
#include "molecule.h"

namespace votca {
namespace csg {

enum class BeadMapType { Spherical, Ellipsoidal };

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap {
 public:
  BeadMap() = default;
  virtual ~BeadMap() = default;
  virtual void Apply(const BoundaryCondition &) = 0;
  virtual void Initialize(const Molecule *in, Bead *out,
                          tools::Property *opts_bead,
                          tools::Property *opts_map) = 0;

 protected:
  const Molecule *_in;
  Bead *_out;
  tools::Property *_opts_map;
  tools::Property *_opts_bead;
};

/*******************************************************
    Mapper class, collection of maps
*******************************************************/
class Map {
 public:
  Map(const Molecule &in, Molecule &out) : _in(in), _out(out) {}
  // Move constructor
  Map(Map &&map);
  // Move assignment
  Map &operator=(Map &&map);
  BeadMap *CreateBeadMap(const BeadMapType type);

  // void AddBeadMap(BeadMap *bmap) { _maps.push_back(bmap); }

  void Apply(const BoundaryCondition &bc);

 protected:
  Molecule _in;
  Molecule _out;
  std::vector<std::unique_ptr<BeadMap>> _maps;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_MAP_H
