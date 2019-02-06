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

#ifndef _VOTCA_CSG_MAP_H
#define _VOTCA_CSG_MAP_H

#include "molecule.h"
#include <vector>
#include <votca/tools/property.h>
#include <votca/tools/vec.h>

namespace votca {
namespace csg {
using namespace votca::tools;

class BeadMap;
/*******************************************************
    Mapper class, collection of maps
*******************************************************/
class Map {
 public:
  Map(Molecule &in, Molecule &out) : _in(in), _out(out) {}
  ~Map();

  void AddBeadMap(BeadMap *bmap) { _maps.push_back(bmap); }

  void Apply();

 protected:
  Molecule _in, _out;
  std::vector<BeadMap *> _maps;
};

/*******************************************************
    Interface for all maps
*******************************************************/
class BeadMap {
 public:
  virtual ~BeadMap(){};
  virtual void Apply() = 0;
  virtual void Initialize(Molecule *in, Bead *out, Property *opts_map,
                          Property *opts_bead);

 protected:
  Molecule *_in;
  Bead *_out;
  Property *_opts_map;
  Property *_opts_bead;
};

inline void BeadMap::Initialize(Molecule *in, Bead *out, Property *opts_bead,
                                Property *opts_map) {
  _in = in;
  _out = out;
  _opts_map = opts_map;
  _opts_bead = opts_bead;
}

/*******************************************************
    Linear map for spherical beads
*******************************************************/
class Map_Sphere : public BeadMap {
 public:
  Map_Sphere() {}
  void Apply();

  void Initialize(Molecule *in, Bead *out, Property *opts_bead,
                  Property *opts_map);

 protected:
  void AddElem(Bead *in, double weight, double force_weight);

  struct element_t {
    Bead *_in;
    double _weight;
    double _force_weight;
  };
  std::vector<element_t> _matrix;
};

inline void Map_Sphere::AddElem(Bead *in, double weight, double force_weight) {
  element_t el;
  el._in = in;
  el._weight = weight;
  el._force_weight = force_weight;
  _matrix.push_back(el);
}

/*******************************************************
    Linear map for ellipsoidal bead
*******************************************************/
class Map_Ellipsoid : public Map_Sphere {
 public:
  Map_Ellipsoid() {}
  void Apply();

 protected:
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_MAP_H */
