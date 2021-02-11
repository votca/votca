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

#ifndef _VOTCA_CSG_TOPOLOGYMAP_H
#define _VOTCA_CSG_TOPOLOGYMAP_H

// Standard includes
#include <memory>
#include <vector>

// Local VOTCA includes
#include "map.h"
#include "topology.h"

namespace votca {
namespace csg {

class TopologyMap {
 public:
  TopologyMap(Topology *in, Topology *out);

  void AddMoleculeMap(std::unique_ptr<Map> map);

  void Apply();

 private:
  Topology *_in;
  Topology *_out;

  using MapContainer = std::vector<Map>;
  MapContainer _maps;
};

inline TopologyMap::TopologyMap(Topology *in, Topology *out)
    : _in(in), _out(out) {}

inline void TopologyMap::AddMoleculeMap(std::unique_ptr<Map> map) {
  _maps.push_back(std::move(map));
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_TOPOLOGYMAP_H */
