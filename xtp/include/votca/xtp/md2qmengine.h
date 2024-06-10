/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#ifndef VOTCA_XTP_MD2QMENGINE_H
#define VOTCA_XTP_MD2QMENGINE_H

// VOTCA includes
#include <votca/csg/topology.h>

// Local VOTCA includes
#include "logger.h"
#include "topology.h"

namespace votca {
namespace xtp {

class Md2QmEngine {
 public:
  Md2QmEngine(std::string mapfile) : mapfile_(mapfile) {};

  Topology map(const csg::Topology& top) const;

 private:
  void CheckMappingFile(tools::Property& topology_map) const;
  template <class T>
  bool SameValueForMultipleEntries(const std::vector<tools::Property*>& props,
                                   std::string tag) const;

  Index DetermineAtomNumOffset(const csg::Molecule* mol,
                               const std::vector<Index>& atom_ids_map) const;

  void MakeSegmentsWholePBC(Topology& top) const;
  bool CheckMolWhole(const Topology& top, const Segment& seg) const;

  std::string mapfile_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MD2QMENGINE_H
