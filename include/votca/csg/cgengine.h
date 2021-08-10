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

#ifndef VOTCA_CSG_CGENGINE_H
#define VOTCA_CSG_CGENGINE_H

// Standard includes
#include <list>
#include <map>
#include <memory>

// Third party includes
#include <boost/program_options.hpp>

// VOTCA includes
#include <votca/tools/datacollection.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "cgengine.h"
#include "cgmoleculedef.h"
#include "cgobserver.h"
#include "molecule.h"
#include "nematicorder.h"
#include "topology.h"
#include "topologymap.h"
#include "topologyreader.h"
#include "trajectoryreader.h"
#include "trajectorywriter.h"

namespace votca {
namespace csg {

/**
    \brief coarse graining engine

    This class manages the coarse graining, at the moment it does the
   measurement stuff

    TODO: split this into an additional VotcaApplication object

*/
class CGEngine {
 public:
  CGEngine();

  /**
      create a coarse grained topolgy based on a given topology
  */
  std::unique_ptr<TopologyMap> CreateCGTopology(const Topology &in,
                                                Topology &out);

  /**
      load molecule type from file
  */
  void LoadMoleculeType(const std::string &filename);

  CGMoleculeDef *getMoleculeDef(const std::string &name);

  /**
   * \brief ignores molecule in mapping process
   * \param pattern glob pattern for molecule ident
   */
  void AddIgnore(const std::string &pattern) { ignores_.push_back(pattern); }

  /**
   * \brief checks whether molecule is ignored
   * \param ident identifyier of molecule
   * \return true if is ignored
   */
  bool IsIgnored(const std::string &ident);

 private:
  std::map<std::string, std::unique_ptr<CGMoleculeDef>> molecule_defs_;

  std::list<std::string> ignores_;
};

inline CGMoleculeDef *CGEngine::getMoleculeDef(const std::string &name) {
  std::map<std::string, std::unique_ptr<CGMoleculeDef>>::iterator iter;

  // if there is only 1 molecule definition, don't care about the name
  if (molecule_defs_.size() == 1 && name == "unnamed") {
    return (*(molecule_defs_.begin())).second.get();
  }

  iter = molecule_defs_.find(name);
  if (iter == molecule_defs_.end()) {
    return nullptr;
  }
  return (*iter).second.get();
}

inline bool CGEngine::IsIgnored(const std::string &ident) {
  for (auto &ignore_ : ignores_) {
    if (tools::wildcmp(ignore_, ident)) {
      return true;
    }
  }
  return false;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CGENGINE_H
