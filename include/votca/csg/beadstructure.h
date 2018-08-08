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

#ifndef _VOTCA_CSG_BEADSTRUCTURE_H
#define	_VOTCA_CSG_BEADSTRUCTURE_H

#include <string>
#include <set>
#include <map>
#include <list>
#include <memory>

#include <votca/csg/basebead.h>

#include <votca/tools/graph.h>

namespace votca { namespace csg {

class BaseBead;
/**
 * \brief Designed to determine if the structure beads passed in
 *
 * Essentially it will have the functionality to determine if the stored beads
 * make up a single molecule. It can also break the stored beads up into
 * molecules. It can compare two bead structures and determine if they are the 
 * same structure. 
 **/

class BeadStructure
{
public:
  BeadStructure() : structureIdUpToDate(false), graphUpToDate(false) {};
  ~BeadStructure() {}

  // This assumes that a bead is never composed of more than a single molecule
  bool isSingleMolecule();

  // Follows same method name as topology class
  int BeadCount() { return beads_.size(); }
  void AddBead(BaseBead * bead); 
  BaseBead * getBead(int id);
  void ConnectBeads(int bead1_id, int bead2_id );

  std::vector<BaseBead *> getNeighBeads(int index);
 
  std::vector<std::shared_ptr<BeadStructure>> breakIntoMolecules();

  bool isStructureEquivalent(BeadStructure &beadstructure); 
private:

  void InitializeGraph_();
  void CalculateStructure_();
  
  bool structureIdUpToDate;
  bool graphUpToDate;
  shared_ptr<votca::tools::Graph> graph_;
  std::set<Edge> connections_;
  std::map<int,BaseBead *> beads_;    
  std::map<int,std::shared_ptr<votca::tools::GraphNode>> graphnodes_;
};

}}

#endif	// _VOTCA_CSG_BEADSTRUCTURE_H 

