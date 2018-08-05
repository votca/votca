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
    BeadStructure() {};
    ~BeadStructure() {}

    // This assumes that a bead is never composed of more than a single molecule
    bool isSingleMolecule();

    // Follows same method name as topology class
    int BeadCount() { return beads_.size(); }
    void AddBead(BaseBead * bead) { beads_[bead->getId()] = bead;}
    BaseBead * getBead(int id) { return beads_[id];}
    void ConnectBeads(int bead1_id, int bead2_id );
/*
    std::vector<BaseBead *> getNeighBeads(int index);
    BaseBead * getBead(int index);

    std::vector<BeadStructure *> breakIntoMolecules();

    bool isStructureEquivalent(const BeadGraph &beadgraph) const; 
*/
private:
  votca::tools::Graph graph_;
  std::map<int,BaseBead *> beads_;    
  std::map<int,std::shared_ptr<votca::tools::GraphNode>> graphnodes_;
};
/*
inline std::vector<bead *> BeadStructure::getNeighBeads(int index){
  auto neighbornodes = graph.getNeighNodes(index);
  std::vector<bead *> neighbeads;
  for( auto node_pair : neighbornodes ) {
    neighbeads.push_back(beads_[node_pair.first]); 
  }
  return neighbeads;
}

inline bool BeadStructure::isStructureEquivalent(const BeadGraph &beadgraph){
  return (graph==beadgraph.graph);
}
*/
inline bool BeadStructure::isSingleMolecule(){
  auto vertices = graph_.getVertices();
  auto isolated_nodes = graph_.getIsolatedNodes();
  if( beads_.size()==0 ) return false;
  if( isolated_nodes.size() !=0 ) return false;
  if( vertices.size()!= beads_.size() ) return false;
  return true; 
}

}}

#endif	// _VOTCA_CSG_BEADSTRUCTURE_H 

