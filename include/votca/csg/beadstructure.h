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
 * same structure. The Ids of a bead will not determine if a structure is
 * equivalent or not. Each bead must have a unique Id.
 *
 * E.g. If a beadstructure is composed of the following:
 *
 * BeadStructure 1
 * Name:"A" Id:1 ---- Name:"B" Id:2
 *
 * BeadStucture 2
 * Name:"A" Id:3 ---- Name:"B" Id:4
 *
 * Here BeadStructure 1 and 2 will compare equal
 *
 **/

class BeadStructure
{
public:
  BeadStructure() : structureIdUpToDate(false), graphUpToDate(false) {};
  ~BeadStructure() {}

  /**
   * \brief Determine if the bead structure consists of a single molecule
   *
   * This function will determine if all beads in the structure are connected
   * somehow to everyother bead. The connection does not have to be direct
   *
   * @return - returns a boolean true if it is a single molecule
   **/
  bool isSingleMolecule();

  /**
   * \brief returns the number of beads in the bead structure
   **/
  int BeadCount() { return beads_.size(); }

  /**
   * \brief add a bead to the bead structure
   *
   * The same bead cannot be added twice. 
   **/
  void AddBead(BaseBead * bead); 

  /**
   * \brief Get the bead with the specified id
   **/
  BaseBead * getBead(int id);

  /**
   * \brief Create a connection between two beads in the structure
   *
   * A bead cannot be connected to itself. It also may not be connected to a 
   * bead that has not yet been added to the structure. 
   **/
  void ConnectBeads(int bead1_id, int bead2_id );

  /**
   * \breif Return a vector of all the beads neighboring the index
   **/
  std::vector<BaseBead *> getNeighBeads(int index);

  /**
   * \breif Bread the beadstructure up into molecular units
   *
   * If a beadstructure is composed of several unconnected networks of beads. 
   * These structures will be broken up into their own bead structures and 
   * returned in a vector. 
   **/ 
  std::vector<std::shared_ptr<BeadStructure>> breakIntoMolecules();

  /**
   * \breif Compare the topology of two bead structures
   *
   * This function looks at how the beads are arranged within the bead structure
   * and determines if the topology is the same.
   *
   * @param[in] - beadstructure to compare with
   * @return - if the same returns true else false
   *
   **/
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

