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
#include <list>
#include "topology.h"

namespace votca { namespace csg {
using namespace votca::tools;

class BaseBead;
/**
 * \brief Designed to determine if the structure beads passed in
 *
 * Essentially it will have the functionality to determine if the stored beads
 * make up a single molecule. It can also break the stored beads up into
 * molecules. It can compare to bead structures and determine if they are the same
 * structure. 
 **/

class BeadStructure
    : protected Graph
{
public:
    BeadGraph() {};
    ~BeadGraph() {}

    // This assumes that a bead is never composed of more than a single molecule
    bool isSingleMolecule();

    // Follows same method name as topology class
    int BeadCount();
    AddBead(BaseBead * bead);

    std::vector<BaseBead *> getNeighBeads(int index);
    BaseBead * getBead(int index);

    std::vector<BeadStructure *> breakIntoMolecules();

    BeadStructure!=(const BeadStructure&) const;
    BeadStructure==(const BeadStructure&) const;
    
       
private:
    
};

}}

#endif	// _VOTCA_CSG_BEADSTRUCTURE_H 

