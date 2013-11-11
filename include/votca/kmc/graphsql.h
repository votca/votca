/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_GRAPHSQL_H_
#define __VOTCA_KMC_GRAPHSQL_H_

#include <vector>
#include <votca/kmc/graph.h>

namespace votca { namespace kmc {
  
class GraphSQL : public Graph {

public:
   
    void Initialize();
    
    void Load_graph_segments(string filename);
    void Load_graph_pairs(string filename);
    void Load_graph_static_event_info(string filename);
    
    
};

void GraphSQL::Initialize(){
    ;
}

void GraphSQL::Load_graph_segments(string filename) {
    
    // Load nodes
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ, UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion, eNeutral, eCation, ucCnNe, ucCnNh FROM segments;");
    
    while (stmt->Step() != SQLITE_DONE) {
        
        DNode *newDNode = new DNode();
        AddNode(newDNode);

        newDNode->node_ID  = stmt->Column<int>(0);
        newDNode->node_type = Normal;
        
        double positionX = stmt->Column<double>(1);
        double positionY = stmt->Column<double>(2);
        double positionZ = stmt->Column<double>(3);
        myvec node_position = myvec (positionX, positionY, positionZ);
        newDNode->node_position = node_position;
        
        newDNode->reorg_intorig_hole= stmt->Column<double>(4); // UnCnNe
        newDNode->reorg_intorig_electron = stmt->Column<double>(5); // UnCnNh
        newDNode->reorg_intdest_hole = stmt->Column<double>(6); // UnNcCe
        newDNode->reorg_intdest_electron = stmt->Column<double>(7); // UcNcCh
        
        double eAnion = stmt->Column<double>(8);
        double eNeutral = stmt->Column<double>(9);
        double eCation = stmt->Column<double>(10);
        
        double internal_energy_electron = stmt->Column<double>(11);
        double internal_energy_hole = stmt->Column<double>(12);
        
        double static_electron_node_energy = eCation + internal_energy_electron;
        double static_hole_node_energy = eAnion + internal_energy_hole;

        newDNode->eAnion = eAnion;
        newDNode->eNeutral = eNeutral;
        newDNode->eCation = eCation;
        
        newDNode->internal_energy_electron = internal_energy_electron;
        newDNode->internal_energy_hole = internal_energy_hole;
        
        newDNode->static_electron_node_energy = static_electron_node_energy;
        newDNode->static_hole_node_energy = static_hole_node_energy;
    }
  
    delete stmt;
    stmt = NULL;
   
}

}}



#endif

