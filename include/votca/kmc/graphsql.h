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
#include <votca/kmc/nodesql.h>
#include <votca/kmc/linksql.h>
#include <votca/tools/database.h>

namespace votca { namespace kmc {
  
class GraphSQL : public Graph<NodeSQL, LinkSQL> {

public:
   
    void Initialize();
    
    void Load_graph_segments(string filename);
    void Load_graph_links(string filename);
    
};

inline void GraphSQL::Initialize(){
    ;
}

inline void GraphSQL::Load_graph_segments(string filename) {
    
    // Load nodes
    votca::tools::Database db;
    db.Open( filename );
    //votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ, UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion, eNeutral, eCation, ucCnNe, ucCnNh FROM segments;");
    
    // only rates are needed if Coulomb interactions are is not calculated
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ FROM segments;");
    
    while (stmt->Step() != SQLITE_DONE) {
        
        int id = stmt->Column<int>(0);

        double X = stmt->Column<double>(1);
        double Y = stmt->Column<double>(2);
        double Z = stmt->Column<double>(3);
        
        NodeSQL *node = AddNode();
         
       // myvec node_position = myvec (X, Y, Z);
        
    }
  
    delete stmt;
    stmt = NULL;
   
}

inline void GraphSQL::Load_graph_links (string filename) {
    
    // Load Node Pairs
    votca::tools::Database db;
    db.Open(filename);
    votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2' FROM pairs UNION "
                                               "SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2' FROM pairs ORDER BY segment1;");

    while (stmt->Step() != SQLITE_DONE) {
        
        int node_ID1 = stmt->Column<int>(0);
        int node_ID2 = stmt->Column<int>(1);
        
        Node* init_node = getnode(node_ID1);
        Node* final_node = getnode(node_ID2);
        
        Link* newLink = new Link();
        init_node->AddLink(newLink);
        
        newLink->setnode1(init_node);
        newLink->setnode2(final_node);

    }
        
    delete stmt;
    stmt = NULL;
    
}

}}



#endif

