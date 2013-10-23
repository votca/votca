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

#ifndef __VOTCA_KMC_STATE_H_
#define __VOTCA_KMC_STATE_H_

#include <vector>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/kmc/carrier.h>
#include <votca/kmc/graph.h>

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;


class State {
public:
    
    // Storage and readout of the node_id's of the nodes on which the carriers are to/from a SQL database
    void Save(string SQL_state_filename);
    void Load(string SQL_state_filename);

    // Buying/Selling of carrier numbers from the reservoir
    unsigned int Buy();
    void Sell(unsigned int carrier_nr);
    
    // Initialization/Growth of the reservoir structure
    void Grow(unsigned int growsize);
    int state_grow_size;
    
    vector<Carrier*> carriers; //Carriers which are either in the simulation box or in the reservoir (as defined by the in_sim_box vector)
    vector<int> carrier_reservoir;
    
    vector<vector<vector<Carrier*> > > carriers_in_mesh;
    
    string SQL_state_filename;

    Graph* graph;

private:
    bool in_sim_box(int carrier_nr) {return carriers[carrier_nr]->is_in_sim_box;}
    
};


void State::Save(string SQL_state_filename){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    db.BeginTransaction();       

    int carrier_ID = 0;
    
    votca::tools::Statement *stmt;
    stmt = db.Prepare( "INSERT INTO carriers ("
                            "carrier_id,    node_id,  distanceX,  "
                            "distanceY, distanceZ)"
                            "VALUES ("
                            "?,     ?,     ?,"
                            "?,     ?)");
    
    for(int carrier_nr = 0;carrier_nr<carriers.size();carrier_nr++) {
        if (in_sim_box(carrier_nr)) {
            stmt->Bind(1, carrier_ID);
            stmt->Bind(2, carriers[carrier_nr]->carrier_node->node_ID);
            myvec carrier_distance = carriers[carrier_nr]->carrier_distance;
            stmt->Bind(3, carrier_distance.x());
            stmt->Bind(4, carrier_distance.y()); 
            stmt->Bind(5, carrier_distance.z());
            stmt->InsertStep();
            stmt->Reset();
        }        
    }
    
    delete stmt;
    stmt = NULL;
    
    db.EndTransaction();
}

void State::Load(string SQL_state_filename){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    
    votca::tools::Statement *stmt; 
    stmt = db.Prepare("SELECT carrier_id, node_id, distanceX, distanceY, distanceZ FROM carriers;");
    
    while (stmt->Step() != SQLITE_DONE)
    {
        if(carrier_reservoir.empty()) {
            Grow(state_grow_size);
        }
        
        int carrier_nr = Buy();

        int carrier_ID = stmt->Column<int>(0);
        int node_ID = stmt->Column<int>(1);
        double distancex = stmt->Column<double>(2);
        double distancey = stmt->Column<double>(3);
        double distancez = stmt->Column<double>(4);
        myvec carrier_distance = myvec(distancex,distancey,distancez);
        
        carriers[carrier_nr]->carrier_ID = carrier_ID;
        carriers[carrier_nr]->carrier_node = graph->nodes[node_ID];
        carriers[carrier_nr]->carrier_distance = carrier_distance;
    }
    delete stmt;
    stmt = NULL;    
}

unsigned int State::Buy() {
  unsigned int carriernr_to_sim_box = carrier_reservoir.back();
  carrier_reservoir.pop_back();
  carriers[carriernr_to_sim_box]->is_in_sim_box = true;
  return carriernr_to_sim_box;
}

void State::Sell(unsigned int remove_from_sim_box) {
  carrier_reservoir.push_back(remove_from_sim_box);
  carriers[remove_from_sim_box]->is_in_sim_box = false;
}

void State::Grow(unsigned int nr_carriers) {
  unsigned int new_nr_carriers = carriers.size() + nr_carriers;
  carriers.resize(new_nr_carriers);
  for (unsigned int i=carriers.size(); i<new_nr_carriers; i++) {
    carrier_reservoir.push_back(i);
    carriers[i]->is_in_sim_box = false;
  }
}

}} 

#endif

