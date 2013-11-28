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
#include <list>
#include <votca/kmc/carrier.h>
#include <votca/kmc/node.h>
#include <votca/kmc/graph.h>

enum CarrierType{ Electron, Hole, Exciton};

namespace votca { namespace kmc {

template <class TGraph>
class State {    

    public:

    State(){}
    ~State(){
        typename std::vector<Carrier*>::iterator it;
        for (it = _carriers.begin(); it != _carriers.end(); it++ ) delete *it;        
    }
    
    /// Add a node to the Graph
    Carrier* AddCarrier(int id, Node* node, int carrier_type) { 
        Carrier* carrier = new Carrier(id, node, carrier_type);
        _carriers.push_back(carrier); 
        return carrier;
    } 

    void Print() {std::cout << "what?" << endl;}
    
    void AddCarrier( Carrier* carrier) { _carriers.push_back(carrier); }
    
    void InitState() {_carriers.clear();}
    
    // Storage and readout of the node_id's of the nodes on which the carriers are to/from a SQL database
    void Save(string SQL_state_filename);
    void Load(string SQL_state_filename, TGraph* graph);
    
private:
    bool In_sim_box(Carrier* carrier) {return carrier->inbox();}
    vector<Carrier*> _carriers;
};

template <class TGraph>
void State<TGraph>::Save(string SQL_state_filename){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    db.BeginTransaction();       

    votca::tools::Statement *stmt;
    stmt = db.Prepare( "INSERT INTO carriers ("
                            "node_id, carrier_type,  distanceX,  "
                            "distanceY, distanceZ)"
                            "VALUES ("
                            "?,     ?,     ?,"
                            "?,     ?)");

    typename std::vector<Carrier*>::iterator it;
    for (it = _carriers.begin(); it != _carriers.end(); it++ ) {
        if(In_sim_box((*it))) {
            stmt->Bind(1, (*it)->node()->id());
            stmt->Bind(2, (*it)->type());
            votca::tools::vec carrier_distance = (*it)->distance();
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

template <class TGraph>
void State<TGraph>::Load(string SQL_state_filename, TGraph* graph){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    
    votca::tools::Statement *stmt; 
    stmt = db.Prepare("SELECT node_id, carrier_type, distanceX, distanceY, distanceZ FROM carriers;");
    
    int carrier_ID = 0;
    
    while (stmt->Step() != SQLITE_DONE)
    {   
        int carrier_nodeID = stmt->Column<int>(0);
        int carrier_type = stmt->Column<int>(1);
        double distancex = stmt->Column<double>(2);
        double distancey = stmt->Column<double>(3);
        double distancez = stmt->Column<double>(4);
        votca::tools::vec carrier_distance(distancex,distancey,distancez);

        Node* carrier_node =  graph->GetNode(stmt->Column<int>(0));
        Carrier* newcarrier = AddCarrier(carrier_ID,carrier_node, carrier_type); 
        carrier_node->AddCarrier(carrier_ID);
        carrier_ID++;
        
        // Set distance travelled
        newcarrier->SetDistance(carrier_distance);
        
        // New carrier is ofcourse in the simulation box
        newcarrier->SetInBox(true);
        
    }

    delete stmt;
    stmt = NULL;    
}

}} 

#endif

