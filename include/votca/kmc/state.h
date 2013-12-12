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

namespace votca { namespace kmc {

template <class TGraph, class TCarrier>
class State {    

    public:

    State(){}
    ~State(){
        typename std::vector<TCarrier*>::iterator it;
        for (it = _carriers.begin(); it != _carriers.end(); it++ ) delete *it;        
    }
    
    /// Add a node to the Graph
    TCarrier* AddCarrier(int id) { 
        TCarrier* carrier = new TCarrier(id);
        _carriers.push_back(carrier); 
        return carrier;
    } 

    void Print(std::ostream& out);
    
    TCarrier* GetCarrier(int id) { return _carriers[id];}
    int GetCarrierSize() {return _carriers.size();}
    
    virtual void AddCarrier( TCarrier* carrier) { _carriers.push_back(carrier); }
    
    void InitState();
    
    // Storage and readout of the node_id's of the nodes on which the carriers are to/from a SQL database
    void Save(string SQL_state_filename);
    void Load(string SQL_state_filename, TGraph* graph);

    bool In_sim_box(TCarrier* carrier) {return carrier->inbox();}    
private:
    vector<TCarrier*> _carriers;
};

template <class TGraph, class TCarrier>
void State<TGraph, TCarrier>::InitState(){
    _carriers.clear();
}

template <class TGraph, class TCarrier>
void State<TGraph, TCarrier>::Save(string SQL_state_filename){
    
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

    typename std::vector<TCarrier*>::iterator it;
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

template <class TGraph, class TCarrier>
void State<TGraph,TCarrier>::Load(string SQL_state_filename, TGraph* graph){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    
    votca::tools::Statement *stmt; 
    stmt = db.Prepare("SELECT node_id, carrier_type, distanceX, distanceY, distanceZ FROM carriers;");
    
    int carrier_ID = 0;
    
    while (stmt->Step() != SQLITE_DONE)
    {   
        int carrier_nodeID = stmt->Column<int>(0);
        Node* carrier_node =  graph->GetNode(stmt->Column<int>(0));
        
        int carrier_type = stmt->Column<int>(1);
        
        double distancex = stmt->Column<double>(2);
        double distancey = stmt->Column<double>(3);
        double distancez = stmt->Column<double>(4);
        votca::tools::vec carrier_distance(distancex,distancey,distancez);

        TCarrier* newcarrier = AddCarrier(carrier_ID);
        newcarrier->SetCarrierNode(carrier_node);
        carrier_node->AddCarrier(carrier_ID);
        newcarrier->SetCarrierType(carrier_type);
        carrier_ID++;
        
        // Set distance travelled
        newcarrier->SetDistance(carrier_distance);
        
        // New carrier is in the simulation box
        newcarrier->SetInBox(true);
        
    }

    delete stmt;
    stmt = NULL;    
}

template <class TGraph, class TCarrier>
void State<TGraph,TCarrier>::Print(std::ostream& out) {
    
    typename std::vector<TCarrier*>::iterator it;    
    for(it = _carriers.begin(); it != _carriers.end(); it++)  out << (*it)->id() << " " << (*it)->node()->position() << " " << 
            (*it)->type() << " " << (*it)->distance() << " " << (*it)->inbox() << endl;
}

}} 

#endif

