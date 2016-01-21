/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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
#include <fstream>
#include <iostream>
#include <votca/xtp/carrier.h>
#include <votca/xtp/node.h>
#include <votca/xtp/graph.h>

namespace votca { namespace xtp {

template <class TGraph, class TCarrier>
class State {    

    public:

    State(){}
    virtual ~State(){
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
    void Save(const char* filename);
    void Load(const char* filename, TGraph* graph);

    bool In_sim_box(TCarrier* carrier) {return carrier->inbox();}    
private:
    vector<TCarrier*> _carriers;
};

template <class TGraph, class TCarrier>
void State<TGraph, TCarrier>::InitState(){
    _carriers.clear();
}

template <class TGraph, class TCarrier>
void State<TGraph, TCarrier>::Save(const char* filename){
    
    ofstream statestore;
    statestore.open(filename);
    
    if(statestore.is_open()) {
        typename std::vector<TCarrier*>::iterator it;
        for (it = _carriers.begin(); it != _carriers.end(); it++ ) {
            if(In_sim_box((*it))) {
                votca::tools::vec c_distance = (*it)->distance();
                statestore << (*it)->node()->id() << " " << (*it)->type() << " " << c_distance.x() << " " << c_distance.y() << " " << c_distance.z() << "\n";
            }
        }
    }
    else {
        std::cout << "WARNING: Can't open state store file" << "\n";
    }
    
    statestore.close();
}

template <class TGraph, class TCarrier>
void State<TGraph,TCarrier>::Load(const char* filename, TGraph* graph){
    
    ifstream statestore;
    statestore.open(filename);
    string line;
    
    int car_nodeID; int car_type;
    double car_distancex; double car_distancey; double car_distancez;
    
    int carrier_ID = 0;
    
    if(statestore.is_open()) {
        while (getline(statestore,line)) {

            //read and parse lines
            std::istringstream linestream(line);
            linestream >> car_nodeID >> car_type >> car_distancex >> car_distancey >> car_distancez;
            
            //add carrier to node
            TCarrier* newcarrier = AddCarrier(carrier_ID);
            Node* carrier_node =  graph->GetNode(car_nodeID);
            newcarrier->SetCarrierNode(carrier_node);
            carrier_node->AddCarrier(carrier_ID);

            //set carrier type
            newcarrier->SetCarrierType(car_type);

            // Set distance travelled
            votca::tools::vec carrier_distance(car_distancex,car_distancey,car_distancez);
            newcarrier->SetDistance(carrier_distance);

            // New carrier is in the simulation box
            newcarrier->SetInBox(true);
            carrier_ID++;
        }
    }
    else {
        std::cout << "WARNING: Can't open state store file" << "\n";
    }

    statestore.close();    
}

template <class TGraph, class TCarrier>
void State<TGraph,TCarrier>::Print(std::ostream& out) {
    
    typename std::vector<TCarrier*>::iterator it;    
    for(it = _carriers.begin(); it != _carriers.end(); it++)  out << (*it)->id() << " " << (*it)->node()->position() << " " << 
            (*it)->type() << " " << (*it)->distance() << " " << (*it)->inbox() << "\n";
}

}} 

#endif

