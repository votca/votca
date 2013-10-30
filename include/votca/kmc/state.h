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
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/kmc/carrier.h>
#include <votca/kmc/graph.h>
#include <votca/kmc/globaleventinfo.h>

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;


class State {
public:
    
    // Storage and readout of the node_id's of the nodes on which the carriers are to/from a SQL database
    void Save(string SQL_state_filename);
    void Load(string SQL_state_filename, Graph* graph, Globaleventinfo* globevent);

    // Buying/Selling of carrier numbers from the reservoir
    unsigned int Buy(vector <Carrier*> carriers, vector <int> carrier_reservoir);
    void Sell(vector <Carrier*> carriers, vector <int> carrier_reservoir, unsigned int remove_from_sim_box);
    void Grow(vector <Carrier*> carriers, vector <int> carrier_reservoir, unsigned int nr_new_carriers, int max_pair_degree);
    
    vector<Carrier*> electrons; //change
    vector<int> electron_reservoir;
    vector<Carrier*> holes; //change
    vector<int> hole_reservoir;

    string SQL_state_filename;

    int meshsizeX; int meshsizeY; int meshsizeZ;
    vector< vector< vector< vector <list<int> > > > > coulomb_mesh;
    void Init_coulomb_mesh(Graph* graph, Globaleventinfo* globevent);
    void Add_to_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent);
    void Remove_from_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent);
  
private:
    bool El_in_sim_box(int electron_nr) {return electrons[electron_nr]->is_in_sim_box;}
    bool Ho_in_sim_box(int hole_nr) {return holes[hole_nr]->is_in_sim_box;}
};

void State::Init_coulomb_mesh(Graph* graph, Globaleventinfo* globevent){
    meshsizeX = ceil(graph->sim_box_size.x()/globevent->coulcut);
    meshsizeY = ceil(graph->sim_box_size.y()/globevent->coulcut);
    meshsizeZ = ceil(graph->sim_box_size.z()/globevent->coulcut);
    
    coulomb_mesh.resize(meshsizeX);
    for(int i = 0;i<meshsizeX;i++) {
        coulomb_mesh[i].resize(meshsizeY);
        for(int j = 0;j<meshsizeY;j++) {
            coulomb_mesh[i][j].resize(meshsizeZ);
            for(int k=0;k<meshsizeZ;k++) {
                coulomb_mesh[i][j][k].resize(2);
            }
        }
    }
    
    for(int ic=0;ic<electrons.size();ic++){
        Add_to_coulomb_mesh(graph, electrons[ic], globevent);
    }

    for(int ic=0;ic<holes.size();ic++){
        Add_to_coulomb_mesh(graph, holes[ic], globevent);
    }
}

void State::Add_to_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent){

    int charge;
    if(carrier->carrier_type == Electron) {
        charge = 0;
    }
    else if(carrier->carrier_type == Hole) {
        charge = 1;
    }
    
    double posx = graph->nodes[carrier->carrier_node_ID]->node_position.x();
    double posy = graph->nodes[carrier->carrier_node_ID]->node_position.y();
    double posz = graph->nodes[carrier->carrier_node_ID]->node_position.z();
        
    int iposx = floor(posx/globevent->coulcut); 
    int iposy = floor(posy/globevent->coulcut); 
    int iposz = floor(posz/globevent->coulcut);
    
    coulomb_mesh[iposx][iposy][iposz][charge].push_back(carrier->carrier_ID);    
}

void State::Remove_from_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent){

    int charge;
    if(carrier->carrier_type == Electron) {
        charge = 0;
    }
    else if(carrier->carrier_type == Hole) {
        charge = 1;
    }
    
    double posx = graph->nodes[carrier->carrier_node_ID]->node_position.x();
    double posy = graph->nodes[carrier->carrier_node_ID]->node_position.y();
    double posz = graph->nodes[carrier->carrier_node_ID]->node_position.z();
        
    int iposx = floor(posx/globevent->coulcut); 
    int iposy = floor(posy/globevent->coulcut); 
    int iposz = floor(posz/globevent->coulcut);
    
    coulomb_mesh[iposx][iposy][iposz][charge].remove(carrier->carrier_ID);    
}    


void State::Save(string SQL_state_filename){
    
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
    
    for(int electron_nr = 0;electron_nr<electrons.size();electron_nr++) {
        if (El_in_sim_box(electron_nr)) {
            stmt->Bind(1, electrons[electron_nr]->carrier_node_ID);
            stmt->Bind(2, 0);
            myvec carrier_distance = electrons[electron_nr]->carrier_distance;
            stmt->Bind(3, carrier_distance.x());
            stmt->Bind(4, carrier_distance.y()); 
            stmt->Bind(5, carrier_distance.z());
            stmt->InsertStep();
            stmt->Reset();
        }        
    }
    
    for(int hole_nr = 0;hole_nr<holes.size();hole_nr++) {
        if (Ho_in_sim_box(hole_nr)) {
            stmt->Bind(1, holes[hole_nr]->carrier_node_ID);
            stmt->Bind(2, 1);
            myvec carrier_distance = holes[hole_nr]->carrier_distance;
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

void State::Load(string SQL_state_filename, Graph* graph, Globaleventinfo* globevent){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    
    votca::tools::Statement *stmt; 
    stmt = db.Prepare("SELECT node_id, carrier_type, distanceX, distanceY, distanceZ FROM carriers;");
    
    while (stmt->Step() != SQLITE_DONE)
    {   
        int cartype = stmt->Column<int>(1);
        if(cartype == 0) { // electron
            if(electron_reservoir.empty()) {Grow(electrons,electron_reservoir,globevent->state_grow_size, graph->max_pair_degree);}
            int electron_nr = Buy(electrons, electron_reservoir);
            int carnode_ID = stmt->Column<int>(0);
            electrons[electron_nr]->carrier_node_ID = carnode_ID;
            electrons[electron_nr]->carrier_type = Electron;
            graph->nodes[carnode_ID]->carriers_on_node.push_back(electrons[electron_nr]);
            double distancex = stmt->Column<double>(2);
            double distancey = stmt->Column<double>(3);
            double distancez = stmt->Column<double>(4);
            electrons[electron_nr]->carrier_distance = myvec(distancex,distancey,distancez);            
        }
        else if(cartype == 1) { // hole
            if(hole_reservoir.empty()) {Grow(holes,hole_reservoir,globevent->state_grow_size, graph->max_pair_degree);}
            int hole_nr = Buy(holes, hole_reservoir);
            int carnode_ID = stmt->Column<int>(0);
            holes[hole_nr]->carrier_node_ID = carnode_ID;
            holes[hole_nr]->carrier_type = Hole;
            graph->nodes[carnode_ID]->carriers_on_node.push_back(holes[hole_nr]);
            double distancex = stmt->Column<double>(2);
            double distancey = stmt->Column<double>(3);
            double distancez = stmt->Column<double>(4);
            holes[hole_nr]->carrier_distance = myvec(distancex,distancey,distancez);                
        }
    }
    delete stmt;
    stmt = NULL;    
}

unsigned int State::Buy(vector <Carrier*> carriers, vector <int> carrier_reservoir) {
    
    unsigned int carriernr_to_sim_box = carrier_reservoir.back();
    carrier_reservoir.pop_back();
    carriers[carriernr_to_sim_box]->is_in_sim_box = true;
    return carriernr_to_sim_box;
}

void State::Sell(vector <Carrier*> carriers, vector <int> carrier_reservoir, unsigned int remove_from_sim_box) {
    
    carrier_reservoir.push_back(remove_from_sim_box);
    carriers[remove_from_sim_box]->is_in_sim_box = false;
}

void State::Grow(vector <Carrier*> carriers, vector <int> carrier_reservoir, unsigned int nr_new_carriers, int max_pair_degree) {
    
    unsigned int new_nr_carriers = carriers.size() + nr_new_carriers;
    for (unsigned int i=carriers.size(); i<new_nr_carriers; i++) {
    
        Carrier *newCarrier = new Carrier();
        carriers.push_back(newCarrier);
        
        carrier_reservoir.push_back(i);
        newCarrier->is_in_sim_box = false;
        newCarrier->carrier_ID = i;
        
        //initialize sr potential storage
        newCarrier->srfrom = 0.0;
        for(int i = 0; i<max_pair_degree; i++) {
            newCarrier->srto.push_back(0.0);
        }
        
    }
}

}} 

#endif

