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
#include <votca/kmc/bsumtree.h>

typedef votca::tools::vec myvec;

namespace votca { namespace kmc {
  
using namespace std;

enum Inject_Type {Equal, Fermi};

class State {
public:
    
    // Storage and readout of the node_id's of the nodes on which the carriers are to/from a SQL database
    void Save(string SQL_state_filename);
    void Load(string SQL_state_filename, Graph* graph, Globaleventinfo* globevent);
    
    // Start with an empty state object
    void Init();

    // Buying/Selling of carrier numbers from the reservoir
    unsigned int Buy();
    void Sell(unsigned int remove_from_sim_box);
    void Grow(unsigned int nr_new_carriers, int max_pair_degree);
    
    vector<Carrier*> carriers; //change
    vector<int> carrier_reservoir;

    string SQL_state_filename;

    // Creation of coulomb mesh
    int meshsizeX; int meshsizeY; int meshsizeZ;
    vector< vector< vector <list<int> > > > coulomb_mesh;
    void Init_coulomb_mesh(Graph* graph, Globaleventinfo* globevent);
    void Add_to_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent);
    void Remove_from_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent);
    
    // Injection and removal of charges (for example in a double carrier bulk setting) (still to be done)
   // Bsumtree* electron_inject;
   // Bsumtree* hole_inject;
   // void Initialize_inject_trees(Graph* graph, Inject_Type injecttype, Globaleventinfo* globevent);
   // void Add_charge_in_box(Node* node, CarrierType carrier_type);
   // void Remove_charge_from_box(Node* node, CarrierType carrier_type);
    
  
private:
    bool Car_in_sim_box(int carrier_nr) {return carriers[carrier_nr]->is_in_sim_box;}
};

/*void State::Initialize_inject_trees(Graph* graph, Inject_Type injecttype, Globaleventinfo* globevent) {
    electron_inject->initialize(graph->nodes.size());
    hole_inject->initialize(graph->nodes.size());
    
    for (int inode = 0; inode < graph->nodes.size(); inode++) {
        if (injecttype == Equal) {
            electron_inject->setrate(inode, 1.0);
            hole_inject->setrate(inode, 1.0);
        }
        else if(injecttype == Fermi) {
            electron_inject->setrate(inode, exp(-1.0*globevent->beta*graph->nodes[inode]->static_electron_node_energy));
            hole_inject->setrate(inode, exp(-1.0*globevent->beta*graph->nodes[inode]->static_hole_node_energy));
        }
    }
}*/

void State::Init(){
    carriers.clear();
    carrier_reservoir.clear();
}

void State::Init_coulomb_mesh(Graph* graph, Globaleventinfo* globevent){
    meshsizeX = ceil(graph->sim_box_size.x()/globevent->coulcut);
    meshsizeY = ceil(graph->sim_box_size.y()/globevent->coulcut);
    meshsizeZ = ceil(graph->sim_box_size.z()/globevent->coulcut);
    
    coulomb_mesh.resize(meshsizeX);
    for(int i = 0;i<meshsizeX;i++) {
        coulomb_mesh[i].resize(meshsizeY);
        for(int j = 0;j<meshsizeY;j++) {
            coulomb_mesh[i][j].resize(meshsizeZ);
        }
    }
    
    for(int ic=0;ic<carriers.size();ic++){
        Add_to_coulomb_mesh(graph, carriers[ic], globevent);
    }
}

void State::Add_to_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent){

    double posx = graph->nodes[carrier->carrier_node_ID]->node_position.x();
    double posy = graph->nodes[carrier->carrier_node_ID]->node_position.y();
    double posz = graph->nodes[carrier->carrier_node_ID]->node_position.z();
        
    int iposx = floor(posx/globevent->coulcut); 
    int iposy = floor(posy/globevent->coulcut); 
    int iposz = floor(posz/globevent->coulcut);
    
    coulomb_mesh[iposx][iposy][iposz].push_back(carrier->carrier_ID);    
}

void State::Remove_from_coulomb_mesh(Graph* graph, Carrier* carrier, Globaleventinfo* globevent){

    double posx = graph->nodes[carrier->carrier_node_ID]->node_position.x();
    double posy = graph->nodes[carrier->carrier_node_ID]->node_position.y();
    double posz = graph->nodes[carrier->carrier_node_ID]->node_position.z();
        
    int iposx = floor(posx/globevent->coulcut); 
    int iposy = floor(posy/globevent->coulcut); 
    int iposz = floor(posz/globevent->coulcut);
    
    coulomb_mesh[iposx][iposy][iposz].remove(carrier->carrier_ID);    
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
    
    for(int carrier_nr = 0;carrier_nr<carriers.size();carrier_nr++) {
        if (Car_in_sim_box(carrier_nr)) {
            stmt->Bind(1, carriers[carrier_nr]->carrier_node_ID);
            
            int cartype;
            if(carriers[carrier_nr]->carrier_type == Electron) {
                cartype = 0;
            }
            else if(carriers[carrier_nr]->carrier_type == Hole) {
                cartype = 1;
            }
            stmt->Bind(2, cartype);
            
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

void State::Load(string SQL_state_filename, Graph* graph, Globaleventinfo* globevent){
    
    votca::tools::Database db;
    db.Open( SQL_state_filename );
    
    votca::tools::Statement *stmt; 
    stmt = db.Prepare("SELECT node_id, carrier_type, distanceX, distanceY, distanceZ FROM carriers;");
    
    while (stmt->Step() != SQLITE_DONE)
    {   
        if(carrier_reservoir.empty()) {Grow(globevent->state_grow_size, graph->max_pair_degree);}
        int carrier_nr = Buy();
        int carnode_ID = stmt->Column<int>(0);
        carriers[carrier_nr]->carrier_node_ID = carnode_ID;
        
        CarrierType cartype;
        if(stmt->Column<int>(1) == 0) {
            cartype == Electron;
        }
        else if(stmt->Column<int>(1) == 1) {
            cartype == Hole;
        }
        carriers[carrier_nr]->carrier_type = cartype;
        
        graph->nodes[carnode_ID]->carriers_on_node.push_back(carriers[carrier_nr]);
        
        double distancex = stmt->Column<double>(2);
        double distancey = stmt->Column<double>(3);
        double distancez = stmt->Column<double>(4);
        carriers[carrier_nr]->carrier_distance = myvec(distancex,distancey,distancez);            
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

void State::Grow(unsigned int nr_new_carriers, int max_pair_degree) {
    
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

