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

#ifndef __VOTCA_KMC_GRAPH_LATTICE_H_
#define __VOTCA_KMC_GRAPH_LATTICE_H_

#include <vector>
#include <list>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/tools/random2.h>
#include <votca/kmc/dnode.h>
#include <votca/kmc/globaleventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

typedef votca::tools::vec myvec;

enum CorrelationType{Uncorrelated, Correlated, Anticorrelated };

class GraphLattice {

public:
     GraphLattice() {};
     
    ~GraphLattice() {
        vector<DNode*>::iterator it;
        for (it = nodes.begin(); it != nodes.end(); it++ ) delete *it;
    };   
    
    void Load_graph(string SQL_graph_filename, double left_electrode_distance, double right_electrode_distance, Globaleventinfo* globevent);    
    void Generate_cubic_graph(  int nx, int ny, int nz, double lattice_constant,
                                double disorder_strength,votca::tools::Random2 *RandomVariable, 
                                double disorder_ratio, CorrelationType correlation_type, double left_electrode_distance, double right_electro_distance,
                                Globaleventinfo* globevent);
    
    vector<DNode*> nodes;
    DNode* left_electrode;
    DNode* right_electrode;
    
    myvec sim_box_size;    
    int max_pair_degree;
    double hopdist;
    
    int nr_left_injector_nodes;
    int nr_right_injector_nodes;
    
    int nodemeshsizeX; int nodemeshsizeY; int nodemeshsizeZ;
    vector< vector< vector <list<DNode*> > > > node_mesh;
    void Init_node_mesh(myvec sim_box_size, double hopdist);
    void Add_to_node_mesh(DNode* node, double hopdist);
    
private:
    
    void Load_graph_nodes(string filename);
    void Load_graph_static_energies(string filename);
    void Load_graph_pairs(string filename);
    void Load_graph_static_event_info(string filename);
    
    void Create_cubic_graph_nodes(int nx, int ny, int nz, double lattice_constant, myvec front, myvec back);
    void Create_static_energies(votca::tools::Random2 *RandomVariable, double disorder_strength, double disorder_ratio, CorrelationType correlation_type);
    
    void Determine_graph_pairs(vector<DNode*> nodes, double hopdist, int nodemeshsizeX, int nodemeshsizeY, int nodemeshsizeZ,
                                         vector< vector< vector <list<DNode*> > > > node_mesh );
    
    void Setup_device_graph(vector<DNode*> nodes, DNode* left_electrode, DNode* right_electrode, double hopdist, double left_electrode_distance, double right_electrode_distance);
    void Break_periodicity(vector<DNode*>nodes , bool x_direction, bool y_direction, bool z_direction);
    
    double Determine_hopping_distance(vector<DNode*> nodes);
    myvec Determine_sim_box_size(vector<DNode*> nodes);
    int Determine_max_pair_degree(vector<DNode*> nodes);

    void Set_all_self_image_potential(vector<DNode*> nodes, myvec sim_box_size, Globaleventinfo* globevent);   
    double Calculate_self_image_potential(double nodeposx, double length, Globaleventinfo* globevent);

    myvec Periodicdistance(myvec init, myvec final, myvec boxsize);    
    
};

void GraphLattice::Init_node_mesh(myvec sim_box_size, double hopdist){
    nodemeshsizeX = ceil(sim_box_size.x()/hopdist);
    nodemeshsizeY = ceil(sim_box_size.y()/hopdist);
    nodemeshsizeZ = ceil(sim_box_size.z()/hopdist);
    
    node_mesh.resize(nodemeshsizeX);
    for(int i = 0;i<nodemeshsizeX;i++) {
        node_mesh[i].resize(nodemeshsizeY);
        for(int j = 0;j<nodemeshsizeY;j++) {
            node_mesh[i][j].resize(nodemeshsizeZ);
        }
    }
    
    for(int inode=0;inode<nodes.size();inode++){
        Add_to_node_mesh(nodes[inode], hopdist);
    }
}

void GraphLattice::Add_to_node_mesh(DNode* node, double hopdist){

    double posx = node->node_position.x();
    double posy = node->node_position.y();
    double posz = node->node_position.z();
        
    int iposx = floor(posx/hopdist); 
    int iposy = floor(posy/hopdist); 
    int iposz = floor(posz/hopdist);
    
    node_mesh[iposx][iposy][iposz].push_back(node);       
}

void GraphLattice::Load_graph(string filename, double left_electrode_distance, double right_electrode_distance, Globaleventinfo* globevent){
    
    Load_graph_nodes(filename);
    Load_graph_static_energies(filename);
    Load_graph_pairs(filename);
    Load_graph_static_event_info(filename);
    
    hopdist = Determine_hopping_distance(nodes);
    sim_box_size = Determine_sim_box_size(nodes);
    max_pair_degree = Determine_max_pair_degree(nodes);
    
    if(globevent->device) {
        Setup_device_graph(nodes,left_electrode,right_electrode,hopdist,left_electrode_distance,right_electrode_distance);
        Set_all_self_image_potential(nodes,sim_box_size,globevent);
        Init_node_mesh(sim_box_size, hopdist);
    }
    
}

void GraphLattice::Generate_cubic_graph(int nx, int ny, int nz, double lattice_constant,
                                double disorder_strength, votca::tools::Random2 *RandomVariable, 
                                double disorder_ratio, CorrelationType correlation_type, double left_electrode_distance, double right_electrode_distance,
                                Globaleventinfo* globevent) {

    Create_cubic_graph_nodes(nx, ny, nz, lattice_constant, myvec(0.0,0.0,0.0), myvec (lattice_constant, lattice_constant, lattice_constant));
    sim_box_size = Determine_sim_box_size(nodes);
    Init_node_mesh(sim_box_size, hopdist);
    Determine_graph_pairs(nodes,hopdist,nodemeshsizeX,nodemeshsizeY,nodemeshsizeZ, node_mesh);
    Create_static_energies(RandomVariable, disorder_strength, disorder_ratio, correlation_type);  

    if(globevent->device){
        Setup_device_graph(nodes,left_electrode,right_electrode,hopdist,left_electrode_distance,right_electrode_distance);
        Set_all_self_image_potential(nodes,sim_box_size,globevent);
    }
    
}

void GraphLattice::Load_graph_nodes(string filename) {
    
    // Load nodes
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ FROM segments;");
    
    while (stmt->Step() != SQLITE_DONE) {
        
        DNode *newNode = new DNode();
        nodes.push_back(newNode);

        newNode->node_ID  = stmt->Column<int>(0);
        newNode->node_type = Normal;
        
        double positionX = stmt->Column<double>(1);
        double positionY = stmt->Column<double>(2);
        double positionZ = stmt->Column<double>(3);
        myvec node_position = myvec (positionX, positionY, positionZ);
        newNode->node_position = node_position;
    }
  
    delete stmt;
    stmt = NULL;

}

void GraphLattice::Load_graph_static_energies(string filename) {
  
    // Load Static energies
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion, eNeutral, eCation, ucCnNe, ucCnNh FROM segments;");

    int read_index = 0;
    while (stmt->Step() != SQLITE_DONE) {
        
        DNode* loadNode = nodes[read_index];
        if(loadNode->node_ID != stmt->Column<int>(0)) {std::cout << "WARNING: mismatch between the node_ID's /n";}
        
        loadNode->reorg_intorig_hole= stmt->Column<double>(1); // UnCnNe
        loadNode->reorg_intorig_electron = stmt->Column<double>(2); // UnCnNh
        loadNode->reorg_intdest_hole = stmt->Column<double>(3); // UnNcCe
        loadNode->reorg_intdest_electron = stmt->Column<double>(4); // UcNcCh
        
        double eAnion = stmt->Column<double>(5);
        double eNeutral = stmt->Column<double>(6);
        double eCation = stmt->Column<double>(7);
        
        double internal_energy_electron = stmt->Column<double>(8);
        double internal_energy_hole = stmt->Column<double>(9);
        
        double static_electron_node_energy = eCation + internal_energy_electron;
        double static_hole_node_energy = eAnion + internal_energy_hole;

        loadNode->eAnion = eAnion;
        loadNode->eNeutral = eNeutral;
        loadNode->eCation = eCation;
        
        loadNode->internal_energy_electron = internal_energy_electron;
        loadNode->internal_energy_hole = internal_energy_hole;
        
        loadNode->static_electron_node_energy = static_electron_node_energy;
        loadNode->static_hole_node_energy = static_hole_node_energy;
        
        read_index++;
    }
    
    delete stmt;
    stmt = NULL;
}

void GraphLattice::Load_graph_pairs(string filename) {
    
    // Load Node Pairs
    votca::tools::Database db;
    db.Open(filename);
    votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2' FROM pairs UNION "
                                               "SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2' FROM pairs ORDER BY segment1;");

    while (stmt->Step() != SQLITE_DONE) {
        
        int node_ID1 = stmt->Column<int>(0);
        int node_ID2 = stmt->Column<int>(1);
        DNode* node1 = nodes[node_ID1];
        DNode* node2 = nodes[node_ID2];
        
        node1->setPair(node2);
    }
        
    delete stmt;
    stmt = NULL;
 
}

void GraphLattice::Load_graph_static_event_info(string filename) {
    
    // Load Node Pairs
    votca::tools::Database db;
    db.Open(filename);
    votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2',drX, drY, drZ, "
                      "rate12e AS 'rate_e', rate12h AS 'rate_h', Jeff2e, Jeff2h, lOe, l0h FROM pairs UNION "
                      "SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2',-drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ', "
                      "rate21e AS 'rate_e', rate21h AS 'rate_h', Jeff2e, Jeff2h, l0e, l0h FROM pairs ORDER BY segment1;");

    while (stmt->Step() != SQLITE_DONE) {
        
      int node_ID1 = stmt->Column<int>(0);
      int node_ID2 = stmt->Column<int>(1);
      DNode* node1 = nodes[node_ID1];
      DNode* node2 = nodes[node_ID2];

      double drX = stmt->Column<double>(2);
      double drY = stmt->Column<double>(3);
      double drZ = stmt->Column<double>(4);
      myvec dr = myvec(drX, drY, drZ); //distance between node2 and node1

      double rate12e = stmt->Column<double>(5);
      double rate12h = stmt->Column<double>(6);
      double Jeff2e = stmt->Column<double>(7);
      double Jeff2h = stmt->Column<double>(8);
      double reorg_oute = stmt->Column<double>(9); 
      double reorg_outh = stmt->Column<double>(10);
    
      nodes[node_ID1]->setStaticeventinfo(node2, dr, rate12e, rate12h, Jeff2e, Jeff2h, reorg_oute, reorg_outh);
    }
    
    delete stmt;
    stmt = NULL;
}

void GraphLattice::Create_cubic_graph_nodes(int NX, int NY, int NZ, double lattice_constant, myvec front, myvec back) {
    
    int node_index = 0;
    
    for(int ix=0; ix<NX; ix++) {
        for(int iy=0; iy<NY; iy++) {
            for(int iz=0; iz<NZ; iz++) {
                DNode *newNode = new DNode();
                nodes.push_back(newNode);

                newNode->node_ID = node_index;
                newNode->node_type = Normal;
                
                myvec nodeposition = myvec(front.x() + ix*lattice_constant,front.y() + iy*lattice_constant,front.z() + iz*lattice_constant);
                newNode->node_position = nodeposition;
                
                node_index++;    
            }
        }
    }
    
    double sim_box_sizeX = front.x() + back.x() + NX*lattice_constant;
    double sim_box_sizeY = front.y() + back.y() + NY*lattice_constant;
    double sim_box_sizeZ = front.z() + back.z() + NZ*lattice_constant;
    
    sim_box_size = myvec(sim_box_sizeX,sim_box_sizeY,sim_box_sizeZ);
}

void GraphLattice::Create_static_energies(votca::tools::Random2 *RandomVariable, double disorder_strength, double disorder_ratio, CorrelationType correlation_type){
      
    for(int inode=0;inode<nodes.size();inode++) {
      
        double el_node_energy = RandomVariable->rand_gaussian(disorder_strength);
        double ho_node_energy;
        nodes[inode]->static_electron_node_energy = el_node_energy;
        
        if(correlation_type == Correlated) {
            ho_node_energy = disorder_ratio*el_node_energy;
            nodes[inode]->static_hole_node_energy = ho_node_energy;
        }
        else if(correlation_type == Anticorrelated) {
            ho_node_energy = -1.0*disorder_ratio*el_node_energy;
            nodes[inode]->static_hole_node_energy = ho_node_energy;
        }
        else {
            ho_node_energy = RandomVariable->rand_gaussian(disorder_ratio*disorder_strength);
            nodes[inode]->static_hole_node_energy = ho_node_energy;
        }
    }
}

void GraphLattice::Setup_device_graph(vector<DNode*> nodes, DNode* left_electrode, DNode* right_electrode, double hopdist, double left_electrode_distance, double right_electrode_distance){

    left_electrode->node_type = LeftElectrode;
    left_electrode->static_electron_node_energy = 0.0;
    left_electrode->static_hole_node_energy = 0.0;
    right_electrode->node_type = RightElectrode;
    right_electrode->static_electron_node_energy = 0.0;
    right_electrode->static_hole_node_energy = 0.0;
    
    // Make sure the periodicity is broken up, i.e. pairs passing over the electrodes should be broken
    Break_periodicity(nodes, true, false, false);
    
    // Translate the graph due to the spatial location of the electrodes and update system box size accordingly, putting the left electrode at x = 0
    // left_electrode_distance is the distance of the left electrode to the node with minimum x-coordinate

    double minX = nodes[0]->node_position.x();    
    
    for(int inode=0; inode<nodes.size(); inode++) {

        if(minX>nodes[inode]->node_position.x()) {minX = nodes[inode]->node_position.x();}

    }
    
    //distance by which the graph should be translated is left_electrode_distance - minX

    double xtranslate = left_electrode_distance - minX;

    for(int inode=0; inode<nodes.size(); inode++) {
        double oldxpos = nodes[inode]->node_position.x();
        double newxpos = oldxpos + xtranslate;
        double ypos = nodes[inode]->node_position.y();
        double zpos = nodes[inode]->node_position.z();
        myvec newposition = myvec(newxpos,ypos,zpos);
        nodes[inode]->node_position = newposition;
    }

    //adjust system box size accordingly

    myvec old_sim_box_size = Determine_sim_box_size(nodes);
    double new_sim_box_sizeX = old_sim_box_size.x() + left_electrode_distance + right_electrode_distance;
    sim_box_size = myvec(new_sim_box_sizeX, old_sim_box_size.y(), old_sim_box_size.z());
    
    
    
    //determine the nodes which are injectable from the left electrode and the nodes which are injectable from the right electrode

    int linjector_ID = 0;
    int rinjector_ID = 0;

    for(int inode=0; inode<nodes.size(); inode++) { 
      
        double left_distance = nodes[inode]->node_position.x();
     
        if(left_distance <= hopdist) {

            myvec dr = myvec(-1.0*left_distance,0.0,0.0);            
            nodes[inode]->setPair(left_electrode);
            nodes[inode]->setStaticeventinfo(left_electrode, dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            left_electrode->setPair(nodes[inode]);
            left_electrode->setStaticeventinfo(nodes[inode], -1.0*dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            nodes[inode]->left_injector_ID = linjector_ID;
            linjector_ID++;
      
        }
      
        double right_distance = sim_box_size.x() - nodes[inode]->node_position.x();
      
        if(right_distance <= hopdist) {
          
            myvec dr = myvec(right_distance,0.0,0.0);            
            nodes[inode]->setPair(right_electrode);
            nodes[inode]->setStaticeventinfo(right_electrode, dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            right_electrode->setPair(nodes[inode]);
            right_electrode->setStaticeventinfo(nodes[inode], -1.0*dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS            
            nodes[inode]->right_injector_ID = rinjector_ID;
            rinjector_ID++;
        }
    }
    nr_left_injector_nodes = linjector_ID;
    nr_right_injector_nodes = rinjector_ID;

    // Recompute the max pair degree, due to the breaking of the periodicity and the connection of nodes to the electrode nodes
    max_pair_degree = Determine_max_pair_degree(nodes);

}

void GraphLattice::Set_all_self_image_potential(vector<DNode*> nodes, myvec sim_box_size, Globaleventinfo* globevent) {
    
    for(int inode=0; inode<nodes.size();inode++){
        myvec nodepos = nodes[inode]->node_position;
        double device_length = sim_box_size.x();
        nodes[inode]->self_image_potential = Calculate_self_image_potential(nodepos.x(),device_length,globevent);
    }
    left_electrode->self_image_potential = 0.0;
    right_electrode->self_image_potential = 0.0;
}

double GraphLattice::Calculate_self_image_potential(double nodeposx, double length, Globaleventinfo* globevent){
    
    double selfimagepot = 0.0;

    double distx_1;
    double distx_2;
    int sign;
    for (int i=0;i<globevent->nr_sr_images; i++) {
        if (div(i,2).rem==0) { // even generation
            sign = -1;
            distx_1 = i*length + 2*nodeposx;
            distx_2 = (i+2)*length - 2*nodeposx; 
        }
        else {
            sign = 1;
            distx_1 = (i+1)*length;
            distx_2 = (i+1)*length;
        }
        selfimagepot += sign*(1.0/distx_1 + 1.0/distx_2);
    }

    return globevent->self_image_prefactor*selfimagepot;        
}

void GraphLattice::Break_periodicity(vector<DNode*> nodes, bool x_direction, bool y_direction, bool z_direction){

    vector<int> remove_pairs;
    
    for(int inode=0; inode<nodes.size();inode++){
        
        remove_pairs.clear();
        
        //  determine which pairs should be removed
        
        for(int ipair=0;ipair<nodes[inode]->static_event_info.size();ipair++) {
            
            bool flagged_for_removal = false;
            
            myvec pnode1 = nodes[inode]->node_position;
            myvec pnode2 = nodes[inode]->static_event_info[ipair].pairnode->node_position;
            myvec dr = nodes[inode]->static_event_info[ipair].distance;
            
            if(x_direction){
                
                if((pnode1.x() + dr.x() > pnode2.x())&&(pnode1.x()>pnode2.x())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                if((pnode1.x() + dr.x() < pnode2.x())&&(pnode1.x()<pnode2.x())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                
            }

            if(y_direction){
                
                if((pnode1.y() + dr.y() > pnode2.y())&&(pnode1.y()>pnode2.y())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                if((pnode1.y() + dr.y() < pnode2.y())&&(pnode1.y()<pnode2.y())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                
            }
            
            if(z_direction){
                
                if((pnode1.z() + dr.z() > pnode2.z())&&(pnode1.z()>pnode2.z())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                if((pnode1.z() + dr.z() < pnode2.z())&&(pnode1.z()<pnode2.z())&&!flagged_for_removal){
                    remove_pairs.push_back(ipair);
                    flagged_for_removal = true;
                }
                
            }
            
        }

        // remove pairs
        int removed = 0;
        for(int iremove = 0; iremove<remove_pairs.size();iremove++){
            int to_be_removed = remove_pairs[iremove] - removed;
            nodes[inode]->removePair(to_be_removed); //removes pairs and static event info objects
            removed++;
        }
        
        // check max pair degree
        
    }
}    

double GraphLattice::Determine_hopping_distance(vector<DNode*> nodes) {
    
    //Determination of hopping distance
    
    double hopdistance = 0.0;
    
    for(int inode=0; inode < nodes.size(); inode++) {
        for(int ipair=0;ipair<nodes[inode]->static_event_info.size();ipair++) {
            
            myvec pairdistancevec = nodes[inode]->static_event_info[ipair].distance;
            double pairdistance = abs(pairdistancevec);
            if(pairdistance>hopdistance) {hopdistance = pairdistance;}
        }
    }
    
    return hopdistance;
}


myvec GraphLattice::Determine_sim_box_size(vector<DNode*> nodes) {
    
    //Determination of simulation box size
    //To do this, we first need to find a node with position vector a and pairing node with position vector b, such that
    //|a-b|>hopping distance
    //Note that it is possible that none of the pairs pass the simulation box boundaries
    //In this special case, we must determine the node with max x/y/z coordinate and min x/y/z coordinate
    
    bool bndcrosspairXfound = false;
    bool bndcrosspairYfound = false;
    bool bndcrosspairZfound = false;
    
    double sim_box_sizeX;
    double sim_box_sizeY;
    double sim_box_sizeZ;
    
    double maxX = nodes[0]->node_position.x(); //initial values
    double maxY = nodes[0]->node_position.y();
    double maxZ = nodes[0]->node_position.z();
    double minX = nodes[0]->node_position.x(); //initial values
    double minY = nodes[0]->node_position.y();
    double minZ = nodes[0]->node_position.z();
    
    for(int inode=0;inode<nodes.size();inode++) {
        
        if(bndcrosspairXfound&&bndcrosspairYfound&&bndcrosspairZfound) {break;}
        
        for(int ipair=0;ipair<nodes[inode]->static_event_info.size();ipair++) {
            
            if(bndcrosspairXfound&&bndcrosspairYfound&&bndcrosspairZfound) {break;}
        
            myvec pnode1 = nodes[inode]->node_position;
            myvec pnode2 = nodes[inode]->static_event_info[ipair].pairnode->node_position;
            myvec dr = nodes[inode]->static_event_info[ipair].distance;
            
            if(maxX<pnode1.x()) {maxX = pnode1.x();}
            if(minX>pnode1.x()) {minX = pnode1.x();}
            if(maxY<pnode1.y()) {maxY = pnode1.y();}
            if(minY>pnode1.y()) {minY = pnode1.y();}
            if(maxZ<pnode1.z()) {maxZ = pnode1.z();}
            if(minZ>pnode1.z()) {minZ = pnode1.z();}
            
            if((pnode1.x() + dr.x() > pnode2.x())&&(pnode1.x()>pnode2.x())){
                bndcrosspairXfound = true;
                sim_box_sizeX = pnode1.x() + dr.x() - pnode2.x(); 
            }
            if((pnode1.x() + dr.x() < pnode2.x())&&(pnode1.x()<pnode2.x())){
                bndcrosspairXfound = true;
                sim_box_sizeX = pnode2.x() - dr.x() - pnode1.x(); 
            }            

            if((pnode1.y() + dr.y() > pnode2.y())&&(pnode1.y()>pnode2.y())){
                bndcrosspairYfound = true;
                sim_box_sizeY = pnode1.y() + dr.y() - pnode2.y(); 
            }
            if((pnode1.y() + dr.y() < pnode2.y())&&(pnode1.y()<pnode2.y())){
                bndcrosspairYfound = true;
                sim_box_sizeY = pnode2.y() - dr.y() - pnode1.y(); 
            }

            if((pnode1.z() + dr.z() > pnode2.z())&&(pnode1.z()>pnode2.z())){
                bndcrosspairZfound = true;
                sim_box_sizeZ = pnode1.z() + dr.z() - pnode2.z(); 
            }
            if((pnode1.z() + dr.z() < pnode2.z())&&(pnode1.z()<pnode2.z())){
                bndcrosspairZfound = true;
                sim_box_sizeZ = pnode2.z() - dr.z() - pnode1.z(); 
            }
        }
    }
    
    //for the possible outcome that none of the pairs have crossed the simulation box boundary
    if(!bndcrosspairXfound) {sim_box_sizeX = maxX-minX;}
    if(!bndcrosspairYfound) {sim_box_sizeY = maxY-minY;}
    if(!bndcrosspairZfound) {sim_box_sizeZ = maxZ-minZ;}

    myvec simboxsize = myvec(sim_box_sizeX, sim_box_sizeY, sim_box_sizeZ);
    return simboxsize;
}

int GraphLattice::Determine_max_pair_degree(vector<DNode*> nodes){
    
    //Determination of the maximum degree in the graph
    
    int maxdegree = 0;
    
    for(int inode=0; inode < nodes.size(); inode++) {

        if(nodes[inode]->pairing_nodes.size()>maxdegree) {maxdegree = nodes[inode]->pairing_nodes.size();}

    }
    
    return maxdegree;    
}

void GraphLattice::Determine_graph_pairs(vector<DNode*> nodes, double hopdist, int nodemeshsizeX, int nodemeshsizeY, int nodemeshsizeZ,
    vector< vector< vector <list<DNode*> > > > node_mesh ) {  
  
    for (int inode = 0; inode<nodes.size(); inode++) {
      
        // Define cubic boundaries in non-periodic coordinates
        DNode* initnode = nodes[inode];
        myvec initnodepos = initnode->node_position;
    
        double ix1 = initnodepos.x()-hopdist; double ix2 = initnodepos.x()+hopdist;
        double iy1 = initnodepos.y()-hopdist; double iy2 = initnodepos.y()+hopdist;
        double iz1 = initnodepos.z()-hopdist; double iz2 = initnodepos.z()+hopdist;

        // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
        int sx1 = floor(ix1/hopdist);
        int sx2 = floor(ix2/hopdist);
        int sy1 = floor(iy1/hopdist);
        int sy2 = floor(iy2/hopdist);
        int sz1 = floor(iz1/hopdist);
        int sz2 = floor(iz2/hopdist);      
 
        // Now visit all relevant sublattices
        for (int isz=sz1; isz<=sz2; isz++) {
            int r_isz = isz;
            while (r_isz < 0) r_isz += nodemeshsizeZ;
            while (r_isz >= nodemeshsizeZ) r_isz -= nodemeshsizeZ;
            for (int isy=sy1; isy<=sy2; isy++) {
                int r_isy = isy;
                while (r_isy < 0) r_isy += nodemeshsizeY;
                while (r_isy >= nodemeshsizeY) r_isy -= nodemeshsizeY;
                for (int isx=sx1; isx<=sx2; isx++) {
                    int r_isx = isx;
                    while (r_isx < 0) r_isx += nodemeshsizeX;
                    while (r_isx >= nodemeshsizeX) r_isx -= nodemeshsizeX;
        
                    // Ask a list of all nodes in this sublattice
                    list<DNode*>::iterator li1,li2,li3;
                    list<DNode* > *nodeList = &node_mesh[r_isx][r_isy][r_isz];
                    li1 = nodeList->begin();
                    li2 = nodeList->end();
                    for (li3=li1; li3!=li2; li3++) {
                        DNode* probenode = *li3;
                        if(inode!=probenode->node_ID){ 
                            myvec probenodepos = probenode->node_position;
                            myvec differ = Periodicdistance(initnodepos,probenodepos,sim_box_size);
                            double distance = abs(differ);
                            if(distance <= hopdist) {
                                nodes[inode]->setPair(probenode);
                                nodes[inode]->setStaticeventinfo(probenode, differ, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0);
                            }
                        }
                    }
                }
            }
        }
    }
}

myvec GraphLattice::Periodicdistance(myvec init, myvec final, myvec boxsize) {
    
  myvec pre = final-init;
  
  double prex = pre.x();
  double prey = pre.y();
  double prez = pre.z();
  
  if(prex<-0.5) {prex+=boxsize.x();}
  if(prex>0.5) {prex-=boxsize.x();}
  if(prey<-0.5) {prey+=boxsize.y();}
  if(prey>0.5) {prey-=boxsize.y();}
  if(prez<-0.5) {prez+=boxsize.z();}
  if(prez>0.5) {prez-=boxsize.z();}
  
  myvec perdif = myvec(prex,prey,prez);
  
  return perdif;       
}   

}}



#endif

