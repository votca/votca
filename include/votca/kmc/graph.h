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

#ifndef __VOTCA_KMC_GRAPH_H_
#define __VOTCA_KMC_GRAPH_H_

#include <vector>
#include <list>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/kmc/node.h>

namespace votca { namespace kmc {
  
using namespace std;

typedef votca::tools::vec myvec;

enum CorrelationType{Uncorrelated, Correlated, Anticorrelated };

class Graph {
public:
    void Load_graph(string SQL_graph_filename);    
    void Generate_cubic_graph(  int nx, int ny, int nz, double lattice_constant,
                                double hopping_distance, double disorder_strength, double disorder_ratio, CorrelationType correlation_type);
    
    vector<Node*> nodes;
    
    myvec sim_box_size;    
    int max_pair_degree;
    double hopping_distance;
 
private:
    void Load_graph_nodes(string filename);
    void Load_graph_static_energies(string filename);
    void Load_graph_pairs(string filename);
    void Load_graph_static_event_info(string filename);
    
    void Create_cubic_graph_nodes(int nx, int ny, int nz, double lattice_constant);
    void Create_static_energies(votca::tools::Random2 *RandomVariable, double disorder_strength, double disorder_ratio, CorrelationType correlation_type);
    
    void Determine_graph_pairs(double hopping_distance);
    
    double Determine_hopping_distance(vector<Node*> nodes);
    myvec Determine_sim_box_size(vector<Node*> nodes, double hopping_distance);
    int Determine_max_pair_degree();
    
};

void Graph::Load_graph(string filename){
    Load_graph_nodes(filename);
    Load_graph_static_energies(filename);
    Load_graph_pairs(filename);
    Load_graph_static_event_info(filename);
    
}

void Graph::Load_graph_nodes(string filename) {
    
    // Load nodes
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ FROM segments;");
   
    while (stmt->Step() != SQLITE_DONE) {
        
        Node *newNode = new Node();
        nodes.push_back(newNode);

        newNode->node_ID  = stmt->Column<int>(0);
        
        double positionX = stmt->Column<double>(1);
        double positionY = stmt->Column<double>(2);
        double positionZ = stmt->Column<double>(3);
        myvec node_position = myvec (positionX, positionY, positionZ);
        newNode->node_position = node_position;

    }
  
    delete stmt;
    stmt = NULL;

}

void Graph::Load_graph_static_energies(string filename) {
  
    // Load Static energies
    votca::tools::Database db;
    db.Open( filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion, eNeutral, eCation, ucCnNe, ucCnNh FROM segments;");

    int read_index = 0;
    while (stmt->Step() != SQLITE_DONE) {
        
        Node* loadNode = nodes[read_index];
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

void Graph::Load_graph_pairs(string filename) {
    
    // Load Node Pairs
    votca::tools::Database db;
    db.Open(filename);
    votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2' FROM pairs UNION "
                                               "SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2' FROM pairs ORDER BY segment1;");

    while (stmt->Step() != SQLITE_DONE) {
        
        int node_ID1 = stmt->Column<int>(0);
        int node_ID2 = stmt->Column<int>(1);
        Node* node1 = nodes[node_ID1];
        Node* node2 = nodes[node_ID2];
        
        node1->setPair(node2);
    }
        
    delete stmt;
    stmt = NULL;
 
}

void Graph::Load_graph_static_event_info(string filename) {
    
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
      Node* node1 = nodes[node_ID1];
      Node* node2 = nodes[node_ID2];

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

double Graph::Determine_hopping_distance(vector<Node*> nodes) {
    
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


myvec Graph::Determine_sim_box_size(vector<Node*> nodes, double hopdistance) {
    
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
            
            if(pnode1.x()-pnode2.x() > hopdistance){
                bndcrosspairXfound = true;
                sim_box_sizeX = pnode1.x() + dr.x() - pnode2.x(); 
            }
            if(pnode1.x()-pnode2.x() < -1.0*hopdistance){
                bndcrosspairXfound = true;
                sim_box_sizeX = pnode2.x() - dr.x() - pnode1.x(); 
            }            

            if(pnode1.y()-pnode2.y() > hopdistance){
                bndcrosspairYfound = true;
                sim_box_sizeY = pnode1.y() + dr.y() - pnode2.y(); 
            }
            if(pnode1.y()-pnode2.y() < -1.0*hopdistance){
                bndcrosspairYfound = true;
                sim_box_sizeY = pnode2.y() - dr.y() - pnode1.y(); 
            }

            if(pnode1.z()-pnode2.z() > hopdistance){
                bndcrosspairZfound = true;
                sim_box_sizeZ = pnode1.z() + dr.z() - pnode2.z(); 
            }
            if(pnode1.z()-pnode2.z() < -1.0*hopdistance){
                bndcrosspairZfound = true;
                sim_box_sizeZ = pnode2.z() - dr.z() - pnode1.z(); 
            }
        }
    }
    
    //for the possible outcome that none of the pairs have crossed the simulation box boundary
    if(!bndcrosspairXfound) {sim_box_sizeX = maxX-minX;}
    if(!bndcrosspairYfound) {sim_box_sizeY = maxY-minY;}
    if(!bndcrosspairZfound) {sim_box_sizeZ = maxZ-minZ;}

    myvec sim_box_size = myvec(sim_box_sizeX, sim_box_sizeY, sim_box_sizeZ);
    return sim_box_size;
}





}}



#endif

