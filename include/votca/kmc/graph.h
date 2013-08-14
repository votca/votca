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

#include <votca/kmc/node.h>
#include <votca/tools/database.h>
#include <votca/tools/vec.h>
#include <votca/tools/random.h>


namespace votca { namespace kmc {
  
using namespace std;

class Graph {
 public:
//     void Graph_Load(){};
     void CreateCubicLattice(int NX, int NY, int NZ, double latt_const, double hopping_distance);
     void CreateGaussianEnergyLandscape(double disorder_strength, int numberofnodes);
     //const vector<Node*> &getNeighbours( Node* node, CarrierType type ) { return node->getNeighbours(type); }
   
 private:
     vector< Node* > nodes;
     
};

/*void Graph::Graph_Load() {
    
    votca::tools::Database db;
    db.Open( _filename );
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ FROM segments;");
    
    int node_id;
    vec nodeposition;
    
    int index = 0;
    int numberofnodes = 0;
    
    while (stmt->Step() != SQLITE_DONE) {
        Node *newNode = new Node();
        nodes.push_back(newNode);
   
        node_id = stmt->Column<int>(0);
        nodeposition = vec(stmt->Column<double>(1),stmt->Coulumb<double(2),stmt->Coulomb<double(3)); //positions in nm

        nodes[index]._id = node_id;
        nodes[index]._position = nodeposition;
        
        index++;
        numberofnodes++;
    }
    delete stmt;

    totalnumberofnodes = numberofnodes;

    // Load pairs and distances    
    
    stmt = db.Prepare("SELECT seg1-1, seg2-1, drX, drY, drZ FROM pairs ORDER BY segment1;");

    int numberofpairs = 0;
    
    while (stmt->Step() != SQLITE_DONE)
    {
        int node1 = stmt->Column<int>(0);
        int node2 = stmt->Column<int>(1);

        nodes[node1]->setNeighbours(node2);
        nodes[node2]->setNeighbours(node1);
    }    
    delete stmt;

   

}*/

void Graph:: CreateCubicLattice(int NX, int NY, int NZ, double latt_const, double hopping_distance) {
    
    //specify lattice types (square lattice, fcc lattice, fractal lattice?)
    //and dimensions
    
    int node_id;
    vec nodeposition;
    int numberofnodes = 0;
    int index = 0;
    
    for(int ix = 0; ix<NX; ix++) {
        for(int iy=0; iy<NY; iy++) {
            for(int iz=0; iz<NZ; iz++) {
                Node *newNode = new Node();
                nodes.push_back(newNode);

                node_id = NX*NY*iz + NX*iy + ix-1;
                nodeposition = vec(latt_const*ix,latt_const*iy,latt_const*iz); //positions in nm

                nodes[index]->_id = node_id;
                nodes[index]->_position = nodeposition;
                    
                index++;
                numberofnodes++;
            }
        }
    }
    
    // totalnumberofnodes = numberofnodes;
    
    //specify pairs
    
    vec position1;
    vec position2;
    double distancesqr;
    int numberofpairs = 0;
    
    for (int index1=0;index1<numberofnodes;index1++) {
        for (int index2 = 0;index2<numberofnodes;index++) {
            position1 = nodes[index1]->getPosition();
            position2 = nodes[index2]->getPosition();
            
            distancesqr = (position1.x()-position2.x())*(position1.x()-position2.x())  +  (position1.y()-position2.y())*(position1.y()-position2.y()) + (position1.z()-position2.z())*(position1.z()-position2.z());
            
            if(distancesqr <= hopping_distance*hopping_distance) {
                nodes[index1]->setNeighbours(nodes[index2]);
                nodes[index2]->setNeighbours(nodes[index1]);
                numberofpairs++;
            }
        }
    }
}

void Graph:: CreateGaussianEnergyLandscape(double disorder_strength, int numberofnodes) {
   
    for (int index1 = 0; index1<numberofnodes; index1++) {
        votca::tools::Random *RandomVariable = new votca::tools::Random();
        RandomVariable->init(rand(), rand(), rand(), rand());
        nodes[index1]->_static_energy = RandomVariable->rand_gaussian(disorder_strength);
    }
}

}



#endif

