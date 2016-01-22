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

#include "votca/xtp/graphcubic.h"


using namespace std;

namespace votca {
    namespace xtp {
 /*       
void GraphCubic::Create_cubic_graph_nodes(int NX, int NY, int NZ, double lattice_constant, myvec front, myvec back) {
    
    int node_index = 0;
    
    for(int ix=0; ix<NX; ix++) {
        for(int iy=0; iy<NY; iy++) {
            for(int iz=0; iz<NZ; iz++) {
                DNode *newDNode = new DNode();
                AddNode(newDNode);

                newDNode->node_ID = node_index;
                newDNode->node_type = Normal;
                
                myvec nodeposition = myvec(front.x() + ix*lattice_constant,front.y() + iy*lattice_constant,front.z() + iz*lattice_constant);
                newDNode->node_position = nodeposition;
                
                node_index++;    
            }
        }
    }
    
//    double sim_box_sizeX = front.x() + back.x() + NX*lattice_constant;
//    double sim_box_sizeY = front.y() + back.y() + NY*lattice_constant;
//    double sim_box_sizeZ = front.z() + back.z() + NZ*lattice_constant;
    
//    sim_box_size = myvec(sim_box_sizeX,sim_box_sizeY,sim_box_sizeZ);
}

void GraphCubic::Create_static_energies(votca::tools::Random2 *RandomVariable, double disorder_strength, double disorder_ratio, CorrelationType correlation_type){
      
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

void GraphCubic::Initialize(){;}
  */      
    }
}
