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

#ifndef __VOTCA_KMC_VISUALISATION_H_
#define __VOTCA_KMC_VISUALISATION_H_

#include <votca/xtp/event.h>
#include <votca/xtp/graph.h>
#include <votca/xtp/eventinfo.h>
#include <iostream>
#include <iomanip>


namespace votca { namespace xtp {
  
using namespace std;

class Visualisation {
    
public:
    
    void Init_visualisation(GraphKMC* graph, Eventinfo* eventinfo);
    void Update_visualisation(Event* event);
    void Print_visualisation(); 
 
private:
    
    void Init_layers();
    void Init_meshes();
    void Init_node_numbers(GraphKMC* graph);
    
    ofstream viz_stream;
    char viz_file[100];

    vector< vector< vector <double> > > viz_mesh;
    vector< vector< vector <int> > > num_mesh;
    vector<double> layer_cur;
    vector<int> layer_num;

    int viz_meshnr_x; int viz_meshnr_y; int viz_meshnr_z;
    double viz_size_x; double viz_size_y; double viz_size_z;
    
};

void Visualisation::Init_visualisation(GraphKMC* graph, Eventinfo* eventinfo)
{
    std::cout << "here" << endl;
    strcpy(viz_file, eventinfo->visualisation_filename.c_str());    
    viz_stream.open(viz_file);
      std::cout << "here" << endl;
  
    viz_meshnr_x = eventinfo->visualisation_number_x; viz_meshnr_y = eventinfo->visualisation_number_y; viz_meshnr_z = eventinfo->visualisation_number_z;
    viz_size_x = eventinfo->simboxsize.x()/viz_meshnr_x; viz_size_y = eventinfo->simboxsize.y()/viz_meshnr_y; viz_size_z = eventinfo->simboxsize.z()/viz_meshnr_z;
     std::cout << "here" << endl;
   
    this->Init_layers();
      std::cout << "here" << endl;

    this->Init_meshes();
        std::cout << "here" << endl;

    this->Init_node_numbers(graph);
    std::cout << "here" << endl;

}

void Visualisation::Update_visualisation(Event* event)
{
    Link* viz_link = event->link();
    Node* node1 = viz_link->node1(); Node* node2 = viz_link->node2();
    votca::tools::vec node1_pos = node1->position(); votca::tools::vec node2_pos = node2->position();

    if(node1->type() == (int) NormalNode && node2->type() == (int) NormalNode) {
        int mesh1_pos_x = floor(node1_pos.x()/viz_size_x); int mesh1_pos_y = floor(node1_pos.y()/viz_size_y); int mesh1_pos_z = floor(node1_pos.z()/viz_size_z);
        int mesh2_pos_x = floor(node2_pos.x()/viz_size_x); int mesh2_pos_y = floor(node2_pos.y()/viz_size_y); int mesh2_pos_z = floor(node2_pos.z()/viz_size_z);

        if(node1_pos.z() < node2_pos.z()) {
            if(!(((mesh1_pos_x == mesh2_pos_x)&&(mesh1_pos_y == mesh2_pos_y))&&(mesh1_pos_z == mesh2_pos_z))) {
                viz_mesh[mesh1_pos_x][mesh1_pos_y][mesh1_pos_z]+=1.0;
                viz_mesh[mesh2_pos_x][mesh2_pos_y][mesh2_pos_z]+=1.0;
                layer_cur[mesh1_pos_z] += 1.0;
                layer_cur[mesh2_pos_z] += 1.0;
            }
        }
        if(node1_pos.z() > node2_pos.z()) {
            if(!(((mesh1_pos_x == mesh2_pos_x)&&(mesh1_pos_y == mesh2_pos_y))&&(mesh1_pos_z == mesh2_pos_z))) {
                viz_mesh[mesh1_pos_x][mesh1_pos_y][mesh1_pos_z]-=1.0;
                viz_mesh[mesh2_pos_x][mesh2_pos_y][mesh2_pos_z]-=1.0;
                layer_cur[mesh1_pos_z] -= 1.0;
                layer_cur[mesh2_pos_z] -= 1.0;
            }
        }
    }
}

void Visualisation::Print_visualisation()
{
    for(int ix = 0; ix < viz_meshnr_x; ix++) {
        for(int iy=0; iy < viz_meshnr_y; iy++) {
            for(int iz=0; iz < viz_meshnr_z; iz++) {
               
               if(num_mesh[ix][iy][iz] == 0) num_mesh[ix][iy][iz] = 1;
               if(layer_num[iz] == 0) layer_num[iz] = 1;
               
               double value;
               double layer_value;

               value = viz_mesh[ix][iy][iz]/(1.0*num_mesh[ix][iy][iz]);
               layer_value = layer_cur[iz]/(1.0*layer_num[iz]);
               if(layer_value == 0) layer_value = 1;
               
               double current = value/layer_value;
               viz_stream << ix << " " << iy << " " << iz << " " << current << "\n";
            }
        }
    } 
    viz_stream.flush();
}

void Visualisation::Init_layers()
{
    layer_cur.resize(viz_meshnr_z);
    layer_num.resize(viz_meshnr_z);
    for(int i=0;i<viz_meshnr_z;i++){
        layer_cur[i]=0.0;
        layer_num[i]=0;
    }    
}

void Visualisation::Init_meshes()
{
    viz_mesh.resize(viz_meshnr_x);
    num_mesh.resize(viz_meshnr_x);
    for(int i = 0;i<viz_meshnr_x;i++) {
        viz_mesh[i].resize(viz_meshnr_y);
        num_mesh[i].resize(viz_meshnr_y);
        for(int j = 0;j<viz_meshnr_y;j++) {
            viz_mesh[i][j].resize(viz_meshnr_z);
            num_mesh[i][j].resize(viz_meshnr_z);
            for(int k = 0; k<viz_meshnr_z;k++) {
                viz_mesh[i][j][k]=0.0;
                num_mesh[i][j][k]=0;
            }
        }
    }    
}

void Visualisation::Init_node_numbers(GraphKMC* graph)
{
    for(int it = 0; it != graph->Numberofnodes(); it++) 
    {
        Node* node = graph->GetNode(it);
        votca::tools::vec position = node->position();
        if(node->type() == (int) NormalNode) {
            int mesh_pos_x = floor(position.x()/viz_size_x);
            int mesh_pos_y = floor(position.y()/viz_size_y);
            int mesh_pos_z = floor(position.z()/viz_size_z);
            if (mesh_pos_x == viz_meshnr_x) mesh_pos_x--; 
            if (mesh_pos_y == viz_meshnr_y) mesh_pos_y--; 
            if (mesh_pos_z == viz_meshnr_z) mesh_pos_z--; 
            
            num_mesh[mesh_pos_x][mesh_pos_y][mesh_pos_z]++;
            layer_num[mesh_pos_z]++;
        }
    }   
}

}}

#endif
