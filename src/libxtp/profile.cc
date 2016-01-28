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

#include "votca/xtp/profile.h"


using namespace std;

namespace votca {
    namespace xtp {

void Profile::Initialize_storage_arrays(Eventinfo* eventinfo) {
    
    _layersize = (eventinfo->simboxsize.z()-eventinfo->left_electrode_distance - eventinfo->right_electrode_distance)/eventinfo->number_of_layers;
 
    _positional_sum.resize(eventinfo->number_of_layers);
    _number_of_nodes.resize(eventinfo->number_of_layers); 

}

void Profile::Add_nodes_to_profile(GraphKMC* graph){
    for(int i =0; i<graph->Numberofnodes(); i++){
        NodeDevice* node = graph->GetNode(i);
        if(node->type() == (int) NormalNode) {
            votca::tools::vec position = node->position();
            Add_node_to_layer(position.z(), node->layer());
        }
    }
}

void Profile::Add_node_to_layer(double posz, int layer){
    _number_of_nodes[layer]++;
    _positional_sum[layer] += posz;
}

void Profile::Calculate_positional_average(Eventinfo* eventinfo){

    _positional_average.clear();
    
    for(int i = 0; i<eventinfo->number_of_layers;i++){
        double position = eventinfo->left_electrode_distance + 1.0*i*eventinfo->lat_const; 
        _positional_average.push_back(position);
        if(_number_of_nodes[i] != 0) _empty_layer.push_back(false);
        else _empty_layer.push_back(true);
    }
}

void Profile::Calculate_layer_boundaries(Eventinfo* eventinfo){
    double boundary = eventinfo->left_electrode_distance;
    for(int i = 0;i<eventinfo->number_of_layers+1;i++) { // one more boundary then there are slabs
        _layer_boundaries.push_back(boundary);
        boundary += _layersize;
    }
}


        
    }
}
