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

#ifndef __VOTCA_KMC_PROFILE_H_
#define __VOTCA_KMC_PROFILE_H_

#include <votca/tools/vec.h>
#include <votca/kmc/graphdevice.h>
#include <votca/kmc/eventinfo.h>

namespace votca { namespace kmc {
  
using namespace std;

class Profile {

public:
    Profile(GraphDevice* graph, Eventinfo* eventinfo) {
        Initialize_storage_arrays(eventinfo);
        Add_nodes_to_profile(graph);
        Calculate_positional_average();
        Calculate_layer_boundaries(eventinfo);
    }
     
    ~Profile(){}     
    
    const int &number_of_layers() const { return _number_of_layers; }
    const double &position(int i) const {return _positional_average[i]; }
    
    inline void Initialize_storage_arrays(Eventinfo* eventinfo);
    inline void Add_nodes_to_profile(GraphDevice* graph);
    inline void Add_node_to_layer(double posx, int layer);
    inline void Calculate_positional_average();
    inline void Calculate_layer_boundaries(Eventinfo* eventinfo);

private:

    int _number_of_layers;

    vector<double> _positional_sum;
    vector<int> _number_of_nodes;
    
    vector<bool> _empty_layer;
    vector<double> _positional_average;
    vector<double> _layer_boundaries;
  
};

inline void Profile::Initialize_storage_arrays(Eventinfo* eventinfo) {
    
    _number_of_layers = ceil(eventinfo->simboxsize.x()/eventinfo->layersize);
 
    _positional_sum.resize(_number_of_layers);
    _number_of_nodes.resize(_number_of_layers); 

}

inline void Profile::Add_nodes_to_profile(GraphDevice* graph){
    for(int i =0; i<graph->Numberofnodes(); i++){
        NodeDevice* node = graph->GetNode(i);
        votca::tools::vec position = node->position();
        Add_node_to_layer(position.x(), node->layer());
    }
}

inline void Profile::Add_node_to_layer(double posx, int layer){
    _number_of_nodes[layer]++;
    _positional_sum[layer] += posx;
};

inline void Profile::Calculate_positional_average(){

    double average = 0.0;
    for(int i = 0;i<_number_of_layers;i++) {
        if(_number_of_nodes[i] != 0) {
            average = _positional_sum[i]/(1.0*_number_of_nodes[i]);
            _positional_average.push_back(average);
            _empty_layer.push_back(false);
        }
        else {
            _positional_average.push_back(average);
            _empty_layer.push_back(true);
        }
    }
}

inline void Profile::Calculate_layer_boundaries(Eventinfo* eventinfo){
    double boundary = 0.0;
    for(int i = 0;i<_number_of_layers;i++) {
        _layer_boundaries.push_back(boundary);
        boundary += eventinfo->layersize;
    }
    _layer_boundaries.push_back(eventinfo->simboxsize.x());
}

}}

#endif /* _VOTCA_KMC_PROFILE_H */
