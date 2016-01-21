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

#ifndef __VOTCA_KMC_PROFILE_H_
#define __VOTCA_KMC_PROFILE_H_

#include <votca/tools/vec.h>
#include <votca/xtp/graphkmc.h>
#include <votca/xtp/eventinfo.h>

namespace votca { namespace xtp {
  
using namespace std;

class Profile {

public:
    Profile(GraphKMC* graph, Eventinfo* eventinfo) {
        Initialize_storage_arrays(eventinfo);
        Add_nodes_to_profile(graph);
        Calculate_positional_average(eventinfo);
        Calculate_layer_boundaries(eventinfo);
    }
    
    Profile() {}
     
    ~Profile(){}     
    
    /// number of profile layers
    const double &layersize() const { return _layersize; }
    
    /// number of nodes in a layer
    const int &number_of_nodes(int i) const { return _number_of_nodes[i];}
    
    /// positional average of the nodes in the layer
    const double &position(int i) const {return _positional_average[i]; }

    /// left boundary of slab i (number_of_layers is right boundary of last slab)
    const double &boundary(int i) const {return _layer_boundaries[i]; }
    
    bool emptylayer(int i) const {return _empty_layer[i]; }
    
    /// Initialize arrays
    void Initialize_storage_arrays(Eventinfo* eventinfo);
    
    /// Read out positions of nodes from graph object
    void Add_nodes_to_profile(GraphKMC* graph);
    
    /// Read out positions per layer
    void Add_node_to_layer(double posz, int layer);
    
    /// Calculate positional average of all layers
    void Calculate_positional_average(Eventinfo* eventinfo);
    
    /// Calculate boundaries between layers
    void Calculate_layer_boundaries(Eventinfo* eventinfo);

private:

    double _layersize;

    vector<double> _positional_sum;
    vector<int> _number_of_nodes;
    
    vector<bool> _empty_layer;
    vector<double> _positional_average;
    vector<double> _layer_boundaries;
  
};


}}

#endif /* _VOTCA_KMC_PROFILE_H */
