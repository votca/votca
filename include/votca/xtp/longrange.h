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

#ifndef __VOTCA_KMC_LONGRANGE_H_
#define __VOTCA_KMC_LONGRANGE_H_

#include <votca/tools/vec.h>
#include <votca/xtp/eventinfo.h>
#include <votca/xtp/profile.h>
#include <votca/xtp/statereservoir.h>

namespace votca { namespace xtp {
  
using namespace std;

class Longrange : public Profile 
{
    
public:

    Longrange(GraphKMC* graph, Eventinfo* eventinfo) : Profile(graph, eventinfo){};

    Longrange() : Profile(){};    

    /// Initialize longrange profile after state is read in
    void Init_Load_State(StateReservoir* state, Eventinfo* eventinfo);    
    
    /// Add charge to longrange object
    void Add_charge(double charge, int layer) {_layercharge[layer] += charge;}

    /// Update longrange coulomb potential cache
    void Update_cache(Eventinfo* eventinfo);
    void Update_cache_slab(GraphKMC* graph, Eventinfo* eventinfo);
    
    /// Get longrange coulomb potential cache
    double Get_cached_longrange(int layer);
    double Get_cached_longrange_slab(int node_index);
    double Get_cached_density(int layer,  Eventinfo* eventinfo);
    double Get_layer_averaged_cached_longrange_slab(int layer);
    
    /// Reser longrange coulomb potential cache
    void Reset(Eventinfo* eventinfo);
    void Reset_slab(GraphKMC* graph, Eventinfo* eventinfo); 
    
    /// Initialize the longrange class: -determine which layers are contributing to which layers -precalculate all cut-out disc contributions
    void Initialize(Eventinfo* eventinfo);
    void Initialize_slab_node(NodeDevice* node, Eventinfo* eventinfo);
    void Initialize_slab(GraphKMC* graph, Eventinfo* eventinfo);
    
    //note that the number of images for the calculation of the long range potential should be considerably larger 
    //than the number for the short range potential
    double Calculate_longrange(int layer, bool cut_out_discs, Eventinfo* eventinfo); // Calculate long-range part of Coulomb interaction
    double Calculate_longrange_slab(Node* node, double left_node_distance, double right_node_distance, bool cut_out_discs,Eventinfo* eventinfo);

    ///precalculate the coulombic contributions from the cut-out discs
    double Calculate_disc_contrib(double calcpos, int contrib_layer, Eventinfo* eventinfo);
    double Calculate_disc_contrib_slab_node(NodeDevice* node, int contrib_layer, Eventinfo* eventinfo);

private:

    vector<double> _layercharge;
    vector<double> _longrange_cache;
    vector<double> _average_longrange_cache;
    
    vector< vector <double> > _precalculate_disc_contrib; // Precalculated disc contributions
  
    vector<int> _first_contributing_layer; // What is the first layer that contributes to the relevant layer?
    vector<int> _final_contributing_layer; // What is the last layer that contributes to the relevant layer?*/
  
};


}}

#endif
