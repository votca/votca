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

#ifndef __VOTCA_KMC_LONGRANGE_H_
#define __VOTCA_KMC_LONGRANGE_H_

#include <votca/tools/vec.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/profile.h>

namespace votca { namespace kmc {
  
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
    inline double Calculate_disc_contrib(int calculate_layer, int contrib_layer, Eventinfo* eventinfo);
    inline double Calculate_disc_contrib_slab_node(NodeDevice* node, int contrib_layer, Eventinfo* eventinfo);

private:

    vector<double> _layercharge;
    vector<double> _longrange_cache;
    vector<double> _average_longrange_cache;
    
    vector< vector <double> > _precalculate_disc_contrib; // Precalculated disc contributions
  
    vector<int> _first_contributing_layer; // What is the first layer that contributes to the relevant layer?
    vector<int> _final_contributing_layer; // What is the last layer that contributes to the relevant layer?*/
  
};

void Longrange::Init_Load_State(StateReservoir* state, Eventinfo* eventinfo) {
    for(int icar =0; icar<state->GetCarrierSize(); icar++) {
        CarrierBulk* carrier = state->GetCarrier(icar);
        Node* node = carrier->node();
        if(carrier->type() == (int) Electron) {
            Add_charge(-1.0,dynamic_cast<NodeDevice*>(node)->layer());
        }
        else if(carrier->type() == (int) Hole) {
            Add_charge(1.0,dynamic_cast<NodeDevice*>(node)->layer());
        }
    }
    Update_cache(eventinfo);
}

void Longrange::Update_cache(Eventinfo* eventinfo) {
    for (int i=0; i<eventinfo->number_of_layers; i++) {
        if(!this->emptylayer(i)) {_longrange_cache[i] = Calculate_longrange(i,true, eventinfo);}
    }
}

void Longrange::Update_cache_slab(GraphKMC* graph, Eventinfo* eventinfo) {

    for(int ilayer=0; ilayer<eventinfo->number_of_layers; ilayer++) {
        _average_longrange_cache[ilayer] = 0.0;
    }    

    for (int i=0; i<graph->Numberofnodes(); i++) {
        Node* node = graph->GetNode(i);
        votca::tools::vec node_pos = node->position();
        if(node->type() == (int) NormalNode) {
            _longrange_cache[i] = Calculate_longrange_slab(node, node_pos.z(), eventinfo->simboxsize.z()-node_pos.z(),true, eventinfo);
        }
        
        int node_layer = dynamic_cast<NodeDevice*>(node)->layer();
        _average_longrange_cache[node_layer] += _longrange_cache[i];
    }
    
    for(int ilayer=0; ilayer<eventinfo->number_of_layers; ilayer++) {
        _average_longrange_cache[ilayer] /= this->number_of_nodes(ilayer);
    }    
}

double Longrange::Get_cached_longrange(int layer) {
    return _longrange_cache[layer];
}

double Longrange::Get_cached_longrange_slab(int node_index) {
    return _longrange_cache[node_index];
}

double Longrange::Get_layer_averaged_cached_longrange_slab(int layer) {
    return _average_longrange_cache[layer];
}

double Longrange::Get_cached_density(int layer, Eventinfo* eventinfo) {
//    return _layercharge[layer]/(this->number_of_nodes(layer));
      return _layercharge[layer]/(this->layersize()*eventinfo->simboxsize.x()*eventinfo->simboxsize.y());
}

void Longrange::Reset(Eventinfo* eventinfo) {
    for (int i=0; i<eventinfo->number_of_layers; i++) {
        _layercharge[i] = 0.0;
        _longrange_cache[i] = 0.0;
    }
}

void Longrange::Reset_slab(GraphKMC* graph, Eventinfo* eventinfo) {
    for (int i=0; i<eventinfo->number_of_layers; i++) {
        _layercharge[i] = 0.0;
    }
    for (int i=0; i<graph->Numberofnodes(); i++){
        _longrange_cache[i] = 0.0;
    }
}

void Longrange::Initialize (Eventinfo* eventinfo) {

    // the properties of the profile object are initialized in the constructor of the profile object itself
    
    _first_contributing_layer.clear();
    _final_contributing_layer.clear();
    
    
    for (int ilayer=0;ilayer<eventinfo->number_of_layers;ilayer++) {

        // define for every layer, how many other layers are within the coulomb cut off radius from this layer
        double define_layerpos = this->position(ilayer);

        int start_index = 0;
        bool startfound = false;
        while (!startfound) {
            while(this->emptylayer(start_index)) start_index++;
            double start_layerpos = this->position(start_index);
            if((define_layerpos-start_layerpos)<=eventinfo->coulomb_cut_off_radius) {
                startfound = true;
                _first_contributing_layer.push_back(start_index);
            }
            start_index++;
        }

        int final_index = eventinfo->number_of_layers-1;
        bool finalfound = false;
        while (!finalfound) {
            while(this->emptylayer(final_index)) final_index--;
            double final_layerpos = this->position(final_index);
            if((final_layerpos-define_layerpos)<=eventinfo->coulomb_cut_off_radius) {
                finalfound = true;
                _final_contributing_layer.push_back(final_index);
            }
            final_index--;
        }

    }

    _precalculate_disc_contrib.resize(eventinfo->number_of_layers);

    for(int ilayer=0; ilayer<eventinfo->number_of_layers; ilayer++) {
        _layercharge.push_back(0.0);
        _longrange_cache.push_back(0.0);   

        if(this->emptylayer(ilayer)) {_precalculate_disc_contrib[ilayer].clear(); }
        else {
            int first_layer = _first_contributing_layer[ilayer];
            int final_layer = _final_contributing_layer[ilayer];

            int number_of_contributing_layers = final_layer - first_layer + 1;        
            _precalculate_disc_contrib[ilayer].resize(number_of_contributing_layers);

            for (int j=first_layer; j<=final_layer; j++) {
                if(this->emptylayer(j)) { _precalculate_disc_contrib[ilayer][j-first_layer] = 0.0;} // no nodes in this layer
                else {_precalculate_disc_contrib[ilayer][j-first_layer] = Calculate_disc_contrib(ilayer,j,eventinfo);}
            }
        }
    }
 
}

void Longrange::Initialize_slab (GraphKMC* graph, Eventinfo* eventinfo) {

    for(int it = 0; it < graph->Numberofnodes(); it++) {
        if(graph->GetNode(it)->type() == (int) NormalNode) {
            this->Initialize_slab_node(graph->GetNode(it),eventinfo);
        }
        _longrange_cache.push_back(0.0);
    }    

    for(int ilayer=0; ilayer<eventinfo->number_of_layers; ilayer++) {
        _layercharge.push_back(0.0);
        _average_longrange_cache.push_back(0.0);
    }
    
}

void Longrange::Initialize_slab_node (NodeDevice* node, Eventinfo* eventinfo) {

    // the properties of the profile object are initialized in the constructor of the profile object itself
    
    // define for every layer, how many other layers are within the coulomb cut off radius from this layer
    votca::tools::vec define_nodepos = node->position();
    double nodeposz = define_nodepos.z();
    int start_index;
    int final_index;
    
    if(node->layer() == 0) {
       node->setfirstcontriblayer(0.0);
       start_index = 0;
    }
    else {   
        start_index = 0;
        bool startfound = false;
        while (!startfound) {
            while(this->emptylayer(start_index)) start_index++;
            double start_layerpos = this->position(start_index) + 0.5*this->layersize();
            if((nodeposz-start_layerpos)<=eventinfo->coulomb_cut_off_radius) {
                startfound = true;
                node->setfirstcontriblayer(start_index);
            }
            else {
                start_index++;
            }
        }
    }

    if(node->layer() == eventinfo->number_of_layers-1) {
       node->setfinalcontriblayer(eventinfo->number_of_layers-1);
       final_index = eventinfo->number_of_layers-1;
    }
    else {   
        final_index = eventinfo->number_of_layers-1;
        bool finalfound = false;
        while (!finalfound) {
            while(this->emptylayer(final_index)) final_index--;
            double final_layerpos = this->position(final_index)-0.5*this->layersize();
            if((final_layerpos-nodeposz)<=eventinfo->coulomb_cut_off_radius) {
                finalfound = true;
                node->setfinalcontriblayer(final_index);
            }
            else {
                final_index--;
            }
        }
    }

    node->disc_coul_clear();
    
    for (int j=start_index; j<=final_index; j++) {
        if(this->emptylayer(j)) {node->disc_coul_set(0.0);} // no nodes in this layer
        else {
            double disc_contrib = Calculate_disc_contrib_slab_node(node,j,eventinfo);
            node->disc_coul_set(disc_contrib);
        }
    }

}

inline double Longrange::Calculate_disc_contrib_slab_node(NodeDevice* node, int contrib_layer, Eventinfo* eventinfo) {

    votca::tools::vec nodepos = node->position();
    double calcpos = nodepos.z();
    int node_layer = node->layer();
    
    double RC = eventinfo->coulomb_cut_off_radius;
    
    double first_contrib_pos;
    first_contrib_pos = this->position(contrib_layer) - 0.5*this->layersize();
    if(calcpos-first_contrib_pos > RC) first_contrib_pos = calcpos - RC;
    double second_contrib_pos;
    second_contrib_pos = this->position(contrib_layer) + 0.5*this->layersize();
    if(second_contrib_pos - calcpos > RC) second_contrib_pos = calcpos + RC;
 
    // determination of relative distances
    
    bool firstleft = true;
    bool secondleft = true;
    
    double firstrdist = first_contrib_pos-calcpos;
    double secondrdist = second_contrib_pos-calcpos;
     
    double direct_contrib;
    direct_contrib = RC*(secondrdist-firstrdist) - 0.5*secondrdist*fabs(secondrdist) + 0.5*firstrdist*fabs(firstrdist);
    
    double L = eventinfo->simboxsize.z();
    double mirror_contrib = 0.0;
   
    for (long i=0; i<eventinfo->number_long_range_images ; i++) {
        
        // Calculate contribution from images
        double dist1;
        double dist2;
        
        long sign;
        double mirror1_firstrdist;
        double mirror1_secondrdist;
        double mirror2_firstrdist;
        double mirror2_secondrdist;

        if (ldiv(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = -1;
            dist1 = 1.0*i*L + 2.0*calcpos;         mirror1_firstrdist = firstrdist;       mirror1_secondrdist = secondrdist;
            dist2 = 1.0*i*L + 2.0*L - 2.0*calcpos; mirror2_firstrdist = -1.0*secondrdist; mirror2_secondrdist = -1.0*firstrdist;
        }
        else { // odd generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = 1; 
            dist1 = (i+1)*L;  mirror1_firstrdist = firstrdist;       mirror1_secondrdist = secondrdist;
            dist2 = (i+1)*L;  mirror2_firstrdist = -1.0*secondrdist; mirror2_secondrdist = -1.0*firstrdist;
        }
        mirror_contrib += sign*(pow((RC*RC + 2.0*dist1*mirror1_secondrdist+dist1*dist1),(3.0/2.0))/(3.0*dist1)
                               -pow((RC*RC + 2.0*dist1*mirror1_firstrdist+dist1*dist1),(3.0/2.0))/(3.0*dist1)
                                - 0.5*(mirror1_secondrdist*mirror1_secondrdist) - dist1*mirror1_secondrdist 
                                + 0.5*(mirror1_firstrdist*mirror1_firstrdist)   + dist1*mirror1_firstrdist);
        mirror_contrib += sign*(pow((RC*RC + 2.0*dist2*mirror2_secondrdist+dist2*dist2),(3.0/2.0))/(3.0*dist2)
                               -pow((RC*RC + 2.0*dist2*mirror2_firstrdist+dist2*dist2),(3.0/2.0))/(3.0*dist2)
                                - 0.5*(mirror2_secondrdist*mirror2_secondrdist) - dist2*mirror2_secondrdist 
                                + 0.5*(mirror2_firstrdist*mirror2_firstrdist)   + dist2*mirror2_firstrdist);
    }
    
    double contrib = direct_contrib + mirror_contrib;
    return contrib;
}

inline double Longrange::Calculate_disc_contrib(int calculate_layer, int contrib_layer, Eventinfo* eventinfo) {
 
    double calcpos = this->position(calculate_layer);
    double contribpos = this->position(contrib_layer);
    double rdist = contribpos-calcpos;
  
    double contrib = eventinfo->coulomb_cut_off_radius-fabs(rdist); // Direct contribution (no image), factor 2 pi is missing, included in calculate_longrange   
    double radiussqr = eventinfo->coulomb_cut_off_radius*eventinfo->coulomb_cut_off_radius-rdist*rdist; //radius of contributing disc
    
    double L = eventinfo->simboxsize.z();
    
    for (long i=0; i<eventinfo->number_long_range_images ; i++) {
   
        // Calculate contribution from images
        double dist1;
        double dist2;
        double sign;
        if (ldiv(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = -1.0;
            dist1 = i*L + 2*contribpos - rdist;
            dist2 = i*L + 2*L - 2*contribpos + rdist;
        }
        else { // odd generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = 1; 
            dist1 = (i+1)*L + rdist;
            dist2 = (i+1)*L - rdist;
        }
        double diag1 = sqrt(radiussqr+dist1*dist1);
        double diag2 = sqrt(radiussqr+dist2*dist2);
        contrib += sign*(diag1-fabs(dist1)+diag2-fabs(dist2));
    }
    return contrib;
}

double Longrange::Calculate_longrange_slab(Node* node, double left_node_distance, double right_node_distance, bool cut_out_discs,Eventinfo* eventinfo) {
    // Potential is expressed in multiples of e/(4*pi*epsilon) (with e the elementary charge>0)
    const double Pi = 3.14159265358979323846264338327950288419716939937510;       

    double slab_contrib1 = 0.0;
    double cut_out_contrib = 0.0;
    double longrangeslab = 0;

    int start_index = dynamic_cast<NodeDevice*>(node)->firstcontriblayer();     
    int layer = dynamic_cast<NodeDevice*>(node)->layer();
    
    for(int i=0; i<layer; i++) {
        if(!this->emptylayer(i)) {
//            double charge_i = 1.0*(_layercharge[i])/(this->number_of_nodes(i));
            double charge_i = 1.0*(_layercharge[i])/(this->layersize()*eventinfo->simboxsize.x()*eventinfo->simboxsize.y());

            double position_i = 1.0*this->position(i);
            slab_contrib1 += position_i*charge_i*this->layersize(); // potential of a charged plate between two electrodes

            // calculation of contribution of disc
            if (i>=start_index) {
                // Cut out short-range sphere (relative to first contributing layer)
                cut_out_contrib += charge_i*dynamic_cast<NodeDevice*>(node)->contrib(i-start_index);
            }
        }
    }

    double slab_contrib2 = 0.0;
    
    for(int i=layer+1; i<eventinfo->number_of_layers; i++) {
        if(!this->emptylayer(i)) {
//            double charge_i = 1.0*_layercharge[i]/(this->number_of_nodes(i));
            double charge_i = 1.0*_layercharge[i]/(this->layersize()*eventinfo->simboxsize.x()*eventinfo->simboxsize.y());
            double rel_position_i = 1.0*(eventinfo->simboxsize.z()-this->position(i));
            slab_contrib2 += rel_position_i*charge_i*this->layersize(); // potential of a charged plate between two electrodes

            int final_index = dynamic_cast<NodeDevice*>(node)->finalcontriblayer();
            
            if (final_index >= i) {
                // Cut out short-range sphere (relative to first contributing layer)
                cut_out_contrib += charge_i*dynamic_cast<NodeDevice*>(node)->contrib(i-start_index);
            }
        }
    }
    
    double slab_contrib3 = 0.0;
//    double charge_i = 1.0*_layercharge[layer]/(this->number_of_nodes(layer));
    double charge_i = 1.0*_layercharge[layer]/(this->layersize()*eventinfo->simboxsize.x()*eventinfo->simboxsize.y());
    double position_i = 1.0*this->position(layer);
    slab_contrib3 += 0.5*position_i*charge_i*this->layersize(); // potential of a charged plate between two electrodes

    double slab_contrib4 = 0.0;
    double rel_position_i = 1.0*(eventinfo->simboxsize.z()-this->position(layer));
    slab_contrib4 += 0.5*rel_position_i*charge_i*this->layersize(); // potential of a charged plate between two electrodes
    
    cut_out_contrib += charge_i*dynamic_cast<NodeDevice*>(node)->contrib(layer-start_index);    
    if (!cut_out_discs) { cut_out_contrib = 0.0; }    
    
    //note that local positioning in the slab itself is calculated on the fly
   
    longrangeslab += 4.0*Pi*(slab_contrib1*(right_node_distance/eventinfo->simboxsize.z()) + slab_contrib2*(left_node_distance/eventinfo->simboxsize.z()));
    longrangeslab += 4.0*Pi*(slab_contrib3*(right_node_distance/eventinfo->simboxsize.z()) + slab_contrib4*(left_node_distance/eventinfo->simboxsize.z()));
    longrangeslab += -2.0*Pi*charge_i*(left_node_distance - this->position(layer))*(left_node_distance - this->position(layer));
    longrangeslab += -2.0*Pi*charge_i*0.25*this->layersize()*this->layersize(); // d = 0.5 layersize, therefore the 0.25
    longrangeslab += -2.0*Pi*cut_out_contrib;
    
    return longrangeslab;

}

double Longrange::Calculate_longrange(int layer, bool cut_out_discs,Eventinfo* eventinfo) {
    // Potential is expressed in multiples of e/(4*pi*epsilon) (with e the elementary charge>0)
    const double Pi = 3.14159265358979323846264338327950288419716939937510;       

    double plate_contrib1 = 0.0;
    double disc_contrib = 0.0;
    double layerpos = this->position(layer);
    
    for(int i=0; i<layer; i++) {
        if(!this->emptylayer(i)) {
            double charge_i = 1.0*_layercharge[i]/(this->number_of_nodes(i));
            double position_i = 1.0*this->position(i);
            plate_contrib1 += position_i*charge_i; // potential of a charged plate between two electrodes

            // calculation of contribution of disc
            double distance = layerpos -position_i;
            if (distance<=eventinfo->coulomb_cut_off_radius) {
                // Cut out short-range sphere (relative to first contributing layer)
                int first_layer = _first_contributing_layer[layer];
                disc_contrib -= charge_i*_precalculate_disc_contrib[layer][i-first_layer];
            }
        }
    }
    
    double plate_contrib2 = 0.0;
    
    for(int i=layer; i<eventinfo->number_of_layers; i++) {
        if(!this->emptylayer(i)) {
            double charge_i = 1.0*_layercharge[i]/(this->number_of_nodes(i));
            double rel_position_i = 1.0*(eventinfo->simboxsize.z()-this->position(i));
            plate_contrib2 += rel_position_i*charge_i; // potential of a charged plate between two electrodes

            // calculation of contribution of disc        
            double distance = this->position(i)-layerpos;

            if (distance<=eventinfo->coulomb_cut_off_radius) {
                // Cut out short-range sphere (relative to first contributing layer)
                int first_layer = _first_contributing_layer[layer];
                disc_contrib -= charge_i*_precalculate_disc_contrib[layer][i-first_layer];
            }
        }
    }
    if (!cut_out_discs) { disc_contrib = 0.0; }
    
    return 4.0*Pi*(plate_contrib1*(1-layerpos/eventinfo->simboxsize.z()) + plate_contrib2*(layerpos/eventinfo->simboxsize.z())) + 2.0*Pi*disc_contrib;

}

}}

#endif
