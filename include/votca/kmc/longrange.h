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
#include <votca/kmc/graphlattice.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/profile.h>

namespace votca { namespace kmc {
  
using namespace std;

class Longrange : public Profile 
{
    
public:

    Longrange(GraphDevice* graph, Eventinfo* eventinfo) : Profile(graph, eventinfo){
    };    

    /// Initialize longrange profile after state is read in
    void Init_Load_State(StateDevice* state, Eventinfo* eventinfo);    
    
    /// Add charge to longrange object
    void Add_charge(double charge, int layer) {_layercharge[layer] += charge;}

    /// Update longrange coulomb potential cache
    void Update_cache(Eventinfo* eventinfo);
    
    /// Get longrange coulomb potential cache
    double Get_cached_longrange(int layer);
    
    /// Reser longrange coulomb potential cache
    void Reset();
    
    /// Initialize the longrange class: -determine which layers are contributing to which layers -precalculate all cut-out disc contributions
    void Initialize(Eventinfo* eventinfo);
    //note that the number of images for the calculation of the long range potential should be considerably larger 
    //than the number for the short range potential
    double Calculate_longrange(int layer, bool cut_out_discs, Eventinfo* eventinfo); // Calculate long-range part of Coulomb interaction

    ///precalculate the coulombic contributions from the cut-out discs
    inline double Calculate_disc_contrib(int calculate_layer, int contrib_layer, Eventinfo* eventinfo);

private:

    vector<double> _layercharge;
    vector<double> _longrange_cache;
    
    vector< vector <double> > _precalculate_disc_contrib; // Precalculated disc contributions
  
    vector<int> _first_contributing_layer; // What is the first layer that contributes to the relevant layer?
    vector<int> _final_contributing_layer; // What is the last layer that contributes to the relevant layer?*/
  
};

void Longrange::Init_Load_State(StateDevice* state, Eventinfo* eventinfo) {
    for(int icar =0; icar<state->GetCarrierSize(); icar++) {
        CarrierDevice* carrier = state->GetCarrier(icar);
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
    for (int i=0; i<this->number_of_layers(); i++) {
        _longrange_cache[i] = Calculate_longrange(i,true, eventinfo);
    }
}

double Longrange::Get_cached_longrange(int layer) {
    return _longrange_cache[layer];
}

void Longrange::Reset() {
    for (int i=0; i<this->number_of_layers(); i++) {
        _layercharge[i] = 0.0;
        _longrange_cache[i] = 0.0;
    }
}

void Longrange::Initialize (Eventinfo* eventinfo) {

    for (int ilayer=0;ilayer<this->number_of_layers();ilayer++) {

        // define for every layer, how many other layers are within the coulomb cut off radius from this layer
        double define_layer = this->position(ilayer);

        int start_index = 0;
        bool startfound = false;
        while (!startfound) {
            double start_layer = this->position(start_index);
            if((define_layer-start_layer)<=eventinfo->coulcut) {
                startfound = true;
                _first_contributing_layer.push_back(start_index);
            }
            start_index++;
        }

        int final_index = this->number_of_layers()-1;
        bool finalfound = false;
        while (!finalfound) {
            double final_layer = this->position(final_index);
            if((final_layer-define_layer)<=eventinfo->coulcut) {
                finalfound = true;
                _final_contributing_layer.push_back(final_index);
            }
            final_index--;
        }

    }
  
    _precalculate_disc_contrib.resize(this->number_of_layers());
    
    for(int ilayer=0; ilayer<this->number_of_layers(); ilayer++) {
        _layercharge.push_back(0.0);
        _longrange_cache.push_back(0.0);
        int first_layer = _first_contributing_layer[ilayer];
        int final_layer = _final_contributing_layer[ilayer];

        int number_of_contributing_layers = final_layer - first_layer + 1;        
        _precalculate_disc_contrib[ilayer].resize(number_of_contributing_layers);

        for (int j=first_layer; j<=final_layer; j++) {
            _precalculate_disc_contrib[ilayer][j-first_layer] = Calculate_disc_contrib(ilayer,j,eventinfo);
        }
    }

}

inline double Longrange::Calculate_disc_contrib(int calculate_layer, int contrib_layer, Eventinfo* eventinfo) {
 
    double calcpos = this->position(calculate_layer);
    double contribpos = this->position(contrib_layer);
    double rdist = contribpos-calcpos;
  
    double contrib = eventinfo->coulcut-fabs(rdist); // Direct contribution (no image), factor 2 pi is missing, included in compute_longrange   
    double radiussqr = eventinfo->coulcut*eventinfo->coulcut-rdist*rdist; //radius of contributing disc
    
    double L = eventinfo->simboxsize.x();
    
    for (long i=0; i<eventinfo->nr_lr_images; i++) {
   
        // Calculate contribution from images
        long dist1;
        long dist2;
        long sign;
        if (ldiv(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
            sign = -1;
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


double Longrange::Calculate_longrange(int layer, bool cut_out_discs,Eventinfo* eventinfo) {
    // Potential is expressed in multiples of e/(4*pi*epsilon) (with e the elementary charge>0)
    const double Pi = 3.14159265358979323846264338327950288419716939937510;       

    double plate_contrib1 = 0.0;
    double disc_contrib = 0.0;

    for(int i=0; i<layer; i++) {
        double charge_i = 1.0*_layercharge[i];
        double position_i = 1.0*this->position(i);
        plate_contrib1 += position_i*charge_i; // potential of a charged plate between two electrodes
        double distance = layer-position_i;
        int first_layer = _first_contributing_layer[layer];
        if (distance<=eventinfo->coulcut) {
            // Cut out short-range sphere
            disc_contrib -= charge_i*_precalculate_disc_contrib[layer][i-first_layer];
        }
    }
    double plate_contrib2 = 0.0;
    for(int i=layer; i<this->number_of_layers(); i++) {
        double charge_i = 1.0*_layercharge[i];
        double rel_position_i = 1.0*(eventinfo->simboxsize.x()-this->position(i));
        plate_contrib2 += rel_position_i*charge_i; // potential of a charged plate between two electrodes
        double distance = this->position(i)-layer;
        int first_layer = _first_contributing_layer[layer];
        if (distance<=eventinfo->coulcut) {
            // Cut out short-range sphere
            disc_contrib -= charge_i*_precalculate_disc_contrib[layer][i-first_layer];
        }
    }
    if (!cut_out_discs) { disc_contrib = 0.0; }
    double layerpos = this->position(layer);
    return 4.0*Pi*(plate_contrib1*(1-layerpos/eventinfo->simboxsize.x()) + plate_contrib2*(layerpos/eventinfo->simboxsize.x()) + 0.5*disc_contrib)/(eventinfo->simboxsize.y()*eventinfo->simboxsize.z());
}

}}

#endif
