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

#ifndef __VOTCA_KMC_LONGRANGE_PLANE_H_
#define __VOTCA_KMC_LONGRANGE_PLANE_H_

#include <votca/tools/vec.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/profile.h>

namespace votca { namespace kmc {
  
using namespace std;

class Longrange_plane : public Longraneg 
{
    
public:

    Longrange_plane(GraphKMC* graph, Eventinfo* eventinfo) : Longrange(graph, eventinfo){};

    Longrange_plane() : Longrange(){};    

    /// Initialize the longrange class: -determine which layers are contributing to which layers -precalculate all cut-out disc contributions
    void Initialize(Eventinfo* eventinfo);

private:

  
};

void Longrange_plane::Initialize (GraphKMC* graph,Eventinfo* eventinfo) {

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
                else {_precalculate_disc_contrib[ilayer][j-first_layer] = Calculate_disc_contrib(this->position(ilayer),j,eventinfo);}
            }
        }
    }
 
}

}}

#endif
