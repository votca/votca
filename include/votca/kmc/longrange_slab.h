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

#ifndef __VOTCA_KMC_LONGRANGE_SLAB_H_
#define __VOTCA_KMC_LONGRANGE_SLAB_H_

#include <votca/tools/vec.h>
#include <votca/kmc/eventinfo.h>
#include <votca/kmc/profile.h>

namespace votca { namespace kmc {
  
using namespace std;

class Longrange_slab : public Longrange 
{
    
public:

    Longrange_slab(GraphKMC* graph, Eventinfo* eventinfo) : Longrange(graph, eventinfo){};

    Longrange_slab() : Longrange(){};    
    /// Initialize the longrange class: -determine which layers are contributing to which layers -precalculate all cut-out disc contributions
    void Initialize(GraphKMC* graph, Eventinfo* eventinfo);
    void Initialize_node(NodeDevice* node, Eventinfo* eventinfo);

private:
  
};

void Longrange_slab::Initialize (GraphKMC* graph, Eventinfo* eventinfo) {

    for(int it = 0; it < graph->Numberofnodes(); it++) {
        if(graph->GetNode(it)->type() == (int) NormalNode) {
            this->Initialize_node(graph->GetNode(it),eventinfo);
        }
        _longrange_cache.push_back(0.0);
    }    

    for(int ilayer=0; ilayer<eventinfo->number_of_layers; ilayer++) {
        _layercharge.push_back(0.0);
        _average_longrange_cache.push_back(0.0);
    }
    
}

void Longrange_slab::Initialize_node (NodeDevice* node, Eventinfo* eventinfo) {

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

    // the correction factor on the double counting correction can also be precalculated    
    
    start_index = 0;
    bool startfound = false;
    while (!startfound) {
        double start_layerbound = this->boundary(start_index);
        if((nodeposz-start_layerbound)<=eventinfo->coulomb_cut_off_radius) {
            startfound = true;
            node->setfirstcontribboundary(start_index);
        }
        else {
            start_index++;
        }
    }
    
    final_index = eventinfo->number_of_layers;
    bool finalfound = false;
    while (!finalfound) {
        double final_layerbound = this->boundary(final_index);
        if((final_layerbound-nodeposz)<=eventinfo->coulomb_cut_off_radius) {
            finalfound = true;
            node->setfinalcontribboundary(final_index);
        }
        else {
            final_index--;
        }
    }    

    node->disc_correct_clear();

    for (int j=start_index; j<=final_index; j++) {
        double disc_correct = Calculate_disc_contrib(nodeposz,j,eventinfo);
        node->disc_correct_set(disc_correct);
    }
    
}


}}

#endif
