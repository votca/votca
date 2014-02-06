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

#ifndef __VOTCA_KMC_STATEDEVICE_H_
#define __VOTCA_KMC_STATEDEVICE_H_

#include <vector>
#include <list>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/kmc/state.h>
#include <votca/kmc/carrierdevice.h>
//include <votca/kmc/carrier.h>
//#include <votca/kmc/graphlattice.h>
//#include <votca/kmc/globaleventinfo.h>
#include <votca/kmc/bsumtree.h>
#include <votca/tools/random2.h>

namespace votca { namespace kmc {

class StateDevice : public State<GraphDevice, CarrierDevice> {
    
public:

    /// clear the carriers vector and the carrier reservoir and set length of both vectors equal to zero
    void InitStateDevice();
    
    /// load a state store file
    void LoadStateDevice(const char* filename, GraphDevice* graph,Eventinfo* eventinfo);
    
    /// save a state store file
    void SaveStateDevice();
    
    /// initialize interactions (set to 0)
    void InitInteractions(int pre_size, int post_size,Eventinfo* eventinfo);
    
    /// clear reservoir
    void InitReservoir() {carrier_reservoir.clear();}
    
    /// Buying/Selling of carrier numbers from the reservoir
    unsigned int Buy();
    void Sell(unsigned int remove_from_sim_box);
    
    /// Growing of carrier and reservoir vector
    void Grow(unsigned int nr_new_carriers, int maxpairdegree);
    
    /// Print carrier list (for debugging)
    void PrintDevice(std::ostream& out);
    
    /// Is the carrier_reservoir empty?
    bool ReservoirEmpty(){return carrier_reservoir.empty();}
    
    /// Random injection of charges (both carrier types)
    void Random_init_double_injection(Bsumtree* ho_site_inject_probs, Bsumtree* el_site_inject_probs,  GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
    
    /// Random injection of charges (one carrier type))
    void Random_init_injection(int carrier_type, Bsumtree* site_inject_probs, GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
  
private:

    /// Random injection of one carrier type
    void Random_injection_one_type(double density, int carrier_type, Bsumtree* site_inject_probs, GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
    
    /// Add carrier of given carrier type to chosen node
    void Add_carrier_to_chosen_node(NodeDevice* chosen_node, int carrier_type, Eventinfo* eventinfo);    
    
    vector<int> carrier_reservoir;
};

void StateDevice::InitStateDevice(){
    this->InitState();
    InitReservoir();
}

void StateDevice::LoadStateDevice(const char* filename, GraphDevice* graph,Eventinfo* eventinfo){
    
    // Initializes the Coulomb interactions
    
    int pre_carrier_size = this->GetCarrierSize();
    this->Load(filename, graph);
    int post_carrier_size = this->GetCarrierSize();

    this->InitInteractions(pre_carrier_size,post_carrier_size,eventinfo);
}

void StateDevice::InitInteractions(int pre_carrier_size,int post_carrier_size,Eventinfo* eventinfo){

    for (unsigned int i=pre_carrier_size; i < post_carrier_size; i++) {
        CarrierDevice* probecarrier = GetCarrier(i);
        probecarrier->Init_to_Coulomb(eventinfo->maxpairdegree);
        probecarrier->Set_from_Coulomb(0.0);        
    }
}


unsigned int StateDevice::Buy() {
    
    unsigned int carriernr_to_sim_box = carrier_reservoir.back();
    carrier_reservoir.pop_back();
    CarrierDevice* newcarrier = this->GetCarrier(carriernr_to_sim_box);
    newcarrier->SetInBox(true);
    return carriernr_to_sim_box;
}

void StateDevice::Sell(unsigned int remove_from_sim_box) {
    
    carrier_reservoir.push_back(remove_from_sim_box);
    this->GetCarrier(remove_from_sim_box)->SetInBox(false);
    this->GetCarrier(remove_from_sim_box)->SetCarrierType((int) Reservoir);
    this->GetCarrier(remove_from_sim_box)->Reset_to_Coulomb();
    this->GetCarrier(remove_from_sim_box)->Set_from_Coulomb(0.0);
}

void StateDevice::Grow(unsigned int nr_new_carriers, int maxpairdegree) {
    
    unsigned int new_nr_carriers = this->GetCarrierSize() + nr_new_carriers;
    for (unsigned int i=this->GetCarrierSize(); i<new_nr_carriers; i++) {

        CarrierDevice* newcarrier = this->AddCarrier(i);         
        carrier_reservoir.push_back(i);
        newcarrier->SetInBox(false);
        newcarrier->SetDistance(votca::tools::vec(0.0,0.0,0.0)); //initialize the travelled distance vector
        newcarrier->SetCarrierType((int) Reservoir); //set carrier in reservoir
        newcarrier->Init_to_Coulomb(maxpairdegree);
        newcarrier->Set_from_Coulomb(0.0);
    }
}

void StateDevice::PrintDevice(std::ostream& out) {
    
    this->Print(out);
    std::cout << endl;
    std::cout << "reservoir indices: ";
    typename std::vector<int>::iterator it;   
    for(it = carrier_reservoir.begin(); it != carrier_reservoir.end(); it++) { 

        std::cout << (*it) << " ";
    }
    std::cout << endl;
}

void StateDevice::Random_init_injection(int carrier_type, Bsumtree* site_inject_probs, GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2* RandomVariable){
   
    double density;
    
    if(carrier_type == (int) Hole)          {density = eventinfo->ho_density;}
    else if(carrier_type == (int) Electron) {density = eventinfo->el_density;}
    
    this->Random_injection_one_type(density,carrier_type, site_inject_probs,graph,eventinfo,RandomVariable);
    
}

void StateDevice::Random_init_double_injection(Bsumtree* ho_site_inject_probs, Bsumtree* el_site_inject_probs, GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable){
   
    double electron_density = eventinfo->el_density;
    double hole_density = eventinfo->ho_density;
    
    this->Random_injection_one_type(electron_density,(int) Electron, el_site_inject_probs,graph,eventinfo,RandomVariable);
    this->Random_injection_one_type(hole_density,(int) Hole, ho_site_inject_probs,graph,eventinfo,RandomVariable);
    
}

void StateDevice::Random_injection_one_type(double density, int carrier_type, Bsumtree* site_inject_probs, GraphDevice* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable){
    
    int number_of_carriers = density*graph->Numberofnodes();
    
    for (int it = 0; it < graph->Numberofnodes() ; it++) {
        if(graph->GetNode(it)->occ() == -1) {
            site_inject_probs->setrate(it,1.0); // 
        }
        else {
            site_inject_probs->setrate(it,0.0); // no double occup
        }
    }
        
    double tot_probsum;
    int carcounter = 0;
    
    while(carcounter!=number_of_carriers){    
        
        tot_probsum = site_inject_probs->compute_sum();
        
        double randn = tot_probsum*RandomVariable->rand_uniform();
        long node_ID;
        NodeDevice* chosen_node;
        
        node_ID = site_inject_probs->search(randn);
        chosen_node = graph->GetNode(node_ID);
  
        if(this->ReservoirEmpty()) this->Grow(eventinfo->growsize, eventinfo->maxpairdegree);
        this->Add_carrier_to_chosen_node(chosen_node, carrier_type, eventinfo);

        site_inject_probs->setrate(node_ID,0.0); // no double occupation
        carcounter++;        
    }
}

void StateDevice::Add_carrier_to_chosen_node(NodeDevice* chosen_node, int carrier_type, Eventinfo* eventinfo){
    
    //add carrier to node
    int carrier_ID = this->Buy();
    CarrierDevice* newcarrier = this->GetCarrier(carrier_ID);
    newcarrier->SetCarrierNode(chosen_node);  
    chosen_node->AddCarrier(carrier_ID);        

    //set carrier type
    newcarrier->SetCarrierType(carrier_type);

    //initialize distance
    newcarrier->SetDistance(votca::tools::vec(0.0,0.0,0.0));

    // New carrier is in the simulation box
    newcarrier->SetInBox(true);

    //initialize coulomb interactions
    newcarrier->Init_to_Coulomb(eventinfo->maxpairdegree);
    newcarrier->Set_from_Coulomb(0.0);
}

}} 

#endif

