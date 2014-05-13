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

#ifndef __VOTCA_KMC_STATERESERVOIR_H_
#define __VOTCA_KMC_STATERESERVOIR_H_

#include <vector>
#include <list>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/kmc/state.h>
#include <votca/kmc/carrierbulk.h>
#include <votca/kmc/bsumtree.h>
#include <votca/tools/random2.h>

namespace votca { namespace kmc {

class StateReservoir : public State<GraphKMC, CarrierBulk> {
    
public:

    /// clear the carriers vector and the carrier reservoir and set length of both vectors equal to zero
    void InitStateReservoir();
    
    /// load a state store file
    void LoadStateReservoir(const char* filename, GraphKMC* graph,Eventinfo* eventinfo);
    
    /// save a state store file
    void SaveStateReservoir();
    
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
    void Random_init_injection(int nrelectrons, int nrholes, Bsumtree* site_inject_probs, GraphKMC* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
  
private:

    /// Random injection of one carrier type
    void Random_injection_one_type(int nr_charges, int carrier_type, Bsumtree* site_inject_probs, GraphKMC* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
    
    /// Add carrier of given carrier type to chosen node
    void Add_carrier_to_chosen_node(NodeDevice* chosen_node, int carrier_type, Eventinfo* eventinfo);    
    
    vector<int> carrier_reservoir;
};

void StateReservoir::InitStateReservoir(){
    this->InitState();
    InitReservoir();
}

void StateReservoir::LoadStateReservoir(const char* filename, GraphKMC* graph,Eventinfo* eventinfo){
    
    // Initializes the Coulomb interactions
    
    int pre_carrier_size = this->GetCarrierSize();
    this->Load(filename, graph);
    int post_carrier_size = this->GetCarrierSize();

    this->InitInteractions(pre_carrier_size,post_carrier_size,eventinfo);
}

void StateReservoir::InitInteractions(int pre_carrier_size,int post_carrier_size,Eventinfo* eventinfo){

    for (unsigned int i=pre_carrier_size; i < post_carrier_size; i++) {
        CarrierBulk* probecarrier = GetCarrier(i);
        probecarrier->Init_to_Coulomb(eventinfo->maxpairdegree);
        probecarrier->Set_from_Coulomb(0.0);        
    }
}


unsigned int StateReservoir::Buy() {
    
    unsigned int carriernr_to_sim_box = carrier_reservoir.back();
    carrier_reservoir.pop_back();
    CarrierBulk* newcarrier = this->GetCarrier(carriernr_to_sim_box);
    newcarrier->SetInBox(true);
    return carriernr_to_sim_box;
}

void StateReservoir::Sell(unsigned int remove_from_sim_box) {
    
    carrier_reservoir.push_back(remove_from_sim_box);
    this->GetCarrier(remove_from_sim_box)->SetInBox(false);
    this->GetCarrier(remove_from_sim_box)->SetCarrierType((int) Reservoir);
    this->GetCarrier(remove_from_sim_box)->Reset_to_Coulomb();
    this->GetCarrier(remove_from_sim_box)->Set_from_Coulomb(0.0);
}

void StateReservoir::Grow(unsigned int nr_new_carriers, int maxpairdegree) {
    
    unsigned int new_nr_carriers = this->GetCarrierSize() + nr_new_carriers;
    for (unsigned int i=this->GetCarrierSize(); i<new_nr_carriers; i++) {

        CarrierBulk* newcarrier = this->AddCarrier(i);         
        carrier_reservoir.push_back(i);
        newcarrier->SetInBox(false);
        newcarrier->SetDistance(votca::tools::vec(0.0,0.0,0.0)); //initialize the travelled distance vector
        newcarrier->SetCarrierType((int) Reservoir); //set carrier in reservoir
        newcarrier->Init_to_Coulomb(maxpairdegree);
        newcarrier->Set_from_Coulomb(0.0);
    }
}

void StateReservoir::PrintDevice(std::ostream& out) {
    
    this->Print(out);
    std::cout << endl;
    std::cout << "reservoir indices: ";
    typename std::vector<int>::iterator it;   
    for(it = carrier_reservoir.begin(); it != carrier_reservoir.end(); it++) { 

        std::cout << (*it) << " ";
    }
    std::cout << endl;
}

void StateReservoir::Random_init_injection(int nr_electrons, int nr_holes, Bsumtree* site_inject_probs, GraphKMC* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable){

    if(nr_holes != 0 ) this->Random_injection_one_type(nr_holes,(int) Hole, site_inject_probs,graph,eventinfo,RandomVariable);
    if(nr_electrons != 0) this->Random_injection_one_type(nr_electrons,(int) Electron, site_inject_probs,graph,eventinfo,RandomVariable);
    
}

void StateReservoir::Random_injection_one_type(int nr_charges, int carrier_type, Bsumtree* site_inject_probs, GraphKMC* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable){
    
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
    
    while(carcounter!=nr_charges){    
        
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

void StateReservoir::Add_carrier_to_chosen_node(NodeDevice* chosen_node, int carrier_type, Eventinfo* eventinfo){
    
    //add carrier to node
    int carrier_ID = this->Buy();
    CarrierBulk* newcarrier = this->GetCarrier(carrier_ID);
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

