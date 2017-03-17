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

#ifndef __VOTCA_KMC_STATERESERVOIR_H_
#define __VOTCA_KMC_STATERESERVOIR_H_

#include <vector>
#include <list>
#include <votca/xtp/graphkmc.h>
#include <votca/xtp/nodedevice.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/vec.h>
#include <votca/xtp/state.h>
#include <votca/xtp/carrierbulk.h>
#include <votca/xtp/bsumtree.h>
#include <votca/tools/random2.h>
#include <votca/xtp/state.h>

namespace votca { namespace xtp {

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
    
    void Init_trajectory(string filename);
    
    void Print_trajectory(double simtime);
  
private:

    /// Random injection of one carrier type
    void Random_injection_one_type(int nr_charges, int carrier_type, Bsumtree* site_inject_probs, GraphKMC* graph, Eventinfo* eventinfo, votca::tools::Random2 *RandomVariable);
    
    /// Add carrier of given carrier type to chosen node
    void Add_carrier_to_chosen_node(NodeDevice* chosen_node, int carrier_type, Eventinfo* eventinfo);    

    void Print_header();
    
    vector<int> carrier_reservoir;
    
    ofstream traj_stream;
    char traj_file[100];
};



}} 

#endif

