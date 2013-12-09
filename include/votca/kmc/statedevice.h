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
//include <votca/kmc/carrier.h>
//#include <votca/kmc/graphlattice.h>
//#include <votca/kmc/globaleventinfo.h>
//#include <votca/kmc/bsumtree.h>

namespace votca { namespace kmc {

template<class TGraph>    
class StateDevice : public State<TGraph> {
    
public:

    /// clear the carriers vector and the carrier reservoir and set length of both vectors equal to zero
    void InitStateDevice();
    
    /// clear reservoir
    void InitReservoir() {carrier_reservoir.clear();}
    
    /// Buying/Selling of carrier numbers from the reservoir
    unsigned int Buy(int growsize);
    void Sell(unsigned int remove_from_sim_box);
    
    /// Growing of carrier and reservoir vector
    void Grow(unsigned int nr_new_carriers);
    
    /// Print carrier list (for debugging)
    void PrintDevice(std::ostream& out);
    
private:

    vector<int> carrier_reservoir;
    
};

template <class TGraph>
void StateDevice<TGraph>::InitStateDevice(){
    this->InitState();
    carrier_reservoir.clear();
}

template <class TGraph>
unsigned int StateDevice<TGraph>::Buy(int growsize) {
    
    if(carrier_reservoir.size()==0) {Grow(growsize);}
    unsigned int carriernr_to_sim_box = carrier_reservoir.back();
    carrier_reservoir.pop_back();
    Carrier* newcarrier = this->GetCarrier(carriernr_to_sim_box);
    newcarrier->SetInBox(true);
    std::cout << newcarrier->inbox() << " a" << endl;
    std::cout << (this->GetCarrier(8))->inbox() << " a" << endl;
    std::cout << newcarrier->type() << endl;
    return carriernr_to_sim_box;
}

template <class TGraph>
void StateDevice<TGraph>::Sell(unsigned int remove_from_sim_box) {
    
    carrier_reservoir.push_back(remove_from_sim_box);
    this->GetCarrier(remove_from_sim_box)->SetInBox(false);
}

template <class TGraph>
void StateDevice<TGraph>::Grow(unsigned int nr_new_carriers) {
    
    unsigned int new_nr_carriers = this->GetCarrierSize() + nr_new_carriers;
    for (unsigned int i=this->GetCarrierSize(); i<new_nr_carriers; i++) {

        Carrier* newcarrier = this->AddCarrier(i);         

        carrier_reservoir.push_back(i);
        newcarrier->SetInBox(false);
        newcarrier->SetDistance(votca::tools::vec(0.0,0.0,0.0)); //initialize the travelled distance vector
        newcarrier->SetCarrierType((int) Reservoir); //set carrier in reservoir
    }
}

template <class TGraph>
void StateDevice<TGraph>::PrintDevice(std::ostream& out) {
    
    this->Print(out);
    std::cout << endl;
    std::cout << "reservoir indices: ";
    typename std::vector<int>::iterator it;   
    for(it = carrier_reservoir.begin(); it != carrier_reservoir.end(); it++) { 

        std::cout << (*it) << " ";
    }
    std::cout << endl;
}

}} 

#endif

