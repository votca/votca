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

#ifndef __VOTCA_KMC_CARRIERDEVICE_H_
#define __VOTCA_KMC_CARRIERDEVICE_H_

#include <votca/kmc/node.h>

namespace votca { namespace kmc {
    
class CarrierDevice : public Carrier {
public:

    CarrierDevice(int id) : Carrier(id){
    };
    
    /// is carrier in box or not?
    const bool &inbox() const { return _in_sim_box; }
    const double &on_site_coulomb() const {return _on_site_coulomb; }
    
    /// set "inbox/outbox" status
    void SetInBox(bool in_sim_box) {_in_sim_box = in_sim_box;}
    
    void SetCoulomb(double coulomb) {_on_site_coulomb = coulomb;}
    
private:
    
    bool _in_sim_box;
    double _on_site_coulomb;
};

}} 

#endif

