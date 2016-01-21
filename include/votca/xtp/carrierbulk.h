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

#ifndef __VOTCA_KMC_CARRIERBULK_H_
#define __VOTCA_KMC_CARRIERBULK_H_

#include <votca/xtp/node.h>

namespace votca { namespace xtp {
    
class CarrierBulk : public Carrier {
public:

    CarrierBulk(int id) : Carrier(id){
    };
    
    /// is carrier in box or not?
    const bool &inbox() const { return _in_sim_box; }
    const double &on_site_coulomb() const {return _from_coulomb; }
    const double &to_site_coulomb(int linkid) const {return _to_coulomb[linkid];}
    
    /// set "inbox/outbox" status
    void SetInBox(bool in_sim_box) {_in_sim_box = in_sim_box;}
    
    void Add_from_Coulomb(double coulomb) {_from_coulomb += coulomb;}
    void Set_from_Coulomb(double coulomb) {_from_coulomb = coulomb;}
    
    void Set_on_node(double sim_time) {_on_node = sim_time;}
    
    void Init_to_Coulomb(int maxpairdegree) {_to_coulomb.resize(maxpairdegree); Reset_to_Coulomb();}
    void Reset_to_Coulomb() { for (unsigned it = 0; it < _to_coulomb.size(); it++ ) { _to_coulomb[it] = 0.0;} }
    void Add_to_Coulomb(double coulomb, int linkID) {_to_coulomb[linkID] += coulomb;}
    void Set_to_Coulomb(double coulomb, int linkID) {_to_coulomb[linkID] = coulomb;}
    
   
private:
    
    bool _in_sim_box;
    double _from_coulomb;
    vector<double> _to_coulomb;
    double _on_node;
};

}} 

#endif

