
 /* Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_KMC_NUMOUTPUT_H
#define	_VOTCA_KMC_NUMOUTPUT_H

#include <votca/kmc/event.h>

namespace votca { namespace kmc {

class Numoutput
{

public:
    
    void Initialize();
    void Initialize_equilibrate();
    void Init_convergence_check(double simtime);
    void Convergence_check(double simtime, Eventinfo* eventinfo);
    void Update(Event* event, double simtime, double timestep);
    void Write(double simtime);
    
    const bool &iv_conv() const {return _direct_iv_convergence;}
    const bool &reco_conv() const {return _direct_reco_convergence;}
    
    const int &iv_count() const {return _direct_iv_counter;}
    const int &reco_count() const {return _direct_reco_counter;}
    
private:
    
    int _nelectrons;
    int _nholes;
    int _ncarriers;
    
    int _ninjections;
    int _ncollections;
    int _nleftinjections;
    int _nrightinjections;
    int _nleftcollections;
    int _nrightcollections;
    
    int _nplaintransfer;
    int _nrecombinations;
    int _ninject_to_recombination;
    
    double _ninjectionrate;
    double _ncollectionrate;
    double _nrecombinationrate;
    
    double _vel_x;
    double _vel_y;
    double _vel_z;
    
    double _electron_vel_x;
    double _electron_vel_y;
    double _electron_vel_z;
    
    double _hole_vel_x;
    double _hole_vel_y;
    double _hole_vel_z;
    
    double _vx_old;
    double _reco_old;
    
    int _direct_iv_counter;
    int _direct_reco_counter;
    
    bool _direct_iv_convergence;
    bool _direct_reco_convergence;
};

void Numoutput::Initialize() {
    _nelectrons = 0; _nholes = 0; _ncarriers = 0;
    
    _direct_iv_counter = 0;
    _direct_reco_counter = 0;
    
    _direct_iv_convergence = false;
    _direct_reco_convergence = false;
    
    Initialize_equilibrate();
}

void Numoutput::Initialize_equilibrate() {
    
    _ninjections = 0; _ncollections = 0; 
    _nleftinjections = 0; _nrightinjections = 0;
    _nleftcollections = 0; _nrightcollections = 0;
    
    _nplaintransfer = 0; _nrecombinations = 0; _ninject_to_recombination = 0;
    _ninjectionrate = 0.0; _ncollectionrate = 0.0; _nrecombinationrate = 0.0;
    
    _vel_x = 0.0; _vel_y = 0.0; _vel_z = 0.0;
    _electron_vel_x = 0.0; _electron_vel_y = 0.0; _electron_vel_z = 0.0;
    _hole_vel_x = 0.0; _hole_vel_y = 0.0; _hole_vel_z = 0.0;       
}

void Numoutput::Init_convergence_check(double simtime) {
    _vx_old = _vel_x/simtime;
    _reco_old = _nrecombinations/simtime;
}

void Numoutput::Convergence_check(double simtime, Eventinfo* eventinfo) {
    
    if (fabs(_vel_x/simtime-_vx_old)/_vx_old < 0.05)               { _direct_iv_counter++;  } else { _direct_iv_counter = 0;  }
    if (fabs(_nrecombinations/simtime-_reco_old)/_reco_old < 0.05) { _direct_reco_counter++;} else { _direct_reco_counter = 0;}

    if(_direct_iv_counter >= eventinfo->number_direct_conv_iv) {_direct_iv_convergence = true;}
    if(_direct_reco_counter >= eventinfo->number_direct_conv_reco) {_direct_reco_convergence = true;}

    _vx_old   = _vel_x/simtime;
    _reco_old = _nrecombinations/simtime;
}

void Numoutput::Update(Event* event, double simtime, double timestep) {
    
    if(event->init_type() == (int) Injection)   { // Injection events
        _ninjections++;
        if(event->link()->node1()->type() == LeftElectrodeNode) { _nleftinjections++; } else { _nrightinjections++; }
    
        if(event->final_type() == (int) TransferTo) {
            _ncarriers++;
            if(event->carrier_type() == (int) Electron) { _nelectrons++; } else { _nholes++;}
        }
        else if(event->final_type() == (int) Recombination) {
            _nrecombinations++;
            _ninject_to_recombination++;
            _ncarriers--; 
            if(event->carrier_type() == (int) Electron) {
                // hole gets recombined away
                _nholes--;
            }
            else {
                // electron gets recombined away
                _nelectrons--;
            }
        }
        else if(event->final_type() == (int) Collection) {
            _ncollections++;
            if(event->link()->node2()->type() == LeftElectrodeNode) { _nleftcollections++;} else { _nrightcollections++;}
        }
    }
    else if(event->init_type() == (int) TransferFrom) { // Normal transfer
        if(event->final_type() == (int) Collection) {
            _ncollections++;
            if(event->link()->node2()->type() == LeftElectrodeNode) { _nleftcollections++;} else { _nrightcollections++;}

            _ncarriers--;
            if(event->carrier_type() == (int) Electron) { _nelectrons--; } else { _nholes--;}    
        }
        if(event->final_type() == (int) Recombination) { //injection to recombination (sink of charges here)
            _nrecombinations++;
            _nelectrons--;    _ncarriers--;
            _nholes--;        _ncarriers--;  
        }
        else if(event->final_type() == (int) TransferTo) {
            _nplaintransfer++;
        }
    }
        
    _ninjectionrate = _ninjections/simtime;
    _ncollectionrate = _ncollections/simtime;
    _nrecombinationrate = _nrecombinations/simtime;
    
    votca::tools::vec travelvec = event->link()->r12();
    double direction;
    if(event->carrier_type() == (int) Electron) { direction = -1.0; } else { direction = 1.0;}
    
    _vel_x += direction*travelvec.x();
    _vel_y += direction*travelvec.y();
    _vel_z += direction*travelvec.z();
    
    if(event->carrier_type() == (int) Electron) {
        _electron_vel_x += direction*travelvec.x();
        _electron_vel_y += direction*travelvec.y();
        _electron_vel_z += direction*travelvec.z();        
    }
    
    if(event->carrier_type() == (int) Hole) {
        _hole_vel_x += direction*travelvec.x();
        _hole_vel_y += direction*travelvec.y();
        _hole_vel_z += direction*travelvec.z();        
    }
}

void Numoutput::Write(double simtime) {
    std::cout << " el " << _nelectrons << " ho " << _nholes  << " ca " << _ncarriers << 
            " tr " << _nplaintransfer << " rec " << _nrecombinations << " irec " << _ninject_to_recombination <<
            " in " << _ninjections << " co " << _ncollections <<
//            " lin " << _nleftinjections << " rin " << _nrightinjections << " lco " << _nleftcollections << " rco " << _nrightcollections <<
            " ira " << _ninjectionrate << " cra " << _ncollectionrate << " rra " << _nrecombinationrate << 
            " vx " << _vel_x/simtime << " vy " << _vel_y/simtime << " vz " << _vel_z/simtime <<
//            " evx " << _electron_vel_x << " evy " << _electron_vel_y << " evz " << _electron_vel_z <<
//            " hvx " << _hole_vel_x << " hvy " << _hole_vel_y << " hvz " << _hole_vel_z << 
            endl;
}

}} 

#endif // _VOTCA_KMC_NUMOUTPUT_H
