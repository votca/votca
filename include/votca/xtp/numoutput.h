
 /* Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/xtp/event.h>
#include <votca/xtp/graph.h>
#include <votca/xtp/visualisation.h>

namespace votca { namespace xtp {

class Numoutput
{

public:
    
    Numoutput() {
        visualisation = new Visualisation();
    }
    
    ~Numoutput() {
        delete visualisation;
    }
    
    void Initialize(Eventinfo* eventinfo);
    void Initialize_equilibrate(Eventinfo* eventinfo);
    void Init_convergence_check(double simtime);
    void Convergence_check(double simtime, Eventinfo* eventinfo);
    void Update(Event* event, double simtime, double timestep);
    void Write(int it, double simtime, double timestep, Eventinfo* eventinfo);
    
    void Repeat_count_init();
    void Repeat_count_update(Event* chosenevent);
    
    void Init_visualisation(GraphKMC* graph, Eventinfo* eventinfo) {visualisation->Init_visualisation(graph, eventinfo);}
    void Update_visualisation(Event* event) {visualisation->Update_visualisation(event);}
    void Print_visualisation() {visualisation->Print_visualisation();}
    
    const bool &iv_conv() const {return _direct_iv_convergence;}
    const bool &reco_conv() const {return _direct_reco_convergence;}
    
    const int &iv_count() const {return _direct_iv_counter;}
    const int &reco_count() const {return _direct_reco_counter;}
    
    const int &holes() const {return _nholes;}
    const int &nr_repeats() const {return _repeat_counter;}
    
private:
    
    Visualisation* visualisation;
    
    void Write_header_one(Eventinfo* eventinfo);
    void Write_header_two();
    void Write_header_three();
    void Write_header_four(Eventinfo* eventinfo);

    vector<double> layercurrent;
    
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
    
    double _intvel_x;
    double _intvel_y;
    double _intvel_z;
    
    double _electron_vel_x;
    double _electron_vel_y;
    double _electron_vel_z;
    
    double _hole_vel_x;
    double _hole_vel_y;
    double _hole_vel_z;
    
    double _vz_old;
    double _reco_old;
    
    int _direct_iv_counter;
    int _direct_reco_counter;
    
    bool _direct_iv_convergence;
    bool _direct_reco_convergence;
    
    int _repeat_counter; 
    int _old_from_node_id;
    int _old_to_node_id;    
    
};

void Numoutput::Initialize(Eventinfo* eventinfo) {
    _nelectrons = 0; _nholes = 0; _ncarriers = 0;
    
    _direct_iv_counter = 0;
    _direct_reco_counter = 0;
    
    _direct_iv_convergence = false;
    _direct_reco_convergence = false;


    
    Initialize_equilibrate(eventinfo);
}

void Numoutput::Initialize_equilibrate(Eventinfo* eventinfo) {
    
    _ninjections = 0; _ncollections = 0; 
    _nleftinjections = 0; _nrightinjections = 0;
    _nleftcollections = 0; _nrightcollections = 0;
    
    _nplaintransfer = 0; _nrecombinations = 0; _ninject_to_recombination = 0;
    _ninjectionrate = 0.0; _ncollectionrate = 0.0; _nrecombinationrate = 0.0;
    
    _vel_x = 0.0; _vel_y = 0.0; _vel_z = 0.0;
    _intvel_x = 0.0; _intvel_y = 0.0; _intvel_z = 0.0;
    _electron_vel_x = 0.0; _electron_vel_y = 0.0; _electron_vel_z = 0.0;
    _hole_vel_x = 0.0; _hole_vel_y = 0.0; _hole_vel_z = 0.0;

    layercurrent.clear();
    for (int i =0; i < eventinfo->number_of_layers; i++) {
        layercurrent.push_back(0.0);
    }    

}

void Numoutput::Init_convergence_check(double simtime) {
    _vz_old = _vel_z/simtime;
    _reco_old = _nrecombinations/simtime;
}

void Numoutput::Convergence_check(double simtime, Eventinfo* eventinfo) {
    
    if (fabs(_vel_z/simtime-_vz_old)/_vz_old < 0.08 && _vel_z > 0.0)               { _direct_iv_counter++;  } else { _direct_iv_counter = 0;  }
    if (fabs(_nrecombinations/simtime-_reco_old)/_reco_old < 0.08) { _direct_reco_counter++;} else { _direct_reco_counter = 0;}

    if(_direct_iv_counter >= eventinfo->number_direct_conv_iv) {_direct_iv_convergence = true;}
    if(_direct_reco_counter >= eventinfo->number_direct_conv_reco) {_direct_reco_convergence = true;}

    _vz_old   = _vel_z/simtime;
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
    
    Node* node1 = event->link()->node1();
    Node* node2 = event->link()->node2();
    //votca::tools::vec nodepos1 = node1->position();
    //votca::tools::vec nodepos2 = node2->position();
    int node1_layer = dynamic_cast<NodeDevice*>(node1)->layer();
    int node2_layer = dynamic_cast<NodeDevice*>(node2)->layer();
    
    votca::tools::vec travelvec = event->link()->r12();
    double direction;
    double dirx; double diry; double dirz;

    if(event->carrier_type() == (int) Electron) { direction = -1.0; } else { direction = 1.0;}
    if(travelvec.x() > 0) {dirx = 1.0;} else {if(travelvec.x() == 0) {dirx = 0.0;} else {dirx = -1.0;}}
    if(travelvec.y() > 0) {diry = 1.0;} else {if(travelvec.y() == 0) {diry = 0.0;} else {diry = -1.0;}}
    if(travelvec.z() > 0) {dirz = 1.0;} else {if(travelvec.z() == 0) {dirz = 0.0;} else {dirz = -1.0;}}
    
    if(node1->type() == (int) NormalNode) {
        if(node1_layer != node2_layer) layercurrent[node1_layer] += 0.5*dirz;
        _intvel_x += 0.5*dirx;        
        _intvel_y += 0.5*diry;
        _intvel_z += 0.5*dirz;              
    }
    
    if(node2->type() == (int) NormalNode) {
        if(node1_layer != node2_layer) layercurrent[node2_layer] += 0.5*dirz;
        _intvel_x += 0.5*dirx;        
        _intvel_y += 0.5*diry;
        _intvel_z += 0.5*dirz;              
    }    
    
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
 
void Numoutput::Write(int it, double simtime, double timestep, Eventinfo* eventinfo) {
    
    this->Write_header_one(eventinfo);
    std::cout << setw(15) << it;
    if(eventinfo->repeat_counting) std::cout << setw(15) << this->nr_repeats();
    std::cout << setw(15) << this->iv_conv();
    std::cout << setw(15) << this->iv_count();
    std::cout << setw(15) << this->reco_conv();
    std::cout << setw(15) << this->reco_count();
    std::cout << setw(20) << simtime;
    std::cout << setw(20) << timestep;
    std::cout << "\n";
    std::cout << "\n";

    this->Write_header_two();
    std::cout << setw(15) << _nelectrons;
    std::cout << setw(15) << _nholes;    
    std::cout << setw(15) << _ncarriers;
    std::cout << setw(15) << _nplaintransfer;
    std::cout << setw(15) << _nrecombinations;
    std::cout << setw(15) << _nrecombinationrate;
    std::cout << "\n";
    std::cout << "\n";
    
    if(eventinfo->device){
        this->Write_header_three();
        std::cout << setw(20) << _ninject_to_recombination;
        std::cout << setw(15) << _ninjections;
        std::cout << setw(15) << _ncollections;
        std::cout << setw(20) << _nleftinjections;
        std::cout << setw(20) << _nleftcollections;
        std::cout << setw(20) << _nrightinjections;
        std::cout << setw(20) << _nrightcollections;
        std::cout << "\n";
        std::cout << "\n";
    }
    
    this->Write_header_four(eventinfo);
    if(eventinfo->device) {
        std::cout << setw(15) << _ninjectionrate;
        std::cout << setw(15) << _ncollectionrate;
    }
    double layercur = 0.0;
    for (int i = 0; i<eventinfo->number_of_layers; i++) {
        layercur += layercurrent[i];
    }
    layercur /= eventinfo->number_of_layers;
     std::cout << setw(15) << _vel_x/simtime;
    std::cout << setw(15) << _vel_y/simtime;
    std::cout << setw(15) << _vel_z/simtime;
    std::cout << setw(15) << _intvel_x/simtime;
    std::cout << setw(15) << _intvel_y/simtime;
    std::cout << setw(15) << _intvel_z/simtime;
    std::cout << setw(15) << layercur/simtime;
    std::cout << "\n";      
    std::cout << "\n";
}

void Numoutput::Write_header_one(Eventinfo* eventinfo) {
    std::cout << setw(15) << "event nr";
    if(eventinfo->repeat_counting) std::cout << setw(15) << "nr repeats";
    std::cout << setw(15) << "iv conv";
    std::cout << setw(15) << "iv count";
    std::cout << setw(15) << "reco conv";
    std::cout << setw(15) << "reco count";
    std::cout << setw(20) << "sim_time";
    std::cout << setw(20) << "timestep";
    std::cout << "\n";
}

void Numoutput::Write_header_two() {
    std::cout << setw(15) << "nr electrons";
    std::cout << setw(15) << "nr holes";    
    std::cout << setw(15) << "nr carriers";
    std::cout << setw(15) << "nr transfers";
    std::cout << setw(15) << "nr recombins";
    std::cout << setw(15) << "rec rate";
    std::cout << "\n";
}

void Numoutput::Write_header_three() {
    std::cout << setw(20) << "nr_rec_to_inject";
    std::cout << setw(15) << "nr injects";
    std::cout << setw(15) << "nr collects";
    std::cout << setw(20) << "nr left injects";
    std::cout << setw(20) << "nr left collects";
    std::cout << setw(20) << "nr right injects";
    std::cout << setw(20) << "nr right collects";
    std::cout << "\n";
}

void Numoutput::Write_header_four(Eventinfo* eventinfo) {
    if(eventinfo->device) {
        std::cout << setw(15) << "inject rate";
        std::cout << setw(15) << "collect rate";
    }
    std::cout << setw(15) << "av vel x";
    std::cout << setw(15) << "av vel y";
    std::cout << setw(15) << "av vel z";
    std::cout << "\n";        
}

void Numoutput::Repeat_count_init(){
    _repeat_counter = 0; 
    _old_from_node_id = -10;
    _old_to_node_id = 10;          
}

void Numoutput::Repeat_count_update(Event* chosenevent){
    int goto_node_id = chosenevent->link()->node2()->id();
    int from_node_id = chosenevent->link()->node1()->id();
    if(goto_node_id == _old_from_node_id && from_node_id == _old_to_node_id) _repeat_counter++;
    _old_from_node_id = from_node_id;
    _old_to_node_id = goto_node_id;
}

}} 

#endif // _VOTCA_KMC_NUMOUTPUT_H
