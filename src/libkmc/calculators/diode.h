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

#ifndef __VOTCA_KMC_DIODE_H
#define	__VOTCA_KMC_DIODE_H

#include <iostream>
#include <votca/kmc/state.h>
#include <votca/kmc/graphsql.h>
#include <votca/kmc/graphcubic.h>
#include <votca/kmc/globaleventinfo.h>
#include <votca/kmc/events.h>
#include <votca/kmc/vssmgroup.h>

using namespace std;

namespace votca { namespace kmc {
    
//typedef votca::tools::vec myvec;

   
class Diode : public KMCCalculator 
{
public:
    
    Graph* graph;
    State* state;
    Events* events;
    Vssmgroup* vssmgroup;
    Globaleventinfo* globevent;
    
    Diode() {};
   ~Diode() {};

    void Initialize(const char *filename, Property *options, const char *outputfile);
    bool EvaluateFrame();

    double sim_time;
    
    // input parameters (put in globaleventinfo?)
    int seed; long nr_equilsteps; long nr_timesteps; long steps_update_longrange;
    int nx; int ny; int nz; double lattice_constant; double hopdist; double disorder_strength; 
    double disorder_ratio; CorrelationType correlation_type; double left_electrode_distance; double right_electrode_distance;
    myvec graph_front; myvec graph_back;
    

protected:
   void RunKMC(void); 
            
private:
    static const double kB   = 8.617332478E-5; // eV/K
    static const double hbar = 6.5821192815E-16; // eV*s
    static const double eps0 = 8.85418781762E-12/1.602176565E-19; // e**2/eV/m = 8.85418781762E-12 As/Vm
    static const double epsr = 3.0; // relative material permittivity
    static const double Pi   = 3.14159265358979323846;
   
};


void Diode::Initialize(const char *filename, Property *options, const char *outputfile) {
    
    cout << "Initializing" << endl;
    
    nx = options->get("options.diode.nx").as<int>();
    ny = options->get("options.diode.ny").as<int>();
    nz = options->get("options.diode.nz").as<int>();
    lattice_constant = options->get("options.diode.lattice_constant").as<double>();
    graph_front = myvec(options->get("options.diode.graph_front_x").as<double>(),
                        options->get("options.diode.graph_front_y").as<double>(),
                        options->get("options.diode.graph_front_z").as<double>());
    graph_back = myvec(options->get("options.diode.graph_back_x").as<double>(),
                       options->get("options.diode.graph_back_y").as<double>(),
                       options->get("options.diode.graph_back_z").as<double>());    
    
    /*
    graph = new GraphCubic();
    graph->Create_cubic_graph_nodes(nx, ny, nz, lattice_constant, graph_front, graph_back);
    graph->Print(std::cout);
    delete graph;    
    */

    graph = new GraphSQL();
    graph->Initialize();
    graph->Print(std::cout);
    delete graph;    
    
    
    exit(0);
    
    state = new State();
    events = new Events();
    vssmgroup = new Vssmgroup();
}

bool Diode::EvaluateFrame() {
    
    // register all QM packages (Gaussian, turbomole, etc))
    // EventFactory::RegisterAll(); 
        
    RunKMC();
}

void Diode::RunKMC() {
 
    exit(0);
    /*    
    //Setup random number generator
    seed = 1;
    srand(seed); // srand expects any integer in order to initialise the random number generator
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());    
    
    //Initialize all structures
    graph->hopdist = hopdist;
    graph->Generate_cubic_graph(nx, ny, nz, lattice_constant, disorder_strength,RandomVariable, disorder_ratio, 
                                correlation_type, left_electrode_distance, right_electro_distance,globevent);   
    state->Init();    
    events->Initialize_eventvector(graph, state, globevent);
    events->Initialize_longrange (graph, globevent);
    events->Recompute_all_injection_events(graph, globevent);
    events->Recompute_all_non_injection_events(graph, state, globevent); 
    
    
    sim_time = 0.0;
    for (long it = 0; it < 2*nr_equilsteps + nr_timesteps; it++) {
        
        // Update longrange cache (expensive, so not done at every timestep)
        if(ldiv(it, steps_update_longrange).rem == 0 && it>0){
            events->longrange->Update_cache(graph->sim_box_size, globevent);
            events->Recompute_all_injection_events(graph, globevent);
            events->Recompute_all_non_injection_events(graph,state,globevent);
        }
        
        vssmgroup->Recompute_in_device(events);
        sim_time += vssmgroup->Timestep(RandomVariable);
        vssmgroup->Perform_one_step_in_device(events,graph,state,globevent,RandomVariable);
    }
     */
}

}}


#endif	/* __VOTCA_KMC_DIODE_H */