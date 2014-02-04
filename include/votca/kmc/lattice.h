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

#ifndef __VOTCA_KMC_LATTICE_H_
#define __VOTCA_KMC_LATTICE_H_


#include <votca/tools/database.h>
#include <votca/tools/vec.h>
#include <votca/tools/random.h>
#include "carrier.h"
#include "node.h"
#include "state.h"
#include "graph.h"
#include "rates.h"
#include "longrange.h"

namespace votca { namespace kmc {
  
using namespace std;

typedef votca::tools::vec myvec;

enum CarrierType{Electron, Hole, Exciton};
enum Action{Add, Remove};
enum GraphType{Square, Load};
enum CorrelationType {Correlated, Anticorrelated, Uncorrelated};

class Lattice {
  public:
    void Initialize(votca::tools::Random2 *RandomVariable);
    void Lattice_Grow(int growsize,CarrierType carriers);
    void Recompute_Injection_Rates();
    void Recompute_Jump_Rates();
    double Compute_Jump_Rate(Node* finalnode, Node* initnode, CarrierType cartype1);
    double Compute_Coulomb_Potential(double startx, myvec dif);
    void Add_Remove_Carrier(Node* carvertex, CarrierType cartype, Action add_remove);

    // aiding methods
    myvec periodic(myvec v); // returns periodic vector
    myvec perdifference(myvec v, myvec w); // returns difference vector (final node - initial node);
   
    int _nholes;
    int _nelectrons;
    int _ncarriers;
    
    long _init_maxholes;
    long _init_maxelectrons;
    long _init_maxexcitons;
    
    double _alpha;
    double _beta;
    double _efield;
    double _binding_energy;
    double _coulomb_strength;
    
    double _electron_prefactor;
    double _hole_prefactor;
    double _injection_prefactor;
    double _collection_prefactor;
    double _recombination_prefactor;
    double _transfer_prefactor;
    
    double _injection_bar_left;
    double _injection_bar_right;
    
    string _formalism;
     
  private:

    Graph* _graph;
    Longrange* _longrange;
    Rates _electron_rates;
    Rates _hole_rates;
    Rates _exciton_rates;
    State* _electron_state;
    State* _hole_state;
    State* _exciton_state;
    
    int _NCHARGE_X; //number of subboxes for coulomb interaction calculations
    int _NCHARGE_Y;
    int _NCHARGE_Z;
    double _coul_box_dimX; //dimensions of these subboxes
    double _coul_box_dimY;
    double _coul_box_dimZ;     

    int _NSITE_X; //number of subboxes for neighbour calculations
    int _NSITE_Y;
    int _NSITE_Z;
    double _site_box_dimX; //dimensions of these subboxes
    double _site_box_dimY;
    double _site_box_dimZ;
     
    double _NX; //simulation box dimensions
    double _NY;
    double _NZ;
     
    double _left_electrode_distance; //distance between left electrode and first leftmost node
    double _right_electrode_distance; //distance between right electrode and first rightmost node
     
    bool _device; //device setting (false means bulk setting)
    double _hopping_distance; //maximum hopping distance
    bool _dual_injection; //injection of both holes and electrons from the same electrode
    double _self_image_prefactor;
     
    CorrelationType _correlation;
    double _disorder_strength;
    double _disorder_factor; //sigma_el/sigma_ho
 
    GraphType _graphtype;
    string _graph_filename;
    bool _load_pairs; //load pairs of nodes (note that the hopping distance as set to determine the pairs can differ from hopping distance as set in program)
    bool _load_static_energies; //load static energies (calculate from random gaussian otherwise)
    bool _load_transfers;
    
    vector<double> _positional_average; //positional average of node positions in every layer    
     
    long _maxholes; // maximum number for state and rates classes
    long _maxelectrons;
    long _maxexcitons; 
    long _maxcarriers;

    double _coulomb_cutoff;
    int _number_of_shortrange_images;
    int _number_of_longrange_images;
    
    std::list<int> _chargelist[_NCHARGE_X][_NCHARGE_Y][_NCHARGE_Z]; // List of charges in each subbox
    std::list<int> _sitelist[_NSITE_X][_NSITE_Y][_NSITE_Z]; // List of charges in each subbox
    
};

void Lattice:: Initialize(votca::tools::Random2 *RandomVariable ) {

  // define subboxing parameters
  _NSITE_X = floor(_NX/_site_box_dimX); // number of subboxes in mesh, size defined by hopping distance
  _NSITE_Y = floor(_NY/_site_box_dimY);
  _NSITE_Z = floor(_NZ/_site_box_dimZ);
    
  _NCHARGE_X = floor(_NX/_coul_box_dimX); // number of subboxes in mesh, size defined by coulomb cut off radius
  _NCHARGE_Y = floor(_NY/_coul_box_dimY);
  _NCHARGE_Z = floor(_NZ/_coul_box_dimZ);

  // create graph
  _graph->_NX = _NX;
  _graph->_NY = _NY;
  _graph->_NZ = _NZ;
  
  _graph->_NSITE_X = _NSITE_X;
  _graph->_NSITE_Y = _NSITE_Y;
  _graph->_NSITE_Z = _NSITE_Z;
  _graph->_site_box_dimX = _site_box_dimX;
  _graph->_site_box_dimY = _site_box_dimY;
  _graph->_site_box_dimZ = _site_box_dimZ;
  
  _graph->_device = _device;
  _graph->_hopping_distance = _hopping_distance;
  _graph->_lattice_constant = _lattice_constant;
  _graph->_left_electrode_distance = _left_electrode_distance;
  _graph->_right_electrode_distance = _right_electrode_distance;
  
  _graph->_correlation = _correlation;
  _graph->_disorder_strength = _disorder_strength;
  _graph->_disorder_factor = _disorder_factor;
  
  
  //define what kind of graph should be created
  
  if(_graphtype == Load) {
    _graph->_graph_filename = _graph_filename;
    _graph->LoadGraphNodes();
    
    if(_load_static_energies) {
      _graph->LoadStaticEnergies();     
    }
    else {
      _graph->CreateRandomStaticEnergies(RandomVariable);
    }

    if(_load_pairs) {
      _graph->LoadGraphPairs();
    }
    else {
      _graph->CalculateGraphPairs();
    }
    if(_load_transfers) {
      _graph->LoadGraphTransfers();        
    }
  }
  else if(_graphtype == Square) {
    _graph->CreateCubicGraphNodes();
    _graph->CreateRandomStaticEnergies(RandomVariable);
    _graph->CalculateGraphPairs();    
  }
  
  _graph->SetupInjectors();
  _graph->InitializeLayers();
  
  _sitelist = _graph->_sitelist;
  
  // initialize longrange class

  _longrange->_positional_average = _graph->_positional_average;
  _longrange->_number_of_layers = _NSITE_X;
  _longrange->_device_length = _NX;
  _longrange->_box_dimension_Y = _NY;
  _longrange->_box_dimension_Z = _NZ;
  _longrange->_number_of_longrange_images = _number_of_longrange_images;
  _longrange->Initialize();

  // initialize the binary sum trees
  
  _maxelectrons = _init_maxelectrons;
  _maxholes = _init_maxholes;
  _maxexcitons = _init_maxexcitons;
  _maxcarriers = _init_maxelectrons + _init_maxholes + _init_maxexcitons;
  
  if(_dual_injection) {
    _electron_rates.initialize(_graph->_number_of_injector_nodes,_maxelectrons,_graph->_max_number_of_el_pairs);
    _hole_rates.initialize(_graph->_number_of_injector_nodes,_maxholes,_graph->_max_number_of_ho_pairs);
  }
  else {
    _electron_rates.initialize(_graph->_number_of_right_injector_nodes,_maxelectrons,_graph->_max_number_of_el_pairs);
    _hole_rates.initialize(_graph->_number_of_left_injector_nodes,_maxholes,_graph->_max_number_of_ho_pairs);      
  }

  // initial growth of structures to accomodate carriers
  _electron_state->grow(_maxelectrons);
  _hole_state->grow(_maxholes);
  _exciton_state->grow(_maxexcitons);
  
  _ncarriers = 0;
  _nholes = 0;
  _nelectrons = 0;
  _nexcitons = 0;
  
  for(int inode=0; inode<_graph->_number_of_nodes;inode++) { // pre-empty the device
    _graph->_nodes[inode]->_hole_number=-1;
    _graph->_nodes[inode]->_electron_number=-1;
  }
  
}

void Lattice:: Lattice_Grow(int growsize, Carriertype carriers){
    if(carriers==Electron) {
        _maxelectrons += growsize;
        _electron_state->grow(growsize);
        _electron_rates.resize(_maxelectrons);
    }
    else if(carriers==Hole) {
        _maxholes += growsize;
        _hole_state->grow(growsize);
        _hole_rates.resize(_maxholes);
    }
    else if(carriers==Exciton) {
        _maxexcitons += growsize;
        _exciton_state->grow(growsize);
        _exciton_rates.resize(_maxexcitons);
    }
    _maxcarriers += growth;
}

void Lattice ::Recompute_Injection_Rates() {
  for(int count=0;count<_graph->_number_of_left_injector_nodes;count++) {
    // Injection rate for layer x=0
    Node* finalnode = _graph->_nodes[_left_injector_nodes[count]];
    Node* initnode = _graph->_left_electrode;
    double hole_injection_rate = compute_jump_rate(finalnode,initnode,Hole);
    
    if(!_dual_injection) {
      int left_injector_index = finalnode->_left_injectorindex;
      _hole_rates.set_injection_rate(left_injector_index,hole_injection_rate);
    }
    else{ //if both electrons and holes can be injected from the same electrode
      int injector_index = finalnode->_injectorindex;
      _hole_rates.set_injection_rate(injector_index,hole_injection_rate);
      double electron_injection_rate = compute_jump_rate(finalnode,initnode,Electron);
      _electron_rates.set_injection_rate(injector_index,electron_injection_rate);
    }
  }

  for(int count=0;count<_graph->_number_of_right_injector_nodes;count++) {
    // Injection rate for layer x=0
    Node* finalnode = _graph->_nodes[_right_injector_nodes[count]];
    Node* initnode = _graph->_right_electrode;
    double electron_injection_rate = compute_jump_rate(finalnode,initnode,Electron);
    
    if(!_dual_injection) {
      int right_injector_index = finalnode->_right_injectorindex;
      _electron_rates.set_injection_rate(right_injector_index,electron_injection_rate);
    }
    else{ //if both electrons and holes can be injected from the same electrode
      int injector_index = finalnode->_injectorindex;
      _electron_rates.set_injection_rate(injector_index,electron_injection_rate);
      double hole_injection_rate = compute_jump_rate(finalnode,initnode,Hole);
      _hole_rates.set_injection_rate(injector_index,hole_injection_rate);
    }    
  }  
}

void Lattice::Recompute_Jump_Rates() {
  for(long ic=0; ic<_nholes; ic++) {
    int itemnr = _hole_state->get_itemnr_by_index(ic);
    pcarrier1 = _hole_state->Get_item(itemnr);
    int node_id = pcarrier1->_node_index;
    Node* initnode = _graph->_nodes[node_id];
    
    int number_of_pairs = initnode->_number_of_ho_pairs;
    
    for(int jump=0; jump<number_of_pairs; jump++) {
      Node* finalnode = initnode->_ho_pair_nodes[jump];
      double rate = compute_jump_rate(finalnode,initnode,Hole);
      _hole_rates.set_jump_rate(itemnr,jump,rate);          
    }
  }
  
  for(long ic=0; ic<_nelectrons; ic++) {
    int itemnr = _electron_state->get_itemnr_by_index(ic);
    pcarrier1 = _electron_state->Get_item(itemnr);
    int node_id = pcarrier1->_node_index;
    Node* initnode = _graph->_nodes[node_id];
    
    int number_of_pairs = initnode->_number_of_el_pairs;
    
    for(int jump=0; jump<number_of_pairs; jump++) {
      Node* finalnode = initnode->_el_pair_nodes[jump];
      double rate = compute_jump_rate(finalnode,initnode,Electron);
      _electron_rates.set_jump_rate(itemnr,jump,rate);                  
    }
  }
}

void Lattice:: Add_Remove_Carrier(Node* carnode, CarrierType cartype, enum action) {

  double RC = _coulomb_cutoff;
  double RCSQR = RC*RC;
  double hopdist = _hopping_distance;
  int arsign;  // arsign = 1 for add, arsign = -1 for remove
  
  int carrier1;
  Carrier* pcarrier;
  int charge1;
  Myvec carpos = carnode->_position;
  int carnumberofpairs;

  if(action == Add) {
    // Add new carrier to the lattice
    arsign = 1;
    if(cartype==Hole) {
      carrier1 = _hole_state.Buy();
      pcarrier = _hole_state.Get_item(carrier1);
      carnode->_hole_number = carrier1;
    }
    else if(cartype==Electron) {
      carrier1 = _electron_state.Buy();
      pcarrier = _electron_state.Get_item(carrier1);        
      carnode->_electron_number = carrier1;
    }
    pcarrier->_carriertype = cartype;
    pcarrier->_node_index = carnode->_node_id;
    carnode->_carrier = pcarrier1;
    
    if (cartype == Hole) {
      if(_nholes == _maxholes){
        int growsize = 10;
        lattice_grow(growsize);
      } 
      _nholes++;
      charge1 = 1;
      pcarrier->_shortrange_coulomb.resize(1+carnode->_number_of_ho_pairs);
      carnumberofpairs = _number_of_ho_pairs;
    }
    else if (cartype == Electron) {
      if(_nelectrons == _maxelectrons){
        int growsize = 10;
        lattice_grow(growsize);
      }
      _nelectrons++;
      charge1 = -1;
      pcarrier->_shortrange_coulomb.resize(1+carnode->_number_of_el_pairs);
      carnumberofpairs = _number_of_el_pairs;
    } 
    _ncarriers++;

    // add charge to charge list
    int iposx = floor(carpos.x()/_coul_box_dimX); 
    int iposy = floor(carpos.y()/_coul_box_dimY); 
    int iposz = floor(carpos.z()/_coul_box_dimZ);
    chargelist[iposx][iposy][iposz].push_back(carrier1);
  }
  else if(action == Remove) {
    // Remove existing carrier from lattice
    arsign = -1;
    Carriertype rem_cartype = carnode->_carrier->_carriertype;
    if (rem_cartype == Hole) {
      _nholes--;
      charge1 = 1;
      carnode->_hole_number = -1; //empty (necessary to set to -1)
    }
    else if (rem_cartype == Electron) {
      _nelectrons--;
      charge1 = -1; 
      carnode->_electron_number = -1; //empty (necessary to set to -1)
    }
    carnode->_carrier = NULL;
    _ncarriers--;
  }

  int layerindex = carnode->_layer_index;
  _longrange->_layercharge[layerindex] += arsign*charge1; 
  

  // Define cubic boundaries in non-periodic coordinates
  double ix1 = carpos.x-RC-hopdist; double ix2 = carpos.x+RC+hopdist;
  double iy1 = carpos.y-RC-hopdist; double iy2 = carpos.y+RC+hopdist;
  double iz1 = carpos.z-RC-hopdist; double iz2 = carpos.z+RC+hopdist;

  // Break periodicity in x-direction
  if(_device) {
    if (ix1<0.0) ix1 = 0.0;
    if (ix2>=_NX) ix2 = _NX;
  }
  
  // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
  int sx1 = floor(ix1,_coul_box_dimX);
  int sx2 = floor(ix2,_coul_box_dimX);
  int sy1 = floor(iy1,_coul_box_dimY);
  int sy2 = floor(iy2,_coul_box_dimY);
  int sz1 = floor(iz1,_coul_box_dimZ);
  int sz2 = floor(iz2,_coul_box_dimZ);
  
  // Now visit all relevant sublattices
  for (int isz=sz1; isz<=sz2; isz++) {
    int r_isz = isz;
    while (r_isz < 0) r_isz += _NCHARGE_Z;
    while (r_isz >= _NCHARGE_Z) r_isz -= _NCHARGE_Z;
    for (int isy=sy1; isy<=sy2; isy++) {
      int r_isy = isy;
      while (r_isy < 0) r_isy += _NCHARGE_Y;
      while (r_isy >= _NCHARGE_Y) r_isy -= _NCHARGE_Y;
      for (int isx=sx1; isx<=sx2; isx++) {
        int r_isx = isx;
        while (r_isx < 0) r_isx += _NCHARGE_X;
        while (r_isx >= _NCHARGE_X) r_isx -= _NCHARGE_X;
        
        // Ask a list of all charges in this sublattice
        std::list<int>::iterator li1,li2,li3;
        std::list<int> *pList = &chargelist[r_isx][r_isy][r_isz];
        li1 = pList->begin();
        li2 = pList->end();
        for (li3=li1; li3!=li2; li3++) {
          int carrier2 = *li3;
          Carrier* probecarrier = carriers.get_item(carrier2);
          Node* probenode = _graph->_nodes[probercarrier->_node_index];
          myvec probepos = probenode->_position;
          Carriertype probecartype = probecarrier->_carriertype;
          int probecharge;
          if(probecartype==Electron) {
            probecharge = -1;
          }
          else {
            probecharge = 1;
          }
          
          int interaction_sign = arsign*charge1*probecharge;
          
       
          // Compute coordinates in non-periodic lattice
          myvec periodic_convert = myvec((isx-r_isx)*_coul_box_dimX,(isy-r_isy)*_coul_box_dimY,(isz-r_isz)*_coul_box_dimZ);
          myvec np_probepos = probepos + periodic_convert;
          myvec difference = np_probepos-carpos;

          double distancesqr = difference.x()*difference.x() + difference.y()*difference.y()+difference.z()*difference.z();

          if (distancesqr < pow(0.5*_graph->_minimum_distance,2.0)) { // Charge interacting with its own images, special case
            if (action==ADD) {
              pcarrier ->_shortrange_coulomb[0] += _self_image_prefactor*compute_Coulomb_potential(carpos.x(),myvec(0.0,0.0,0.0));
 
              // Adjust Coulomb potential for neighbours of carrier1
              for (int jump=0; jump < carnumberofpairs; jump++) {
                myvec jumpdistance;
                if(cartype==Hole) {
                  jumpdistance = carnode->_ho_pair_distances[jump];
                }
                else {
                  jumpdistance = carnode->_el_pair_distances[jump];
                }
                pcarrier->_shortrange_coulomb[jump+1] += self_image_prefactor*compute_Coulomb_potential(carpos.x(),jumpdistance);
              }
            }
          }
          else { // dists>0, Charge interacting with another charge and its images
            if (distancesqr <= RCSQR) {
                if (action==ADD) pcarrier->_shortrange_coulomb[0] +=interaction_sign*compute_Coulomb_potential(np_probepos.x(),-1.0*difference);
                probecarrier->_shortrange_coulomb[0] += interaction_sign*compute_Coulomb_potential(carpos.x(),difference);
            }

            if (action==ADD) {
              // Adjust Coulomb potential for neighbours of carrier1
              for (int jump=0; jump < carnumberofpairs; jump++) {
                myvec jumpdistance;
                if(cartype == Hole){         
                  jumpdistance = carnode->_ho_pair_distances[jump];
                }
                else {
                  jumpdistance = carnode->_el_pair_distances[jump];
                }
                myvec jumppos = carnode->_position + jumpdistance;
                myvec jumpdifference = np_probepos - jumppos;
                double distsqr = jumpdifference.x()*jumpdifference.x() + jumpdifference.y()*jumpdifference.y()+jumpdifference.z()*jumpdifference.z*();

                if(distsqr <= RCSQR) {
                  pcarrier->_shortrange_coulomb[jump+1] += interaction_sign*compute_Coulomb_potential(np_probepos.x(),jumpdifference);
                }
              }
            }
           
            // Adjust Coulomb potential for neighbours of carrier2
            if (probecartype == Hole) {
                probenumberofpairs = probenode->_number_of_ho_pairs;
            }
            else if(probecartype == Electron) {
                probenumberofpairs = probenode->_number_of_el_pairs;
            }
            for (int jump=0; jump < probenumberofpairs; jump++) {
              myvec jumpdistance;
              Node* jumpnode;
              if(probecartype == Hole) {
                jumpdistance = probenode->_ho_pair_distances[jump];
                jumpnode = probenode->_ho_pair_nodes[jump];
              }
              else {
                jumpdistance = probenode->_el_pair_distances[jump];
                jumpnode = probenode->_ho_pair_nodes[jump];
              }
              myvec jumppos = np_probepos + jumpdistance;
              myvec jumpdifference = carpos - jumppos;
              double distsqr = jumpdifference.x()*jumpdifference.x() + jumpdifference.y()*jumpdifference.y()+jumpdifference.z()*jumpdifference.z*();
              if(distsqr <= RCSQR) {
                probecarrier->_shortrange_coulomb[jump+1] += interaction_sign*compute_Coulomb_potential(carpos.x(),jumpdifference);
              }

              if(probecartype == Electron) {
                  
                double rate = compute_jump_rate(jumpnode,probenode,Electron);
                _electron_rates.set_jump_rate(carrier2,jump,rate);
              }
              else {
                double rate = compute_jump_rate(jumpnode,probenode,Hole);
                _hole_rates.set_jump_rate(carrier2,jump,rate);                  
              }
            }
          }
          if (action==REMOVE) {
            // Reset Coulomb potential for carrier1 and its neighbours
            for (int jump=0; jump < carnumberofpairs+1; jump++) {
              pcarrier->_shortrange_coulomb[jump] = 0.0;  
            }
          }
        }
      }
    }
  }  

  
  // update the injection rates
  
  if (carpos.x() <= RC) { // Distance to left electrode
  
    //maximum entry inside device is hopdist (make this variable!!)

    double bound = sqrt(double(RCSQR - carpos.x()*carpos.x()));

    // Define cubic boundaries in non-periodic coordinates
    double iy1 = carpos.y-bound-hopdist; double iy2 = carpos.y+bound+hopdist;
    double iz1 = carpos.z-bound-hopdist; double iz2 = carpos.z+bound+hopdist;

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sy1 = floor(iy1,_site_box_dimY);
    int sy2 = floor(iy2,_site_box_dimY);
    int sz1 = floor(iz1,_site_box_dimZ);
    int sz2 = floor(iz2,_site_box_dimZ);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
      int r_isz = isz;
      while (r_isz < 0) r_isz += _NSITE_Z;
      while (r_isz >= _NSITE_Z) r_isz -= _NSITE_Z;
      for (int isy=sy1; isy<=sy2; isy++) {
        int r_isy = isy;
        while (r_isy < 0) r_isy += _NSITE_Y;
        while (r_isy >= _NSITE_Y) r_isy -= _NSITE_Y;
       
        // Ask a list of all charges in this sublattice
        std::list<int>::iterator li1,li2,li3;
        std::list<int> *pList = &sitelist[0][r_isy][r_isz];
        li1 = pList->begin();
        li2 = pList->end();
        for (li3=li1; li3!=li2; li3++) {
          int node_index = *li3;
          Node* probenode = _graph->_nodes[node_index];
          myvec probepos = probenode->_position;
          
          int interaction_sign = arsign*charge1;
          
          // Compute coordinates in non-periodic lattice
          
          myvec periodic_convert = myvec(0.0,(isy-r_isy)*_site_box_dimY,(isz-r_isz)*_site_box_dimZ);
          myvec np_probepos = probepos + periodic_convert;
          myvec difference = np_probepos-carpos;

          double distancesqr = difference.x()*difference.x()+difference.y*difference.y+difference.z*difference.z;

          if (distancesqr <= RCSQR) {
            probenode->_ho_injection_potential +=interaction_sign*compute_Coulomb_potential(carpos.x(),difference);
            double injector_index = probenode->_injectorindex;
            double hole_injection_rate = compute_jump_rate(probenode,_graph->_left_electrode,Hole);
            _hole_rates.set_injection_rate(injector_index,hole_injection_rate);
            if(_dual_injection) { //if both electrons and holes can be injected from the same electrode
              interaction_sign *= -1.0;
              probenode->_el_injection_potential +=interaction_sign*compute_Coulomb_potential(carpos.x(),difference);
              double electron_injection_rate = compute_jump_rate(probenode,_graph->_left_electrode,Electron);
              _electron_rates.set_injection_rate(injector_index,electron_injection_rate);
            }
          }
        }
      }
    }
  }  

  if ((NX-carpos.x()) <= RC) { // Distance to right electrode
  
    //maximum entry inside device is hopdist (make this variable!!)

    double bound = sqrt(double(RCSQR - (NX-carpos.x())*(NX-carpos.x())));

    // Define cubic boundaries in non-periodic coordinates
    double iy1 = carpos.y-bound-hopdist; double iy2 = carpos.y+bound+hopdist;
    double iz1 = carpos.z-bound-hopdist; double iz2 = carpos.z+bound+hopdist;

    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sy1 = floor(iy1,_site_box_dimY);
    int sy2 = floor(iy2,_site_box_dimY);
    int sz1 = floor(iz1,_site_box_dimZ);
    int sz2 = floor(iz2,_site_box_dimZ);
  
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
      int r_isz = isz;
      while (r_isz < 0) r_isz += _NSITE_Z;
      while (r_isz >= _NSITE_Z) r_isz -= _NSITE_Z;
      for (int isy=sy1; isy<=sy2; isy++) {
        int r_isy = isy;
        while (r_isy < 0) r_isy += _NSITE_Y;
        while (r_isy >= _NSITE_Y) r_isy -= _NSITE_Y;
       
        // Ask a list of all charges in this sublattice
        std::list<int>::iterator li1,li2,li3;
        std::list<int> *pList = &sitelist[0][r_isy][r_isz];
        li1 = pList->begin();
        li2 = pList->end();
        for (li3=li1; li3!=li2; li3++) {
          int node_index = *li3;
          Node* probenode = _graph->_nodes[node_index];
          myvec probepos = probenode->_position;
          
          int interaction_sign = -1.0*arsign*charge1; // electrons are injected at the right electrode
          
          // Compute coordinates in non-periodic lattice
          
          myvec periodic_convert = myvec(0.0,(isy-r_isy)*_site_box_dimY,(isz-r_isz)*_site_box_dimZ);
          myvec np_probepos = probepos + periodic_convert;
          myvec difference = np_probepos-carpos;

          double distancesqr = difference.x()*difference.x()+difference.y*difference.y+difference.z*difference.z;

          if (distancesqr <= RCSQR) {
            probenode->_el_injection_potential +=interaction_sign*compute_Coulomb_potential(carpos.x(),difference);
            double injector_index = probenode->_injectorindex;
            double electron_injection_rate = compute_jump_rate(probenode,_graph->_right_electrode,Electron);
            _electron_rates.set_injection_rate(injector_index,electron_injection_rate);
            if(_dual_injection) { //if both electrons and holes can be injected from the same electrode
              interaction_sign *= -1.0;
              probenode->_ho_injection_potential +=interaction_sign*compute_Coulomb_potential(carpos.x(),difference);
              double hole_injection_rate = compute_jump_rate(probenode,_graph->_right_electrode,Hole);
              _hole_rates.set_injection_rate(injector_index,hole_injection_rate);
            }
          }
        }
      }
    }
  }  
  
  // update hopping rates for carrier 1
  for (int jump=0; jump < carnode->_numberofneighbours; jump++) {
    Node* jumpnode;
    if(cartype == Electron) {
      jumpnode = probenode->_ho_pair_nodes[jump];
      double electron_rate;
      if(action== ADD) {
        electron_rate = compute_jump_rate(jumpnode,carnode,Electron);
      }
      else {
        electron_rate = 0.0;
      }
      _electron_rates.set_jump_rate(carrier1,jump,electron_rate);
    }
    else if(cartype==Hole) {
      jumpnode = probenode->_ho_pair_nodes[jump];    
      double hole_rate; 
      if(action==ADD) {
        hole_rate = compute_jump_rate(jumpnode,carnode,Hole);
      }
      else {
        hole_rate = 0.0;
      }
      _hole_rates.set_jump_rate(carrier1,jump,rate);
    }
  }
  
  if (action==REMOVE) {
      if(cartype == Hole) {
        _hole_state.Sell(carrier1);
      }
      else if(cartype == Electron) {
        _electron_state.Sell(carrier1);
      }
    // Remove charge from chargelist
    int iposx = floor(carpos.x()/_coul_box_dimX); 
    int iposy = floor(carpos.y()/_coul_box_dimY); 
    int iposz = floor(carpos.z()/_coul_box_dimZ);
    chargelist[iposx][iposy][iposz].remove(carrier1);
  }
}

double Lattice::Compute_Jump_Rate(Node* finalnode, Node* initnode, CarrierType cartype1){

  enum InjectionType{ LeftInjection, ElectronInjection, NoInjection};
  InjectionType Injection_Type;
  enum EventType{ ElectronCollection, HoleCollection, ElectronTransfer, HoleTransfer, Blocking, Recombination};
  EventType Event_type;
  int charge;
    
  double prefactor = 1.0;
  double injection_penalty = 0.0;
    
  if(initnode->_node_type  == LeftElectrode) { 
    Injection_type = LeftInjection;
    prefactor *= _injection_prefactor;
    injection_penalty = _left_injection_bar;
  }
  else if(initnode->_node_type == RightElectrode){
    Injection_type = RightInjection;
    prefactor *= _injection_prefactor;
    injection_penalty = _right_injection_bar;
  }
  else {
    Injection_type = NoInjection;
  }
    
  if(finalnode->_node_type  != Normal) { 
    if(cartype1==Electron) {
      Event_type = ElectronCollection;
      charge = -1;
    }
    else if(cartype1==Hole) {
      Event_type = HoleCollection;
      charge = 1;
    }
  }
  else {
    int elcarrier2 = finalnode->_electron_number;
    int hocarrier2 = finalnode->_hole_number;
    if ((elcarrier2 == -1)&&(hocarrier2 == -1)) { //second site is empty
      if(cartype1==Electron) {
        Event_type = ElectronTransfer;
        charge = -1;
      }
      else if(cartype1==Hole) {
        Event_type = HoleTransfer;
        charge = 1;
      }
    }
    else { //second site is non-empty
        Carrier* pcarrier2;
        if(elcarrier2 != -1) {
          pcarrier2 = _electron_state.Get_item(elcarrier2);
        }
        else if(hocarrier2 != -1) {
            pcarrier2 = _hole_state.Get_item(hocarrier2);
        }
      CarrierType cartype2 = pcarrier2->_carriertype;
      if(((cartype1 == Electron) && (cartype2 == Electron))||((cartype1 == Hole) && (cartype2 == Hole))) {
        Event_type = Blocking; // Equal carrier type encounter
      }
      else if(((cartype1 == Electron) && (cartype2 == Hole))||((cartype1 == Hole) && (cartype2 == Electron))) {
        Event_type = Recombination; // Unequal carrier type encounter
      }
    }
  }
    
  //first calculate quantum mechanical wavefunction overlap
    
  double prefactor = 1.0;
    
  if(initnode->_node_type == Normal) {//no injection
    myvec initpos = initnode->_position;
    if (finalnode->_node_type == Normal) {
      myvec finalpos = finalnode->_position;
      myvec distancevec = perdifference(initpos,finalpos);
    }
    else if (finalnode->_node_type == LeftElectrode) { //collection
      myvec distancevec = myvec(-1.0*initpos.x(),0.0,0.0);
    }
    else { //collection to right electrode
      myvec distancevec = myvec(_NX-1.0*initpos.x(),0.0,0.0); 
    }
  }
  else { //injection
    if (finalnode->_node_type == Normal) { //no collection
      myvec finalpos = finalnode->_position;
      if(initnode->_node_type == LeftElectrode) {
        myvec distancevec = myvec(finalpos.x(),0.0,0.0);
      }
      else {
        myvec distancevec = myvec(-1.0*(_NX-finalpos.x()),0.0,0.0);
      }
    }
    else { //collection
      if(initnode->_node_type == LeftElectrode) { // from left to right electrode
        myvec distancevec = myvec(_NX,0.0,0.0);
      }
      else { // from right to left electrode
        myvec distancevec = myvec(-1.0*_NX,0.0,0.0);
      }
    }
  }

  double distance = abs(distancevec);    
  double distancefactor = exp(-1.0*_alpha*distance);
  
  //second, calculate boltzmann factor (Coulomb interaction still have to be implemented)
   
  double init_energy;
  double final_energy;
    
  double inject_in_image_pot = 0.0;
  myvec finalpos = finalnode->_position;
  inject_in_image_pot = Calculate_Coulomb_Potential(finalpos.x(),myvec(0.0,0.0,0.0));
    
  if(cartype1 == Electron) {
    int carrier1_number = _initnode->_electron_number;
    Carrier* pcarrier1 = _electron_state.Get_item(carrier1_number);
    if (Injection_Type == NoInjection && Event_Type != ElectronCollection) {
      init_energy = initnode->_static_electron_energy + _coulomb_strength*(pcarrier1->_shortrange_coulomb[0] + charge*_Longrange->get_cached_longrange(initnode->_layer_index));
      final_energy = finalnode->_static_electron_energy + _coulomb_strength*(pcarrier1->_shortrange_coulomb[jump+1] + charge*_Longrange->get_cached_longrange(finalnode->_layer_index));
    }
    else if (Injection_Type == NoInjection && Event_Type == ElectronCollection) {
      init_energy = initnode->_static_electron_energy + injection_penalty + _coulomb_strength*(pcarrier1->_shortrange_coulomb[0] + charge*_Longrange->get_cached_longrange(initnode->_layer_index));
      final_energy = finalnode->_static_electron_energy;
    }
    else {
      init_energy = initnode->_static_electron_energy;
      final_energy = finalnode->_static_electron_energy + injection_penalty + _coulomb_strength*(finalnode->_el_injection_potential + charge*_Longrange->get_cached_longrange(finalnode->_layer_index) + inject_in_image_pot);
    }
    prefactor *= _electron_prefactor;
  }
  else if(cartype1 == Hole) {
    int carrier1_number = _initnode->_hole_number;
    Carrier* pcarrier1 = _hole_state.Get_item(carrier1_number);
    if (Injection_Type == NoInjection && Event_Type != HoleCollection) {
      init_energy = initnode->_static_hole_energy + _coulomb_strength*(pcarrier1->_shortrange_coulomb[0] + charge*_Longrange->get_cached_longrange(initnode->_layer_index));
      final_energy = finalnode->_static_hole_energy + _coulomb_strength*(pcarrier1->_shortrange_coulomb[jump+1] + charge*_Longrange->get_cached_longrange(finalnode->_layer_index));
    }
    else if (Injection_Type == NoInjection && Event_Type == HoleCollection) {
      init_energy = initnode->_static_hole_energy + injection_penalty + _coulomb_strength*(pcarrier1->_shortrange_coulomb[0] + charge*_Longrange->get_cached_longrange(initnode->_layer_index));
      final_energy = finalnode->_static_hole_energy;
    }
    else {
      init_energy = initnode->_static_hole_energy;
      final_energy = finalnode->_static_hole_energy + injection_penalty + _coulomb_strength*(finalnode->_ho_injection_potential + charge*_Longrange->get_cached_longrange(finalnode->_layer_index) + inject_in_image_pot);
    }
    prefactor *= _hole_prefactor;
  }
    
  double energycontrib;
  double energyfactor;
    
  if (_formalism == "Miller") {
    if(Event_type == Blocking) {
      energyfactor = 0.0; // Keep this construct here for eventual simulation of bipolaron formation for example
    }
    else if(Event_type == Recombination) {
      energycontrib = -1.0*_binding_energy + final_energy - init_energy;
      if (energycontrib>0.0) {
        energyfactor = exp(-1.0*_beta*energycontrib);
      }
      else {
        energyfactor = 1.0;
      }
      prefactor *= _recombination_prefactor;
    }
    else if((Event_type == ElectronTransfer)||(Event_type==HoleTransfer)) {
      energycontrib = final_energy - init_energy;
      if (energycontrib>0.0) {
        energyfactor = exp(-1.0*_beta*energycontrib);
      }
      else {
        energyfactor = 1.0;
      }
      prefactor *= _transfer_prefactor;
    }
    else if((Event_type == ElectronCollection)||(Event_type==HoleCollection)) {
      energycontrib = final_energy - init_energy;
      if (energycontrib>0.0) {
        energyfactor = exp(-1.0*_beta*energycontrib);
      }
      else {
        energyfactor = 1.0;
      }
      prefactor *= _collection_prefactor;
    }
  }
    
  double jump_rate = prefactor*distancefactor*energyfactor;
  return jump_rate;    
}

double Lattice::Compute_Coulomb_potential(double startx, myvec dif) {

  double coulpot;
  double RC = _coulomb_cutoff
  double RCSQR = RC*RC;
  int nimages = _number_of_shortrange_images;
   
  if(!_device) {
    coulpot = 1.0/abs(dif)-1.0/RC;
  }
  else {
    coulpot = 0.0;
    double distsqr_planar = dif.y()*dif.y() + dif.z()*dif.z();
    double distsqr = dif.x()*dif.x() + distsqr_planar;
    if (distsqr >= 0.5*_graph->_minimum_distance) {
      coulpot += 1.0/abs(dif)-1.0/RC;
    }
      
    int sign;
    double distx_1;
    double distx_2;
    double distancesqr_1;
    double distancesqr_2;
    bool outside_cut_off1 = false;
    bool outside_cut_off2 = false;
      
    while(!(outside_cut_off1&&outside_cut_off2)) {
      for (int i=0;i<nimages; i++) {
        if (div(i,2).rem==0) { // even generation
          sign = -1;
          distx_1 = i*_NX + 2*startx + dif.x();
          distx_2 = (i+2)*_NX - 2*startx - dif.x(); 
        }
        else {
          sign = 1;
          distx_1 = (i+1)*_NX + dif.x();
          distx_2 = (i+1)*_NX - dif.x();
        }
        distancesqr_1 = distx_1*distx_1 + distsqr_planar;
        if (distancesqr_1<=RCSQR) {
          coulpot += sign*1.0/sqrt(distancesqr_1)-1.0/(RC);
        }
        else {
          outside_cut_off1 = true;
        }
        distancesqr_2 = distx_2*distx_2 + distsqr_planar;
        if (distancesqr_2<=RCSQR) {
          coulpot += sign*1.0/sqrt(distancesqr_2)-1.0/(RC);
        }
        else {
          outside_cut_off2 = true;
        }
      }
    }
  }
      
  return coulpot;
}

myvec Lattice::periodic(myvec v) {

  while (v.y() < 0) v.y() += _NY;
  while (v.y() >= _NY) v.y() -= _NY;
  while (v.z() < 0) v.z() += NZ;
  while (v.z() >= _NZ) v.z() -= _NZ; 

  if(!_device) {
    while (v.x() < 0) v.x() += _NX;
    while (v.x() >= _NX) v.x() -= _NX;                
  }
      
  return &v;
}
   
myvec Lattice::perdifference(myvec init, myvec final) { // 

  myvec pre = final-init;
  double prex = pre.x();
  double prey = pre.y();
  double prez = pre.z();
  if(!_device) {
    if(prex<-0.5) {prex+=_NX;}
    if(prex>0.5) {prex-=_NX;}
  }
  if(prey<-0.5) {prey+=_NY;}
  if(prey>0.5) {prey-=_NY;}
  if(prez<-0.5) {prez+=_NZ;}
  if(prez>0.5) {prez-=NZ;}       
  myvec perdif = myvec(prex,prey,prez);

  return perdif;       
}

}}



#endif

