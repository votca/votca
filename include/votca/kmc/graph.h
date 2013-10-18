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

#ifndef __VOTCA_KMC_GRAPH_H_
#define __VOTCA_KMC_GRAPH_H_


#include <votca/tools/database.h>
#include <votca/tools/vec.h>
#include <votca/tools/random.h>
#include "node.h"
#include <list>

namespace votca { namespace kmc {
  
using namespace std;

typedef votca::tools::vec myvec;

enum CorrelationType {Correlated, Anticorrelated, Uncorrelated};
//enum NodeType {Normal, LeftElectrode, RightElectrode};

class Graph {
public:
  void LoadGraphNodes();
  void LoadGraphStaticEnergies();
  void LoadGraphPairs(); // Load only node pairs
  void LoadGraphTransfers();
  
  void CreateCubicGraphNodes();
  void CreateRandomStaticEnergies(votca::tools::Random2 *RandomVariable);
  void CalculateGraphPairs();
  
  void SetupInjectors();
   
  void InitializeLayers();

  vector<Node*> _nodes;   

  string _graph_filename;
  bool _load_energies;
  bool _load_pairs;

  double _NX;
  double _NY;
  double _NZ;

  double _hopping_distance;
  double _lattice_constant;
  
  double _site_box_dimX; //dimensions of boxes for pairfinding algorithm
  double _site_box_dimY;
  double _site_box_dimZ;
  int _NSITE_X; // number of boxes
  int _NSITE_Y;
  int _NSITE_Z;
  
  CorrelationType _correlation;
  double _disorder_strength;  
  double _disorder_factor;
  
  int _max_number_of_el_pairs;
  int _max_number_of_ho_pairs;
  int _max_number_of_exciton_pairs;
  
  double _minimum_distance; //minimum distance for a pair of nodes (needed for program)

  bool _device;
  double _left_electrode_distance;
  double _right_electrode_distance;
  
  vector<int> _left_injector_nodes;
  vector<int> _right_injector_nodes;
  int _number_of_left_injector_nodes;
  int _number_of_right_injector_nodes;
  int _number_of_injector_nodes;
  
  Node* _left_electrode;
  Node* _right_electrode;  

  vector<double> _positional_average; // positions of layers  

  int _number_of_nodes;
  
private:
  
  int _number_of_pairs;

  std::list<int> sitelist[_NSITE_X][_NSITE_Y][_NSITE_Z]; // List of nodes in each subbox  
  std::list<int> layerlist[_NSITE_X]; // list of nodes in each layer
  vector<double> _positional_sum;

  myvec perdifference(myvec v, myvec w); // returns difference vector (final node - initial node);
};

void Graph::LoadGraphNodes() {
    
  // Load nodes
  votca::tools::Database db;
  db.Open( _graph_filename );
  votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, posX, posY, posZ FROM segments;");
   
  int read_index = 0;
  while (stmt->Step() != SQLITE_DONE) {
    Node *newNode = new Node();
    _nodes.push_back(newNode);

    int node_id = stmt->Column<int>(0);
    _nodes[read_index]->_node_id = node_id;
    _nodes[read_index]->_node_type = Normal;
    myvec nodeposition;
    if(_device) { 
      nodeposition = myvec(stmt->Column<double>(1) + _left_electrode_distance, stmt->Column<double>(2), stmt->Column<double>(3)); // in nm
    }
    else {
      nodeposition = myvec(stmt->Column<double>(1), stmt->Column<double>(2), stmt->Column<double>(3)); // in nm
    }
    _nodes[read_index]->_position = nodeposition;

    read_index++;
  }
  delete stmt;
  _number_of_nodes = read_index;
}

void Graph::LoadGraphStaticEnergies() {
  
  // Load Static energies
  votca::tools::Database db;
  db.Open(_graph_filename);
  votca::tools::Statement *stmt = db.Prepare("SELECT UnCnNe, UnCnNh, UcNcCe, UcNcCh, eAnion,eNeutral,eCation,ucCnNe, ucCnNh FROM segments;");

  int read_index = 0;
  while (stmt->Step() != SQLITE_DONE) {
    _nodes[read_index]->_reorg_intorig_el = stmt->Column<double>(0); // UnCnNe
    _nodes[read_index]->_reorg_intorig_ho = stmt->Column<double>(1); // UnCnNh
    _nodes[read_index]->_reorg_intdest_el = stmt->Column<double>(2); // UcNcCe
    _nodes[read_index]->_reorg_intdest_ho = stmt->Column<double>(3); // UcNcCh
    
    double eAnion = stmt->Column<double>(4);
    _nodes[read_index]->_eAnion = eAnion;
    double eNeutral = stmt->Column<double>(5);
    _nodes[read_index]->_eNeutral = eNeutral;
    double eCation = stmt->Column<double>(6);
    _nodes[read_index]->_eCation = eCation;
    double internalenergy_el = stmt->Column<double>(7);
    _nodes[read_index]->_internalenergy_el = internalenergy_el;
    double internalenergy_ho = stmt->Column<double>(8);
    _nodes[read_index]->_internalenergy_ho = internalenergy_ho;
    
    _nodes[read_index]->_static_electron_energy = eAnion + internalenergy_el;
    _nodes[read_index]->_static_hole_energy = eAnion + internalenergy_ho;

    _nodes[read_index]->_ho_injection_potential = 0.0;
    _nodes[read_index]->_el_injection_potential = 0.0;    
    
    read_index++;
  }
  delete stmt;
}

void Graph::LoadGraphPairs() {
    
  // Load Node Pairs
  votca::tools::Database db;
  db.Open(_graph_filename);
  votca::tools::Statement *stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', drX, drY, drZ FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ' FROM pairs ORDER BY segment1;");

  _minimum_distance = sqrt(_NX*_NX+_NY*_NY+_NZ*_NZ);
  
  int read_index = 0;
  while (stmt->Step() != SQLITE_DONE) {
    int seg1 = stmt->Column<int>(0);
    int seg2 = stmt->Column<int>(1);

    
    // which nodes has these indices?
    
    bool seg1_found = false;
    int seg1_search_index = 0;
    bool seg2_found = false;
    int seg2_search_index = 0;
    while(!seg1_found) {
        seg1_search_index++;
        if(_nodes[seg1_search_index]->_node_id == seg1) { seg1_found = true; }
    }
    while(!seg2_found) {
        seg2_search_index++;
        if(_nodes[seg2_search_index]->_node_id == seg2) { seg2_found = true; }
    }

    myvec dr = myvec(stmt->Column<double>(2), stmt->Column<double>(3),stmt->Column<double>(4)); // distance is seg2 - seg1
    myvec pos1 = _nodes[seg1_search_index]->_position;

    if(!(_device && ((pos1.x + dr.x > _NX) || (pos1.x - dr.x < 0.0)))) {
      if(abs(dr)<_minimum_distance) {_minimum_distance = abs(dr);}

      _nodes[seg1_search_index]->setElPair(_nodes[seg2_search_index],dr);
      _nodes[seg1_search_index]->setHoPair(_nodes[seg2_search_index],dr);
      _nodes[seg1_search_index]->setExcitonPair(_nodes[seg2_search_index],dr);
    
      read_index++;
    }
  }    
  delete stmt;
  
  _number_of_pairs = read_index;
  
  // Determine number of pairs for each nodes and maximum number of pairs
  int max_el_number = 0;
  int max_ho_number = 0;
  int max_exciton_number = 0;
  
  for(node_index=0;node_index<_number_of_nodes;node_index++) {
    int el_number_of_pairs = _nodes[node_index]->_el_pair_nodes.size();
    int ho_number_of_pairs = _nodes[node_index]->_ho_pair_nodes.size();
    int exciton_number_of_pairs = _nodes[node_index]->_exciton_pair_nodes.size();
    
    _nodes[node_index]->_number_of_el_pairs = el_number_of_pairs;
    _nodes[node_index]->_number_of_ho_pairs = ho_number_of_pairs;
    _nodes[node_index]->_number_of_exciton_pairs = exciton_number_of_pairs;
    
    if(el_number_of_pairs>max_el_number) {max_el_number = el_number_of_neigbours;}
    if(ho_number_of_pairs>max_ho_number) {max_ho_number = ho_number_of_neigbours;}
    if(exciton_number_of_pairs>max_exciton_number) {max_exciton_number = exciton_number_of_neigbours;}    
  }
  
  _max_number_of_el_pairs = max_el_number;
  _max_number_of_ho_pairs = max_ho_number;
  _max_number_of_exciton_pairs = max_exciton_number;  
}

void Graph::LoadGraphTransfers() {
    
  // Load Node Pairs
  votca::tools::Database db;
  db.Open(_graph_filename);
  stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', rate12e AS 'rate_e', rate12h AS 'rate_h', Jeff2e, Jeff2h, lOe, l0h FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', rate21e AS 'rate_e', rate21h AS 'rate_h', Jeff2e, Jeff2h, l0e, l0h FROM pairs ORDER BY segment1;");

  while (stmt->Step() != SQLITE_DONE) {
    int seg1 = stmt->Column<int>(0);
    int seg2 = stmt->Column<int>(1);

    // which nodes has these indices?
    
    bool seg1_found = false;
    int seg1_search_index = 0;
    bool seg2_found = false;
    int seg2_search_index = 0;
    while(!seg1_found) {
        seg1_search_index++;
        if(_nodes[seg1_search_index]->_node_id == seg1) { seg1_found = true; }
    }
    while(!seg2_found) {
        seg2_search_index++;
        if(_nodes[seg2_search_index]->_node_id == seg2) { seg2_found = true; }
    }

    double rate12e = stmt->Column<double>(2);
    double rate12h = stmt->Column<double>(3);
    double Jeff2e = stmt->Column<double>(4);
    double Jeff2h = stmt->Column<double>(5);
    double reorg_oute = stmt->Column<double>(6); 
    double reorg_outh = stmt->Column<double>(7);
    
    _nodes[seg1_search_index]->setElTransfers(rate12e,Jeff2e,reorg_oute);
    _nodes[seg1_search_index]->setHoTransfers(rate12h,Jeff2h,reorg_outh);
  }    
  delete stmt;
}
          
void Graph::CreateCubicGraphNodes() {

  int node_index = 0;
  
  if(device) {
    _cubeNX = floor((_NX - _left_electrode_distance - _right_electrode_distance)/_lattice_constant)+1;
  }
  else {
    _cubeNX = floor(_NX/_lattice_constant);
  }
  
  _cubeNY = floor(_NY/_lattice_constant);
  _cubeNZ = floor(_NZ/_lattice_constant);
  
    
  for(int ix=0; ix<_cubeNX; ix++) {
    for(int iy=0; iy<_cubeNY; iy++) {
      for(int iz=0; iz<_cubeNZ; iz++) {
        Node *newNode = new Node();
        _nodes.push_back(newNode);

        _nodes[node_index]->_node_id = node_index;
        _nodes[read_index]->_nodetype = Normal;

        if(device) {
          nodeposition = myvec(ix*_lattice_constant+_left_electrode_distance,iy*_lattice_constant,iz*_lattice_constant); //positions in nm
        }
        else {
          nodeposition = myvec(ix*_lattice_constant,iy*_lattice_constant,iz*_lattice_constant);
        }
        _nodes[node_index]->_position = nodeposition;
        node_index++;    
      }
    }
  }
  _number_of_nodes = node_index;
}
    
void Graph::CreateRandomGaussianStaticEnergies(votca::tools::Random2 *RandomVariable) {
           
  for(int node_index=0;node_index<_number_of_nodes;node_index++) {
      double el_site_energy = RandomVariable->rand_gaussian(_disorder_strength);
      double ho_site_energy;
    _nodes[node_index]->_static_electron_energy = el_site_energy;
    if(_correlation == Correlated) {
        ho_site_energy = _disorder_factor*el_site_energy;
        _nodes[node_index]->_static_hole_energy = ho_site_energy;
    }
    else if(_correlation == Anticorrelated} {
      ho_site_energy = -1.0*_disorder_factor*el_site_energy;
      _nodes[node_index]->_static_hole_energy = ho_site_energy;
    }
    else {
       ho_site_energy = RandomVariable->rand_gaussian(_disorder_factor*_disorder_strength);
      _nodes[node_index]->_static_hole_energy = ho_site_energy);
    }
    _nodes[node_index]->_ho_injection_potential = 0.0;
    _nodes[node_index]->_el_injection_potential = 0.0;
  }
}
      
  
void Graph::CalculateGraphPairs() {  
  
  // define a sitelist for pairfinding algorithm
  
  for (int inode = 0; inode<_number_of_nodes; inode++) {
    double posx = _nodes[inode].x();
    double posy = _nodes[inode].y();
    double posz = _nodes[inode].z();
        
    int iposx = floor(posx/(_site_box_dimX)); 
    int iposy = floor(posy/(_site_box_dimY)); 
    int iposz = floor(posz/(_site_box_dimZ));
      
    sitelist[iposx][iposy][iposz].push_back(inode);
  }
  
  _minimum_distance = _hopping_distance;
  _number_of_pairs = 0;
  
  for (int inode = 0; inode<_number_of_nodes; inode++) {
    // Define cubic boundaries in non-periodic coordinates
    Node* initnode = _nodes[inode];
    myvec initnodepos = initnode->_position
    double ix1 = initnodepos.x()-hopdist; double ix2 = initnodepos.x()+hopdist;
    double iy1 = initnodepos.y()-hopdist; double iy2 = initnodepos.y()+hopdist;
    double iz1 = initnodepos.z()-hopdist; double iz2 = initnodepos.z()+hopdist;

    // Break periodicity in x-direction in case of device setting  
      
    if(device) {
      if(ix1<0) ix1=0.0;
      if(ix2>=NX) ix2 = _NX;
    }
      
    // Translate cubic boundaries to sublattice boundaries in non-periodic coordinates
    int sx1 = floor(ix1/_site_box_dimX);
    int sx2 = floor(ix2/_site_box_dimX);
    int sy1 = floor(iy1/_site_box_dimY);
    int sy2 = floor(iy2/_site_box_dimY);
    int sz1 = floor(iz1/_site_box_dimZ);
    int sz2 = floor(iz2/_site_box_dimZ);      
 
    // Now visit all relevant sublattices
    for (int isz=sz1; isz<=sz2; isz++) {
      int r_isz = isz;
      while (r_isz < 0) r_isz += _NSITE_Z;
      while (r_isz >= _NSITE_Z) r_isz -= _NSITE_Z;
      for (int isy=sy1; isy<=sy2; isy++) {
        int r_isy = isy;
        while (r_isy < 0) r_isy += _NSITE_Y;
        while (r_isy >= _NSITE_Y) r_isy -= _NSITE_Y;
        for (int isx=sx1; isx<=sx2; isx++) {
          int r_isx = isx;
          while (r_isx < 0) r_isx += _NSITE_X;
          while (r_isx >= _NSITE_X) r_isx -= _NSITE_X;
        
          // Ask a list of all sites in this sublattice
          std::list<int>::iterator li1,li2,li3;
          std::list<int> *pList = &sitelist[r_isx][r_isy][r_isz];
          li1 = pList->begin();
          li2 = pList->end();
          for (li3=li1; li3!=li2; li3++) {
            int probenode_index = *li3;
            if(inode!=probenode_index){ 
              Node* probenode = _nodes[probenode_index];
              myvec probenodepos = probenode->_position;
              myvec differ = perdifference(initnodepos,probenodepos);
              double distance = abs(differ);
              if(distance < _minimum_distance) {_minimum_distance = distance;}
              if(distance <=_hopping_distance) {
                _nodes[inode]->setElPair(probenode,differ);
                _nodes[inode]->setHoPair(probenode,differ);
                _nodes[inode]->setExcitonPair(probenode,differ);
                _number_of_pairs++;
              }
            }
          }
        }
      }
    }
  }
  
  // Determine number of pairs for each nodes and maximum number of pairs
  int max_el_number = 0;
  int max_ho_number = 0;
  int max_exciton_number = 0;
  
  for(node_index=0;node_index<_number_of_nodes;node_index++) {
    int el_number_of_pairs = _nodes[node_index]->_el_pair_nodes.size();
    int ho_number_of_pairs = _nodes[node_index]->_ho_pair_nodes.size();
    int exciton_number_of_pairs = _nodes[node_index]->_exciton_pair_nodes.size();
    
    _nodes[node_index]->_number_of_el_pairs = el_number_of_pairs;
    _nodes[node_index]->_number_of_ho_pairs = ho_number_of_pairs;
    _nodes[node_index]->_number_of_exciton_pairs = exciton_number_of_pairs;
    
    if(el_number_of_pairs>max_el_number) {max_el_number = el_number_of_neigbours;}
    if(ho_number_of_pairs>max_ho_number) {max_ho_number = ho_number_of_neigbours;}
    if(exciton_number_of_pairs>max_exciton_number) {max_exciton_number = exciton_number_of_neigbours;}    
  }
  
  _max_number_of_el_pairs = max_el_number;
  _max_number_of_ho_pairs = max_ho_number;
  _max_number_of_exciton_pairs = max_exciton_number;  
}
    

void Graph::SetupInjectors() {  

  _left_electrode->_node_type = LeftElectrode;
  _right_electrode->_node_type = RightElectrode;
  
  int read_injector_index = 0;
  int read_left_injector_index = 0;
  int read_right_injector_index = 0;
  int number_of_left_nodes = 0;
  int number_of_right_nodes = 0;
  
  for(inode=0;inode<_number_of_nodes;inode++) { 
      
    myvec nodeposition = _nodes[inode]->_position;
     
    if(nodeposition.x()<=_hopping_distance) {

      myvec dr = myvec(-1.0*nodeposition.x(),0.0,0.0);          
         
      _nodes[inode]->_left_injectable = true;
      _nodes[inode]->setElPair(_left_electrode,dr); //collection event
      _nodes[inode]->setHoPair(_left_electrode,dr); //collection event
      _nodes[inode]->setExcitonPair(_left_electrode,dr); //collection event

      int number_of_el_pairs = _nodes[inode]->_number_of_el_pairs;        
      int number_of_ho_pairs = _nodes[inode]->_number_of_ho_pairs;  
      int number_of_exciton_pairs = _nodes[inode]->_number_of_exciton_pairs;          
        
      _nodes[inode]->_number_of_el_pairs = number_of_el_pairs+1;
      _nodes[inode]->_number_of_ho_pairs = number_of_ho_pairs+1;
      _nodes[inode]->_number_of_exciton_pairs = number_of_exciton_pairs+1;
        
      if((number_of_el_pairs+1)>_max_number_of_el_pairs) {_max_number_of_el_pairs = number_of_el_pairs+1;}
      if((number_of_ho_pairs+1)>_max_number_of_ho_pairs) {_max_number_of_ho_pairs = number_of_ho_pairs+1;}
      if((number_of_exciton_pairs+1)>_max_number_of_exciton_pairs) {_max_number_of_exciton_pairs = number_of_exciton_pairs+1;}
        
      _number_of_pairs++;

      _nodes[inode]->_injectorindex = read_injector_index;
      _nodes[inode]->_left_injectorindex = read_left_injector_index;
      _left_electrode->setElPair(_nodes[inode],-1.0*dr);
      _left_electrode->setHoPair(_nodes[inode],-1.0*dr);
      read_left_injector_index++;
      read_injector_index++;

      number_of_left_nodes++;        
      _left_injector_nodes.pushback(inode);
    }
      
    if((NX-nodeposition.x())<=_hopping_distance) {
          
      myvec dr = myvec((NX-nodeposition.x()),0.0,0.0);
         
      _nodes[inode]->_right_injectable = true;
      _nodes[inode]->setElPair(_right_electrode,dr); //collection event
      _nodes[inode]->setHoPair(_right_electrode,dr); //collection event
      _nodes[inode]->setExcitonPair(_right_electrode,dr); //collection event

      int number_of_el_pairs = _nodes[inode]->_number_of_el_pairs;        
      int number_of_ho_pairs = _nodes[inode]->_number_of_ho_pairs;  
      int number_of_exciton_pairs = _nodes[inode]->_number_of_exciton_pairs;          
        
      _nodes[inode]->_number_of_el_pairs = number_of_el_pairs+1;
      _nodes[inode]->_number_of_ho_pairs = number_of_ho_pairs+1;
      _nodes[inode]->_number_of_exciton_pairs = number_of_exciton_pairs+1;

      if((number_of_el_pairs+1)>_max_number_of_el_pairs) {_max_number_of_el_pairs = number_of_el_pairs+1;}
      if((number_of_ho_pairs+1)>_max_number_of_ho_pairs) {_max_number_of_ho_pairs = number_of_ho_pairs+1;}
      if((number_of_exciton_pairs+1)>_max_number_of_exciton_pairs) {_max_number_of_exciton_pairs = number_of_exciton_pairs+1;}
        
      _number_of_pairs++;

      _nodes[inode]->_injectorindex = read_injector_index;
      _nodes[inode]->_right_injectorindex = read_right_injector_index;
      _right_electrode->setElPair(_nodes[inode],-1.0*dr);
      _right_electrode->setHoPair(_nodes[inode],-1.0*dr);
      read_right_injector_index++;
      read_injector_index++;

      number_of_right_nodes++;        
      _right_injector_nodes.pushback(inode);
    }
  }
  _number_of_left_injector_nodes = number_of_left_nodes;
  _number_of_right_injector_nodes = number_of_right_nodes;
  _number_of_injection_events = number_of_left_nodes + number_of_right_nodes;
}
 
void Graph::InitializeLayers() {

  _positional_sum.resize(_NSITE_X);
  _positional_average.resize(_NSITE_X);
    
  for (int ix = 0;ix<_NSITE_X;ix++){
    for (int iy =0;iy<_NSITE_Y;iy++){
      for (int iz = 0;iz<_NSITE_Z;iz++){
        std::list<int>::iterator li1,li2,li3;
        std::list<int> *pList = &sitelist[ix][iy][iz];
        li1 = pList->begin();
        li2 = pList->end();
        for (li3=li1; li3!=li2; li3++) {
          int node_index = *li3;
          layerlist[ix].push_back(node_index);
        }
      }
    }
  }
    
  for (int ix =0 ;ix<_NSITE_X; ix++){
    std::list<int>::iterator li1,li2,li3;
    std::list<int> *pList = &layerlist[ix];
    li1 = pList->begin();
    li2 = pList->end();
    double x_position = 0.0;
    int site_counter = 0;
    for (li3=li1; li3!=li2; li3++) {
      int node_index = *li3;
      Node* probenode = _graph->_nodes[node_index];
      probenode->_layer_index = ix;
      x_position +=probenode->_position.x();
      site_counter++;
    }
    _positional_average[ix]=x_position/(1.0*site_counter);
  }
}

myvec Graph::perdifference(myvec init, myvec final) { // 
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
  if(prez>0.5) {prez-=_NZ;}       
  myvec perdif = myvec(prex,prey,prez);
  return perdif;       
}          

}}



#endif

