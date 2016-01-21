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

#ifndef __VOTCA_KMC_GRAPHKMC_H_
#define __VOTCA_KMC_GRAPHKMC_H_

#include <vector>
#include <votca/xtp/graphsql.h>
#include <votca/xtp/nodedevice.h>
#include <votca/xtp/linkdevice.h>
#include <votca/xtp/eventinfo.h>

namespace votca { namespace xtp {

enum NodeType{NormalNode, LeftElectrodeNode, RightElectrodeNode};
    
class GraphKMC : public GraphSQL<NodeDevice,LinkDevice> {

public:

    ///setup graph for bulk simulations
    void Setup_bulk_graph(Eventinfo* eventinfo);
    
    ///setup graph for device simulations, define electrode nodes and form links between those nodes and neighbouring nodes and set maxpairdegree/hopping_distance/sim_box_size
    void Setup_device_graph(Eventinfo* eventinfo);       
  
    /// max_pair_degree
    const int &maxpairdegree() const { return _max_pair_degree; }

    /// hopping distance
    const double &hopdist() const { return _hop_distance; }
    
    /// minimum distance
    const double &mindist() const { return _min_distance; }

    /// simulation box size
    const votca::tools::vec &simboxsize() const { return _sim_box_size; }
    
    /// left electrode node
    NodeDevice* &left() { return _left_electrode; }
    
    /// right electrode node
    NodeDevice* &right() { return _right_electrode; }

    /// total sum of r12.x
    const double &total_link_distance_z() const {return _total_link_distance_z;}    
    
    /// average hole node energy
    double Average_hole_node_energy();
    
    double Average_hole_left_electrode_energy(); 
    double Average_hole_right_electrode_energy(); 
    
    /// variance hole node energy
    double stddev_hole_node_energy();
    
    double stddev_hole_left_electrode_energy();
    double stddev_hole_right_electrode_energy();
    
    /// average electron node energy
    double Average_electron_node_energy();

    /// variance hole node energy
    double stddev_electron_node_energy();
    
    /// initialize output values (current for example)
    void Initialize_output_values();

    double Electron_inject_reorg();
    double Hole_inject_reorg();

    ///calculate the simulation box size (also determines crossing types of links)
    votca::tools::vec Determine_Sim_Box_Size();     
    
private:

    ///calculate hopping_distance (maximum distance between a pair of nodes) ... needed for injection and coulomb potential calculations
    double Determine_Hopping_Distance();    
    
    ///calculate minimum distance (minimum distance between a pair of nodes)
    double Determine_Minimum_Distance();    

    ///calculate the maximum of all degrees in the graph
    unsigned Determine_Max_Pair_Degree();

    ///gives the total sum of r12.x() for all links (including injection and collection links)
    double Sum_of_link_distances_z();
    
    ///determine the crossing types of the links    
    void Determine_cross_types();    
    
    ///translate graph by given vector
    void Translate_graph(double translate_x, double translate_y, double translate_z);

    ///translate graph in such a way that the minimum coordinates are at 0
    void Put_at_zero_graph();    

    ///reperiodicize graph (put all nodes in box)
    void Push_in_box();

    ///resize the graph by periodically repeating the morphology and depending on the bool flags, breaking of (deleting nodes and boundary crossing links)
    void Resize(double dimX, bool breakX, double dimY, bool breakY, double dimZ, bool breakZ);

    ///break the periodicity of the graph (breaking boundary crossing pairs) .. (run before linksort)
    void Break_periodicity(bool break_x, bool break_y, bool break_z);    

    /// initialize node types to normal type
    void Initialize_node_types();

    ///add electrode nodes
    void Add_electrodes();
    
    ///associate all links in links vector to the corresponding nodes
    void LinkSort();

    /// set self-image coulomb potential on all nodes and links
    void Set_Self_Image_Coulomb_Potential(double device_length,Eventinfo* eventinfo);    

    /// attach layer indices to nodes
    void Set_Layer_indices(Eventinfo* eventinfo);

    /// renumber the Id's of all links
    void Renumber_id(); 

    /// 'roll' morphology, so injection is in a different part of the morphology
    void Roll_morph(double translate_x, double translate_y, double translate_z);     
    
    int _max_pair_degree;
    
    double _hop_distance;
    double _min_distance;
    double _av_el_distance;
    
    votca::tools::vec _sim_box_size;

    NodeDevice* _left_electrode;
    NodeDevice* _right_electrode;
    
    double _total_link_distance_z;
    
};



}}

#endif

