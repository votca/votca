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

#ifndef __VOTCA_KMC_GRAPHDEVICE_H_
#define __VOTCA_KMC_GRAPHDEVICE_H_

#include <vector>
#include <votca/kmc/graphsql.h>
#include <votca/kmc/nodedevice.h>
#include <votca/kmc/linkdevice.h>
#include <votca/kmc/eventinfo.h>

namespace votca { namespace kmc {

enum NodeType{NormalNode, LeftElectrodeNode, RightElectrodeNode};
    
class GraphDevice : public GraphSQL<NodeDevice,LinkDevice> {

public:

    ///setup graph for bulk simulations
    void Setup_bulk_graph(bool resize, Eventinfo* eventinfo);
    
    ///define electrode nodes and form links between those nodes and neighbouring nodes and set maxpairdegree/hopping_distance/sim_box_size
    void Setup_device_graph(double left_distance, double right_distance, bool resize, Eventinfo* eventinfo);       

    ///determine the crossing types of the links    
    void Determine_cross_types();
    
    ///determine injection nodes for Ohmic contacts
    //void Determine_source_nodes();
    
    ///translate graph in such a way that the minimum coordinates are at 0
    void Put_at_zero_graph();
    
    ///reperiodicize graph
    void Push_in_box();
    
    ///resize the graph (by periodically repeating the morphology)
    void Resize(int dimX, int dimY, int dimZ);
    
    ///translate graph to accomodate device geometry
    void Translate_graph(double translate_x, double translate_y, double translate_z);
    
    ///add electrode nodes
    void Add_electrodes();
    
    ///associate all links in links vector to the corresponding nodes
    void LinkSort();
    
    ///break the periodicity of the graph (breaking boundary crossing pairs) .. (run before linksort)
    void Break_periodicity(bool break_x, double dimX, bool break_y, double dimY, bool break_z, double dimZ);

    ///break the periodicity of the graph (breaking boundary crossing pairs) .. (run before linksort)
    void Break_periodicity(bool break_x, bool break_y, bool break_z);    
    
    ///calculate the maximum of all degrees in the graph
    int Determine_Max_Pair_Degree();
    
    ///calculate hopping_distance (maximum distance between a pair of nodes) ... needed for injection and coulomb potential calculations
    double Determine_Hopping_Distance();
    
    ///calculate the simulation box size
    votca::tools::vec Determine_Sim_Box_Size();   
    
    /// max_pair_degree
    const int &maxpairdegree() const { return _max_pair_degree; }

    /// hopping distance
    const double &hopdist() const { return _hop_distance; }

    /// simulation box size
    const votca::tools::vec &simboxsize() const { return _sim_box_size; }
    
    /// left electrode node
    NodeDevice* &left() { return _left_electrode; }
    
    /// right electrode node
    NodeDevice* &right() { return _right_electrode; }

    /// set self-image coulomb potential on all nodes and links
    void Set_Self_Image_Coulomb_Potential(double device_length,Eventinfo* eventinfo);
    
    /// renumber the Id's of all links
    void RenumberId();      
     
    /// attach layer indices to nodes
    void Set_Layer_indices(Eventinfo* eventinfo);
    
    /// initialize node types to normal type
    void Init_node_types();
    
    /// 'roll' morphology, so injection is in a different part of the morphology
    void Roll_morph(double translate_x, double translate_y, double translate_z); 
    
    /// rotate morphology (if rotate_axis = 0, rotation around z axis, if rotate_axis = 1, rotation around y axis)
    // void Rotate_morph(double rotate_in_rads, int rotate_axis);
    
    double Av_hole_node_energy();
    double Av_electron_node_energy();
    
private:

    int _max_pair_degree;  
    double _hop_distance;
    votca::tools::vec _sim_box_size;

    NodeDevice* _left_electrode;
    NodeDevice* _right_electrode;
    
    
    
};

double GraphDevice::Av_hole_node_energy() {
    
    double ho_energy = 0.0;
    double av_ho_energy = 0.0;
    
    typename std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        ho_energy += (*it)->eAnion() + (*it)->UcCnNh();
    }
    
    av_ho_energy = ho_energy/this->Numberofnodes();
    return av_ho_energy;
}

double GraphDevice::Av_electron_node_energy() {
    
    double elect_energy = 0.0;
    double av_elect_energy = 0.0;
    
    typename std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        elect_energy += (*it)->eCation() + (*it)->UcCnNe();
    }
    
    av_elect_energy = elect_energy/this->Numberofnodes();
    return av_elect_energy;
}

void GraphDevice::Translate_graph(double translate_x, double translate_y, double translate_z) {
    
    typename std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        votca::tools::vec oldpos = (*it)->position();
        double newxpos = oldpos.x() + translate_x;
        double newypos = oldpos.y() + translate_y; 
        double newzpos = oldpos.z() + translate_z;
        votca::tools::vec newpos = votca::tools::vec(newxpos,newypos,newzpos);
        (*it)->SetPosition(newpos);
    }
}

void GraphDevice::Setup_bulk_graph( bool resize, Eventinfo* eventinfo){

    // Determine hopping distance
    _hop_distance = this->Determine_Hopping_Distance();

    // Determine boundary crossing types            
    this->Determine_cross_types();
    
    // Detemine simulation box size
    _sim_box_size = this->Determine_Sim_Box_Size();            

    // Push nodes back in box (crossing types are changed by this operation, so re-evaluate)    
    this->Push_in_box();
    this->Determine_cross_types();    

    // Resize by copying the box (crossing types are changed by this operation, so re-evaluate)    
    if(resize) {this->Resize(eventinfo->size_x, eventinfo->size_y, eventinfo->size_z); this->Determine_cross_types(); _sim_box_size = this->Determine_Sim_Box_Size();}
    
    //set node types for existing nodes as Normal
    this->Init_node_types();

    // associate links in links vector with the corresponding nodes
    this->LinkSort();
    
    // determine maximum degree of graph
    _max_pair_degree = this->Determine_Max_Pair_Degree();    
    
    // Set the layer indices on every node
    this->Set_Layer_indices(eventinfo);
  
    // clear the occupation of the graph
    this->Clear();

    //renumber link id's
    this->RenumberId();   

}    

void GraphDevice::Setup_device_graph(double left_distance, double right_distance, bool resize, Eventinfo* eventinfo){
    
    // Determine hopping distance before breaking periodicity
    _hop_distance = this->Determine_Hopping_Distance();
    
    // Make sure injection hops are possible
    if(left_distance > _hop_distance) {_hop_distance = left_distance;}
    if(right_distance > _hop_distance) {_hop_distance = right_distance;}

    // Determine boundary crossing types
    this->Determine_cross_types();
  
    // Determine simulation box size
    _sim_box_size = this->Determine_Sim_Box_Size();
  
    // Translate the graph due to the spatial location of the electrodes and update system box size accordingly, putting the left electrode at x = 0
    // left_electrode_distance is the distance of the left electrode to the node with minimum x-coordinate
    // distance by which the graph should be translated is left_electrode_distance - minX

    // Translate graph so that minimum coordinates are at 0
    this->Put_at_zero_graph();

    // Push nodes back in box (crossing types are changed by this operation, so re-evaluate)
    this->Push_in_box();
    this->Determine_cross_types();
    
    // Resize by copying the box (crossing types are changed by this operation, so re-evaluate)
    if(resize) {this->Resize(eventinfo->size_x, eventinfo->size_y, eventinfo->size_z); this->Determine_cross_types(); _sim_box_size = this->Determine_Sim_Box_Size();}
  
    // Break periodicity
    if(resize) {
        this->Break_periodicity(true,eventinfo->size_x,false,eventinfo->size_y,false,eventinfo->size_z);
    }
    else {
        this->Break_periodicity(true,false,false);
    }
    
    // Recalculate simulation box size after periodicity breaking
    _sim_box_size = this->Determine_Sim_Box_Size();

    // Translate graph to accomodate for device geometry
    this->Translate_graph(left_distance, 0.0, 0.0);    

    //adjust simulation box size accordingly to given electrode distances
    
    votca::tools::vec old_sim_box_size = _sim_box_size;
    double new_sim_box_sizeX = old_sim_box_size.x() + left_distance + right_distance;
     _sim_box_size =  votca::tools::vec(new_sim_box_sizeX, old_sim_box_size.y(), old_sim_box_size.z());

    //set node types for existing nodes as Normal
    this->Init_node_types();
 
    this->Add_electrodes();
 
    // associate links in links vector with the corresponding nodes
    this->LinkSort();
    
    // determine maximum degree of graph
    _max_pair_degree = this->Determine_Max_Pair_Degree();
   
    // Set the self-image coulomb potential on every node
    this->Set_Self_Image_Coulomb_Potential(_sim_box_size.x(), eventinfo);

    // Set the layer indices on every node
    this->Set_Layer_indices(eventinfo);
  
    // clear the occupation of the graph
    this->Clear();

    //renumber link id's
    this->RenumberId(); 
   
}

void GraphDevice::Init_node_types() {
    typename std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) (*it)->SetType((int) NormalNode);    
}

/* void GraphDevice::Determine_source_nodes() {

    //determine the nodes which are injectable from the left side and the nodes which are injectable from the right side

    typename std::vector<NodeDevice*>::iterator it;    
    for(it  = this->_nodes.begin(); it != this->_nodes.end(); it++) { 
   
        //default
        (*it)->SetInjectable(false);
        
        votca::tools::vec nodepos = (*it)->position();
        double left_distance = nodepos.x();
        int linkID = this->_links.size();

        if(left_distance <= _hop_distance) (*it)->SetInjectable(true);
        double right_distance = _sim_box_size.x() - nodepos.x();
        if(right_distance <= _hop_distance) (*it)->SetInjectable(true);     

    }    
    
}*/

void GraphDevice::Add_electrodes() {
    
    //introduce the electrode nodes (might make a special electrode node header file for this)
    _left_electrode = new NodeDevice(-1, tools::vec(0.0,0.0,0.0));
    _right_electrode = new NodeDevice(-2, tools::vec(_sim_box_size.x(),0.0,0.0));
    _left_electrode->SetType((int) LeftElectrodeNode);
    _right_electrode->SetType((int) RightElectrodeNode);

    //determine the nodes which are injectable from the left electrode and the nodes which are injectable from the right electrode

    typename std::vector<NodeDevice*>::iterator it;    
    for(it  = this->_nodes.begin(); it != this->_nodes.end(); it++) { 
      
        // when electrodes are present physically, no need to keep track of injectability of the nodes
        //(*it)->SetInjectable(false);
        
        votca::tools::vec nodepos = (*it)->position();
        double left_distance = nodepos.x();
        int linkID = this->_links.size();

        if(left_distance <= _hop_distance) {
            votca::tools::vec dr = votca::tools::vec(-1.0*left_distance,0.0,0.0);   
            LinkSQL* newLinkCollect = this->AddLink(linkID,(*it), _left_electrode, dr); 
            linkID++;
            LinkSQL* newLinkInject = new LinkSQL(linkID, _left_electrode, (*it), -1.0*dr);
            _left_electrode->AddLink(newLinkInject);
        }
      
        double right_distance = _sim_box_size.x() - nodepos.x();
        if(right_distance <= _hop_distance) {
            votca::tools::vec dr = votca::tools::vec(right_distance,0.0,0.0);   
            LinkSQL* newLinkCollect = this->AddLink(linkID,(*it), _right_electrode, dr); 
            linkID++;
            LinkSQL* newLinkInject = new LinkSQL(linkID, _right_electrode, (*it), -1.0*dr);
            _right_electrode->AddLink(newLinkInject);
        }
    }

    this->AddNode(_left_electrode); //in this way the electrode nodes are caught by destructor
    this->AddNode(_right_electrode);
    _left_electrode->RemoveCarrier();
    _right_electrode->RemoveCarrier();    

}

void GraphDevice::LinkSort(){
    
    typename std::vector<LinkDevice*>::iterator it;
    for (it = this->_links.begin(); it != this->_links.end(); it++ ) {
        NodeDevice* node1 = dynamic_cast<NodeDevice*>((*it)->node1());
        if(node1->type() == NormalNode) node1->AddLink((*it));
    }
    
}

int GraphDevice::Determine_Max_Pair_Degree(){
    
    int max_pair_degree = 0;
    typename std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        if((*it)->type() == NormalNode) {
            if((*it)->links().size()>max_pair_degree) max_pair_degree = (*it)->links().size();
        }
    }
    return max_pair_degree;
    
}

double GraphDevice::Determine_Hopping_Distance(){
    
    double hop_distance = 0.0;
    typename std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) {
        votca::tools::vec dR = (*it)->r12();
        double distance = abs(dR);
        if(distance>hop_distance) hop_distance = distance;
    }
    return hop_distance;
}


votca::tools::vec GraphDevice::Determine_Sim_Box_Size(){ 
    
    votca::tools::vec sim_box_size;
    //Determination of simulation box size
    //Note that it is possible that none of the pairs pass the simulation box boundaries
    //In this special case, we must determine the node with max x/y/z coordinate and min x/y/z coordinate
    
    bool pairXfound = false; bool pairYfound = false; bool pairZfound = false;
    
    //for determination of initial value of simboxsize
    bool initXfound = false; bool initYfound = false; bool initZfound = false;
    double newX; double newY; double newZ;
    
    //dimensions of the simulation box
    double simX; double simY; double simZ;
    
    //initial values for maximum and minimum coordinates
    votca::tools::vec initpos = this->_nodes[0]->position();
    double maxX = initpos.x(); double maxY = initpos.y(); double maxZ = initpos.z();
    double minX = initpos.x(); double minY = initpos.y(); double minZ = initpos.z();
    
    typename std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) {
        
        votca::tools::vec pos1 = (*it)->node1()->position();
        votca::tools::vec pos2 = (*it)->node2()->position();
        votca::tools::vec dr = (*it)->r12();

        if(maxX<pos1.x()) {maxX = pos1.x();}   if(minX>pos1.x()) {minX = pos1.x();}
        if(maxY<pos1.y()) {maxY = pos1.y();}   if(minY>pos1.y()) {minY = pos1.y();}
        if(maxZ<pos1.z()) {maxZ = pos1.z();}   if(minZ>pos1.z()) {minZ = pos1.z();}

        //theoretical possibility that hopping distance is larger than the actual simulation box size (pair crosses boundaries more than once), in which the value of simboxsize is to large
        //check for minimum to catch these exceptional cases
        
        if((*it)->crossxtype() == (int) PosxCross ) {pairXfound = true; newX = pos1.x() + dr.x() - pos2.x(); if(!initXfound) {initXfound = true; simX = newX;} else if(simX<newX) { simX = newX;} }
        if((*it)->crossxtype() == (int) NegxCross ) {pairXfound = true; newX = pos2.x() - dr.x() - pos1.x(); if(!initXfound) {initXfound = true; simX = newX;} else if(simX<newX) { simX = newX;} }
        if((*it)->crossytype() == (int) PosyCross ) {pairYfound = true; newY = pos1.y() + dr.y() - pos2.y(); if(!initYfound) {initYfound = true; simY = newY;} else if(simY<newY) { simY = newY;} }
        if((*it)->crossytype() == (int) NegyCross ) {pairYfound = true; newY = pos2.y() - dr.y() - pos1.y(); if(!initYfound) {initYfound = true; simY = newY;} else if(simY<newY) { simY = newY;} }
        if((*it)->crossztype() == (int) PoszCross ) {pairZfound = true; newZ = pos1.z() + dr.z() - pos2.z(); if(!initZfound) {initZfound = true; simZ = newZ;} else if(simZ<newZ) { simZ = newZ;} }
        if((*it)->crossztype() == (int) NegzCross ) {pairZfound = true; newZ = pos2.z() - dr.z() - pos1.z(); if(!initZfound) {initZfound = true; simZ = newZ;} else if(simZ<newZ) { simZ = newZ;} }
        
    }
    
    //for the possible outcome that none of the pairs are crossing the simulation box boundary
    if(!pairXfound) {simX = maxX-minX;}
    if(!pairYfound) {simY = maxY-minY;}
    if(!pairZfound) {simZ = maxZ-minZ;}

    sim_box_size = votca::tools::vec(simX,simY,simZ);
    return sim_box_size;
}

void GraphDevice::Determine_cross_types(){

    int number_of_x_crossing_links = 0;
    int number_of_y_crossing_links = 0;
    int number_of_z_crossing_links = 0;
    
    typename std::vector<LinkDevice*>::iterator it;
    for(it = this->_links.begin(); it != this->_links.end(); it++) {

        votca::tools::vec pos1 = (*it)->node1()->position();
        votca::tools::vec pos2 = (*it)->node2()->position();
        votca::tools::vec dr = (*it)->r12();

        int crossx_flag = (int) NoxCross;
        int crossy_flag = (int) NoyCross;
        int crossz_flag = (int) NozCross;

        if(pos2.x()-pos1.x() < 0.0  && dr.x() > 0.0 ) { crossx_flag = (int) PosxCross; number_of_x_crossing_links++; }
        if(pos2.x()-pos1.x() > 0.0  && dr.x() < 0.0 ) { crossx_flag = (int) NegxCross; number_of_x_crossing_links++; }  
        if(pos2.y()-pos1.y() < 0.0  && dr.y() > 0.0 ) { crossy_flag = (int) PosyCross; number_of_y_crossing_links++; }  
        if(pos2.y()-pos1.y() > 0.0  && dr.y() < 0.0 ) { crossy_flag = (int) NegyCross; number_of_y_crossing_links++; }
        if(pos2.z()-pos1.z() < 0.0  && dr.z() > 0.0 ) { crossz_flag = (int) PoszCross; number_of_z_crossing_links++; }  
        if(pos2.z()-pos1.z() > 0.0  && dr.z() < 0.0 ) { crossz_flag = (int) NegzCross; number_of_z_crossing_links++; }
        
        (*it)->setCrossxType(crossx_flag);
        (*it)->setCrossyType(crossy_flag);
        (*it)->setCrosszType(crossz_flag);
        
    }
}

void GraphDevice::Put_at_zero_graph(){

    // Put min x at 0
    
    votca::tools::vec pos = this->_nodes[0]->position(); // initial value 
    double minX = pos.x();
    double minY = pos.y();
    double minZ = pos.z();

    typename std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        pos = (*it)->position(); 
        if(pos.x() < minX) {minX = pos.x();}
        if(pos.y() < minY) {minY = pos.y();}
        if(pos.z() < minZ) {minZ = pos.z();}
    }

   //distance by which the graph should be translated is left_electrode_distance - minX
 
    double xtranslate = -1.0*minX;
    double ytranslate = -1.0*minY;
    double ztranslate = -1.0*minZ;

    this->Translate_graph(xtranslate,ytranslate,ztranslate);
}

void GraphDevice::Push_in_box(){
   
    //checks which nodes are falling outside the simboxsize
    //periodically determine where they are in the simulation box

    int push_in_box_ops = 0;
    
    typename std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        
        votca::tools::vec pos = (*it)->position();
        
        double newxpos = pos.x(); double newypos = pos.y(); double newzpos = pos.z();
        
        while(newxpos>_sim_box_size.x()) { newxpos -= _sim_box_size.x(); push_in_box_ops++; }
        while(newxpos<0.0)               { newxpos += _sim_box_size.x(); push_in_box_ops++; }
        while(newypos>_sim_box_size.y()) { newypos -= _sim_box_size.y(); push_in_box_ops++; }
        while(newypos<0.0)               { newypos += _sim_box_size.y(); push_in_box_ops++; }
        while(newzpos>_sim_box_size.z()) { newzpos -= _sim_box_size.z(); push_in_box_ops++; }        
        while(newzpos<0.0)               { newzpos += _sim_box_size.z(); push_in_box_ops++; }
        
        votca::tools::vec newpos = votca::tools::vec(newxpos,newypos,newzpos);
        (*it)->SetPosition(newpos);

    }
    
}

void GraphDevice::Break_periodicity(bool break_x, double dimX, bool break_y, double dimY, bool break_z, double dimZ){

    // Break crossing links
    for(int it = this->_links.size()-1; it >= 0; it--) {

        LinkDevice* ilink = this->_links[it];
        
        votca::tools::vec pos1 = ilink->node1()->position();
        votca::tools::vec pos2 = ilink->node2()->position();
        votca::tools::vec dr = ilink->r12();

        bool remove_flag = false;
        
        // remove crossing links and links which are connected to nodes which will fall out of the simulation box
        
        if(break_x){ if(ilink->crossxtype() != (int) NoxCross || (pos1.x() > dimX || pos2.x() > dimX)) { remove_flag = true; } }        
        if(break_y){ if(ilink->crossytype() != (int) NoyCross || (pos1.y() > dimY || pos2.x() > dimX)) { remove_flag = true; } } 
        if(break_z){ if(ilink->crossztype() != (int) NozCross || (pos1.z() > dimZ || pos2.x() > dimX)) { remove_flag = true; } }  
        
        if(remove_flag) {
            this->RemoveLink(it);
            delete ilink;   
        }
        //std::cout << it << endl;
    }
  
    // Remove nodes
    
    for(int it = this->_nodes.size()-1; it >= 0; it--) {
        
        NodeDevice* inode = this->_nodes[it];
        
        votca::tools::vec pos = inode->position();
        double xpos = pos.x(); double ypos = pos.y(); double zpos = pos.z();
        
        bool remove_flag = false;        
        
        if(break_x && xpos > dimX) { remove_flag = true; }
        if(break_y && ypos > dimY) { remove_flag = true; }
        if(break_z && zpos > dimZ) { remove_flag = true; }        
 
        if(remove_flag) {
            this->RemoveNode(it);
            delete inode;
        }

    }    
    
}

void GraphDevice::Break_periodicity(bool break_x, bool break_y, bool break_z){

    // Break crossing links
    for(int it = this->_links.size()-1; it >= 0; it--) {

        LinkDevice* ilink = this->_links[it];
        
        bool remove_flag = false;
        
        // remove crossing links and links which are connected to nodes which will fall out of the simulation box
        
        if(break_x && ilink->crossxtype() != (int) NoxCross) { remove_flag = true; }        
        if(break_y && ilink->crossytype() != (int) NoyCross) { remove_flag = true; } 
        if(break_z && ilink->crossztype() != (int) NozCross) { remove_flag = true; }  
        
        if(remove_flag) {
            this->RemoveLink(it);
            delete ilink;   
        }
        
    }   
    
}

void GraphDevice::Resize(int dimX, int dimY, int dimZ) {

    
    int repeatX = ceil(dimX/_sim_box_size.x());
    int repeatY = ceil(dimY/_sim_box_size.y());
    int repeatZ = ceil(dimZ/_sim_box_size.z());
 
    int number_of_nodes = this->Numberofnodes();

    // repeat nodes
 

    for(int ix = 0; ix<repeatX; ix++){
        for(int iy = 0; iy<repeatY; iy++){
            for(int iz = 0; iz<repeatZ; iz++){
                if(!((ix==0)&&((iy==0)&&(iz==0)))) {
                    for(int it = 0; it < number_of_nodes; it++) {
                        NodeDevice* probenode = this->GetNode(it);
                        int node_id = probenode->id();
                        votca::tools::vec node_pos = probenode->position();
                        // remapping of the positions

                        double posX = node_pos.x() + ix*_sim_box_size.x();
                        double posY = node_pos.y() + iy*_sim_box_size.y();
                        double posZ = node_pos.z() + iz*_sim_box_size.z();
                        votca::tools::vec new_node_pos = votca::tools::vec(posX,posY,posZ);
                        // remapping of the indices
                        int new_node_id = node_id + (iz+iy*repeatZ+ix*repeatY*repeatZ)*number_of_nodes;
                        // add node to nodes vector
                        NodeDevice* newNodeDevice = this->AddNode(new_node_id,new_node_pos);
                        Node* testnode = this->GetNode(new_node_id);

                        // copy data to the periodically repeated nodes
                        newNodeDevice->setU(probenode->UnCnNe(), probenode->UnCnNh(),   probenode->UcNcCe(), probenode->UcNcCh());
                        newNodeDevice->setE(probenode->eAnion(), probenode->eNeutral(), probenode->eCation());
                        newNodeDevice->setu(probenode->UcCnNe(), probenode->UcCnNh());
                    }
                }
            }
        }
    }
    
    
    // repeat links
    
    long number_of_links = this->Numberoflinks();
    
    for(long ilink = 0; ilink < number_of_links; ilink++) {
        
        LinkDevice* probelink = this->GetLink(ilink);
        long link_id = ilink;
        int node1_id = probelink->node1()->id();
        int node2_id = probelink->node2()->id();
        NodeDevice* node1 = this->GetNode(node1_id);
        NodeDevice* node2 = this->GetNode(node2_id);

        votca::tools::vec link_r12 = probelink->r12();

        int lx; int ly; int lz;
        int crossxtype = probelink->crossxtype();
        int crossytype = probelink->crossytype();
        int crossztype = probelink->crossztype();
        
        // make copies of the links

        for(int ix = 0; ix<repeatX; ix++){
            for(int iy = 0; iy<repeatY; iy++){
                for(int iz = 0; iz<repeatZ; iz++){
                    if(!((ix==0)&&((iy==0)&&(iz==0)))) {
                        
                        // id is in vector format not important
                        long new_link_id = link_id + (iz+iy*repeatZ+ix*repeatY*repeatZ)*number_of_links;  
                         
                        // shift indices to allow for correct crossing over boundaries
                        if(crossxtype == (int) NoxCross) { lx = ix;}
                        if(crossxtype == (int) PosxCross) { lx = ix+1; if(lx == repeatX) {lx = 0;}}
                        if(crossxtype == (int) NegxCross) { lx = ix-1; if(lx == -1) {lx = repeatX-1;}}

                        if(crossytype == (int) NoyCross) { ly = iy;}
                        if(crossytype == (int) PosyCross) { ly = iy+1; if(ly == repeatY) {ly = 0;}}
                        if(crossytype == (int) NegyCross) { ly = iy-1; if(ly == -1) {ly = repeatY-1;}}

                        if(crossztype == (int) NozCross) { lz = iz;}
                        if(crossztype == (int) PoszCross) { lz = iz+1; if(lz == repeatZ) {lz = 0;}}
                        if(crossztype == (int) NegzCross) { lz = iz-1; if(lz == -1) {lz = repeatZ-1;}}                        
                        
                        // obtain the new node pointers via the mapped id's
                        int new_node1_id = node1_id + (iz+iy*repeatZ+ix*repeatY*repeatZ)*number_of_nodes;
                        int new_node2_id = node2_id + (lz+ly*repeatZ+lx*repeatY*repeatZ)*number_of_nodes;

                        NodeDevice* new_node1 = this->GetNode(new_node1_id);
                        NodeDevice* new_node2 = this->GetNode(new_node2_id);
                        
                        // r12 is merely translated, not changed
                        
                        LinkDevice* newLinkDevice = this->AddLink(new_link_id,new_node1,new_node2,link_r12);
                        
                        // copy data to the periodically repeated nodes
                        newLinkDevice->setRate(  probelink->rate12e(),probelink->rate12h(),probelink->rate21e(),probelink->rate21h());
                        newLinkDevice->setJeff2( probelink->Jeff2e(), probelink->Jeff2h());
                        newLinkDevice->setlO(    probelink->lOe(),    probelink->lOh());                        
                    }
                }
            }
        }
        
        // correctly reconnect the original links
        if(crossxtype == (int) NoxCross)  { lx = 0;         }
        if(crossxtype == (int) PosxCross) { lx = 1;         }
        if(crossxtype == (int) NegxCross) { lx = repeatX-1; } 

        if(crossytype == (int) NoyCross)  { ly = 0;         }
        if(crossytype == (int) PosyCross) { ly = 1;         }
        if(crossytype == (int) NegyCross) { ly = repeatY-1; }

        if(crossztype == (int) NozCross)  { lz = 0;         }
        if(crossztype == (int) PoszCross) { lz = 1;         }
        if(crossztype == (int) NegzCross) { lz = repeatY-1; }
        
        int new_node1_id = node1_id;
        int new_node2_id = node2_id + (lz+ly*repeatZ+lx*repeatY*repeatZ)*number_of_nodes;
        NodeDevice* new_node1 = this->GetNode(new_node1_id);
        NodeDevice* new_node2 = this->GetNode(new_node2_id);

        probelink->SetNodes(new_node1, new_node2);        
//        std::cout << probelink->Jeff2e() << " " << probelink->Jeff2h() << " " << probelink->lOe() << " " << probelink->lOh() << endl;                        
    } 
}

void GraphDevice::Set_Self_Image_Coulomb_Potential(double device_length, Eventinfo* eventinfo){
    
    typename std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        votca::tools::vec node_pos = (*it)->position();
        (*it)->Compute_Self_Image_Coulomb_Potential(node_pos.x(),device_length,eventinfo);
    }
    
    _left_electrode->setSelfImage(0.0);
    _right_electrode->setSelfImage(0.0);
}

void GraphDevice::Set_Layer_indices(Eventinfo* eventinfo){
    typename std::vector<NodeDevice*>::iterator it;
    for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        votca::tools::vec pos = (*it)->position();
        double posx = pos.x();
        int iposx = floor(posx/eventinfo->layersize);
        (*it)->setLayer(iposx);
    }
}

void GraphDevice::RenumberId() {

    typename std::vector<NodeDevice*>::iterator it;
    for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        
        int renum_ID = 0;
        for (int ilink = 0 ; ilink < (*it)->links().size(); ilink++) {
            (*it)->links()[ilink]->SetID(renum_ID); renum_ID++;
        }        
    }
}

void GraphDevice::Roll_morph(double translate_x, double translate_y, double translate_z) {
 
    this->Determine_cross_types();

    _sim_box_size = this->Determine_Sim_Box_Size();            

    this->Push_in_box();
    this->Determine_cross_types();
    
    this->Translate_graph(translate_x, translate_y, translate_z);
    this->Push_in_box();
    this->Determine_cross_types();
}

/*void GraphDevice::Rotate_morph(double rotate_in_rads, int rotate_axis) {
    
    this->Determine_cross_types();
    
    _sim_box_size = this->Determine_Sim_Box_Size();
    
    this->Push_in_box();
    this->Determine_cross_types();

    // first resize the morphology, rotate and cut off
    
}*/
    
}}

#endif

