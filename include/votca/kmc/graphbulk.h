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

#ifndef __VOTCA_KMC_GRAPHBULK_H_
#define __VOTCA_KMC_GRAPHBULK_H_

#include <vector>
#include <votca/kmc/graphsql.h>
//#include <votca/kmc/nodedevice.h>
//#include <votca/kmc/graphcubic.h>

namespace votca { namespace kmc {

enum NodeType {NormalNode, LeftElectrodeNode, RightElectrodeNode};    
    
template<class TGraph, class TNode, class TLink>    
class GraphBulk : public TGraph {

public:

    ///determine maxpairdegree, hopping distance and simulation box size
    void Determine_Graph_Parameters(){    
    _max_pair_degree = Determine_Max_Pair_Degree();
    _hop_distance = Determine_Hopping_Distance();
    _sim_box_size = Determine_Sim_Box_Size();

    };

    ///add all links to the nodes
    void LinkSort();
    
    ///calculate the maximum of all degrees in the graph
    int Determine_Max_Pair_Degree();
    
    ///calculate hopping_distance (maximum distance between a pair of nodes) ... needed for injection and coulomb potential calculations
    double Determine_Hopping_Distance();
    
    ///calculate the simulation box size
    votca::tools::vec Determine_Sim_Box_Size();

    ///break the periodicity of the graph (breaking boundary crossing pairs) .. (run before linksort)
    void Break_periodicity(bool break_x, bool break_y, bool break_z);

    ///define electrode nodes and form links between those nodes and neighbouring nodes
    void Setup_device_graph(double left_distance, double right_distance);   
    
private:

    void RenumberId();        
    
    int _max_pair_degree;  
    double _hop_distance;
    votca::tools::vec _sim_box_size;
    
    TNode* left_electrode;
    TNode* right_electrode;
     
};


template<class TGraph, class TNode, class TLink>
void GraphBulk<TGraph, TNode, TLink>::LinkSort(){
    typename std::vector<TLink*>::iterator it;
    for (it = this->_links.begin(); it != this->_links.end(); it++ ) {
        Node* node1 = (*it)->node1();
        votca::tools::vec pos = node1->position();
//        std::cout << pos.x() << " " << pos.y() << " " << pos.z() << endl;
        node1->AddLink((*it));
    }
}

template<class TGraph, class TNode, class TLink>
int GraphBulk<TGraph, TNode, TLink>::Determine_Max_Pair_Degree(){
    
    int max_pair_degree = 0;
    typename std::vector<TNode*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        if((*it)->links().size()>max_pair_degree) max_pair_degree = (*it)->links().size();
    }
    return max_pair_degree; 
}

template<class TGraph, class TNode, class TLink>
double GraphBulk<TGraph, TNode, TLink>::Determine_Hopping_Distance(){
    
    double hop_distance = 0.0;
    typename std::vector<TLink*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) {
        votca::tools::vec dR = (*it)->r12();
        double distance = abs(dR);
        if(distance>hop_distance) hop_distance = distance;
    }
    return hop_distance;
}


template<class TGraph, class TNode, class TLink>
votca::tools::vec GraphBulk<TGraph, TNode, TLink>::Determine_Sim_Box_Size(){ 
    
    //Determination of simulation box size
    //To do this, we first need to find a node with position vector a and pairing node with position vector b, such that
    //a+r12>b or a+r12<b
    //Note that it is possible that none of the pairs pass the simulation box boundaries
    //In this special case, we must determine the node with max x/y/z coordinate and min x/y/z coordinate
    
    //is a boundary crossing pair being found?
    bool pairXfound = false; bool pairYfound = false; bool pairZfound = false;
    
    double sim_box_sizeX; double sim_box_sizeY; double sim_box_sizeZ;
    
    //initial values
    votca::tools::vec initpos = this->_nodes[0]->position();
    double maxX = initpos.x(); double maxY = initpos.y(); double maxZ = initpos.z();
    double minX = initpos.x(); double minY = initpos.y(); double minZ = initpos.z();
    
    typename std::vector<TLink*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) {
        
        if(pairXfound&&pairYfound&&pairZfound) break;
        
        votca::tools::vec pos1 = (*it)->node1()->position();
        votca::tools::vec pos2 = (*it)->node2()->position();
        votca::tools::vec dr = (*it)->r12();

        if(maxX<pos1.x()) {maxX = pos1.x();}
        if(minX>pos1.x()) {minX = pos1.x();}
        if(maxY<pos1.y()) {maxY = pos1.y();}
        if(minY>pos1.y()) {minY = pos1.y();}
        if(maxZ<pos1.z()) {maxZ = pos1.z();}
        if(minZ>pos1.z()) {minZ = pos1.z();}
        
        if(pos1.x() > pos2.x() && dr.x()>0) {   pairXfound = true;  sim_box_sizeX = pos1.x() + dr.x() - pos2.x();}
        if(pos1.x() < pos2.x() && dr.x()<0) {   pairXfound = true;  sim_box_sizeX = pos2.x() - dr.x() - pos1.x();}
        if(pos1.y() > pos2.y() && dr.y()>0) {   pairYfound = true;  sim_box_sizeY = pos1.y() + dr.y() - pos2.y();}
        if(pos1.y() < pos2.y() && dr.y()<0) {   pairYfound = true;  sim_box_sizeY = pos2.y() - dr.y() - pos1.y();}
        if(pos1.z() > pos2.z() && dr.z()>0) {   pairZfound = true;  sim_box_sizeZ = pos1.z() + dr.z() - pos2.z();}
        if(pos1.z() < pos2.z() && dr.z()<0) {   pairZfound = true;  sim_box_sizeZ = pos2.z() - dr.z() - pos1.z();}
        
        //possibility that hopping distance is larger than the simulation box size
                                        
    }
    
    //for the possible outcome that none of the pairs have crossed the simulation box boundary
    if(!pairXfound) {sim_box_sizeX = maxX-minX;}
    if(!pairYfound) {sim_box_sizeY = maxY-minY;}
    if(!pairZfound) {sim_box_sizeZ = maxZ-minZ;}

    votca::tools::vec sim_box_size(sim_box_sizeX, sim_box_sizeY, sim_box_sizeZ);
    return sim_box_size;

}

template<class TGraph, class TNode, class TLink>
void GraphBulk<TGraph, TNode, TLink>::Break_periodicity(bool break_x, bool break_y, bool break_z){

//    typename std::vector<TLink*>::iterator it;
    int debugcount = 0;
    std::cout <<  this->_links.size() << endl;
    for(int it = this->_links.size()-1; it != 0; it--) {

        TLink* ilink = this->_links[it];
        
        votca::tools::vec pos1 = ilink->node1()->position();
        votca::tools::vec pos2 = ilink->node2()->position();
        votca::tools::vec dr = ilink->r12();

        bool remove_flag = false;
        
        if(break_x){
            if(pos1.x() > pos2.x() && dr.x()>0) remove_flag = true;
            if(pos1.x() < pos2.x() && dr.x()<0) remove_flag = true;
        }        
        if(break_y){
            if(pos1.y() > pos2.y() && dr.y()>0) remove_flag = true;
            if(pos1.y() < pos2.y() && dr.y()<0) remove_flag = true;
        } 
        if(break_z){
            if(pos1.z() > pos2.z() && dr.z()>0) remove_flag = true;
            if(pos1.z() < pos2.z() && dr.z()<0) remove_flag = true;
        }
        
        if(remove_flag) {
            std::cout << pos1 << " " << pos2 << "\n";
            //this->_links.erase(this->_links.begin()+it);/* std::cout << debugcount << endl;debugcount++;*/
            this->RemoveLink(it);
        }
    }
    std::cout <<  this->_links.size() << endl;
    this->RenumberId();
}

template<class TGraph, class TNode, class TLink>
void GraphBulk<TGraph, TNode, TLink>::Setup_device_graph(double left_distance, double right_distance){

    // Translate the graph due to the spatial location of the electrodes and update system box size accordingly, putting the left electrode at x = 0
    // left_electrode_distance is the distance of the left electrode to the node with minimum x-coordinate
    
 /* 
    votca::tools::vec pos = this->_nodes[0]->position(); // initial value 
    double minX = pos.x();

    typename std::vector<TNode*>::iterator it;      
     for(int it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        pos = (*it)->position(); if(pos.x() < minX) {minX = pos.x();}
    }
    std::cout << minX << "\n";
   //distance by which the graph should be translated is left_electrode_distance - minX

    double xtranslate = left_distance - minX;

    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        votca::tools::vec oldpos = (*it)->position();
        double newxpos = oldpos.x() + xtranslate;
        votca::tools::vec newpos = votca::tools::vec(newxpos,oldpos.y(),oldpos.z());
        (*it)->SetPosition(newpos);
    }

    //adjust system box size accordingly

    votca::tools::vec old_sim_box_size = _sim_box_size;
    double new_sim_box_sizeX = old_sim_box_size.x() + left_distance + right_distance;
    _sim_box_size =  votca::tools::vec(new_sim_box_sizeX, old_sim_box_size.y(), old_sim_box_size.z());
            
    //introduce electrode nodes
    votca::tools::vec leftpos = votca::tools::vec(0.0,0.0,0.0);
    votca::tools::vec rightpos = votca::tools::vec(_sim_box_size.x(),0.0,0.0);
     left_electrode = this->AddNode(-1, leftpos);
    right_electrode = this->AddNode(-2, rightpos);

    std::cout << (int) LeftElectrode << " " << (int) RightElectrode << "\n";
    

    left_electrode->SetType((int) LeftElectrode);
    right_electrode->SetType((int) RightElectrode);
//    left_electrode->SetType(1);
//    right_electrode->SetType(2);
    std::cout << left_electrode->type() << " " << right_electrode->type() << "\n";    
    
   //determine the nodes which are injectable from the left electrode and the nodes which are injectable from the right electrode

    int linjector_ID = 0;
    int rinjector_ID = 0;

    for(int inode=0; inode<nodes.size(); inode++) { 
      
        double left_distance = nodes[inode]->node_position.x();
     
        if(left_distance <= hopdist) {

            myvec dr = myvec(-1.0*left_distance,0.0,0.0);            
            nodes[inode]->setPair(left_electrode);
            nodes[inode]->setStaticeventinfo(left_electrode, dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            left_electrode->setPair(nodes[inode]);
            left_electrode->setStaticeventinfo(nodes[inode], -1.0*dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            nodes[inode]->left_injector_ID = linjector_ID;
            linjector_ID++;
      
        }
      
        double right_distance = sim_box_size.x() - nodes[inode]->node_position.x();
      
        if(right_distance <= hopdist) {
          
            myvec dr = myvec(right_distance,0.0,0.0);            
            nodes[inode]->setPair(right_electrode);
            nodes[inode]->setStaticeventinfo(right_electrode, dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS
            right_electrode->setPair(nodes[inode]);
            right_electrode->setStaticeventinfo(nodes[inode], -1.0*dr, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0); //NEED TO DETERMINE THIS            
            nodes[inode]->right_injector_ID = rinjector_ID;
            rinjector_ID++;
        }
    }
    nr_left_injector_nodes = linjector_ID;
    nr_right_injector_nodes = rinjector_ID;

    // Recompute the max pair degree, due to the breaking of the periodicity and the connection of nodes to the electrode nodes
    max_pair_degree = Determine_max_pair_degree(nodes);*/

}

template<class TGraph, class TNode, class TLink>
void GraphBulk<TGraph,TNode, TLink>::RenumberId() {
    int resetID = 0;
    typename std::vector<TLink*>::iterator it;
    for (it = this->_links.begin(); it != this->_links.end(); it++) {
        (*it)->SetID(resetID);
        resetID++;
    }
    // are pair IDs necessary?
}
    
}}

/*this->_links.erase(link_id-id_shift)*/

#endif

