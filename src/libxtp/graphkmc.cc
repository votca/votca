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

#include <votca/xtp/graphkmc.h>


using namespace std;

namespace votca {
    namespace xtp {
void GraphKMC::Setup_bulk_graph( Eventinfo* eventinfo)
{

    // Determine hopping distance
    _hop_distance = this->Determine_Hopping_Distance();
 
    // Detemine simulation box size
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
    if(eventinfo->resize_morphology) {
        this->Resize(eventinfo->size_x, false, eventinfo->size_y, false, eventinfo->size_z, false); 
        this->Determine_cross_types(); 
        _sim_box_size = this->Determine_Sim_Box_Size();
    }
    
    //set node types for existing nodes as Normal
    this->Initialize_node_types();

    // determine minimum hopping distance    
    _min_distance = this->Determine_Minimum_Distance();

    // determine sum of r12.x    
    _total_link_distance_z = this->Sum_of_link_distances_z();
    
    // associate links in links vector with the corresponding nodes
    this->LinkSort();
    
    // determine maximum degree of graph
    _max_pair_degree = this->Determine_Max_Pair_Degree();    
    
    // Set the layer indices on every node
    this->Set_Layer_indices(eventinfo);
  
    // clear the occupation of the graph
    this->Clear();

    //renumber link id's
    this->Renumber_id();   
}    

void GraphKMC::Setup_device_graph(Eventinfo* eventinfo)
{
    // Determine hopping distance before breaking periodicity
    _hop_distance = this->Determine_Hopping_Distance();

    // Detemine simulation box size
    _sim_box_size = this->Determine_Sim_Box_Size();    
    
 
    
    
    // Make sure injection hops are possible
    if(eventinfo->left_electrode_distance > _hop_distance) {_hop_distance = eventinfo->left_electrode_distance;}
    if(eventinfo->right_electrode_distance > _hop_distance) {_hop_distance = eventinfo->right_electrode_distance;}
    
    // Determine simulation box size
    _sim_box_size = this->Determine_Sim_Box_Size();
    std::cout << _sim_box_size << endl;
    // Translate the graph due to the spatial location of the electrodes and update system box size accordingly, putting the left electrode at x = 0
    // left_electrode_distance is the distance of the left electrode to the node with minimum x-coordinate
    // distance by which the graph should be translated is left_electrode_distance - minX
    // Translate graph so that minimum coordinates are at 0
    this->Put_at_zero_graph();

    // Push nodes back in box (crossing types are changed by this operation, so re-evaluate)
    this->Push_in_box();
    this->Determine_cross_types();

    _sim_box_size = this->Determine_Sim_Box_Size();
    
    // Resize by copying the box (crossing types are changed by this operation, so re-evaluate, which is done during determination of simulation box size)
    if(eventinfo->resize_morphology) 
    {
        this->Resize(eventinfo->size_x,false, eventinfo->size_y,false, eventinfo->size_z-eventinfo->left_electrode_distance - eventinfo->right_electrode_distance, true); 
    }
    else
    {
        this->Break_periodicity(false,false,true);
    }

    // Recalculate simulation box size after periodicity breaking
    _sim_box_size = this->Determine_Sim_Box_Size();
    
    // Translate graph to accomodate for device geometry
    this->Translate_graph(0.0, 0.0, eventinfo->left_electrode_distance);    

    // adjust simulation box size accordingly to given electrode distances
    votca::tools::vec old_sim_box_size = _sim_box_size;
    double new_sim_box_sizeZ = old_sim_box_size.z() + eventinfo->left_electrode_distance + eventinfo->right_electrode_distance;
     _sim_box_size =  votca::tools::vec(old_sim_box_size.x(), old_sim_box_size.y(), new_sim_box_sizeZ);
   
    // set node types for existing nodes as Normal
    this->Initialize_node_types();
 
    // add electrode nodes and connect those nodes via links to existing nodes
    this->Add_electrodes();

    // determine minimum hopping distance
    _min_distance = this->Determine_Minimum_Distance();
   
    // determine sum of r12.x
    _total_link_distance_z = this->Sum_of_link_distances_z();

    // associate links in links vector with the corresponding nodes
    this->LinkSort();
  
    // determine maximum degree of graph
    _max_pair_degree = this->Determine_Max_Pair_Degree();
 
    // Set the self-image coulomb potential on every node
    this->Set_Self_Image_Coulomb_Potential(_sim_box_size.z(), eventinfo);

    // Set the layer indices on every node
    this->Set_Layer_indices(eventinfo);
    
    // clear the occupation of the graph
    this->Clear();

    // renumber link id's
    this->Renumber_id(); 
   
}

double GraphKMC::Determine_Hopping_Distance()
{
    double hop_distance = 0.0;
    
    std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
        votca::tools::vec dR = (*it)->r12();
        double distance = abs(dR);
        if(distance>hop_distance) hop_distance = distance;
    }
    
    return hop_distance;
}

double GraphKMC::Determine_Minimum_Distance()
{
    double min_distance = abs(_sim_box_size);

    std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
        votca::tools::vec dR = (*it)->r12();
        double distance = abs(dR);
        if(distance<min_distance) min_distance = distance;
    }
    
    return min_distance;
}

unsigned GraphKMC::Determine_Max_Pair_Degree()
{
    unsigned max_pair_degree = 0;

    std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        if((*it)->type() == NormalNode) 
        {
            if((*it)->links().size()>max_pair_degree) max_pair_degree = (*it)->links().size();
        }
    }
    
    return max_pair_degree;    
}

votca::tools::vec GraphKMC::Determine_Sim_Box_Size()
{ 
    votca::tools::vec sim_box_size;

    this->Determine_cross_types();
    
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
    
    std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
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
    if(!pairXfound) {simX = maxX-minX; std::cout << "non-periodic in x direction" << endl;}
    if(!pairYfound) {simY = maxY-minY;}
    if(!pairZfound) {simZ = maxZ-minZ;}

    sim_box_size = votca::tools::vec(simX,simY,simZ);

    return sim_box_size;
}

void GraphKMC::Determine_cross_types()
{
    std::vector<LinkDevice*>::iterator it;
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
        votca::tools::vec pos1 = (*it)->node1()->position();
        votca::tools::vec pos2 = (*it)->node2()->position();
        votca::tools::vec dr = (*it)->r12();

        int crossx_flag = (int) NoxCross;
        int crossy_flag = (int) NoyCross;
        int crossz_flag = (int) NozCross;

        if(pos2.x()-pos1.x() < 0.0  && dr.x() > 0.0 ) { crossx_flag = (int) PosxCross; }
        if(pos2.x()-pos1.x() > 0.0  && dr.x() < 0.0 ) { crossx_flag = (int) NegxCross; }  
        if(pos2.y()-pos1.y() < 0.0  && dr.y() > 0.0 ) { crossy_flag = (int) PosyCross; }  
        if(pos2.y()-pos1.y() > 0.0  && dr.y() < 0.0 ) { crossy_flag = (int) NegyCross; }
        if(pos2.z()-pos1.z() < 0.0  && dr.z() > 0.0 ) { crossz_flag = (int) PoszCross; }  
        if(pos2.z()-pos1.z() > 0.0  && dr.z() < 0.0 ) { crossz_flag = (int) NegzCross; }
        
        (*it)->setCrossxType(crossx_flag);
        (*it)->setCrossyType(crossy_flag);
        (*it)->setCrosszType(crossz_flag);
    }
}

double GraphKMC::Sum_of_link_distances_z()
{
    double totdistz = 0.0;
    
    std::vector<LinkDevice*>::iterator it;
    for (it = this->_links.begin(); it != this->_links.end(); it++ )
    {
        votca::tools::vec distvec = (*it)->r12();
        if(distvec.z() > 0.0) { totdistz += distvec.z();}
    }
    
    return totdistz;
}

void GraphKMC::Translate_graph(double translate_x, double translate_y, double translate_z) 
{
    std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        votca::tools::vec oldpos = (*it)->position();
        double newxpos = oldpos.x() + translate_x;
        double newypos = oldpos.y() + translate_y; 
        double newzpos = oldpos.z() + translate_z;
        votca::tools::vec newpos = votca::tools::vec(newxpos,newypos,newzpos);
        (*it)->SetPosition(newpos);
    }
}

void GraphKMC::Put_at_zero_graph()
{
    votca::tools::vec pos = this->_nodes[0]->position(); // initial value to search for minimum x-, y- and z- coordinate 
    double minX = pos.x();
    double minY = pos.y();
    double minZ = pos.z();

    std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        pos = (*it)->position(); 
        if(pos.x() < minX) {minX = pos.x();}
        if(pos.y() < minY) {minY = pos.y();}
        if(pos.z() < minZ) {minZ = pos.z();}
    }

    double xtranslate = -1.0*minX;
    double ytranslate = -1.0*minY;
    double ztranslate = -1.0*minZ;

    std::cout << "transl " << xtranslate << " " << ytranslate << " " << ztranslate << endl; 
    this->Translate_graph(xtranslate,ytranslate,ztranslate);
}

void GraphKMC::Push_in_box()
{  
    //checks which nodes are falling outside the simboxsize
    //periodically determine where they are in the simulation box

    std::vector<NodeDevice*>::iterator it;      
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        votca::tools::vec pos = (*it)->position();
        
        double newxpos = pos.x(); double newypos = pos.y(); double newzpos = pos.z();
        
        while(newxpos>_sim_box_size.x()) { newxpos -= _sim_box_size.x();}
        while(newxpos<0.0)               { newxpos += _sim_box_size.x();}
        while(newypos>_sim_box_size.y()) { newypos -= _sim_box_size.y();}
        while(newypos<0.0)               { newypos += _sim_box_size.y();}
        while(newzpos>_sim_box_size.z()) { newzpos -= _sim_box_size.z();}        
        while(newzpos<0.0)               { newzpos += _sim_box_size.z();}
        
        votca::tools::vec newpos = votca::tools::vec(newxpos,newypos,newzpos);
        (*it)->SetPosition(newpos);
    }
}

void GraphKMC::Resize(double dimX, bool breakX, double dimY, bool breakY, double dimZ, bool breakZ) 
{
    int repeatX = ceil(dimX/_sim_box_size.x());
    int repeatY = ceil(dimY/_sim_box_size.y());
    int repeatZ = ceil(dimZ/_sim_box_size.z());
    
    std::cout << dimX/_sim_box_size.x() << " " << repeatX << endl;
    std::cout << dimY/_sim_box_size.y() << " " << repeatY << endl;
    std::cout << dimZ/_sim_box_size.z() << " " << repeatZ << endl;
    
    int number_of_nodes = this->Numberofnodes();

    // repeat nodes

    for(int ix = 0; ix<repeatX; ix++)
    {
        for(int iy = 0; iy<repeatY; iy++)
        {
            for(int iz = 0; iz<repeatZ; iz++)
            {
                if(!((ix==0)&&((iy==0)&&(iz==0)))) 
                {
                    for(int it = 0; it < number_of_nodes; it++) 
                    {
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
    for(long ilink = number_of_links - 1; ilink >= 0; ilink--) 
    {
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

        for(int ix = 0; ix<repeatX; ix++)
        {
            for(int iy = 0; iy<repeatY; iy++)
            {
                for(int iz = 0; iz<repeatZ; iz++)
                {
                    if(!((ix==0)&&((iy==0)&&(iz==0)))) 
                    {
                        // id is in vector format not important (remapped later on)

                        long new_link_id = link_id + (iz+iy*repeatZ+ix*repeatY*repeatZ)*number_of_links;  

                        // check for link breaking conditions, those links will not be formed
                        
                        bool nocopy_link = false;
                        
                        // shift indices to allow for correct crossing over boundaries

                        if(crossxtype == (int) NoxCross) { lx = ix;}
                        if(crossxtype == (int) PosxCross) { lx = ix+1; if(lx == repeatX) {lx = 0;    if(breakX) {nocopy_link = true;}}}
                        if(crossxtype == (int) NegxCross) { lx = ix-1; if(lx == -1) {lx = repeatX-1; if(breakX) {nocopy_link = true;}}}

                        if(crossytype == (int) NoyCross) { ly = iy;}
                        if(crossytype == (int) PosyCross) { ly = iy+1; if(ly == repeatY) {ly = 0;    if(breakY) {nocopy_link = true;}}}
                        if(crossytype == (int) NegyCross) { ly = iy-1; if(ly == -1) {ly = repeatY-1; if(breakY) {nocopy_link = true;}}}

                        if(crossztype == (int) NozCross) { lz = iz;}
                        if(crossztype == (int) PoszCross) { lz = iz+1; if(lz == repeatZ) {lz = 0;    if(breakZ) {nocopy_link = true;}}}
                        if(crossztype == (int) NegzCross) { lz = iz-1; if(lz == -1) {lz = repeatZ-1; if(breakZ) {nocopy_link = true;}}}                        
                        
                        // obtain the new node pointers via the mapped id's
                        
                        int new_node1_id = node1_id + (iz+iy*repeatZ+ix*repeatY*repeatZ)*number_of_nodes;
                        int new_node2_id = node2_id + (lz+ly*repeatZ+lx*repeatY*repeatZ)*number_of_nodes;

                        NodeDevice* new_node1 = this->GetNode(new_node1_id);
                        votca::tools::vec newnode1pos = new_node1->position();
                        
                        NodeDevice* new_node2 = this->GetNode(new_node2_id);
                        votca::tools::vec newnode2pos = new_node2->position();
                       
                        if ((newnode1pos.x() > dimX || newnode2pos.x() > dimX) && breakX) nocopy_link = true;
                        if ((newnode1pos.y() > dimY || newnode2pos.y() > dimY) && breakY) nocopy_link = true;
                        if ((newnode1pos.z() > dimZ || newnode2pos.z() > dimZ) && breakZ) nocopy_link = true;
                        
                        if(!nocopy_link) 
                        {
                            // add link to links vector
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
        }

        // check for link breaking conditions, those links will be removed
        
        bool break_link = false;
                                

        // correctly reconnect the original links (split up)
        
        if(crossxtype == (int) NoxCross)  { lx = 0;                                                                   }
        if(crossxtype == (int) PosxCross) { lx = 1;         if(repeatX == 1) {lx = 0; if(breakX) {break_link = true; }}}
        if(crossxtype == (int) NegxCross) { lx = repeatX-1; if(breakX) {break_link = true;}} 

        if(crossytype == (int) NoyCross)  { ly = 0;                                                                   }
        if(crossytype == (int) PosyCross) { ly = 1;         if(repeatY == 1) {ly = 0; if(breakY) {break_link = true; }}}
        if(crossytype == (int) NegyCross) { ly = repeatY-1; if(breakY) {break_link = true; }}

        if(crossztype == (int) NozCross)  { lz = 0;                                                                   }
        if(crossztype == (int) PoszCross) { lz = 1;         if(repeatZ == 1) {lz = 0; if(breakZ) {break_link = true;}}}
        if(crossztype == (int) NegzCross) { lz = repeatZ-1; if(breakZ) {break_link = true; }}
     
        if(!break_link){
            votca::tools::vec pos1 = node1->position();
            double xpos1 = pos1.x(); double ypos1 = pos1.y(); double zpos1 = pos1.z();       

            votca::tools::vec pos2 = node2->position();
            double xpos2 = pos2.x(); double ypos2 = pos2.y(); double zpos2 = pos2.z();
        
            if(breakX && (xpos1 > dimX || xpos2 > dimX) ) { break_link = true; }
            if(breakY && (ypos1 > dimY || ypos2 > dimY) ) { break_link = true; }
            if(breakZ && (zpos1 > dimZ || zpos2 > dimZ) ) { break_link = true; }
        }

        int new_node1_id = node1_id;       
        int new_node2_id = node2_id + (lz+ly*repeatZ+lx*repeatY*repeatZ)*number_of_nodes;

        Node* new_node1 = this->GetNode(new_node1_id);
        Node* new_node2 = this->GetNode(new_node2_id);

        if(!break_link){
            votca::tools::vec newpos1 = new_node1->position();
            double newxpos1 = newpos1.x(); double newypos1 = newpos1.y(); double newzpos1 = newpos1.z();       

            votca::tools::vec newpos2 = new_node2->position();
            double newxpos2 = newpos2.x(); double newypos2 = newpos2.y(); double newzpos2 = newpos2.z();

            if(breakX && (newxpos1 > dimX || newxpos2 > dimX) ) { break_link = true; }
            if(breakY && (newypos1 > dimY || newypos2 > dimY) ) { break_link = true; }
            if(breakZ && (newzpos1 > dimZ || newzpos2 > dimZ) ) { break_link = true; } 
        }
        
        if(!break_link) 
        {
            // nodes which are paired in link are remapped, not necessary when no copies are into play
            probelink->SetNodes(new_node1, new_node2);
            probelink->setRemove(false);
        }
        else 
        {
            //this->RemoveLink(ilink);  
            probelink->setRemove(true);
        }
    }
    
    std::vector<LinkDevice*> temp_links;
    
    temp_links.clear();
    std::vector<LinkDevice*>::iterator it;    
    for(it = this->_links.begin(); it != this->_links.end(); it++) {
        if(!((*it)->remove())) {
            temp_links.push_back((*it));
        }
    }
    
    this->ClearLinks();
    
    for(it = temp_links.begin(); it != temp_links.end(); it++) {
        this->PushLink((*it));
    }

    temp_links.clear();
            
    // remove nodes which fall outside of the simulation box
    
    for(int it = this->_nodes.size()-1; it >= 0; it--) 
    {
        NodeDevice* inode = this->_nodes[it];
        
        votca::tools::vec pos = inode->position();
        double xpos = pos.x(); double ypos = pos.y(); double zpos = pos.z();
        bool remove_flag = false;        
        
        if(breakX && xpos > dimX) { remove_flag = true; }
        if(breakY && ypos > dimY) { remove_flag = true; }
        if(breakZ && zpos > dimZ) { remove_flag = true; }        
 
        if(remove_flag) 
        {
            this->RemoveNode(it);
            delete inode;
        }
   }  

}

void GraphKMC::Break_periodicity(bool break_x, bool break_y, bool break_z)
{
    // Break crossing links

    for(int it = this->_links.size()-1; it >= 0; it--) 
    {
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

void GraphKMC::Initialize_node_types() 
{
    std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) (*it)->SetType((int) NormalNode);    
}

void GraphKMC::Add_electrodes() 
{
    //introduce the electrode nodes (might make a special electrode node header file for this)

    _left_electrode = new NodeDevice(-1, tools::vec(0.0,0.0,0.0));
    _right_electrode = new NodeDevice(-2, tools::vec(0.0,0.0,_sim_box_size.z()));
    _left_electrode->SetType((int) LeftElectrodeNode);
    _right_electrode->SetType((int) RightElectrodeNode);

    int linkcount = 0;
    double av_left = 0.0;
    double av_right = 0.0;
    
    int linkID = this->Numberoflinks();
    
    //determine the nodes which are injectable from the left electrode and the nodes which are injectable from the right electrode

    std::vector<NodeDevice*>::iterator it;    
    for(it  = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    { 
        votca::tools::vec nodepos = (*it)->position();
        double left_distance = nodepos.z();
        (*it)->Set_injectable(false);

        // check nodes for left electrode
        
        if(left_distance <= _hop_distance) 
        {
            votca::tools::vec dr = votca::tools::vec(0.0,0.0,-1.0*left_distance);   

            //LinkSQL* newLinkCollect = this->AddLink(linkID,(*it), _left_electrode, dr); 
            this->AddLink(linkID,(*it), _left_electrode, dr); 
            linkID++;
            
            LinkSQL* newLinkInject = new LinkSQL(linkID, _left_electrode, (*it), -1.0*dr);
            _left_electrode->AddLink(newLinkInject);

            (*it)->Set_injectable(true);
            
            av_left += fabs(left_distance);
            linkcount++;
        }
      
        // check nodes for right electrode
        
        double right_distance = _sim_box_size.z() - nodepos.z();
        if(right_distance <= _hop_distance) 
        {
            votca::tools::vec dr = votca::tools::vec(0.0,0.0,right_distance);   

            //LinkSQL* newLinkCollect = this->AddLink(linkID,(*it), _right_electrode, dr); 
            this->AddLink(linkID,(*it), _right_electrode, dr); 
            linkID++;
            
            LinkSQL* newLinkInject = new LinkSQL(linkID, _right_electrode, (*it), -1.0*dr);
            _right_electrode->AddLink(newLinkInject);

            (*it)->Set_injectable(true);

            av_right += fabs(right_distance);
            linkcount++;            
        }
    }

    _av_el_distance = (av_left+av_right)/(1.0*linkcount);
    std::cout << "average electrode distance " << _av_el_distance << endl;
    
    this->AddNode(_left_electrode); //in this way the electrode nodes are caught by graph destructor
    this->AddNode(_right_electrode);
    _left_electrode->RemoveCarrier();
    _right_electrode->RemoveCarrier();    
}

void GraphKMC::LinkSort()
{
    std::vector<LinkDevice*>::iterator it;
    for (it = this->_links.begin(); it != this->_links.end(); it++ ) 
    {
        NodeDevice* node1 = dynamic_cast<NodeDevice*>((*it)->node1());
        if(node1->type() == NormalNode) node1->AddLink((*it));
    }
    
}

void GraphKMC::Set_Self_Image_Coulomb_Potential(double device_length, Eventinfo* eventinfo)
{
    std::vector<NodeDevice*>::iterator it;    
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        votca::tools::vec node_pos = (*it)->position();
        (*it)->Compute_Self_Image_Coulomb_Potential(node_pos.z(),device_length,eventinfo);
    }
    
    _left_electrode->setSelfImage(0.0);
    _right_electrode->setSelfImage(0.0);
}

void GraphKMC::Set_Layer_indices(Eventinfo* eventinfo)
{
    std::vector<NodeDevice*>::iterator it;
    for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        votca::tools::vec pos = (*it)->position();
        double posz = pos.z();
        
        double layersize = (_sim_box_size.z()-eventinfo->left_electrode_distance - eventinfo->right_electrode_distance)/eventinfo->number_of_layers;
        int iposz = floor((posz-eventinfo->left_electrode_distance)/layersize);
        if(iposz == eventinfo->number_of_layers) iposz = eventinfo->number_of_layers - 1;
        
        (*it)->setLayer(iposz);
    }
}

void GraphKMC::Renumber_id() 
{
    std::vector<NodeDevice*>::iterator it;
    for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        int renum_ID = 0;
        for (unsigned ilink = 0 ; ilink < (*it)->links().size(); ilink++) 
        {
            (*it)->links()[ilink]->SetID(renum_ID); renum_ID++;
        }        
    }
    
 /*   for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        int node_id = (*it)->id();
        for (int ilink = 0 ; ilink < (*it)->links().size(); ilink++) 
        {
            Node* probe_node = (*it)->links()[ilink]->node2();
            for (int ireverse = 0 ; ireverse < probe_node->links().size(); ireverse++)
            {
                Link* probe_link = probe_node->links()[ireverse];
                if(probe_link->node2()->id()==node_id) (*it)->links()[ilink]->SetReverseID(ireverse);
            }
        }        
    }    */
    
}


void GraphKMC::Initialize_output_values() 
{
    std::vector<NodeDevice*>::iterator it;
    for (it = this->_nodes.begin(); it != this->_nodes.end(); it++) {
        (*it)->Initialize_output_values();    
    }
}

void GraphKMC::Roll_morph(double translate_x, double translate_y, double translate_z) 
{
    this->Determine_cross_types();

    _sim_box_size = this->Determine_Sim_Box_Size();            

    this->Push_in_box();
    this->Determine_cross_types();
    
    this->Translate_graph(translate_x, translate_y, translate_z);
    this->Push_in_box();
    this->Determine_cross_types();
}

double GraphKMC::Average_hole_node_energy() 
{
    double ho_energy = 0.0;
    double av_ho_energy = 0.0;
    
    std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        if((*it)->type() == (int) NormalNode){
            ho_energy += (*it)->eCation() + (*it)->UcCnNh();
        }
    }
    
    av_ho_energy = ho_energy/this->Numberofnodes();
    return av_ho_energy;
}

double GraphKMC::Average_hole_left_electrode_energy() 
{
    double ho_energy = 0.0;
    double av_ho_energy = 0.0;
    int linker = 0;
    
    for(unsigned ilink = 0; ilink < _left_electrode->links().size(); ilink++)    {
        Node* probe_node = _left_electrode->links()[ilink]->node2();
        ho_energy +=  dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh();
        linker++;
    }
    
    av_ho_energy = ho_energy/linker;
    return av_ho_energy;
}

double GraphKMC::Average_hole_right_electrode_energy() 
{
    double ho_energy = 0.0;
    double av_ho_energy = 0.0;
    int linker = 0;
    
    for(unsigned ilink = 0; ilink < _right_electrode->links().size(); ilink++)    {
        Node* probe_node = _right_electrode->links()[ilink]->node2();
        ho_energy +=  dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh();
        linker++;
    }
    
    av_ho_energy = ho_energy/linker;
    return av_ho_energy;
}

double GraphKMC::stddev_hole_node_energy() 
{
    double temp_energy = 0.0;
    double stddev_energy = 0.0;
    double av_ho_energy = this->Average_hole_node_energy();
    
    std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        if((*it)->type() == (int) NormalNode){
            temp_energy += ((*it)->eCation() + (*it)->UcCnNh()-av_ho_energy)*((*it)->eCation() + (*it)->UcCnNh()-av_ho_energy);
        }
    }
    
    stddev_energy = sqrt(temp_energy/this->Numberofnodes());
    return stddev_energy;
}

double GraphKMC::stddev_hole_left_electrode_energy() 
{
    double temp_energy = 0.0;
    double stddev_energy = 0.0;
    double av_ho_energy = this->Average_hole_left_electrode_energy();
    int linker = 0;
    
    for(unsigned ilink = 0; ilink < _left_electrode->links().size(); ilink++)    {
        Node* probe_node = _left_electrode->links()[ilink]->node2();
        temp_energy +=  (dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh() - av_ho_energy)*(dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh() - av_ho_energy);
        linker++;
    }
    
    stddev_energy = sqrt(temp_energy/linker);
    return stddev_energy;
}

double GraphKMC::stddev_hole_right_electrode_energy() 
{
    double temp_energy = 0.0;
    double stddev_energy = 0.0;
    double av_ho_energy = this->Average_hole_right_electrode_energy();
    int linker = 0;
    
    for(unsigned ilink = 0; ilink < _right_electrode->links().size(); ilink++)    {
        Node* probe_node = _right_electrode->links()[ilink]->node2();
        temp_energy +=  (dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh() - av_ho_energy)*(dynamic_cast<NodeDevice*>(probe_node)->eCation() +  dynamic_cast<NodeDevice*>(probe_node)->UcCnNh() - av_ho_energy);
        linker++;
    }
    
    stddev_energy = sqrt(temp_energy/linker);
    return stddev_energy;
}

double GraphKMC::Average_electron_node_energy() 
{
    
    double elect_energy = 0.0;
    double av_elect_energy = 0.0;
    
    std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        if((*it)->type() == (int) NormalNode){
            elect_energy += (*it)->eAnion() + (*it)->UcCnNe();
        }
    }
    
    av_elect_energy = elect_energy/this->Numberofnodes();
    return av_elect_energy;
}

double GraphKMC::stddev_electron_node_energy() 
{
    double temp_energy = 0.0;
    double stddev_energy = 0.0;
    double av_el_energy = this->Average_electron_node_energy();
    
    std::vector<NodeDevice*>:: iterator it;
    for(it = this->_nodes.begin(); it != this->_nodes.end(); it++) 
    {
        if((*it)->type() == (int) NormalNode){
            temp_energy += ((*it)->eAnion() + (*it)->UcCnNe()-av_el_energy)*((*it)->eAnion() + (*it)->UcCnNe()-av_el_energy);
        }
    }
    
    stddev_energy = sqrt(temp_energy/this->Numberofnodes());
    return stddev_energy;
}

double GraphKMC::Electron_inject_reorg()
{
    double electron_reorg = 0.0;
    double temp_reorg;
    int link_count = 0;
    
    std::vector<LinkDevice*>:: iterator it;
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
        Node* node1 = (*it)->node1();
        Node* node2 = (*it)->node2();
        if(node1->type() == (int) NormalNode && node2->type() == (int) NormalNode) {
            temp_reorg = dynamic_cast<NodeSQL*>((*it)->node1())->UnCnNe() + dynamic_cast<NodeSQL*>((*it)->node2())->UcNcCe() + (*it)->lOe();
            electron_reorg += temp_reorg;
            link_count++;
        }
    }

    return (electron_reorg/link_count);    
}

double GraphKMC::Hole_inject_reorg()
{
    double hole_reorg = 0.0;
    double temp_reorg;
    int link_count = 0;
    
    std::vector<LinkDevice*>:: iterator it;
    for(it = this->_links.begin(); it != this->_links.end(); it++) 
    {
        Node* node1 = (*it)->node1();
        Node* node2 = (*it)->node2();

        if(node1->type() == (int) NormalNode && node2->type() == (int) NormalNode) {
            temp_reorg = dynamic_cast<NodeSQL*>((*it)->node1())->UnCnNh() + dynamic_cast<NodeSQL*>((*it)->node2())->UcNcCh() + (*it)->lOh();
            hole_reorg += temp_reorg;
            link_count++;
        }
    }

    return (hole_reorg/link_count);    
}



        
    }
}
