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

#ifndef _VOTCA_KMC_NODE_H
#define	_VOTCA_KMC_NODE_H

#include <votca/tools/vec.h>
#include <votca/xtp/link.h>

namespace votca { namespace xtp {
class Node
{
public:
    Node( int id, votca::tools::vec &position) {
        _id  = id;
        _position = position;
        _occupation = -1; //empty
    }
    
    virtual ~Node(){
        std::vector<Link*>::iterator it;
        for (it = _links.begin(); it != _links.end(); it++ ) delete *it;        
    }

    /// adds a link to a Node
    virtual void AddLink( Link* link ) { _links.push_back(link); }
    /// link ID - syncing with the pair ID 
    const int &id() const { return _id; } 
    /// type
    const int &type() const { return _type; } 
    /// position (nm))
    const votca::tools::vec &position() const { return _position; } 
    /// print Node info
    
    /// links
    const vector<Link* > &links() const {return _links;}
    
    /// carriers
    const int &occ() const {return _occupation;}
    
    virtual void Print( std::ostream &out ) {
        vector< Link* >::iterator it;
        for (it = _links.begin(); it != _links.end(); it++ ) (*it)->Print( out );
    }
    
    /// Set and remove carrier occupation (carrier index) of node (-1 being empty)
    void AddCarrier(int carrier_ID) {_occupation = carrier_ID;}
    void RemoveCarrier() {_occupation = -1;}
    
    /// Set node type
    void SetType(int type) { _type = type;}
    
    /// Set node position
    void SetPosition(votca::tools::vec position) { _position = position;}

    const int &vel_x() const {return _vel_x;} 
    const int &vel_y() const {return _vel_y;} 
    const int &vel_z() const {return _vel_z;}     
    
    void Add_velx(int dis) { _vel_x += dis;}
    void Add_vely(int dis) { _vel_y += dis;}
    void Add_velz(int dis) { _vel_z += dis;}
    
    void Initialize_output_values() {_vel_x = 0; _vel_y = 0; _vel_z = 0;}

       
    
protected:
    
    int _id;
    int _type;
    votca::tools::vec _position;
    int _occupation;

//    double _vel_x;
    int _vel_x;
    int _vel_y;
    int _vel_z;
    
    vector< Link* > _links;
    
};


}}

#endif	/* _VOTCA_KMC_NODE_H */

