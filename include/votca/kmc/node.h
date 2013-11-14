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

#ifndef _VOTCA_KMC_NODE_H
#define	_VOTCA_KMC_NODE_H

#include <votca/tools/vec.h>
#include <votca/kmc/linksql.h>

namespace votca { namespace kmc {
class Node
{
public:
    Node( int id, votca::tools::vec &position) {
        _id  = id;
        _position = position;
    } 

    /// adds a link to a Node
    void AddLink( LinkSQL* link ) { _links.push_back(link); }
    /// link ID - syncing with the pair ID 
    const int &id() const { return _id; } 
    /// type
    const int &type() const { return _type; } 
    /// position (nm))
    const votca::tools::vec &position() const { return _position; } 
    /// print Node info
    virtual void Print( std::ostream &out ) {
        vector< LinkSQL* >::iterator it;
        for (it = _links.begin(); it != _links.end(); it++ ) (*it)->Print( out );
    }

protected:

    int _id;
    int _type;
    votca::tools::vec _position;
    
    vector< LinkSQL* > _links;
};


}}

#endif	/* _VOTCA_KMC_NODE_H */

