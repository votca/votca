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

#ifndef _VOTCA_KMC_LINKDEVICE_H
#define	_VOTCA_KMC_LINKDEVICE_H

#include <votca/tools/vec.h>
#include <votca/kmc/linksql.h>

namespace votca { namespace kmc {

/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separations 
 * 
 */
class LinkDevice: public LinkSQL
{

public:

    LinkDevice( int id, Node* node1, Node* node2, votca::tools::vec r12) : LinkSQL(id,node1,node2,r12){
    };
    
    void setSelfImage(double selfimage){_self_image = selfimage;}

    const double &self_image() const { return _self_image; }         
    
private:
    double _self_image;
    
};

}} 

#endif // _VOTCA_KMC_LINKDEVICE_H
