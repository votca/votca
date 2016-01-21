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

#ifndef _VOTCA_KMC_LINKSQL_H
#define	_VOTCA_KMC_LINKSQL_H

#include <votca/tools/vec.h>
#include <votca/xtp/link.h>

namespace votca { namespace xtp {

/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separations 
 * 
 */
class LinkSQL: public Link
{

public:

    LinkSQL( long id, Node* node1, Node* node2, votca::tools::vec r12) : Link(id,node1,node2,r12){
    };
    
    void setRate(double rate12e, double rate12h, double rate21e, double rate21h){
        _rate12e  = rate12e;
        _rate12h  = rate12h;
        _rate21e  = rate21e;
        _rate21h  = rate21h;
    }
    
    void setJeff2(double Jeff2e, double Jeff2h) {
        _Jeff2e  = Jeff2e;
        _Jeff2h  = Jeff2h;
    }

    void setlO(double lOe, double lOh) {
        _lOe = lOe;
        _lOh = lOh;
    }

    const double &rate12e() const { return _rate12e; }     
    const double &rate12h() const { return _rate12h; }   
    const double &rate21e() const { return _rate21e; }   
    const double &rate21h() const { return _rate21h; }   

    const double &Jeff2e() const { return _Jeff2e; }       
    const double &Jeff2h() const { return _Jeff2h; }   

    const double &lOe() const { return _lOe; }   
    const double &lOh() const { return _lOh;; }     
    
    
private:

    
    double _rate12e;
    double _rate12h;
    double _rate21e;
    double _rate21h;
    
    double _Jeff2e;
    double _Jeff2h;
    
    double _lOe;
    double _lOh;
    
    
    
};

}} 

#endif // _VOTCA_KMC_LINKSQL_H
