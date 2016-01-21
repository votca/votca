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

#ifndef _VOTCA_KMC_LINKDEVICE_H
#define	_VOTCA_KMC_LINKDEVICE_H

#include <votca/tools/vec.h>
#include <votca/xtp/linksql.h>

namespace votca { namespace xtp {

enum CrossxType{NoxCross, PosxCross, NegxCross };        
enum CrossyType{NoyCross, PosyCross, NegyCross };        
enum CrosszType{NozCross, PoszCross, NegzCross };        
    
/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separations 
 * 
 */
class LinkDevice: public LinkSQL
{

public:

    LinkDevice( long id, Node* node1, Node* node2, votca::tools::vec r12) : LinkSQL(id,node1,node2,r12){
    };
    
    /// Set self image coulomb potential (potential of image charges  of a charge on the charge itself)
    void setSelfImage(double selfimage){_self_image = selfimage;}

    void setCrossxType(int xtype) {_crossxtype = xtype;}
    void setCrossyType(int ytype) {_crossytype = ytype;}
    void setCrosszType(int ztype) {_crossztype = ztype;}
    
    void setRemove(bool remove) {_to_be_removed = remove;}
    
    /// self image coulomb potential (potential of image charges  of a charge on the charge itself)
    const double &self_image() const { return _self_image; }
    
    const bool &remove() const { return _to_be_removed; }
    
    const int &crossxtype() const { return _crossxtype; }
    const int &crossytype() const { return _crossytype; }
    const int &crossztype() const { return _crossztype; }

    

private:
    double _self_image;
    int _crossxtype;
    int _crossytype;
    int _crossztype;
    
    bool _to_be_removed;
    
   
};

}} 

#endif // _VOTCA_KMC_LINKDEVICE_H
