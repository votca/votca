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

#ifndef _VOTCA_KMC_NODEDEVICE_H
#define	_VOTCA_KMC_NODEDEVICE_H

#include <votca/kmc/nodesql.h>

namespace votca { namespace kmc {

class NodeDevice : public NodeSQL
{
public:
    
    NodeDevice(int id, tools::vec position) : NodeSQL(id, position){
    };

    void setLayer(int layer){_layer = layer;}
    void setSelfImage(double self_image){_self_image = self_image;}

    const double &self_image() const { return _self_image; }     
    const int &layer() const {return _layer;}   
    
private:

    int _layer;
    double _self_image;

    ;
    
};


}}

#endif	/* _VOTCA_KMC_NODEDEVICE_H */

