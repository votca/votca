// 
// File:   topologymap.cc
// Author: ruehle
//
// Created on January 16, 2008, 11:24 AM
//

#include "topologymap.h"

TopologyMap::~TopologyMap()
{
    MapContainer::iterator i;
    
    for(i=_maps.begin();i!=_maps.end();++i)
        delete *i;
    _maps.clear();
}
    
void TopologyMap::Apply()
{
    MapContainer::iterator iter;    
    int i=0;
    _out->setStep(_in->getStep());
    _out->setTime(_in->getTime());
    _out->setBox(_in->getBox());

    for(iter=_maps.begin();iter!=_maps.end();++iter) {
        (*iter)->Apply(*_in->getMolecule(i), *_out->getMolecule(i));
        i++;
    }
}
