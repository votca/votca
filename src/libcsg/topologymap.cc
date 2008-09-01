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
    
void TopologyMap::Apply(Configuration &conf_in, Configuration &conf_out)
{
    MapContainer::iterator iter;    
    int i=0;
    conf_out.setStep(conf_in.getStep());
    conf_out.setTime(conf_in.getTime());
    conf_out.setBox(conf_in.getBox());

    conf_out.HasPos(conf_in.HasPos());
    conf_out.HasVel(conf_in.HasVel());
    conf_out.HasF(conf_in.HasF());

    for(iter=_maps.begin();iter!=_maps.end();++iter) {
        Molecule min(conf_in, *_in->getMolecule(i));
        Molecule mout(conf_out, *_out->getMolecule(i));
        (*iter)->Apply(min, mout);
        i++;
    }
}
