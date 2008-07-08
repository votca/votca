// 
// File:   bondedstatistics.cc
// Author: victor
//
// Created on 4. Juni 2008, 17:43
//

#include "bondedstatistics.h"

void BondedStatistics::BeginCG(Topology *top, Topology *top_atom)
{
    InteractionContainer &ic = top->getBondedInteractions();
    InteractionContainer::iterator ia;
    
    _bonded_values.clear();
    for(ia=ic.begin(); ia!=ic.end(); ++ia) {
        _bonded_values.CreateArray((*ia)->getName());
    }
}

void BondedStatistics::EndCG()
{
}

void BondedStatistics::EvalConfiguration(Configuration *conf, Configuration *conv_atom)
{
    InteractionContainer &ic = conf->getTopology()->getBondedInteractions();
    InteractionContainer::iterator ia;
    
    DataCollection<double>::container::iterator is;
    
    for(ia=ic.begin(), is = _bonded_values.begin(); ia != ic.end(); ++ia, ++is) {
//        const string &name = (*ia)->getName();        
        (*is)->push_back((*ia)->EvaluateVar(*conf));
    }
}