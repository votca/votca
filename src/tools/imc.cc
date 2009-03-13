/* 
 * File:   imc.cc
 * Author: ruehle
 *
 * Created on March 10, 2009, 3:43 PM
 */

#include <iostream>
#include "beadlist.h"
#include "nblist.h"
#include "imc.h"

void Imc::BeginCG(Topology *top, Topology *top_atom) {
    cout << "begin to calculate inverse monte carlo parameters\n";
    cout << "# of bonded interactions: " << _bonded.size() << endl;
    cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;
    
   _nframes = 0;
   
   // initialize non-bonded structures
   for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter)
        AddInteraction(*iter);
};

interaction *Imc::AddInteraction(Property *p)
{
    string name = p->get("name").value();
    interaction_t *i = new interaction_t;

    _interactions[name] = i;
    // \todo check weather exists       
    
    i->_step = p->get("step").as<double>();
    i->_min = p->get("min").as<double>();
    i->_max = p->get("max").as<double>();

    int n = (i->_max - i->_min) / i->_step;

    i->_current.Initialize(i->_min, i->_max, n);
    i->_average.Initialize(i->_min, i->_max, n);
    
    return *i;
}

void Imc::EndCG()
{        
    map<string, interaction_t *>::iterator iter;
    for(iter=_interactions.begin(); iter!=_interactions.end(); ++iter) {
        cout << (iter->first) << endl;
        cout << iter->second->_average;        
        delete iter->second;
    }
    _interactions.clear();
}
    
void Imc::LoadOptions(const string &file) {
    load_property_from_xml(_options, file);
    _bonded = _options.Select("cg.bonded");
    _nonbonded = _options.Select("cg.non-bonded");
        
}

void Imc::EvalConfiguration(Topology *top, Topology *top_atom) {
    _nframes++;
    DoNonbonded(top);
}

void Imc::DoNonbonded(Topology *top)
{
    for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter) {
        string name = (*iter)->get("name").value();
        
        interaction_t &i = *_interactions[name];
        
        BeadList beads1, beads2;
        
        beads1.Generate(*top, (*iter)->get("type1").value());
        beads2.Generate(*top, (*iter)->get("type2").value());
        
        NBList nb;
        nb.setCutoff(i._max);
        
        if((*iter)->get("type1").value() == (*iter)->get("type2").value())
            nb.Generate(beads1);
        else
            nb.Generate(beads1, beads2);
        
        i._current.Clear();
        NBList::iterator pair_iter;

        for(pair_iter = nb.begin(); pair_iter!=nb.end();++pair_iter) {
                i._current.Process((*pair_iter)->dist());            
        }
        i._average.data().y() += i._current.data().y();
    }    
}