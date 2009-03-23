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

Imc::interaction_t *Imc::AddInteraction(Property *p)
{
    string name = p->get("name").value();
    string group = p->get("group").value();
    
    interaction_t *i = new interaction_t;

    _interactions[name] = i;
    // \todo check weather exists       
    getGroup(group)->_interactions.push_back(i);
    
    i->_step = p->get("step").as<double>();
    i->_min = p->get("min").as<double>();
    i->_max = p->get("max").as<double>();

    int n = (i->_max - i->_min) / i->_step;

    i->_current.Initialize(i->_min, i->_max, n);
    i->_average.Initialize(i->_min, i->_max, n);
    
    return i;
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

Imc::group_t *Imc::getGroup(const string &name)
{
    map<string, group_t *>::iterator iter;
    iter = _groups.find(name);
    if(iter == _groups.end())
        return _groups[name] = new group_t;
    return (*iter).second;
}

void Imc::InitializeGroups()
{
    map<string, group_t *>::iterator group_iter;

    // iterator over all groups
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;
        grp->_correlations.clear();
        int n = grp->_interactions.size();
        // resize vector for correlation matrices
        grp->_correlations.resize(0.5 * n * (n + 1));
        
        vector< ub::matrix<double> >::iterator corr;
        corr = grp->_correlations.begin();
        
        // iterate over all possible compinations for correlations
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1) {
            for (list<interaction_t*>::iterator i2 = i1;
                    i2 != grp->_interactions.end(); ++i2, ++corr) {
                // allocate memory for correlation matrix
                corr->resize((*i1)->_current.getNBins(), (*i2)->_current.getNBins());
            }
        }
    }
}

void Imc::UpdateCorrelations()
{
        map<string, group_t *>::iterator group_iter;

    // iterator over all groups
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;
        grp->_correlations.clear();
        int n = grp->_interactions.size();
        
        vector< ub::matrix<double> >::iterator corr;
        corr = grp->_correlations.begin();
        
        // iterate over all possible compinations for correlations
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1) {
            for (list<interaction_t*>::iterator i2 = i1;
                    i2 != grp->_interactions.end(); ++i2, ++corr) {
                Correlate(*corr, (*i1)->_current.data().y(),
                    (*i1)->_current.data().y());
            }
        }
    }
}

void Imc::Correlate(ub::matrix<double> &corr, ub::vector<double> &v1,
        ub::vector<double> &v2)
{
    ub::matrix<double>::size_type x;
    ub::matrix<double>::size_type y;
    
    for(x=0; x<v1.size(); ++x)
        for(y=0; y<v2.size(); ++y)
            corr(x)(y)=v1[x]*v2[y];
}
