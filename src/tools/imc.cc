/* 
 * File:   imc.cc
 * Author: ruehle
 *
 * Created on March 10, 2009, 3:43 PM
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include "beadlist.h"
#include "nblist.h"
#include "imc.h"

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void Imc::BeginCG(Topology *top, Topology *top_atom) {
    // do some output
    cout << "begin to calculate inverse monte carlo parameters\n";
    cout << "# of bonded interactions: " << _bonded.size() << endl;
    cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;
    
   // we didn't process any frames so far
    _nframes = 0;
   
// initialize non-bonded structures
   for (list<Property*>::iterator iter = _nonbonded.begin();
        iter != _nonbonded.end(); ++iter) {
            interaction_t *i = AddInteraction(*iter);
            // calculate normalization factor for rdf
            BeadList beads1, beads2;
        
            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());

            i->_norm = top->BoxVolume()/(4.*M_PI* i->_step * beads1.size()*(beads2.size()+1.)/2.);
   }
   
   // initialize the group structures
    InitializeGroups();
};

// create an entry for interactions
Imc::interaction_t *Imc::AddInteraction(Property *p)
{
    string name = p->get("name").value();
    string group = p->get("imc.group").value();
    
    interaction_t *i = new interaction_t;    
    _interactions[name] = i;
    getGroup(group)->_interactions.push_back(i);
    
    i->_step = p->get("step").as<double>();
    i->_min = p->get("min").as<double>();
    i->_max = p->get("max").as<double>();
    i->_norm = 1.0;
    i->_p=p;
    
    // initialize the current and average histogram
    int n = (int)((i->_max - i->_min) / i->_step) + 1;

    i->_current.Initialize(i->_min, i->_max+i->_step, n);
    i->_average.Initialize(i->_min, i->_max+i->_step, n);
    
    return i;
}

// end of trajectory, post processing data
void Imc::EndCG()
{
    CalcDeltaS();
    // transform correlation matrix to update matrix
    CalcMatrix();

    // write out interactions
    {
        map<string, interaction_t *>::iterator iter;

        // for all interactions
        for (iter = _interactions.begin(); iter != _interactions.end(); ++iter) {
            
            // calculate the rdf
            Table &t = iter->second->_average.data();            
            Table rdf(t);
            ub::scalar_vector<double> step(rdf.x().size(), 0.5*iter->second->_step);

            rdf.y() = iter->second->_norm * element_div(rdf.y(),
                    element_prod(rdf.x()+step, rdf.x()+step)
                    );
            rdf.Save((iter->first) + ".dist.new");
            cout << "written " << (iter->first) + ".dist.new\n";
            
            delete iter->second;
        }
    }

    // write matrix
    {
        map<string,group_t *>::iterator iter;

        for (iter = _groups.begin(); iter != _groups.end(); ++iter) {
            string filename = (iter->first) + ".gmc";
            ofstream out;
            out.open(filename.c_str());
            group_matrix &M = iter->second->_corr;
            if (!out)
                throw string("error, cannot open file ") + filename;
            // out << M.size1() << " " << M.size2() << endl;
            out << setprecision(8);
            for (int i = 0; i < M.size1(); ++i) {
                for (int j = 0; j < M.size2(); ++j) {
                    out << M(i, j) << " ";
                }
                out << endl;
            }
            cout << "written " << filename << endl;
            delete (*iter).second;
        }
    }

    // clear interactions and groups
    _interactions.clear();
    _groups.clear();
}

// load options from xml file
void Imc::LoadOptions(const string &file) {
    load_property_from_xml(_options, file);
    _bonded = _options.Select("cg.bonded");
    _nonbonded = _options.Select("cg.non-bonded");
        
}

// evaluate current conformation
void Imc::EvalConfiguration(Topology *top, Topology *top_atom) {
    _nframes++;
    // process non-bonded interactions
    DoNonbonded(top);
    // update correlation matrices
    UpdateCorrelations();
}

// process non-bonded interactions for current frame
void Imc::DoNonbonded(Topology *top)
{
    for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter) {
        string name = (*iter)->get("name").value();
        
        interaction_t &i = *_interactions[name];
        
        // generate the bead lists
        BeadList beads1, beads2;
        
        beads1.Generate(*top, (*iter)->get("type1").value());
        beads2.Generate(*top, (*iter)->get("type2").value());
        
        // generate the neighbour list
        NBList nb;
        nb.setCutoff(i._max);
        
        // is it same types or different types?
        if((*iter)->get("type1").value() == (*iter)->get("type2").value())
            nb.Generate(beads1);
        else
            nb.Generate(beads1, beads2);
        
        // clear the current histogram
        i._current.Clear();
        
        // process all pairs
        NBList::iterator pair_iter;
        for(pair_iter = nb.begin(); pair_iter!=nb.end();++pair_iter) {
                i._current.Process((*pair_iter)->dist());            
        }
        
        // update the average
        i._average.data().y() += i._current.data().y();
    }    
}

// returns a group, creates it if doesn't exist
Imc::group_t *Imc::getGroup(const string &name)
{
    map<string, group_t *>::iterator iter;
    iter = _groups.find(name);
    if(iter == _groups.end()) {
        return _groups[name] = new group_t;
    }
    return (*iter).second;
}

// initialize the groups after interactions are added
void Imc::InitializeGroups()
{
    map<string, group_t *>::iterator group_iter;

    // clear all the pairs
    _pairs.clear();

    // iterator over all groups
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;
        
        int n = 0;
        // count number of bins needed in matrix
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1)
            n+=(*i1)->_current.getNBins();
    
        // handy access to matrix
        group_matrix &M = grp->_corr;
        
        // initialize matrix with zeroes
        M.resize(n);
        M = ub::zero_matrix<double>(n, n);
        
        
        // now create references to the sub matrices
        int i, j;
        i=0;j=0;
        // iterate over all possible compinations of pairs
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1) {
            int n1 = (*i1)->_current.getNBins();
            for (list<interaction_t*>::iterator i2 = i1;
                    i2 != grp->_interactions.end(); ++i2) {
                int n2 = (*i2)->_current.getNBins();
                
                // create matrix proxy with sub-matrix
                pair_matrix corr(M, ub::range(i, i+n1), ub::range(j, j+n2));
                // add the pair
                _pairs.push_back(pair_t(*i1, *i2, corr));
                j+=n2;
            }
            i+=n1;
        }
    }
}

// update the correlation matrix
void Imc::UpdateCorrelations() {
    vector<pair_t>::iterator pair;
    
    // update correlation for all pairs
    for (pair = _pairs.begin(); pair != _pairs.end(); ++pair) {
        ub::vector<double> &a = pair->_i1->_current.data().y();
        ub::vector<double> &b = pair->_i2->_current.data().y();
        pair_matrix &M = pair->_corr;

        // M_ij += a_i*b_j
        M += ub::outer_prod(a, b);
    }
}

// transform correlations to update matrix
void Imc::CalcMatrix() {
    map<string, interaction_t *>::iterator ic_iter;

    double norm = (1. / (double) _nframes);

    vector<pair_t>::iterator pair;

    for (pair = _pairs.begin(); pair != _pairs.end(); ++pair) {
        ub::vector<double> &a = pair->_i1->_average.data().y();
        ub::vector<double> &b = pair->_i2->_average.data().y();
        pair_matrix &M = pair->_corr;

        // divide by number of frames
        M *= norm;
        // M_ij = -(M_ij  - a_ib_j)
        M = -(M - ub::outer_prod(a, b));
    }
}

// calculate deviation from target vectors
void Imc::CalcDeltaS()
{
    map<string, interaction_t *>::iterator iter;
    double norm = (1. / (double) _nframes);
    for(iter = _interactions.begin(); iter != _interactions.end(); ++iter) {
        interaction_t *i = iter->second;
        const string &name = iter->first;
        i->_average.data().y() *= norm;
                
        Table target;
        target.Load(i->_p->get("target").as<string>());
        
        i->_average.data().Save(name + ".S");
        
        ub::scalar_vector<double> step(target.x().size(), 0.5*i->_step);
        
        target.y() = (1.0 / i->_norm)*ub::element_prod(target.y(), 
            (ub::element_prod(target.x()+step, target.x()+step))
            ) ;
        target.Save(name + ".St");
        
        
        if(target.y().size() !=  i->_average.data().y().size())
            throw std::runtime_error("number of grid points in target does not match the grid");
        
        target.y() = i->_average.data().y() - target.y();
        // write S
        target.Save(name + ".imc");
            
        cout << "written " << name + ".imc\n";
    }
}
