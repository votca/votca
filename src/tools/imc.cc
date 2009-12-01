/* 
 * File:   imc.cc
 * Author: ruehle
 *
 * Created on March 10, 2009, 3:43 PM
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <votca/tools/rangeparser.h>
#include "beadlist.h"
#include "nblistgrid.h"
#include "imc.h"
#include "imcio.h"

Imc::Imc()
   : _write_every(0), _do_blocks(false), _do_imc(false)
{
}

Imc::~Imc()
{
}

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void Imc::BeginCG(Topology *top, Topology *top_atom) {
    // do some output
    cout << "begin to calculate inverse monte carlo parameters\n";
    cout << "# of bonded interactions: " << _bonded.size() << endl;
    cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;
    
   // we didn't process any frames so far
    _nframes = 0;
    _nblock = 0;
    
   // initialize non-bonded structures
   for (list<Property*>::iterator iter = _nonbonded.begin();
        iter != _nonbonded.end(); ++iter) {
            interaction_t *i = AddInteraction(*iter);
            // calculate normalization factor for rdf
            BeadList beads1, beads2;
        
            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());

            if((*iter)->get("type1").value() ==  (*iter)->get("type2").value())
                i->_norm = top->BoxVolume()/(4.*M_PI* i->_step * beads1.size()*(beads2.size()-1.)/2.);
            else
                i->_norm = top->BoxVolume()/(4.*M_PI* i->_step * beads1.size()*beads2.size());

            i->_is_bonded = false;
   }
   
    // initialize non-bonded structures
   for (list<Property*>::iterator iter = _bonded.begin();
        iter != _bonded.end(); ++iter) {
            interaction_t *i = AddInteraction(*iter);
            i->_is_bonded = true;
   }
   
   // initialize the group structures
    if(_do_imc)
        InitializeGroups();
};

// create an entry for interactions
Imc::interaction_t *Imc::AddInteraction(Property *p)
{
    string name = p->get("name").value();
    string group;
    if(_do_imc)
        group = p->get("inverse.imc.group").value();
    else
        group = "none";
    
    interaction_t *i = new interaction_t;    
    _interactions[name] = i;
    getGroup(group)->_interactions.push_back(i);
    
    i->_step = p->get("step").as<double>();
    i->_min = p->get("min").as<double>();
    i->_max = p->get("max").as<double>();
    i->_norm = 1.0;
    i->_p=p;
    
    // initialize the current and average histogram
    int n = (int)((i->_max - i->_min) / i->_step + 1.000000001);

    i->_current.Initialize(i->_min, i->_max+i->_step, n);
    i->_average.Initialize(i->_min, i->_max+i->_step, n);
    
    return i;
}

// end of trajectory, post processing data
void Imc::EndCG()
{
    if(!_do_blocks) {
        WriteDist();
        if(_do_imc)
            WriteIMCData();
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
    // process bonded interactions
    DoBonded(top);
    // update correlation matrices
    if(_do_imc)
        DoCorrelations();
    
    if(_write_every != 0) {
        if((_nframes % _write_every)==0) {
            _nblock++;
            string suffix = string("_") + boost::lexical_cast<string>(_nblock);
            WriteDist(suffix);
            WriteIMCData(suffix);
            WriteIMCBlock(suffix);
            if(_do_blocks)
                ClearAverages();
        }
    }
            
}

void Imc::ClearAverages()
{
    map<string, interaction_t *>::iterator ic_iter;
    map<string, group_t *>::iterator group_iter;
    
    _nframes = 0;
    for (ic_iter = _interactions.begin(); ic_iter != _interactions.end(); ++ic_iter)
        ic_iter->second->_average.Clear();    
    
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter)
        group_iter->second->_corr.clear();      
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
        nb.setCutoff(i._max + i._step);
        
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
        i._average.data().y() = (((double)_nframes-1.0)*i._average.data().y() 
                + i._current.data().y())/(double)_nframes;
    }    
}

// process non-bonded interactions for current frame
void Imc::DoBonded(Topology *top)
{
    for (list<Property*>::iterator iter = _bonded.begin();
            iter != _bonded.end(); ++iter) {
        string name = (*iter)->get("name").value();
        
        interaction_t &i = *_interactions[name];

        // clear the current histogram
        i._current.Clear();

        // now fill with new data
        std::list<Interaction *> list = top->InteractionsInGroup(name);

        std::list<Interaction *>::iterator ic_iter;
        for(ic_iter=list.begin(); ic_iter!=list.end(); ++ic_iter) {
            Interaction *ic = *ic_iter;
            double v = ic->EvaluateVar(*top);
            i._current.Process(v);
        }
        // update the average
        i._average.data().y() = (((double)_nframes-1.0)*i._average.data().y()
                + i._current.data().y())/(double)_nframes;
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
    if(!_do_imc) return;
    map<string, group_t *>::iterator group_iter;

    // clear all the pairs
    
    // iterator over all groups
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;      
        grp->_pairs.clear();
    
        int n = 0;
        // count number of bins needed in matrix
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1)
            n+=(*i1)->_current.getNBins();
    
        // handy access to matrix
        group_matrix &M = grp->_corr;
        
        // initialize matrix with zeroes
        M.resize(n,n);
        M = ub::zero_matrix<double>(n, n);
        
        // now create references to the sub matrices
        int i, j;
        i=0;j=0;
        // iterate over all possible compinations of pairs
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1) {
            int n1 = (*i1)->_current.getNBins();
            j = i;
            for (list<interaction_t*>::iterator i2 = i1;
                    i2 != grp->_interactions.end(); ++i2) {
                int n2 = (*i2)->_current.getNBins();
                
                // create matrix proxy with sub-matrix
                pair_matrix corr(M, ub::range(i, i+n1), ub::range(j, j+n2));
                // add the pair
                grp->_pairs.push_back(pair_t(*i1, *i2, i, j, corr));
                j+=n2;
            }
            i+=n1;
        }
    }
}

// update the correlation matrix
void Imc::DoCorrelations() {
    if(!_do_imc) return;
    vector<pair_t>::iterator pair;
    map<string, group_t *>::iterator group_iter;
        
     for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;      
        // update correlation for all pairs
        for (pair = grp->_pairs.begin(); pair != grp->_pairs.end(); ++pair) {
            ub::vector<double> &a = pair->_i1->_current.data().y();
            ub::vector<double> &b = pair->_i2->_current.data().y();
            pair_matrix &M = pair->_corr;

            // M_ij += a_i*b_j
            //for(int i=0; i<M.size1(); ++i)
            //    for(int j=i; j<M.size2(); ++j)
            //        M(i,j) = ((((double)_nframes-1.0)*M(i,j)) + (a(i)*b(j)))/(double)_nframes;
            M = ((((double)_nframes-1.0)*M) + ub::outer_prod(a, b))/(double)_nframes;
        }
    }
}

// write the distribution function
void Imc::WriteDist(const string &suffix)
{
    map<string, interaction_t *>::iterator iter;

    // for all interactions
    for (iter = _interactions.begin(); iter != _interactions.end(); ++iter) {            
        // calculate the rdf
        Table &t = iter->second->_average.data();            
        Table dist(t);

        if(!iter->second->_is_bonded) {
            dist.y() = iter->second->_norm *
                element_div(dist.y(),
                    element_prod(dist.x(), dist.x())
                );
        }
        else {
            dist.y() = iter->second->_norm * dist.y();
        }
        
        dist.Save((iter->first) + suffix + ".dist.new");
        cout << "written " << (iter->first) + suffix + ".dist.new\n";            
    }
}


/**
 *  Here the inverse monte carlo matrix is calculated and written out
 *
 *  steps:
 *      - calculate th
 */
void Imc::WriteIMCData(const string &suffix) {            
    if(!_do_imc) return;
    //map<string, interaction_t *>::iterator ic_iter;
    map<string, group_t *>::iterator group_iter;
    
    // iterate over all groups
    for(group_iter = _groups.begin(); group_iter!=_groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;    
        string grp_name = (*group_iter).first;
        list<interaction_t *>::iterator iter;
        
        // number of total bins for all interactions in group is matrix dimension
        int n=grp->_corr.size1();
                
        // build full set of equations + copy some data to make
        // code better to read
        group_matrix gmc(grp->_corr);
        ub::vector<double> dS(n);
        ub::vector<double> r(n);
        // the next two variables are to later extract the individual parts
        // from the whole data after solving equations
        vector<RangeParser> ranges; // sizes of the individual interactions
        vector<string> names; // names of the interactions
                        
        // copy all averages+r of group to one vector
        n=0;
        int begin=1;
        for(iter=grp->_interactions.begin(); iter != grp->_interactions.end(); ++iter) {
            interaction_t *ic = *iter;
            
            // sub vector for dS
            ub::vector_range< ub::vector<double> > sub_dS(dS, 
                    ub::range(n, n + ic->_average.getNBins()));
            
            // sub vector for r
            ub::vector_range< ub::vector<double> > sub_r(r, 
                    ub::range(n, n + ic->_average.getNBins()));
            
            // read in target and calculate dS
            CalcDeltaS(ic, sub_dS);
            
            // copy r
            sub_r = ic->_average.data().x();
            
            // save size

            RangeParser rp;
            int end = begin  + ic->_average.getNBins() -1;
            rp.Add(begin, end);
            ranges.push_back(rp);
            begin = end+1;
            // save name
            names.push_back(ic->_p->get("name").as<string>());
            
            // shift subrange by size of current
            n+=ic->_average.getNBins();
        }
        
        // now we need to calculate the 
        // A_ij = <S_i*S_j> - <S_i>*<S_j>        
        vector<pair_t>::iterator pair;
        for (pair = grp->_pairs.begin(); pair != grp->_pairs.end(); ++pair) {
            interaction_t *i1 = pair->_i1;
            interaction_t *i2 = pair->_i2;
            
            // make reference to <S_i>
            ub::vector<double> &a = i1->_average.data().y();
            // make reference to <S_j>
            ub::vector<double> &b = i2->_average.data().y();
            
            int i=pair->_offset_i;
            int j=pair->_offset_j;
            int n1=i1->_average.getNBins();
            int n2=i2->_average.getNBins();
            
            // sub matrix for these two interactions
            // we only need to take care about one sub-matrix and not the mirrored
            // one since ublas makes sure the matrix is symmetric
            pair_matrix M(gmc, ub::range(i, i+n1),
                               ub::range(j, j+n2));
            // A_ij = -(<a_i*a_j>  - <a_i>*<b_j>)
            //for(i=0; i<M.size1(); ++i)
            //    for(j=i; j<M.size2(); ++j)
            //        M(i,j) = -(M(i,j) - a(i)*b(j));
            M = -(M - ub::outer_prod(a, b));
        }
        
        imcio_write_dS(grp_name + suffix + ".imc", r, dS);
        imcio_write_matrix(grp_name + suffix + ".gmc", gmc);
        imcio_write_index(grp_name + suffix + ".idx", names, ranges);

    }
}

// calculate deviation from target vectors
void Imc::CalcDeltaS(interaction_t *interaction, ub::vector_range< ub::vector<double> > &dS)
{
    const string &name = interaction->_p->get("name").as<string>();
                
    Table target;
    target.Load(name + ".dist.tgt");
                      
    if(!interaction->_is_bonded) {
        target.y() = (1.0 / interaction->_norm)*ub::element_prod(target.y(),
            (ub::element_prod(target.x(), target.x()))
            ) ;
    }
    else {
            target.y() = (1.0 / interaction->_norm)*target.y();
    }
    if(target.y().size() !=  interaction->_average.data().y().size())
        throw std::runtime_error("number of grid points in target does not match the grid");

    dS = interaction->_average.data().y() - target.y();
}


void Imc::WriteIMCBlock(const string &suffix)
{
    if(!_do_imc) return;
    //map<string, interaction_t *>::iterator ic_iter;
    map<string, group_t *>::iterator group_iter;

    // iterate over all groups
    for(group_iter = _groups.begin(); group_iter!=_groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;
        string grp_name = (*group_iter).first;
        list<interaction_t *>::iterator iter;

        // number of total bins for all interactions in group is matrix dimension
        int n=grp->_corr.size1();

        // build full set of equations + copy some data to make
        // code better to read
        group_matrix gmc(grp->_corr);
        ub::vector<double> dS(n);
        ub::vector<double> r(n);
        // the next two variables are to later extract the individual parts
        // from the whole data after solving equations
        vector<int> sizes; // sizes of the individual interactions
        vector<string> names; // names of the interactions

        // copy all averages+r of group to one vector
        n=0;
        for(iter=grp->_interactions.begin(); iter != grp->_interactions.end(); ++iter) {
            interaction_t *ic = *iter;

            // sub vector for dS
            ub::vector_range< ub::vector<double> > sub_dS(dS,
                    ub::range(n, n + ic->_average.getNBins()));

            // sub vector for r
            ub::vector_range< ub::vector<double> > sub_r(r,
                    ub::range(n, n + ic->_average.getNBins()));

            // read in target and calculate dS
            sub_dS = ic->_average.data().y();
            // copy r
            sub_r = ic->_average.data().x();
            // save size
            sizes.push_back(ic->_average.getNBins());
            // save name
            names.push_back(ic->_p->get("name").as<string>());

            // shift subrange by size of current
            n+=ic->_average.getNBins();
        }

        // write the dS
        ofstream out_dS;
        string name_dS = grp_name + suffix + ".S";
        out_dS.open(name_dS.c_str());
        out_dS << setprecision(8);
        if(!out_dS)
            throw runtime_error(string("error, cannot open file ") + name_dS);

        for(int i=0; i<dS.size(); ++i) {
            out_dS << r[i] << " " << dS[i] << endl;
        }

        out_dS.close();
        cout << "written " << name_dS << endl;

        // write the correlations
        ofstream out_cor;
        string name_cor = grp_name + suffix + ".cor";
        out_cor.open(name_cor.c_str());
        out_cor << setprecision(8);

        if(!out_cor)
            throw runtime_error(string("error, cannot open file ") + name_cor);

        for(group_matrix::size_type i=0; i<grp->_corr.size1(); ++i) {
            for(group_matrix::size_type j=0; j<grp->_corr.size2(); ++j) {
                out_cor << grp->_corr(i, j) << " ";
            }
            out_cor << endl;
        }
        out_cor.close();
        cout << "written " << name_cor << endl;
    }
}

