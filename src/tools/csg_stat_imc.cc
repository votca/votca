/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <votca/tools/rangeparser.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/csg/nblistgrid_3body.h>
#include "csg_stat_imc.h"
#include <votca/csg/imcio.h>


namespace votca { namespace csg {

Imc::Imc()
   : _block_length(0), _do_imc(false), _processed_some_frames(false)
{
}

Imc::~Imc()
{
}

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void Imc::Initialize()
{
    // do some output
    if(_do_imc)
	    cout << "begin to calculate inverse monte carlo parameters\n";
    else
	    cout << "begin to calculate distribution functions\n";
    cout << "# of bonded interactions: " << _bonded.size() << endl;
    cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;

    if ( _bonded.size()+_nonbonded.size() == 0 )
            throw std::runtime_error("No interactions defined in options xml-file - nothing to be done");

    
   // initialize non-bonded structures
   for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter) {
        interaction_t *i = AddInteraction(*iter);
        i->_is_bonded = false;
    }

    // initialize bonded structures
   for (list<Property*>::iterator iter = _bonded.begin();
        iter != _bonded.end(); ++iter) {
            interaction_t *i = AddInteraction(*iter);
            i->_is_bonded = true;
   }
   
   // initialize the group structures
    if(_do_imc)
        InitializeGroups();
};

void Imc::BeginEvaluate(Topology *top, Topology *top_atom)
{
  // we didn't process any frames so far
    _nframes = 0;
    _nblock = 0;
    _processed_some_frames=false;

    // initialize non-bonded structures
   for (list<Property*>::iterator iter = _nonbonded.begin();
            iter != _nonbonded.end(); ++iter) {
        string name = (*iter)->get("name").value();

        interaction_t &i = *_interactions[name];
        
        //Preleminary: Quickest way to incorporate 3 body correlations
        if (i._threebody){
            
            // generate the bead lists
            BeadList beads1, beads2, beads3;

            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());
            beads3.Generate(*top, (*iter)->get("type3").value());

            if(beads1.size() == 0)
                throw std::runtime_error("Topology does not have beads of type \"" + (*iter)->get("type1").value() + "\"\n"
                        "This was specified in type1 of interaction \"" + name+ "\"");
            if(beads2.size() == 0)
                throw std::runtime_error("Topology does not have beads of type \"" + (*iter)->get("type2").value() + "\"\n"
                        "This was specified in type2 of interaction \"" + name + "\"");
            if(beads3.size() == 0)
                throw std::runtime_error("Topology does not have beads of type \"" + (*iter)->get("type3").value() + "\"\n"
                        "This was specified in type3 of interaction \"" + name + "\"");             
        }
        //2body
        if (!i._threebody){

            // generate the bead lists
            BeadList beads1, beads2;

            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());

            if(beads1.size() == 0)
                throw std::runtime_error("Topology does not have beads of type \"" + (*iter)->get("type1").value() + "\"\n"
                        "This was specified in type1 of interaction \"" + name+ "\"");
            if(beads2.size() == 0)
                throw std::runtime_error("Topology does not have beads of type \"" + (*iter)->get("type2").value() + "\"\n"
                        "This was specified in type2 of interaction \"" + name + "\"");

            // calculate normalization factor for rdf
            if ((*iter)->get("type1").value() == (*iter)->get("type2").value())
                i._norm = 1. / (beads1.size()*(beads2.size()) / 2.);
            else
                i._norm = 1. / (beads1.size() * beads2.size());            
        }


    }

    for (list<Property*>::iterator iter = _bonded.begin();
            iter != _bonded.end(); ++iter) {
        string name = (*iter)->get("name").value();

        std::list<Interaction *> list = top->InteractionsInGroup(name);
        if (list.empty() )
            throw std::runtime_error("Bonded interaction '" + name + "' defined in options xml-file, but not in topology - check name definition in the mapping file again");
    }
}

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
    i->_index = _interactions.size();
    _interactions[name] = i;
    getGroup(group)->_interactions.push_back(i);
    
    i->_step = p->get("step").as<double>();
    i->_min = p->get("min").as<double>();
    i->_max = p->get("max").as<double>();
    i->_norm = 1.0;
    i->_p=p;
    i->_threebody=0;
    i->_force=0;    
    //check if option threebody exists
    if (p->exists("threebody")) {
        i->_threebody = p->get("threebody").as<bool>();            
    }
    //check if option force exists
    if (p->exists("force")) {
        i->_force = p->get("force").as<bool>();            
    }    
    //preliminary. should be changed
    i->_cut = 0.37; //(0.37 nm)
    //check if option threebody exists
    if (p->exists("cut")) {
        i->_cut = p->get("cut").as<double>();            
    }    
    
    // initialize the current and average histogram
    int n = static_cast<int>((i->_max - i->_min) / i->_step + 1.000000001);

    i->_average.Initialize(i->_min, i->_max, n);
    if (i->_force){
        i->_average_force.Initialize(i->_min, i->_max, n);
    }
    
    return i;
}

// end of trajectory, post processing data
void Imc::EndEvaluate()
{
    if(_nframes > 0) {
        if(_block_length == 0) {
	    string suffix = string(".") + _extension;
            WriteDist(suffix);
            if(_do_imc)
                WriteIMCData();
        }
    }
    // clear interactions and groups
    _interactions.clear();
    _groups.clear();
    if(!_processed_some_frames)
        throw std::runtime_error("no frames were processed. Please check your input");
}

// load options from xml file
void Imc::LoadOptions(const string &file) {
    load_property_from_xml(_options, file);
    _bonded = _options.Select("cg.bonded");
    _nonbonded = _options.Select("cg.non-bonded");
        
}

// evaluate current conformation
void Imc::Worker::EvalConfiguration(Topology *top, Topology *top_atom) {
    
    _cur_vol = top->BoxVolume();
    // process non-bonded interactions
    DoNonbonded(top);
    // process bonded interactions
    DoBonded(top);            
}

void Imc::ClearAverages()
{
    map<string, interaction_t *>::iterator ic_iter;
    map<string, group_t *>::iterator group_iter;
    
    _nframes = 0;
    for (ic_iter = _interactions.begin(); ic_iter != _interactions.end(); ++ic_iter)
        ic_iter->second->_average.Clear();
        if (ic_iter->second->_force){
            ic_iter->second->_average_force.Clear();
        }
    
    for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter){
        group_iter->second->_corr.setZero(); 
    }
}

class IMCNBSearchHandler {
public:
    IMCNBSearchHandler(HistogramNew *hist)
            : _hist(hist) {}

    HistogramNew *_hist;

    bool FoundPair(Bead *b1, Bead *b2, const vec &r, const double dist) {
        _hist->Process(dist);
        return false;
    }
};


// process non-bonded interactions for current frame
void Imc::Worker::DoNonbonded(Topology *top)
{
    for (list<Property*>::iterator iter = _imc->_nonbonded.begin();
            iter != _imc->_nonbonded.end(); ++iter) {
        string name = (*iter)->get("name").value();
        
        interaction_t &i = *_imc->_interactions[name];

        // clear the current histogram
        _current_hists[i._index].Clear();
        _current_hists_force[i._index].Clear();      
        
        bool gridsearch=true;

        if(_imc->_options.exists("cg.nbsearch")) {
            if(_imc->_options.get("cg.nbsearch").as<string>() == "grid")
                gridsearch=true;
            else if(_imc->_options.get("cg.nbsearch").as<string>() == "simple")
                gridsearch=false;
            else throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
        }    
        
        //Preleminary: Quickest way to incorporate 3 body correlations
        if (i._threebody){
            
            // generate the bead lists
            BeadList beads1, beads2, beads3;
        
            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());            
            beads3.Generate(*top, (*iter)->get("type3").value());            
            
            // generate the neighbour list    
            NBList_3Body *nb; 
 
            if(gridsearch) 
                nb = new NBListGrid_3Body();            
            else
                nb = new NBList_3Body();                        
            
            nb->setCutoff(i._cut); // implement different cutoffs for different interactions!
            //Here, a is the distance between two beads of a triple, where the 3-body interaction is zero              
            
            //check if type1 and type2 are the same
            if ((*iter)->get("type1").value() == (*iter)->get("type2").value()){
                //if all three types are the same
                if ((*iter)->get("type2").value() == (*iter)->get("type3").value()){        
                    nb->Generate(beads1, true);
                }
                //if type2 and type3 are different, use the Generate function for 2 bead types
                if ((*iter)->get("type2").value() != (*iter)->get("type3").value()){           
                    nb->Generate(beads1, beads3, true); 
                }
            }
            //if type1 != type2
            if ((*iter)->get("type1").value() != (*iter)->get("type2").value()){   
                //if the last two types are the same, use Generate function with them as the first two bead types
                //Neighborlist_3body is constructed in a way that the two equal bead types have two be the first 2 types
                if ((*iter)->get("type2").value() == (*iter)->get("type3").value()){
                    nb->Generate(beads1, beads2, true);         
                }
                if ((*iter)->get("type2").value() != (*iter)->get("type3").value()){
                    //type1 = type3 !=type2
                    if ((*iter)->get("type1").value() == (*iter)->get("type3").value()){          
                        nb->Generate(beads2, beads1, true);
                    }
                    //type1 != type2 != type3
                    if ((*iter)->get("type1").value() != (*iter)->get("type3").value()){              
                        nb->Generate(beads1, beads2, beads3, true);
                    }
                }
            }            
            
            NBList_3Body::iterator triple_iter;
            // iterate over all triples
            for (triple_iter = nb->begin(); triple_iter != nb->end(); ++triple_iter){
                vec rij  = (*triple_iter)->r12();
                vec rik  = (*triple_iter)->r13();
                double var = acos(rij*rik/sqrt((rij*rij) * (rik*rik))); 
                _current_hists[i._index].Process(var);               
            }

            delete nb;   
            
        }
        //2body interaction
        if (!i._threebody){
            
            // generate the bead lists
            BeadList beads1, beads2;
        
            beads1.Generate(*top, (*iter)->get("type1").value());
            beads2.Generate(*top, (*iter)->get("type2").value());
               
            // generate the neighbour list
            NBList *nb;

            if(gridsearch)
                nb = new NBListGrid();
            else
                nb = new NBList();

            nb->setCutoff(i._max + i._step);

            IMCNBSearchHandler h(&(_current_hists[i._index]));

            nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

            // is it same types or different types?
            if((*iter)->get("type1").value() == (*iter)->get("type2").value())
                nb->Generate(beads1);
            else
                nb->Generate(beads1, beads2);

            delete nb;

            //if one wants to calculate the mean force
            if (i._force){
                if(gridsearch)
                    nb = new NBListGrid();
                else
                    nb = new NBList();            
         
                nb->setCutoff(i._max + i._step);
            
                // is it same types or different types?
                if((*iter)->get("type1").value() == (*iter)->get("type2").value())
                    nb->Generate(beads1);
                else
                    nb->Generate(beads1, beads2);            

                //process all pairs to calculate the projection of the 
                //mean force on bead 1 on the pair distance: F1 * r12
                NBList::iterator pair_iter;
                for(pair_iter = nb->begin(); pair_iter!=nb->end();++pair_iter) {
                    vec F2 = (*pair_iter)->second->getF();
                    vec F1 = (*pair_iter)->first->getF();                    
                    vec r12 = (*pair_iter)->r();
                    r12.normalize();
                    double var = (*pair_iter)->dist();
                    double scale = F2*r12;
                    _current_hists_force[i._index].Process(var,scale);
                }
                delete nb;                 
            }           
        }
    }    
}

// process non-bonded interactions for current frame
void Imc::Worker::DoBonded(Topology *top)
{
    for (list<Property*>::iterator iter = _imc->_bonded.begin();
            iter != _imc->_bonded.end(); ++iter) {
        string name = (*iter)->get("name").value();
        
        interaction_t &i = *_imc->_interactions[name];

        // clear the current histogram
        _current_hists[i._index].Clear();

        // now fill with new data
        std::list<Interaction *> list = top->InteractionsInGroup(name);

        std::list<Interaction *>::iterator ic_iter;
        for(ic_iter=list.begin(); ic_iter!=list.end(); ++ic_iter) {
            Interaction *ic = *ic_iter;
            double v = ic->EvaluateVar(*top);
            _current_hists[i._index].Process(v);
        }
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
            n+=(*i1)->_average.getNBins();
    
        // handy access to matrix
        group_matrix &M = grp->_corr;
        
        // initialize matrix with zeroes
        M=Eigen::MatrixXd::Zero(n,n);
        
        // now create references to the sub matrices
        int i, j;
        i=0;j=0;
        // iterate over all possible compinations of pairs
        for (list<interaction_t*>::iterator i1 = grp->_interactions.begin();
                i1 != grp->_interactions.end(); ++i1) {
            int n1 = (*i1)->_average.getNBins();
            j = i;
            for (list<interaction_t*>::iterator i2 = i1;
                    i2 != grp->_interactions.end(); ++i2) {
                int n2 = (*i2)->_average.getNBins();
                
                // create matrix proxy with sub-matrix
                pair_matrix corr=M.block(i,j,n1,n2);
                // add the pair
                grp->_pairs.push_back(pair_t(*i1, *i2, i, j, corr));
                j+=n2;
            }
            i+=n1;
        }
    }
}

// update the correlation matrix
void Imc::DoCorrelations(Imc::Worker *worker) {
    if(!_do_imc) return;
    vector<pair_t>::iterator pair;
    map<string, group_t *>::iterator group_iter;
        
     for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter) {
        group_t *grp = (*group_iter).second;      
        // update correlation for all pairs
        for (pair = grp->_pairs.begin(); pair != grp->_pairs.end(); ++pair) {
            Eigen::VectorXd &a = worker->_current_hists[pair->_i1->_index].data().y();
            Eigen::VectorXd &b = worker->_current_hists[pair->_i2->_index].data().y();
            pair_matrix &M = pair->_corr;

            M = ((((double)_nframes-1.0)*M) + a*b.transpose())/(double)_nframes;
        }
    }
}

// write the distribution function
void Imc::WriteDist(const string &suffix)
{
    map<string, interaction_t *>::iterator iter;

    cout << std::endl;  // Cosmetic, put \n before printing names of distribution files.
    // for all interactions
    for (iter = _interactions.begin(); iter != _interactions.end(); ++iter) {            
        // calculate the rdf
        Table &t = iter->second->_average.data();            
        Table dist(t);
        
        //if no average force calculation, dummy table
        Table force;
        //if average force calculation, table force contains force data
        if (iter->second->_force){
            Table &f = iter->second->_average_force.data();
            force = f;
        }    
        
        if(!iter->second->_is_bonded) {
            //Quickest way to incorporate 3 body correlations
            if (iter->second->_threebody){
	        // \TODO normalize bond and angle differently....
                double norm=dist.y().cwiseAbs().sum();
                if ( norm > 0 ) {
                    dist.y() = iter->second->_norm * dist.y() / ( norm * iter->second->_step );
                }           
                
            }
            
            //2body
            if (!iter->second->_threebody){
                //force normalization
                //normalize by number of pairs found at a specific distance
                for(unsigned int i=0; i<force.y().size(); ++i) {
                    //check if any number of pairs has been found at this distance, then normalize
                    if (dist.y()[i]!=0){
                        force.y()[i]/=dist.y()[i];
                    }
                    //else set to zero
                    else{
                        force.y()[i]=0;
                    }
                }                
                
                // normalization is calculated using exact shell volume (difference of spheres)
                for(unsigned int i=0; i<dist.y().size(); ++i) {
                    double x1 = dist.x()[i] - 0.5*iter->second->_step;
                    double x2 = x1 + iter->second->_step;
                    if(x1<0) {
                        dist.y()[i]=0;
                    }
                    else {
                        dist.y()[i] = _avg_vol.getAvg()*iter->second->_norm *
                            dist.y()[i]/(4./3.*M_PI*(x2*x2*x2 - x1*x1*x1));                
                    }
                }                
            }

        }
        else {
	    // \TODO normalize bond and angle differently....
            double norm=dist.y().cwiseAbs().sum();
            if ( norm > 0 ) {
                dist.y() = iter->second->_norm * dist.y() / ( norm * iter->second->_step );
            }
        }
        
        dist.Save((iter->first) + suffix);
        cout << "written " << (iter->first) + suffix << "\n";

        //preliminary
        if (iter->second->_force){
            force.Save((iter->first) + ".force.new");
            cout << "written " << (iter->first) + ".force.new" << "\n";            
        } 
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
        int n=grp->_corr.rows();
                
        // build full set of equations + copy some data to make
        // code better to read
        group_matrix gmc(grp->_corr);
        Eigen::VectorXd dS(n);
        Eigen::VectorXd r(n);
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
            Eigen::VectorBlock< Eigen::VectorXd > sub_dS=dS.segment(n,ic->_average.getNBins());
            
            // sub vector for r
           Eigen::VectorBlock< Eigen::VectorXd > sub_r=r.segment(n,ic->_average.getNBins());
            
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
           Eigen::VectorXd &a = i1->_average.data().y();
            // make reference to <S_j>
            Eigen::VectorXd &b = i2->_average.data().y();
            
            int i=pair->_offset_i;
            int j=pair->_offset_j;
            int n1=i1->_average.getNBins();
            int n2=i2->_average.getNBins();
            
          
            pair_matrix M=gmc.block(i,j,n1,n2);
            M = -(M - a*b.transpose());
            //matrix is symmetric
            gmc.block(j,i,n2,n1)=M.transpose().eval();
        }
        
        imcio_write_dS(grp_name + suffix + ".imc", r, dS);
        imcio_write_matrix(grp_name + suffix + ".gmc", gmc);
        imcio_write_index(grp_name + suffix + ".idx", names, ranges);

    }
}

// calculate deviation from target vectors
void Imc::CalcDeltaS(interaction_t *interaction, Eigen::VectorBlock< Eigen::VectorXd > &dS)
{
    const string &name = interaction->_p->get("name").as<string>();
                
    Table target;
    target.Load(name + ".dist.tgt");

    if(!interaction->_is_bonded) {
        for(unsigned int i=0; i<target.y().size(); ++i) {
            double x1 = target.x()[i] - 0.5*interaction->_step;
            double x2 = x1 + interaction->_step;
            if(x1<0)
                x1=x2=0;
            target.y()[i] = 1./(_avg_vol.getAvg()*interaction->_norm) *
                target.y()[i] * (4./3.*M_PI*(x2*x2*x2 - x1*x1*x1));                
        }        
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
        int n=grp->_corr.rows();

        // build full set of equations + copy some data to make
        // code better to read
        group_matrix gmc(grp->_corr);
        Eigen::VectorXd dS(n);
        Eigen::VectorXd r(n);
        // the next two variables are to later extract the individual parts
        // from the whole data after solving equations
        vector<int> sizes; // sizes of the individual interactions
        vector<string> names; // names of the interactions

        // copy all averages+r of group to one vector
        n=0;
        for(iter=grp->_interactions.begin(); iter != grp->_interactions.end(); ++iter) {
            interaction_t *ic = *iter;

            // sub vector for dS
            Eigen::VectorBlock< Eigen::VectorXd > sub_dS=dS.segment(n,ic->_average.getNBins());
            
            // sub vector for r
           Eigen::VectorBlock< Eigen::VectorXd > sub_r=r.segment(n,ic->_average.getNBins());

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

        for(int i=0; i<grp->_corr.rows(); ++i) {
            for(int j=0; j<grp->_corr.cols(); ++j) {
                out_cor << grp->_corr(i, j) << " ";
            }
            out_cor << endl;
        }
        out_cor.close();
        cout << "written " << name_cor << endl;
    }
}

CsgApplication::Worker *Imc::ForkWorker()
{
    Imc::Worker *worker;
    worker = new Imc::Worker;
    map<string, interaction_t *>::iterator ic_iter;

    worker->_current_hists.resize(_interactions.size());
    worker->_current_hists_force.resize(_interactions.size());
    worker->_imc = this;

    for (ic_iter = _interactions.begin(); ic_iter != _interactions.end(); ++ic_iter) {
        interaction_t *i = ic_iter->second;
        worker->_current_hists[i->_index].Initialize(
        i->_average.getMin(),
        i->_average.getMax(),
        i->_average.getNBins());
        //preliminary
        if (ic_iter->second->_force){
            worker->_current_hists_force[i->_index].Initialize(
            i->_average_force.getMin(),
            i->_average_force.getMax(),
            i->_average_force.getNBins());          
        }
    }
    return worker;
}

void Imc::MergeWorker(CsgApplication::Worker* worker_)
{
    _processed_some_frames = true;
    Imc::Worker *worker = dynamic_cast<Imc::Worker *>(worker_);
        // update the average
    map<string, interaction_t *>::iterator ic_iter;
    //map<string, group_t *>::iterator group_iter;

    ++_nframes;
    _avg_vol.Process(worker->_cur_vol);
    for (ic_iter = _interactions.begin(); ic_iter != _interactions.end(); ++ic_iter) {
        interaction_t *i=ic_iter->second;
        i->_average.data().y() = (((double)_nframes-1.0)*i->_average.data().y()
            + worker->_current_hists[i->_index].data().y())/(double)_nframes;
        //preliminary
        if (i->_force){
            i->_average_force.data().y() = (((double)_nframes-1.0)*i->_average_force.data().y()
                + worker->_current_hists_force[i->_index].data().y())/(double)_nframes;
        }       
    }

        // update correlation matrices
    if(_do_imc)
        DoCorrelations(worker);

    if(_block_length != 0) {
        if((_nframes % _block_length)==0) {
            _nblock++;
            string suffix = string("_") + boost::lexical_cast<string>(_nblock) + string(".") + _extension;
            WriteDist(suffix);
            WriteIMCData(suffix);
            WriteIMCBlock(suffix);
            ClearAverages();
        }
    }

}

}}
