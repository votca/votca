/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#ifndef __VOTCA_XTP_NEIGHBORLIST_H
#define __VOTCA_XTP_NEIGHBORLIST_H

#include <votca/tools/globals.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/qmpair.h>
#include <boost/progress.hpp>
#include <boost/format.hpp>
#include <votca/tools/random2.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace votca { namespace xtp {

namespace TOOLS = votca::tools;
namespace CTP = votca::ctp;

using namespace std;

class XNeighborlist : public XQMCalculator
{

public:

    XNeighborlist() { };
   ~XNeighborlist() {
       // cleanup the list of superexchange pair types
       for ( std::list<CTP::QMNBList::SuperExchangeType*>::iterator it = _superexchange.begin() ; it != _superexchange.end(); it++  ) {
           delete *it;
       }
    };

    std::string Identify() { return "xneighborlist"; }
    
    void Initialize(CTP::Property *options);
    bool EvaluateFrame(CTP::Topology *top);
    void GenerateFromFile(CTP::Topology *top, std::string filename);
    void StochasticConnectivity(CTP::Topology *top, std::string filename);
    bool StochasticConnectOrNot(double thisdistance, std::vector<double> distances, std::vector<double> probabilities, votca::tools::Random2 *RandomVariable);

private:

    std::map< std::string, std::map<std::string,double> > _cutoffs;
    bool                              _useConstantCutoff;
    double                            _constantCutoff;
    bool                              _useExcitonCutoff;
    double                            _excitonqmCutoff;
    std::string                            _generate_from;
    bool                              _generate_from_file;
    std::string                            _probabilityfile;
    bool                              _stochastic_connectivity;
    bool                              _generate_unsafe;
    bool                              _do_bridging;
    std::list<CTP::QMNBList::SuperExchangeType*>        _superexchange;

};
    

void XNeighborlist::Initialize(CTP::Property *options) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();
    
    if (options->exists(key + ".probabilityfile")) {
        _probabilityfile = options->get(key+".probabilityfile").as< string >();
        _stochastic_connectivity = true;
        cout << endl << "... ... generating stochastic connectivity based on file " << _probabilityfile << endl;
    }
    else{
        _stochastic_connectivity = false;
    }
    
    
      
    list< CTP::Property* > segs = options->Select(key+".segments");
    list< CTP::Property* > ::iterator segsIt;

    for (segsIt = segs.begin();
         segsIt != segs.end();
         segsIt++) {

        string types = (*segsIt)->get("type").as<string>();
        double cutoff = (*segsIt)->get("cutoff").as<double>();

        CTP::Tokenizer tok(types, " ");
        vector< string > names;
        tok.ToVector(names);

        if (names.size() != 2) {
            cout << "ERROR: Faulty pair definition for cut-off's: "
                 << "Need two segment names separated by a space" << endl;
            throw std::runtime_error("Error in options file.");
        }

        _cutoffs[names[0]][names[1]] = cutoff;
        _cutoffs[names[1]][names[0]] = cutoff;

    }

    if (options->exists(key+".constant")) {
        _useConstantCutoff = true;
        _constantCutoff = options->get(key+".constant").as< double >();
    }
    else {
        _useConstantCutoff = false;
    }
    if (options->exists(key+".exciton_cutoff")) {
        _useExcitonCutoff = true;
        _excitonqmCutoff = options->get(key+".exciton_cutoff").as< double >();
    }
    else {
        _useExcitonCutoff = false;
    }
    if (options->exists(key+".generate_from")) {
        _generate_from_file = true;
        _generate_from = options->get(key+".generate_from").as< string >();
    }
    else {
        _generate_from_file = false;
        _generate_from = "nofile";
    }
    
  
    
    // if superexchange is given
    if (options->exists(key + ".superexchange")) {
        _do_bridging = true;
    
    
    
        list< CTP::Property* > _se = options->Select(key + ".superexchange");
        list< CTP::Property* > ::iterator seIt;

        for (seIt = _se.begin(); seIt != _se.end(); seIt++) {
            string types = (*seIt)->get("type").as<string>();
            CTP::QMNBList::SuperExchangeType* _su = new CTP::QMNBList::SuperExchangeType(types);
            _superexchange.push_back(_su); 
        }
    }
    else{
       _do_bridging=false; 
    }
            
}

bool XNeighborlist::EvaluateFrame(CTP::Topology *top) {
  

    top->NBList().Cleanup();

    if (_generate_from_file) {        
        this->GenerateFromFile(top, _generate_from);        
    }

    else if (_stochastic_connectivity) {        
        this->StochasticConnectivity(top, _probabilityfile);        
    }

    
    else {        

       
        
     
        
        if (TOOLS::globals::verbose) {
            cout << endl <<  "... ..." << flush;
        }
        boost::progress_display show_progress( top->Segments().size() ); 
       
        for (vector< CTP::Segment* > ::iterator segit1 = top->Segments().begin();              
                segit1 < top->Segments().end();
                segit1++) {
                
               CTP::QMNBList templist=CTP::QMNBList();
                vector< CTP::Segment* > ::iterator segit2;
                vector< CTP::Fragment* > ::iterator fragit1;
                vector< CTP::Fragment* > ::iterator fragit2;
                double cutoff;
                TOOLS::vec r1;
                TOOLS::vec r2;

                CTP::Segment *seg1 = *segit1;
                if (TOOLS::globals::verbose) {
                    cout << "\r ... ... NB List Seg " << seg1->getId() << flush;
                    ++show_progress;
                }

            for (segit2 = segit1 + 1;
                    segit2 < top->Segments().end();
                    segit2++) {

                CTP::Segment *seg2 = *segit2;

                if (!_useConstantCutoff) {
                    // Find cut-off
                    try {
                        cutoff = _cutoffs.at(seg1->getName())
                                         .at(seg2->getName());
                    }
                    catch (out_of_range) {
                        cout << "ERROR: No cut-off specified for segment pair "
                             << seg1->getName() << " | " << seg2->getName() 
                             << ". " << endl;
                        throw std::runtime_error("Missing input in options.");
                    }
                }

                else { cutoff = _constantCutoff; }

                if (cutoff>(0.5*top->getBox().get(0,0))|| cutoff>(0.5*top->getBox().get(1,1)) ||cutoff>(0.5*top->getBox().get(1,1))){
                    throw std::runtime_error("Specified cutoff is larger than half the boxlength.");
                }
                bool stopLoop = false;
                for (fragit1 = seg1->Fragments().begin();
                        fragit1 < seg1->Fragments().end();
                        fragit1 ++) {

                    if (stopLoop) { break; }

                    for (fragit2 = seg2->Fragments().begin();
                            fragit2 < seg2->Fragments().end();
                            fragit2++) {


                        r1 = (*fragit1)->getPos();
                        r2 = (*fragit2)->getPos();
                        if( abs( top->PbShortestConnect(r1, r2) ) > cutoff ) {
                            continue;
                        }
                        else {
                            seg1->calcPos();
                            seg2->calcPos();
                            
                          
                            
                            {
                            top->NBList().Add(seg1, seg2);
                            }
                            stopLoop = true;
                            break;
                        }                

                    } /* exit loop frag2 */
               } /* exit loop frag1 */
            } /* exit loop seg2 */
                
               // break;
#
            
        } /* exit loop seg1 */       

    }

    cout << endl << " ... ... Created " << top->NBList().size() << " direct pairs.";
    if(_useExcitonCutoff){
        cout << endl << " ... ... Determining classical pairs "<<endl;
        CTP::QMNBList &nblist = top->NBList();
        for (CTP::QMNBList::iterator pit = nblist.begin(); pit != nblist.end(); ++pit) {
            TOOLS::vec r1;
            TOOLS::vec r2;
            vector< CTP::Fragment* > ::iterator fragit1;
            vector< CTP::Fragment* > ::iterator fragit2;
            
            bool stopLoop = false;
                for (fragit1 =  (*pit)->Seg1()->Fragments().begin();
                        fragit1 <  (*pit)->Seg1()->Fragments().end();
                        fragit1 ++) {
                    if (stopLoop) { break; }

                    for (fragit2 =  (*pit)->Seg2()->Fragments().begin();
                            fragit2 <  (*pit)->Seg2()->Fragments().end();
                            fragit2++) {


                        r1 = (*fragit1)->getPos();
                        r2 = (*fragit2)->getPos();
                        if( abs( top->PbShortestConnect(r1, r2) ) > _excitonqmCutoff ) {
                            (*pit)->setType(3);
                            continue;
                        }
                        else {   
                            (*pit)->setType(0);
                            stopLoop = true;
                            break;
                        }                

                    } /* exit loop frag2 */
                } /* exit loop frag1 */          
            } //Type 3 Exciton_classical approx
        }        
    if (_do_bridging){
    // add superexchange pairs
    top->NBList().setSuperExchangeTypes(_superexchange);
    top->NBList().GenerateSuperExchange();
    
    // DEBUG output
    if (votca::tools::globals::verbose) {

	CTP::Property bridges_summary;
        CTP::Property *_bridges = &bridges_summary.add("bridges","");

        cout << "Bridged Pairs \n [idA:idB] com distance" << endl;
        for (CTP::QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {
                CTP::QMPair *pair = *ipair;
                CTP::Segment* segment1 = pair->Seg1PbCopy();
                CTP::Segment* segment2 = pair->Seg2PbCopy();
                
                cout << " [" << segment1->getId() << ":" << segment2->getId()<< "] " 
                             << pair->Dist()<< " bridges: " 
                             << (pair->getBridgingSegments()).size() 
                             << " type: " 
                             << pair->getType() 
                             << " | " << flush;
                
                vector<CTP::Segment*> bsegments = pair->getBridgingSegments();
 
                CTP::Property *_pair_property = &_bridges->add("pair","");
                                   
                _pair_property->setAttribute("id1", segment1->getId());
                _pair_property->setAttribute("id2", segment2->getId());
                _pair_property->setAttribute("name1", segment1->getName());
                _pair_property->setAttribute("name2", segment2->getName());
                _pair_property->setAttribute("r12", pair->Dist());
                                    
                CTP::Property *_bridge_property = &_pair_property->add("bridge","");

                for ( vector<CTP::Segment*>::iterator itb = bsegments.begin(); itb != bsegments.end(); itb++ ) {
                    cout << (*itb)->getId() << " " ;
                    _bridge_property->setAttribute("id", (*itb)->getId());
                }        
                
                cout << endl;
        }
        //cout << bridges_summary;
    }
    }
    return true;        
}

bool XNeighborlist::StochasticConnectOrNot(double thisdistance, vector<double> distances, vector<double> probabilities, votca::tools::Random2 *RandomVariable){
    double MINR = distances[0];
    double MAXR = distances[distances.size()-1];
    if(thisdistance == 0){
        // no connection with itself
        return false;
    }
    else if(thisdistance < MINR) {
        return true;
    }
    else if(thisdistance > MAXR){
        return false;
    }
    else{
        //int numberofpoints = distances.size();
        double thisprobability=0.0;
        for(unsigned i = 0; i<distances.size()-1; i++){
            if(distances[i] < thisdistance && thisdistance < distances[i+1]){
                // linear interpolations
                thisprobability = (probabilities[i+1]-probabilities[i])/(distances[i+1]-distances[i])*(thisdistance-distances[i])+probabilities[i];
                break;
            }
        }
        double uniformnumber = RandomVariable->rand_uniform();
        if(uniformnumber < thisprobability) {
            // accept connection
            return true;
        }
        else{
            // don't accept connection
            return false;
        }
    }
}

void XNeighborlist::StochasticConnectivity(CTP::Topology *top, string filename) {
    cout << endl << "... ... Creating connectivity based on the provided probability function that" << endl;
    cout << "... ... describes the centre of mass-dependent probability of two sites being" << endl;
    cout << "... ... conncected (coarse grained model). This file can be generated using the" << endl;
    cout << "... ... panalyze calculator." << endl;
    
    // read in probability function
    cout << endl << "... ... probability function:" << endl;
    vector<double> distances;
    vector<double> probabilities;
    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());
    int linenumber = 0;
    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            if (linenumber > 0){
                vector<string> split;
                CTP::Tokenizer toker(line, " \t");
                toker.ToVector(split);

                if ( !split.size()      ||
                      split[0] == "!"   ||
                      split[0].substr(0,1) == "!" ) { continue; }
             
                double distance          = boost::lexical_cast<double>(split[0]);
                double probability       = boost::lexical_cast<double>(split[1]);
                cout << "        " << distance << "    " << probability << endl;
                if(linenumber > 1){
                    if (probability == 1 && probabilities[probabilities.size()-1] == 1){
                        // first entry should be largest distance with probability 1
                        probabilities.pop_back();
                        distances.pop_back();
                    }
                }
                distances.push_back(distance);
                probabilities.push_back(probability);
            }
            linenumber++;
        }
    }
    else { cout << endl << "ERROR: No such file " << filename << endl;
           throw std::runtime_error("Supply input file."); }
    
    // Initialise random number generator
    if(votca::tools::globals::verbose) { cout << endl << "Initialising random number generator" << endl; }
    srand(12345); 
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());

    
    // get segments
    vector <CTP::Segment*> segments = top->Segments();
    vector <CTP::Segment*>::iterator seg1;
    vector <CTP::Segment*>::iterator seg2;
    
    for (seg1 = segments.begin(); seg1!= segments.end(); seg1++){
        cout << "\r... ... ..." << " evaluating segment ID = "
             << (*seg1)->getId() << flush;

        for (seg2 = seg1; seg2!= segments.end(); seg2++){ // for (seg2 = segments.begin(); seg2!= segments.end(); seg2++):q
            TOOLS::vec r1 = (*seg1)->getPos();
            TOOLS::vec r2 = (*seg2)->getPos();
            double distance = abs( top->PbShortestConnect(r1, r2));
            bool accept = StochasticConnectOrNot(distance, distances, probabilities, RandomVariable);
            if(accept == true){
                // add pair to neighbor list
                top->NBList().Add((*seg1),(*seg2));
                //QMPair* pair12 = top->NBList().Add((*seg1),(*seg2));
                //cout << "add" << endl;
            }
        }
    }

   
}


void XNeighborlist::GenerateFromFile(CTP::Topology *top, string filename) {
    
    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());
    
    if (_generate_unsafe) {
        cout << endl << "... ... Generate unsafe = true ..." << flush;
    }

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            CTP::Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }
            
            int seg1id          = boost::lexical_cast<int>(split[1]);
            int seg2id          = boost::lexical_cast<int>(split[2]);
            
            CTP::Segment* seg1 = top->getSegment(seg1id);
            CTP::Segment* seg2 = top->getSegment(seg2id);
            
            if (not _generate_unsafe) {
                string seg1name     = boost::lexical_cast<string>(split[7]);
                string seg2name     = boost::lexical_cast<string>(split[8]);
                assert(seg1->getName() == seg1name);
                assert(seg2->getName() == seg2name);
            }        
            
            //QMPair* pair12 = 
	    (void)top->NBList().Add(seg1,seg2);
            
    
     /*       
     1  1000 1010 2.4292699e-03 1.61313482160154 -1.05173043628102 0.759048038980236 DCV DCV
     2  1000 1020 1.0551418e-03 1.4977782788484 -0.466982612402543 0.876438986736797 DCV DCV
     3  1000 1023 5.9645622e-03 1.51684342052626 0.189056522949882 0.763447935684869 DCV DCV
     4  1000 1027 2.1161184e-02 -0.121730375289917 0.483095637611721 0.078926185939622 DCV DCV
     5  1000 1034 1.5198626e-03 0.586534707442574 -1.59841490776642 0.695082730832308 DCV DCV
     6  1000 1048 1.0121481e-03 -0.296308693678482 -1.02535652660805 0.347373638982358 DCV DCV
     7  1000 1050 9.3073820e-04 1.34660870303278 -1.49037826725322 0.571647867949114 DCV DCV
     8  1000 1052 1.0803526e-03 -0.337469581935717 -0.853313051455695 0.592304403885553 DCV DCV
     9  1000 1065 4.6567327e-04 0.45481307817542 -1.44727391982856 1.05151722120202 DCV DCV
    10  1000 1073 5.7739082e-03 -0.388582683646161 -0.221439142589984 0.731973764170771 DCV DCV
    */            


        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << filename << endl;
           throw std::runtime_error("Supply input file."); }
    
}



}}

#endif  /* __NEIGHBORLIST2_H */
