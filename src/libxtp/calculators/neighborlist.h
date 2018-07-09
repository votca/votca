/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/qmnblist.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/atom.h>
#include <votca/tools/property.h>
#include <boost/progress.hpp>
#include <boost/format.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace votca { namespace xtp {

class Neighborlist : public xtp::QMCalculator
{

public:

    Neighborlist() { };
   ~Neighborlist() {
       // cleanup the list of superexchange pair types
       for ( std::list<xtp::QMNBList::SuperExchangeType*>::iterator it = _superexchange.begin() ; it != _superexchange.end(); it++  ) {
           delete *it;
       }
    };

     std::string Identify() { return "neighborlist"; }
    
    void Initialize(tools::Property *options);
    bool EvaluateFrame(xtp::Topology *top);
    void GenerateFromFile(xtp::Topology *top,  std::string filename);
  

private:

    
    std::vector<std::string> _included_segments;
    std::map<  std::string,  std::map< std::string,double> > _cutoffs;
    bool                              _useConstantCutoff;
    double                            _constantCutoff;
    bool                              _useExcitonCutoff;
    double                            _excitonqmCutoff;
     std::string                            _generate_from;
    bool                              _generate_from_file;
    bool                              _generate_unsafe;
    bool                              _do_bridging;
    std::list<xtp::QMNBList::SuperExchangeType*>        _superexchange;

};
    

void Neighborlist::Initialize(tools::Property *options) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();
    
      
     std::list< tools::Property* > segs = options->Select(key+".segments");
     std::list< tools::Property* > ::iterator segsIt;

    for (segsIt = segs.begin();
         segsIt != segs.end();
         segsIt++) {

         std::string types = (*segsIt)->get("type").as< std::string>();
        double cutoff = (*segsIt)->get("cutoff").as<double>();

        tools::Tokenizer tok(types, " ");
         std::vector<  std::string > names;
        tok.ToVector(names);
        

        if (names.size() != 2) {
            std::cout << "ERROR: Faulty pair definition for cut-off's: "
                 << "Need two segment names separated by a space" << std::endl;
            throw std::runtime_error("Error in options file.");
        }

        _cutoffs[names[0]][names[1]] = cutoff;
        _cutoffs[names[1]][names[0]] = cutoff;
        if(std::find(_included_segments.begin(), _included_segments.end(), names[0]) == _included_segments.end()){
                             _included_segments.push_back(names[0]);
                        }
        if(std::find(_included_segments.begin(), _included_segments.end(), names[1]) == _included_segments.end()){
                             _included_segments.push_back(names[1]);
                        }
     

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
        _generate_from = options->get(key+".generate_from").as<  std::string >();
    }
    else {
        _generate_from_file = false;
        _generate_from = "nofile";
    }
    
  
    
    // if superexchange is given
    if (options->exists(key + ".superexchange")) {
        _do_bridging = true;
    
    
    
         std::list< tools::Property* > _se = options->Select(key + ".superexchange");
         std::list< tools::Property* > ::iterator seIt;

        for (seIt = _se.begin(); seIt != _se.end(); seIt++) {
             std::string types = (*seIt)->get("type").as< std::string>();
            xtp::QMNBList::SuperExchangeType* _su = new xtp::QMNBList::SuperExchangeType(types);
            _superexchange.push_back(_su); 
        }
    }
    else{
       _do_bridging=false; 
    }
            
}

bool Neighborlist::EvaluateFrame(xtp::Topology *top) {
  

    top->NBList().Cleanup();

    if (_generate_from_file) {        
        GenerateFromFile(top, _generate_from);        
    }    
    else {        

       
        if (tools::globals::verbose) {
            std::cout << std::endl <<  "... ..." << std::flush;
        }
        
        const tools::matrix box=top->getBox();
        double min=box.get(0,0);
        if(min>box.get(1,1)){min=box.get(1,1);}
        if(min>box.get(2,2)){min=box.get(2,2);}
    
        std::vector< xtp::Segment* > segs;
    
        for (std::vector< xtp::Segment* > ::iterator segit1 = top->Segments().begin();              
                    segit1 < top->Segments().end();++segit1) {
            xtp::Segment *seg1 = *segit1;
            if(_useConstantCutoff || std::find(_included_segments.begin(), _included_segments.end(), seg1->getName()) != _included_segments.end()){
                segs.push_back(seg1);
                seg1->calcPos();
                seg1->calcApproxSize();
            }    
        }
        std::cout<<std::endl;
        std::cout <<"Evaluating "<<segs.size()<<" segments for neighborlist. "<< top->Segments().size()-segs.size()<<" segments are not taken into account as specified"<< std::endl;
        if(!_useConstantCutoff){
        std::cout << "The following segments are used in the neigborlist creation"<<std::endl;      
        std::cout<<"\t"<<std::flush;        
        for(std::vector< std::string >::iterator st=_included_segments.begin();st!=_included_segments.end();++st){
            std::cout<<" "<<(*st)<<std::flush;
        }
        std::cout <<std::endl;
        }
        
        std::cout << "\r ... ... Evaluating " <<std::flush; 
        std::vector<std::string> skippedpairs;
       
        for (std::vector< xtp::Segment* > ::iterator segit1 = segs.begin();segit1 < segs.end();++segit1) {
                xtp::Segment *seg1 = *segit1;
                
                std::vector< xtp::Segment* > ::iterator segit2;
                std::vector< xtp::Fragment* > ::iterator fragit1;
                std::vector< xtp::Fragment* > ::iterator fragit2;
                double cutoff;
                tools::vec r1;
                tools::vec r2;

                
                if (tools::globals::verbose) {
                    std::cout << "\r ... ... NB List Seg " << seg1->getId() << std::flush;
                    
                }

            for (segit2 = segit1 + 1;segit2 < segs.end();++segit2) {

                xtp::Segment *seg2 = *segit2;

                if (!_useConstantCutoff) {
                    // Find cut-off
                    try {
                        cutoff = _cutoffs.at(seg1->getName())
                                         .at(seg2->getName());
                    }
                    catch (const std::exception& out_of_range) {
                        std::string pairstring=seg1->getName()+"/"+seg2->getName();
                        if(std::find(skippedpairs.begin(), skippedpairs.end(), pairstring) == skippedpairs.end()){
                            skippedpairs.push_back(pairstring);
                        }
                        continue;
                    }
                }

                else { cutoff = _constantCutoff; }

                if (cutoff>0.5*min){             
                    throw std::runtime_error((boost::format("Cutoff is larger than half the box size. Maximum allowed cutoff is %1$1.1f") % (0.5*min)).str());
                }
                double cutoff2=cutoff*cutoff;
                tools::vec segdistance=top->PbShortestConnect(seg1->getPos(),seg2->getPos());
                double segdistance2=segdistance*segdistance;
                double outside=cutoff+seg1->getApproxSize()+seg2->getApproxSize();
                
                if(segdistance2<cutoff2){
                    top->NBList().Add(seg1, seg2);
                }
                else if(segdistance2>(outside*outside)){
                    continue;
                }
                else {
                    bool stopLoop = false;
                    for (fragit1 = seg1->Fragments().begin();
                            fragit1 < seg1->Fragments().end();
                            fragit1++) {

                        if (stopLoop) {
                            break; }
                        r1 = (*fragit1)->getPos();
                        for (fragit2 = seg2->Fragments().begin();
                                fragit2 < seg2->Fragments().end();
                                fragit2++) {



                            r2 = (*fragit2)->getPos();
                                    tools::vec distance = top->PbShortestConnect(r1, r2);
                                    double dist2 = distance*distance;
                            if (dist2 > cutoff2) {
                                continue;
                            } else {
                                top->NBList().Add(seg1, seg2);
                                        stopLoop = true;
                                break;
                            }

                        } /* exit loop frag2 */
                    } /* exit loop frag1 */
                }
                
                
            } /* exit loop seg2 */
                
               
           
        } /* exit loop seg1 */   
        
       
        
        
        if(skippedpairs.size()>0){
        std::cout << "WARNING: No cut-off specified for segment pairs of type "<<std::endl;              
        for(std::vector< std::string >::iterator st=skippedpairs.begin();st!=skippedpairs.end();++st){
            std::cout<<(*st)<<std::endl;
        }
        std::cout << "pairs were skipped"<<std::endl;
        }
        
    
    std::cout << std::endl << " ... ... Created " << top->NBList().size() << " direct pairs.";
    if(_useExcitonCutoff){
        std::cout << std::endl << " ... ... Determining classical pairs "<<std::endl;
        xtp::QMNBList &nblist = top->NBList();
        for (xtp::QMNBList::iterator pit = nblist.begin(); pit != nblist.end(); ++pit) {
            tools::vec r1;
            tools::vec r2;
             std::vector< xtp::Fragment* > ::iterator fragit1;
             std::vector< xtp::Fragment* > ::iterator fragit2;
            
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

	tools::Property bridges_summary;
        tools::Property *_bridges = &bridges_summary.add("bridges","");

        std::cout << "Bridged Pairs \n [idA:idB] com distance" << std::endl;
        for (xtp::QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {
                xtp::QMPair *pair = *ipair;
                xtp::Segment* segment1 = pair->Seg1PbCopy();
                xtp::Segment* segment2 = pair->Seg2PbCopy();
                
                std::cout << " [" << segment1->getId() << ":" << segment2->getId()<< "] " 
                             << pair->Dist()<< " bridges: " 
                             << (pair->getBridgingSegments()).size() 
                             << " type: " 
                             << pair->getType() 
                             << " | " << std::flush;
                
                std::vector<xtp::Segment*> bsegments = pair->getBridgingSegments();
 
                tools::Property *_pair_property = &_bridges->add("pair","");
                                   
                _pair_property->setAttribute("id1", segment1->getId());
                _pair_property->setAttribute("id2", segment2->getId());
                _pair_property->setAttribute("name1", segment1->getName());
                _pair_property->setAttribute("name2", segment2->getName());
                _pair_property->setAttribute("r12", pair->Dist());
                                    
                tools::Property *_bridge_property = &_pair_property->add("bridge","");

                for (  std::vector<xtp::Segment*>::iterator itb = bsegments.begin(); itb != bsegments.end(); itb++ ) {
                    std::cout << (*itb)->getId() << " " ;
                    _bridge_property->setAttribute("id", (*itb)->getId());
                }        
                
                std::cout << std::endl;
        }
        //std::cout << bridges_summary;
    }
    }
    }
    return true;        
}




void Neighborlist::GenerateFromFile(xtp::Topology *top,  std::string filename) {
    
    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());
    
    if (_generate_unsafe) {
        std::cout << std::endl << "... ... Generate unsafe = true ..." << std::flush;
    }

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
             std::vector< std::string> split;
            tools::Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }
            
            int seg1id          = boost::lexical_cast<int>(split[1]);
            int seg2id          = boost::lexical_cast<int>(split[2]);
            
            xtp::Segment* seg1 = top->getSegment(seg1id);
            xtp::Segment* seg2 = top->getSegment(seg2id);
            
            if (not _generate_unsafe) {
                 std::string seg1name     = boost::lexical_cast<std::string>(split[7]);
                 std::string seg2name     = boost::lexical_cast<std::string>(split[8]);
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
    else { std::cout << std::endl << "ERROR: No such file " << filename << std::endl;
           throw std::runtime_error("Supply input file."); }
    
}



}}

#endif  /* __NEIGHBORLIST2_H */
