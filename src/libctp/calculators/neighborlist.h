/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#ifndef __NEIGHBORLIST2_H
#define __NEIGHBORLIST2_H

#include <votca/tools/globals.h>
#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/qmpair.h>

namespace TOOLS = votca::tools;

namespace votca { namespace ctp {


class Neighborlist : public QMCalculator
{
public:

    Neighborlist() { };
   ~Neighborlist() { };

    string Identify() { return "Neighborlist"; }
    
    void Initialize(Topology *top, Property *options);
    bool EvaluateFrame(Topology *top);
    void GenerateFromFile(Topology *top, string filename);

private:

    map< string, map<string,double> > _cutoffs;
    bool                              _useConstantCutoff;
    double                            _constantCutoff;
    string                            _generate_from;
    bool                              _generate_from_file;

};
    

void Neighborlist::Initialize(Topology* top, Property *options) {

    string key = "options.neighborlist";

    list< Property* > segs = options->Select(key+".segments");
    list< Property* > ::iterator segsIt;

    for (segsIt = segs.begin();
         segsIt != segs.end();
         segsIt++) {

        string types = (*segsIt)->get("type").as<string>();
        double cutoff = (*segsIt)->get("cutoff").as<double>();

        Tokenizer tok(types, " ");
        vector< string > names;
        tok.ToVector(names);

        if (names.size() != 2) {
            cout << "ERROR: Faulty pair definition for cut-off's: "
                 << "Need two segment names separated by a ' '" << endl;
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
    if (options->exists(key+".generate_from")) {
        _generate_from_file = true;
        _generate_from = options->get(key+".generate_from").as< string >();
    }
    else {
        _generate_from_file = false;
        _generate_from = "nofile";
    }
    
}

bool Neighborlist::EvaluateFrame(Topology *top) {

    top->NBList().Cleanup();

    if (_generate_from_file) {        
        this->GenerateFromFile(top, _generate_from);        
    }
    
    else {        

        vector< Segment* > ::iterator segit1;
        vector< Segment* > ::iterator segit2;
        vector< Fragment* > ::iterator fragit1;
        vector< Fragment* > ::iterator fragit2;

        double cutoff;
        vec r1;
        vec r2;

        for (segit1 = top->Segments().begin();
                segit1 < top->Segments().end();
                segit1++) {

                Segment *seg1 = *segit1;

            for (segit2 = segit1 + 1;
                    segit2 < top->Segments().end();
                    segit2++) {

                Segment *seg2 = *segit2;

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
                            top->NBList().Add(seg1, seg2);
                            stopLoop = true;
                            break;
                        }                

                    } /* exit loop frag2 */
                } /* exit loop frag1 */
            } /* exit loop seg2 */
        } /* exit loop seg1 */       

    }
    
    cout << endl << "... ... Created " << top->NBList().size() << " pairs.";

    if (TOOLS::globals::verbose) {
        cout << "[idA:idB] com distance" << endl;
        QMNBList& nblist = top->NBList();
        for (QMNBList::iterator ipair = nblist.begin(); ipair != nblist.end(); ++ipair) {
                QMPair *pair = *ipair;
                Segment* segment1 = pair->Seg1PbCopy();
                Segment* segment2 = pair->Seg2PbCopy();
                cout << " [" << segment1->getId() << ":" << segment2->getId()<< "] " << pair->Dist()<< endl;
        }
    }
    
    return true;        
}


void Neighborlist::GenerateFromFile(Topology *top, string filename) {
    
    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }
            
            int seg1id          = boost::lexical_cast<int>(split[1]);
            int seg2id          = boost::lexical_cast<int>(split[2]);
            string seg1name     = boost::lexical_cast<string>(split[7]);
            string seg2name     = boost::lexical_cast<string>(split[8]);
            double j2           = boost::lexical_cast<double>(split[3]);
            
            Segment* seg1 = top->getSegment(seg1id);
            Segment* seg2 = top->getSegment(seg2id);
            
            
            assert(seg1->getName() == seg1name);
            assert(seg2->getName() == seg2name);
            
            QMPair* pair12 = top->NBList().Add(seg1,seg2);
            
    
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
