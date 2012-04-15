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

private:

    map< string, map<string,double> > _cutoffs;

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
}

bool Neighborlist::EvaluateFrame(Topology *top) {

    top->NBList().Cleanup();

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

            // Find cut-off
            try {
                cutoff = _cutoffs.at(seg1->getName()).at(seg2->getName());
            }
            catch (out_of_range) {
                cout << "ERROR: No cut-off specified for segment pair "
                     << seg1->getName() << " | " << seg2->getName() << ". "
                     << endl;
                throw std::runtime_error("Missing input in options file.");
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
                        top->NBList().Add(seg1, seg2);
                        stopLoop = true;
                        break;
                    }                

                } /* exit loop frag2 */
            } /* exit loop frag1 */
        } /* exit loop seg2 */
    } /* exit loop seg1 */

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

    
}




}}

#endif  /* __NEIGHBORLIST2_H */
