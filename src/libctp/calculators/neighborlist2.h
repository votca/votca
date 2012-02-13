#ifndef __NEIGHBORLIST2_H
#define __NEIGHBORLIST2_H

#include <votca/tools/globals.h>
#include <votca/ctp/qmcalculator2.h>
#include <votca/ctp/qmpair2.h>

namespace TOOLS = votca::tools;

namespace votca { namespace ctp {


class Neighborlist2 : public QMCalculator2
{
public:

    Neighborlist2() { };
   ~Neighborlist2() { };

    string Identify() { return "Neighborlist"; }
    
    void Initialize(Topology *top, Property *options);
    bool EvaluateFrame(Topology *top);

private:

    map< string, map<string,double> > _cutoffs;

};
    

void Neighborlist2::Initialize(Topology* top, Property *options) {

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

bool Neighborlist2::EvaluateFrame(Topology *top) {

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
        QMNBList2& nblist = top->NBList();       
        for (QMNBList2::iterator ipair = nblist.begin(); ipair != nblist.end(); ++ipair) {
                QMPair2 *pair = *ipair;
                Segment* segment1 = pair->Seg1PbCopy();
                Segment* segment2 = pair->Seg2PbCopy();
                cout << " [" << segment1->getId() << ":" << segment2->getId()<< "] " << pair->Dist()<< endl;
        }
    }

    
}




}}

#endif  /* __NEIGHBORLIST2_H */