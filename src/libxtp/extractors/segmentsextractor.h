#ifndef VOTCA_XTP_SEGMENTSEXTRACTOR_H
#define VOTCA_XTP_SEGMENTSEXTRACTOR_H

#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/qmcalculator.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {


class SegmentsExtractor : public QMCalculator
{
public:

    SegmentsExtractor() { };
   ~SegmentsExtractor() { };

    string Identify() { return "extract.segments"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

};


void SegmentsExtractor::Initialize(Property *options) {
    return;
}


bool SegmentsExtractor::EvaluateFrame(Topology *top) {
    
    // Rigidify system (if possible)
    if (!top->Rigidify()) return 0;
    
    string xmlfile = Identify() + ".xml";    
    
    Property state("state", "", "");
    Property &segs = state.add("segments","");
    Property *next = NULL;    
    
    using boost::format;
    
    // SEGMENTS
    vector<Segment*> ::iterator sit;
    next = &segs;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
        Segment *seg = *sit;
        Property &segprop = next->add("segment", "");
        segprop.add("id", (format("%1$d") % seg->getId()).str());
        segprop.add("name", seg->getName());
        segprop.add("molecule", seg->getMolecule()->getName());
        segprop.add("xyz", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
            % seg->getPos().getX() % seg->getPos().getY() % seg->getPos().getZ()).str());
        if (seg->hasChrgState(+1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "hole");
            channel.add("energy_hn", (format("%1$+1.7f") % seg->getSiteEnergy(+1)).str());
            channel.add("lambda_hn", (format("%1$+1.7f") % seg->getU_nC_nN(+1)).str());
            channel.add("lambda_nh", (format("%1$+1.7f") % seg->getU_cN_cC(+1)).str());
            channel.add("occupation_h", (format("%1$+1.7e") % seg->getOcc(+1)).str());
        }
        if (seg->hasChrgState(-1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "electron");
            channel.add("energy_en", (format("%1$+1.7f") % seg->getSiteEnergy(-1)).str());
            channel.add("lambda_en", (format("%1$+1.7f") % seg->getU_nC_nN(-1)).str());
            channel.add("lambda_ne", (format("%1$+1.7f") % seg->getU_cN_cC(-1)).str());
            channel.add("occupation_e", (format("%1$+1.7e") % seg->getOcc(-1)).str());
        }
        
        if (tools::globals::verbose) {
            
            // FRAGMENTS
            Property &fragprop = segprop.add("fragment", "");
            vector< Fragment* > ::iterator fit;
            for (fit = seg->Fragments().begin(); fit < seg->Fragments().end(); ++fit) {
                Fragment *frag = *fit;
                fragprop.add("id", (format("%1$d") % frag->getId()).str());
                fragprop.add("name", frag->getName());
                fragprop.add("xyz", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
                        % frag->getPos().getX() % frag->getPos().getY() % frag->getPos().getZ()).str());                
                
                // ATOMS
                vector<Atom*>::iterator ait;
                for (ait = frag->Atoms().begin(); ait < frag->Atoms().end(); ++ait) {
                    Property &atomprop = fragprop.add("atom", "");
                    Atom *atm = *ait;
                    atomprop.add("id", (format("%1$d") % atm->getId()).str());
                    atomprop.add("element", atm->getElement());
                    atomprop.add("name", atm->getName());
                    atomprop.add("weight", (format("%1$1.2f") % atm->getWeight()).str());
                    atomprop.add("pos", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
                        % atm->getPos().getX() % atm->getPos().getY() % atm->getPos().getZ()).str());
                    if (atm->HasQMPart()) {
                        atomprop.add("qmid", (format("%1$d") % atm->getQMId()).str());
                        atomprop.add("qmpos", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
                            % atm->getQMPos().getX() % atm->getQMPos().getY() % atm->getQMPos().getZ()).str());
                    }
                }
            }
        }
    }
    
    ofstream ofs;    
    ofs.open(xmlfile.c_str(), ofstream::out);
    if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + xmlfile);
    }
    ofs << tools::XML << state;
    ofs.close();
    
    return true;
}


}}

#endif // VOTCA_XTP_SEGMENTSEXTRACTOR_H
