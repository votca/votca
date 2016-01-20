#ifndef VOTCA_XTP_OCCUPATIONSEXTRACTOR_H
#define VOTCA_XTP_OCCUPATIONSEXTRACTOR_H

#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/qmcalculator.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {


class OccupationsExtractor : public QMCalculator
{
public:

    OccupationsExtractor() { };
   ~OccupationsExtractor() { };

    string Identify() { return "extract.occupations"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

};


void OccupationsExtractor::Initialize(Property *options) {
    return;
}


bool OccupationsExtractor::EvaluateFrame(Topology *top) {
    
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
        segprop.add("xyz", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
            % seg->getPos().getX() % seg->getPos().getY() % seg->getPos().getZ()).str());
        if (seg->hasChrgState(+1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "hole");
            channel.add("occupation_h", (format("%1$+1.7e") % seg->getOcc(+1)).str());
        }
        if (seg->hasChrgState(-1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "electron");
            channel.add("occupation_e", (format("%1$+1.7e") % seg->getOcc(-1)).str());
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

#endif // VOTCA_XTP_OCCUPATIONSEXTRACTOR_H
