#ifndef VOTCA_XTP_ENERGYEXTRACTOR_H
#define VOTCA_XTP_ENERGYEXTRACTOR_H

#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/qmcalculator.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {


class EnergyExtractor : public QMCalculator
{
public:

    EnergyExtractor() { };
   ~EnergyExtractor() { };

    string Identify() { return "extract.energy"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

};


void EnergyExtractor::Initialize(Property *options) {
    return;
}


bool EnergyExtractor::EvaluateFrame(Topology *top) {
    
    string xmlfile = Identify() + ".xml";    
    
    Property state("state", "", "");
    Property &segs = state.add("segments","");
    Property &pairs = state.add("pairs","");
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
        if (seg->hasChrgState(+1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "hole");
            channel.add("energy_hn", (format("%1$+1.7f") % seg->getSiteEnergy(+1)).str());
            channel.add("lambda_hn", (format("%1$+1.7f") % seg->getU_nC_nN(+1)).str());
            channel.add("lambda_nh", (format("%1$+1.7f") % seg->getU_cN_cC(+1)).str());
        }
        if (seg->hasChrgState(-1)) {
            Property &channel = segprop.add("channel","");
            channel.setAttribute("type", "electron");
            channel.add("energy_en", (format("%1$+1.7f") % seg->getSiteEnergy(-1)).str());
            channel.add("lambda_en", (format("%1$+1.7f") % seg->getU_nC_nN(-1)).str());
            channel.add("lambda_ne", (format("%1$+1.7f") % seg->getU_cN_cC(-1)).str());
        }
    }
    
    if (tools::globals::verbose) {
        // PAIRS
        QMNBList::iterator pit;
        QMNBList &nb = top->NBList();
        next = &pairs;
        for (pit = nb.begin(); pit != nb.end(); ++pit) {
            QMPair *qmp = *pit;
            Property &pairprop = next->add("pair", "");
            pairprop.add("id1",   (format("%1$d")   % qmp->Seg1()->getId()).str());
            pairprop.add("name1", (format("%1$s")   % qmp->Seg1()->getName()).str());
            pairprop.add("id2",   (format("%1$d")   % qmp->Seg2()->getId()).str());
            pairprop.add("name2", (format("%1$s")   % qmp->Seg2()->getName()).str());

            if (qmp->isPathCarrier(+1)) {
                Property &channel = pairprop.add("channel", "");
                channel.setAttribute("type","hole");
                channel.add("deltaE_h12", (format("%1$1.7f") % qmp->getdE12(+1)).str());
                channel.add("lambdaI_h12", (format("%1$1.7f") % qmp->getReorg12(+1)).str());
                channel.add("lambdaI_h21", (format("%1$1.7f") % qmp->getReorg21(+1)).str());
                channel.add("lambdaO_h", (format("%1$1.7f") % qmp->getLambdaO(+1)).str());
            }
            if (qmp->isPathCarrier(-1)) {
                Property &channel = pairprop.add("channel", "");
                channel.setAttribute("type","electron");
                channel.add("deltaE_e12", (format("%1$1.7f") % qmp->getdE12(-1)).str());
                channel.add("lambdaI_e12", (format("%1$1.7f") % qmp->getReorg12(-1)).str());
                channel.add("lambdaI_e21", (format("%1$1.7f") % qmp->getReorg21(-1)).str());
                channel.add("lambdaO_e", (format("%1$1.7f") % qmp->getLambdaO(-1)).str());
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

#endif // VOTCA_XTP_ENERGYEXTRACTOR_H
