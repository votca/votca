#ifndef VOTCA_XTP_PAIRSEXTRACTOR_H
#define VOTCA_XTP_PAIRSEXTRACTOR_H

#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/qmcalculator.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {


class PairsExtractor : public QMCalculator
{
public:

    PairsExtractor() { };
   ~PairsExtractor() { };

    string Identify() { return "extract.pairs"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

};


void PairsExtractor::Initialize(Property *options) {
    return;
}


bool PairsExtractor::EvaluateFrame(Topology *top) {
    
    // Rigidify std::system (if possible)
    if (!top->Rigidify()) return 0;
    
    string xmlfile = Identify() + ".xml";    
    
    Property state("state", "", "");
    Property &pairs = state.add("pairs","");    
    
    using boost::format;
    
    // PAIRS
    QMNBList::iterator pit;
    QMNBList &nb = top->NBList();
    for (pit = nb.begin(); pit != nb.end(); ++pit) {
        QMPair *qmp = *pit;
        Property &pairprop = pairs.add("pair", "");
        pairprop.add("id1",   (format("%1$d")   % qmp->Seg1()->getId()).str());
        pairprop.add("name1", (format("%1$s")   % qmp->Seg1()->getName()).str());
        pairprop.add("id2",   (format("%1$d")   % qmp->Seg2()->getId()).str());
        pairprop.add("name2", (format("%1$s")   % qmp->Seg2()->getName()).str());
        
        pairprop.add("pairvector", (format("%1$+1.4f %2$+1.4f %3$+1.4f")
            % qmp->R().getX() % qmp->R().getY() % qmp->R().getZ()).str());

        if (qmp->isPathCarrier(+1)) {
            Property &channel = pairprop.add("channel", "");
            channel.setAttribute("type","hole");
            channel.add("jeff2_h", (format("%1$1.7e") % qmp->getJeff2(+1)).str());
            channel.add("deltaE_h12", (format("%1$1.7f") % qmp->getdE12(+1)).str());
            channel.add("lambdaI_h12", (format("%1$1.7f") % qmp->getReorg12(+1)).str());
            channel.add("lambdaI_h21", (format("%1$1.7f") % qmp->getReorg21(+1)).str());
            channel.add("lambdaO_h", (format("%1$1.7f") % qmp->getLambdaO(+1)).str());
            channel.add("rate_h12", (format("%1$1.7e") % qmp->getRate12(+1)).str());
            channel.add("rate_h21", (format("%1$1.7e") % qmp->getRate12(+1)).str());
        }
        if (qmp->isPathCarrier(-1)) {
            Property &channel = pairprop.add("channel", "");
            channel.setAttribute("type","electron");
            channel.add("jeff2_e", (format("%1$1.7e") % qmp->getJeff2(-1)).str());
            channel.add("deltaE_e12", (format("%1$1.7f") % qmp->getdE12(-1)).str());
            channel.add("lambdaI_e12", (format("%1$1.7f") % qmp->getReorg12(-1)).str());
            channel.add("lambdaI_e21", (format("%1$1.7f") % qmp->getReorg21(-1)).str());
            channel.add("lambdaO_e", (format("%1$1.7f") % qmp->getLambdaO(-1)).str());
            channel.add("rate_e12", (format("%1$1.7e") % qmp->getRate12(-1)).str());
            channel.add("rate_e21", (format("%1$1.7e") % qmp->getRate12(-1)).str());
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

#endif // VOTCA_XTP_PAIRSEXTRACTOR_H
