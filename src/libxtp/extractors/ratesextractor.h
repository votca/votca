#ifndef VOTCA_XTP_RATESEXTRACTOR_H
#define VOTCA_XTP_RATESEXTRACTOR_H

#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/qmcalculator.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {


class RatesExtractor : public QMCalculator
{
public:

    RatesExtractor() { };
   ~RatesExtractor() { };

    string Identify() { return "extract.rates"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

};


void RatesExtractor::Initialize(Property *options) {
    return;
}


bool RatesExtractor::EvaluateFrame(Topology *top) {
    
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

        if (qmp->isPathCarrier(+1)) {
            Property &channel = pairprop.add("channel", "");
            channel.setAttribute("type","hole");
            channel.add("rate_h12", (format("%1$1.7e") % qmp->getRate12(+1)).str());
            channel.add("rate_h21", (format("%1$1.7e") % qmp->getRate21(+1)).str());
        }
        if (qmp->isPathCarrier(-1)) {
            Property &channel = pairprop.add("channel", "");
            channel.setAttribute("type","electron");
            channel.add("rate_e12", (format("%1$1.7e") % qmp->getRate12(-1)).str());
            channel.add("rate_e21", (format("%1$1.7e") % qmp->getRate21(-1)).str());
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

#endif // VOTCA_XTP_RATESEXTRACTOR_H
