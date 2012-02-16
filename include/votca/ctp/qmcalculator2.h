#ifndef _QMCALCULATOR2_H
#define _QMCALCULATOR2_H


#include <votca/tools/property.h>
#include <votca/ctp/topology.h>

namespace CTP = votca::ctp;

namespace votca { namespace ctp {

class QMCalculator2
{
public:

                    QMCalculator2() { };
    virtual        ~QMCalculator2() { };

    virtual string  Identify() { return "Generic calculator"; }

    virtual void    Initialize(CTP::Topology *top, Property *options) { }
    virtual bool    EvaluateFrame(CTP::Topology *top) { return true; }
    virtual void    EndEvaluate(CTP::Topology *top) { }

    void            setnThreads(int nThreads) { _nThreads = nThreads; }

protected:

    int _nThreads;

};

}}

#endif /* _QMCALCULATOR2_H */