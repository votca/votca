#ifndef SANDBOX2_H
#define SANDBOX2_H

#include <votca/ctp/qmcalculator2.h>
#include <votca/tools/property.h>

namespace TOOLS = votca::tools;
namespace votca { namespace ctp {

class Sandbox2 : public QMCalculator2
{
public:
    Sandbox2() { };
   ~Sandbox2() { };

    string  Identify() { return "Sandbox"; }
    void    Initialize(Topology *top, TOOLS::Property *opt);
    bool    EvaluateFrame(Topology *top);

private:
    int     _ID;
    double  _p1;
    double  _p2;

};

}} /* exit namespace votca::ctp */

#endif /* SANDBOX2_H */
