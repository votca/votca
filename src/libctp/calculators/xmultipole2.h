#ifndef XMULTIPOLE2_H
#define XMULTIPOLE2_H


#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

class XMP : public QMCalculator
{

    XMP() {};
   ~XMP() {};

    string          Identify() { return "XMP"; }

    void            Initialize(Topology *, Property *);
    
    void            Collect_EMP(string emp_file);
    void            Collect_XMP(string xmp_file);




    bool            EvaluateFrame(Topology *top);




};




}}

#endif

