#ifndef VOTCA_CTP_JOBWRITER_H
#define VOTCA_CTP_JOBWRITER_H

#include<votca/ctp/topology.h>
#include<votca/ctp/qmcalculator.h>


namespace votca { namespace ctp {
    
  
class JobWriter : public QMCalculator
{

public:

    typedef void (JobWriter::*WriteFunct)(Topology*);
    
    string Identify() { return "jobwriter"; }
    void Initialize(Topology *top, Property *options);
    bool EvaluateFrame(Topology *top);    
    
    // NEED TO REGISTER ALL WRITE MEMBERS IN ::Initialize
    void xqmultipole_ct(Topology *top);
    void xqmultipole_chrg(Topology *top);
    void xqmultipole_kmc(Topology *top);
    
    void edft(Topology *top);
    void idft(Topology *top);
    

private:

    Property *_options;
    vector<string> _keys;
    map<string,WriteFunct> _key_funct;
};




    
    
    
}}

#endif