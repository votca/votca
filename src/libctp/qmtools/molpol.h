#ifndef _MOLPOL_H
#define _MOLPOL_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>


namespace votca { namespace ctp {


class MolPol : public QMTool
{
public:

    MolPol() { };
   ~MolPol() { };

    string Identify() { return "Molecular Polarizability Calculator"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    //matrix CalculateMolPol(vector<APolarSite*> &poles, bool verb = true);
    //int    SCF_Induce(vector<APolarSite*> &poles);


private:

    string          _mps_input;
    string          _mps_output;
    string          _pol_output;
    
    double          _aDamp;
    double          _wSOR;
    int             _maxIter;
    double          _epsTol;

    //BasicInteractor _actor;

};


}}

#endif