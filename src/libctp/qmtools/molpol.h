#ifndef MOLPOLTOOL_H
#define MOLPOLTOOL_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>


namespace votca { namespace ctp {


class MolPolTool : public QMTool
{
public:

    MolPolTool() { };
   ~MolPolTool() { };

    string Identify() { return "Molecular Polarizability Calculator"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    matrix CalculateMolPol(vector<APolarSite*> &poles, bool verb = true);
    int    SCF_Induce(vector<APolarSite*> &poles);


private:

    string          _mps_input;
    string          _mps_output;
    string          _pol_output;
    
    double          _aDamp;
    double          _wSOR;
    int             _maxIter;
    double          _epsTol;

    BasicInteractor _actor;

};


}}

#endif