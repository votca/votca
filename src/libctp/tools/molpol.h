#ifndef MOLPOLTOOL_H
#define MOLPOLTOOL_H


#include <votca/ctp/qmtool.h>
#include <votca/ctp/topology.h>
#include <votca/ctp/molpolengine.h>


namespace votca { namespace ctp {


class MolPolTool : public QMTool
{
public:

    MolPolTool() { };
   ~MolPolTool() { };

    string Identify() { return "molpol"; }

    void   Initialize(Property *options);
    bool   Evaluate();

private:

    // FILE I/O OPTIONS
    string          _mps_input;
    string          _mps_output;
    string          _pol_output;
    // MOLPOL OPTIMIZATION
    bool            _do_optimize;
    matrix          _target;
    double          _tolerance;
    // SCF MOLPOL OPTIONS
    double          _aDamp;
    double          _wSOR;
    int             _maxIter;
    double          _epsTol;
    // MOLPOL ENGINE
    MolPolEngine    _molpolengine;

};


}}

#endif