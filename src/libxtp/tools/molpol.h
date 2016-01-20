#ifndef MOLPOLTOOL_H
#define MOLPOLTOOL_H


#include <votca/xtp/qmtool.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/molpolengine.h>


namespace votca { namespace xtp {


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
    vector<string>  _scaling_pattern;
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