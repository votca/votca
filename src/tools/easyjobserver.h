/* 
 * File:   tofobserver.h
 * Author: ruehle
 *
 * Created on December 9, 2009, 9:39 AM
 */

#ifndef _EASYJOBSERVER_H
#define	_EASYJTOFOBSERVER_H

#include <votca/csg/cgobserver.h>
#include <moo/jcalc.h>
#include "qmtopology.h"
#include <votca/tools/property.h>
#include "statesaver.h"
class EasyJObserver
    : public CGObserver
{
public:
    EasyJObserver();
    ~EasyJObserver();

    
    void Initialize(QMTopology &qmtop, Property &opts);

    /// begin coarse graining a trajectory
    void BeginCG(Topology *top, Topology *top_atom);

    /// end coarse graining a trajectory
    void EndCG();

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);

    void setCutoff(const double & cutoff){
        _cutoff = cutoff;
    }
    void setNNnames(string  nnnames);

    void CalcRates(QMNBList &nblist);

    void MakeRatesSIUnits(QMNBList &nblist);

    void print_nbs_to_file(QMNBList &nblist);

protected:
    QMTopology *_qmtop;
    double _cutoff;
    vector <string> _nnnames;
    vector <double> _Js;
    /// electric field [V/m]
    vec _E;
    /// thermal energy [eV]
    double _kT;

    bool MatchNNnames(CrgUnit * crg1, CrgUnit * crg2);
};


#endif	/* _TOFOBSERVER_H */

