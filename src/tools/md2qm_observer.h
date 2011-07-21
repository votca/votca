/* 
 * File:   md2qm_observer.h
 * Author: vehoff
 *
 * Created on March 5, 2010, 3:20 PM
 */

#ifndef _MD2QM_OBSERVER_H
#define	_MD2QM_OBSERVER_H

#include <votca/csg/cgobserver.h>
#include <votca/ctp/qmtopology.h>
#include <votca/tools/property.h>
#include <votca/ctp/statesaversqlite.h>
#include <votca/moo/units.h>

class MD2QMObserver
    : public CGObserver
{
public:
    MD2QMObserver();
    ~MD2QMObserver();


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
    void setOut(const string &  out){
        _out=out;
    }
    void setNNnames(string  nnnames);

    void print_nbs_to_file(QMNBList &nblist);

protected:
    QMTopology *_qmtop;
    /// nearest neighbor cut-off radius
    double _cutoff;
    ///  output streams for velocity averaging & diffusion
    ofstream _out_cont;
    ofstream _out_diff;
    string _out;
    StateSaverSQLite _save;
};

#endif	/* _MD2QM_OBSERVER_H */

