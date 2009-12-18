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
#include <libxml/parser.h>
#include "qmtopology.h"
#include <sys/stat.h>

class ProJObserver
    : public CGObserver
{
public:
    ProJObserver();
    ~ProJObserver();

    
    void setQMTopology(QMTopology &qmtop) {_qmtop = &qmtop; }

    /// begin coarse graining a trajectory
    void BeginCG(Topology *top, Topology *top_atom);

    /// end coarse graining a trajectory
    void EndCG();

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom = 0);

    void setCutoff(const double & cutoff){
        _cutoff = cutoff;
    }
protected:
    QMTopology *_qmtop;
    double _cutoff;
};


#endif	/* _TOFOBSERVER_H */

