/* 
 * File:   qmbeadpair.h
 * Author: james
 *
 * Created on 14 December 2009, 08:47
 */

#ifndef _QMPAIR_H
#define	_QMPAIR_H

#include <moo/crgunit.h>
#include <utility>
#include "qmtopology.h"


class QMPair :
    public std::pair<CrgUnit *, CrgUnit *>
{
public:
    QMPair() {}
    QMPair(CrgUnit *crg1, CrgUnit *crg2, QMTopology * top)
    {
        _r = top->BCShortestConnection(crg1->GetCom(), crg2->GetCom());
        _dist = abs(_r);

        // check if PBC:
        vec d = crg2->GetCom() - crg1->GetCom();
        if (d !=  _r){
            _ghost = new CrgUnit(*crg2);
            vec displ = _r - d;
            _ghost->shift(displ);
            std::pair<CrgUnit *, CrgUnit *>(crg1, _ghost);
        }
        else {
            std::pair<CrgUnit *, CrgUnit *>(crg1,crg2);
            _ghost=NULL;
        }
    }

    ~QMPair(){
        if(_ghost != NULL)
            delete _ghost;
    }
    /// \brief the vector connecting two beads
    vec &r() { return _r; }
    /// \brief the distance of the beads
    double &dist() { return _dist; }
protected:
    vec _r;
    double _dist;
    CrgUnit * _ghost;
};

#endif	/* _QMBEADPAIR_H */

