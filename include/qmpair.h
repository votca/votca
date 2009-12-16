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

class QMTopology;

class QMPair :
    public std::pair<CrgUnit *, CrgUnit *>
{
public:
    QMPair() {}
    QMPair(CrgUnit *crg1, CrgUnit *crg2, QMTopology * top);

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

