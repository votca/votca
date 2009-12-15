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

class QMPair :
    public std::pair<CrgUnit *, CrgUnit *>
{
public:
    QMPair() {}
    QMPair(CrgUnit *crg1, CrgUnit *crg2, vec r)
      : std::pair<CrgUnit *, CrgUnit *>(crg1, crg2), _r(r), _dist(abs(r)) {}

    /// \brief the vector connecting two beads
    vec &r() { return _r; }
    /// \brief the distance of the beads
    double &dist() { return _dist; }
protected:
    vec _r;
    double _dist;
};

#endif	/* _QMBEADPAIR_H */

