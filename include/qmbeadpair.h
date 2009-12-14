/* 
 * File:   qmbeadpair.h
 * Author: james
 *
 * Created on 14 December 2009, 08:47
 */

#ifndef _QMBEADPAIR_H
#define	_QMBEADPAIR_H

#include <votca/tools/vec.h>

class QMBeadPair:: public std::pair<CrgUnit *, CrgUnit *>{
public:
    QMBeadPair(){}
    QMBeadPair(QMBead *bead1, QMBead *bead2, vec r)
    :std::pair<CrgUnit *, CrgUnit *> (bead1->GetCrgUnit(),bead2->GetCrgUnit()), _r(r), _d(abs(r){}


    /// \brief the vector connecting two beads
    vec &r() { return _r; }
    /// \brief the distance of the beads
    double &dist() { return _dist; }
    
private:
    vec _r;
    double _d;
};

#endif	/* _QMBEADPAIR_H */

