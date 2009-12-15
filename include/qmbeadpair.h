/* 
 * File:   qmbeadpair.h
 * Author: james
 *
 * Created on 14 December 2009, 08:47
 */

#ifndef _QMBEADPAIR_H
#define	_QMBEADPAIR_H

#include <votca/csg/beadpair.h>
#include "qmbead.h"

class QMBeadPair : public BeadPair
{
public:
    QMBeadPair() {}
    QMBeadPair(Bead *bead1, Bead *bead2, vec r);

    CrgUnit *CrgUnit1() { return dynamic_cast<QMBead*>(first)->GetCrgUnit(); }
    CrgUnit *CrgUnit2() { return dynamic_cast<QMBead*>(second)->GetCrgUnit(); }
};

inline QMBeadPair::QMBeadPair(Bead *bead1, Bead *bead2, vec r)
    : BeadPair(bead1, bead2, r)
{
}

#endif	/* _QMBEADPAIR_H */

