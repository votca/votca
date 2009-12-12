/* 
 * File:   qmbeads.h
 * Author: james
 *
 * Created on 12 December 2009, 16:46
 */

#ifndef _QMBEAD_H
#define	_QMBEAD_H

#include <votca/csg/bead.h>
#include <moo/crgunit.h>

/**
    \brief contains the bead associated to a crg unit.

    It contains the usual bead info + a pointer to a crg unit + position which
    it should update

*/

class QMBead : public Bead
{
public:
    
    ///at each frame read update the QM bead and the associated crgunit
    void UpdateCrgUnit();
private:

    ///the charge unit
    CrgUnit * _crg;
    int _pos;
    
};


#endif	/* _QMBEAD_H */

