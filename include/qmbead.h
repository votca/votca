/* 
 * File:   qmbeads.h
 * Author: james
 *
 * Created on 12 December 2009, 16:46
 */

#ifndef _QMBEAD_H
#define	_QMBEAD_H
/**
    \brief contains the bead associated to a crg unit.

 * It contains the usual bead info + a vector of all the CG beads that make up a bead
*/

#include <votca/csg/bead.h>
#include <moo/crgunit.h>

class QMBead:Bead{
public:
    Bead();
    ~Bead();

    CrgUnit* GetCrgUnit{
        return _crg;
    }
    ///at each frame read update the QM bead and the associated crgunit
    void UpdateQMBead();
private:

    ///the charge unit
    CrgUnit * _crg;
    /// the
    int _pos;
    /// parents contains the information necessary to update the CrgUnit
    vector <Bead *> _parents;

};


#endif	/* _QMBEAD_H */

