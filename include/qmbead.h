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
#include "qmtopology.h"

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
    /// constructor if bead is created with no rerference bead (e.g. when reading from file)
    QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry, string name,
            int resnr, double m, double q);

    ///the charge unit
    CrgUnit * _crg;
    int _pos;

    friend class QMTopology;
};


inline QMBead::QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry,
    string name, int resnr, double m, double q)
    : Bead(owner, id, *type, symmetry, name, resnr, m, q)
{
}

#endif	/* _QMBEAD_H */

