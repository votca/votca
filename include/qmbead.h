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

class QMTopology;

class QMBead : public Bead{
public:
    QMBead();
    ~QMBead();

    void UpdateCrg(){
        if (_crg != NULL){
            _crg->SetPos(_ipos, Bead::GetPos());
            _crg->SetNorm(_ipos, Bead::GetU());
            _crg->SetPlane(_ipos, Bead::GetV());
        }
    }
    CrgUnit* GetCrgUnit(){
        return _crg;
    }
private:
    /// constructor if bead is created with no rerference bead (e.g. when reading from file)
    QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry, string name,
            int resnr, double m, double q);

    ///the charge unit
    CrgUnit * _crg;
    /// the integer describing the position in the Crgunit
    int _ipos;

    friend class QMTopology;
};


inline QMBead::QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry,
    string name, int resnr, double m, double q)
    : Bead(owner, id, type, symmetry, name, resnr, m, q)
{
}

#endif	/* _QMBEAD_H */

