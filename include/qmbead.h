/* 
 * File:   qmbeads.h
 * Author: james
 *
 * Created on 12 December 2009, 16:46
 */

#ifndef _QMBEAD_H
#define	_QMBEAD_H


#include <votca/csg/bead.h>
#include <votca/moo/crgunit.h>
#include "qmtopology.h"

using namespace votca::csg;

/**
    \brief contains the bead associated to a crg unit.

 * It contains the usual bead info + a vector of all the CG beads that make up a bead
*/

class QMBead : public Bead{
public:
    ~QMBead(){}

    void UpdateCrg(){
        if (_crg != NULL){
            _crg->SetPos(_ipos, Bead::getPos());
            _crg->SetNorm(_ipos, Bead::getU());
            _crg->SetPlane(_ipos, Bead::getV());
        }
    }
    void setCrg(QMCrgUnit * a){
        _crg=a;
    }

    void setiPos(const int & a){
        _ipos=a;
    }

    int getiPos(){
        return _ipos;
    }

    QMCrgUnit* GetCrgUnit(){
        return _crg;
    }

private:
    /// constructor if bead is created with no rerference bead (e.g. when reading from file)
    QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry, string name,
            int resnr, double m, double q);

    ///the charge unit
    QMCrgUnit * _crg;
    /// the integer describing the position in the Crgunit
    int _ipos;

    friend class QMTopology;
};


inline QMBead::QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry,
    string name, int resnr, double m, double q)
    : Bead(owner, id, type, symmetry, name, resnr, m, q), _crg(NULL), _ipos(-1)
{
}

#endif	/* _QMBEAD_H */

