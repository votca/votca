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

    void setPosCrg(){
        _crg->SetPos(_ipos, r);
    }

    void setNormCrg(){
        _crg->SetNorm(_ipos, Bead::GetU());
    }

    void setPlaneCrg(){
        _crg->SetPlane(_ipos, Bead::GetV());
    }
    CrgUnit* GetCrgUnit(){
        return _crg;
    }
    ///at each frame read update the QM bead and the associated crgunit
    ///this is not necessary. in order to update the _crg it is sufficient to overload
    ///the functions setU/setV/setW.
    ///void UpdateQMBead();
private:
    /// constructor if bead is created with no rerference bead (e.g. when reading from file)
    QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry, string name,
            int resnr, double m, double q);

    ///the charge unit
    CrgUnit * _crg;
    /// the
    int _ipos;
    /// parents contains the information necessary to update the CrgUnit
    ///vector <Bead *> _parents;

    friend class QMTopology;
};


inline QMBead::QMBead(Topology *owner, int id, BeadType *type, byte_t symmetry,
    string name, int resnr, double m, double q)
    : Bead(owner, id, type, symmetry, name, resnr, m, q)
{
}

#endif	/* _QMBEAD_H */

