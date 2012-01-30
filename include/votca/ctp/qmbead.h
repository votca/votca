/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _QMBEAD_H
#define	_QMBEAD_H


#include <votca/csg/bead.h>
#include <votca/moo/crgunit.h>
#include "qmtopology.h"

namespace votca { namespace ctp {

using namespace votca::csg;
namespace CSG = votca::csg;
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
    QMBead(CSG::Topology *owner, int id, BeadType *type, byte_t symmetry, string name,
            int resnr, double m, double q);

    ///the charge unit
    QMCrgUnit * _crg;
    /// the integer describing the position in the Crgunit
    int _ipos;

    friend class QMTopology;
};


inline QMBead::QMBead(CSG::Topology *owner, int id, BeadType *type, byte_t symmetry,
    string name, int resnr, double m, double q)
    : Bead(owner, id, type, symmetry, name, resnr, m, q), _crg(NULL), _ipos(-1)
{
}

}}

#endif	/* _QMBEAD_H */

