/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_BOUNDARYCONDITION_H
#define	_VOTCA_CSG_BOUNDARYCONDITION_H

#include <votca/tools/matrix.h>

namespace votca { namespace csg {
using namespace votca::tools;

class BoundaryCondition {

public:
    virtual ~BoundaryCondition() {};

    /**
     * set the simulation box
     * \param box triclinic box matrix
     */
    void setBox(const matrix &box) { _box = box; };

    /**
     * get the simulation box
     * \return triclinic box matrix
     */
    const matrix &getBox() { return _box; };

    /**
     * get the volume of the box
     * \return box volume as double
     */
    virtual double BoxVolume();

    /**
     * get shortest connection vector between r_i and r_j with respect to the (periodic) box
     * \return shortest distance vector
     */
    virtual vec BCShortestConnection(const vec &r_i, const vec &r_j) const = 0;

    enum eBoxtype {
        typeAuto = 0,
        typeTriclinic,
        typeOrthorhombic,
        typeOpen
    };
    virtual eBoxtype getBoxType() = 0;

protected:
    matrix _box;
};

}}

#endif	/* _VOTCA_CSG_BOUNDARYCONDITION_H */

