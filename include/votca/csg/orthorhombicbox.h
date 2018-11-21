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

#ifndef _VOTCA_CSG_ORTHORHOMBICBOX_H
#define	_VOTCA_CSG_ORTHORHOMBICBOX_H

#include "boundarycondition.h"

namespace votca { namespace csg {
using namespace std;
using namespace votca::tools;

class OrthorhombicBox : public BoundaryCondition {

public:
    vec BCShortestConnection(const vec &r_i, const vec &r_j) const;

    eBoxtype getBoxType() {
        return typeOrthorhombic;
    }

protected:
};

}}

#endif	/* _VOTCA_CSG_ORTHORHOMBICBOX_H */

