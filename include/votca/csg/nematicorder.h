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

#ifndef _VOTCA_CSG_NEMATICORDER_H
#define	_VOTCA_CSG_NEMATICORDER_H

#include "topology.h"
#include "topology.h"
#include <votca/tools/matrix.h>

namespace votca { namespace csg {
using namespace votca::tools;

class NematicOrder
{
public:
    NematicOrder() {}
    ~NematicOrder() {}
    
    void Process(Topology &top, const string &filter = "*");
    
    matrix::eigensystem_t &NematicU() {return _nemat_u; }
    matrix::eigensystem_t &NematicV() {return _nemat_v; }
    matrix::eigensystem_t &NematicW() {return _nemat_w; }

private:
    matrix _mu,_mv,_mw;
    matrix::eigensystem_t _nemat_u, _nemat_v, _nemat_w;
};

}}

#endif	/* _VOTCA_CSG_NEMATICORDER_H */

