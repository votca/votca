/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#define _VOTCA_CSG_ORTHORHOMBICBOX_H

#include "boundarycondition.h"

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

class OrthorhombicBox : public BoundaryCondition {

 public:
  TOOLS::vec BCShortestConnection(const TOOLS::vec &r_i,
                                  const TOOLS::vec &r_j) const;

  eBoxtype getBoxType() { return typeOrthorhombic; }

 protected:
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_ORTHORHOMBICBOX_H */
