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

#ifndef _VOTCA_CSG_NBLISTINTRA_H
#define _VOTCA_CSG_NBLISTINTRA_H

// Standard includes
#include <vector>

// VOTCA includes
#include <votca/tools/eigen.h>

// Local VOTCA includes
#include "nblist.h"
#include "votca/tools/NDimVector.h"

namespace votca {
namespace csg {

class NBListIntra : public NBList {
 public:
  void Generate(BeadList &list1, BeadList &list2,
          bool do_exclusions = true) override;
  void Generate(BeadList &list, bool do_exclusions = true) override;

 protected:
  struct cell_t {
    BeadList _beads;
    std::vector<cell_t *> _neighbours;
  };
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_NBLISTINTRA */
