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

#ifndef _VOTCA_CSG_CGOBSERVER_H
#define _VOTCA_CSG_CGOBSERVER_H

#include "topology.h"

namespace votca {
namespace csg {

/**
   \brief Observer class for analysis hook

   Each application which performs analysis operations should use CGEngine. It
   offers a hook (callback class) during the coarse-graining process to evaluate
   each frame. The user does not have to take care about mapping and other
   stoff. Just oberload this class and analyze properties of interest.

 */

class CGObserver {
 public:
  /// \brief called before the first frame
  virtual void BeginCG(Topology *top, Topology *top_atom = 0) = 0;
  /// \brief called after the last frame
  virtual void EndCG() = 0;
  // \brief called for each frame which is mapped
  virtual void EvalConfiguration(Topology *top, Topology *top_atom = 0) = 0;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_CGOBSERVER_H */
