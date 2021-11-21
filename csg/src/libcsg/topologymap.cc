/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// Local VOTCA includes
#include "votca/csg/topologymap.h"
#include "votca/csg/boundarycondition.h"

namespace votca {
namespace csg {

void TopologyMap::Apply() {
  out_->setStep(in_->getStep());
  out_->setTime(in_->getTime());
  out_->setBox(in_->getBox());

  for (auto& map_ : maps_) {
    map_.Apply(out_->getBoundary());
  }
}

}  // namespace csg
}  // namespace votca
