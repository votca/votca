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

#include "bondedstatistics.h"

using namespace votca::tools;

namespace votca {
namespace csg {

void BondedStatistics::BeginCG(Topology *top, Topology *top_atom) {
  InteractionContainer &ic = top->BondedInteractions();
  InteractionContainer::iterator ia;

  _bonded_values.clear();
  for (ia = ic.begin(); ia != ic.end(); ++ia) {
    _bonded_values.CreateArray((*ia)->getName());
  }
}

void BondedStatistics::EndCG() {}

void BondedStatistics::EvalConfiguration(Topology *conf, Topology *conv_atom) {
  InteractionContainer &ic = conf->BondedInteractions();
  InteractionContainer::iterator ia;

  DataCollection<double>::container::iterator is;
  for (ia = ic.begin(), is = _bonded_values.begin(); ia != ic.end();
       ++ia, ++is) {
    (*is)->push_back((*ia)->EvaluateVar(*conf));
  }
}

} // namespace csg
} // namespace votca
