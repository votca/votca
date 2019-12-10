/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "../statefilters/DeltaQ_filter.h"
#include "../statefilters/Density_filter.h"
#include "../statefilters/Localisation_filter.h"
#include "../statefilters/OscillatorStrength_filter.h"
#include "../statefilters/Overlap_filter.h"
#include <votca/xtp/filterfactory.h>

namespace votca {
namespace xtp {

void FilterFactory::RegisterAll(void) {
  Filter().Register<DeltaQ_filter>("chargeTransfer");
  Filter().Register<Density_filter>("density");
  Filter().Register<Localisation_filter>("localisation");
  Filter().Register<OscillatorStrength_filter>("oscillatorstrength");
  Filter().Register<Overlap_filter>("overlap");
}
}  // namespace xtp
}  // namespace votca
