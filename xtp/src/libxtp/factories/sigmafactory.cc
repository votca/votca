
/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

// Local VOTCA includes
#include <votca/xtp/sigmafactory.h>

// Local private VOTCA includes
#include "self_energy_evaluators/sigma_cda.h"
#include "self_energy_evaluators/sigma_exact.h"
#include "self_energy_evaluators/sigma_ppm.h"

namespace votca {
namespace xtp {

void SigmaFactory::RegisterAll(void) {
  this->Register<Sigma_CDA>("cda");
  this->Register<Sigma_Exact>("exact");
  this->Register<Sigma_PPM>("ppm");
}
}  // namespace xtp
}  // namespace votca
