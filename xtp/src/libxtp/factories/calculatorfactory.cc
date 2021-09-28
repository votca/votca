/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#include "votca/xtp/calculatorfactory.h"

// Local private VOTCA includes
#include "calculators/eanalyze.h"
#include "calculators/einternal.h"
#include "calculators/ianalyze.h"
#include "calculators/kmclifetime.h"
#include "calculators/kmcmultiple.h"
#include "calculators/mapchecker.h"
#include "calculators/neighborlist.h"
#include "calculators/vaverage.h"

namespace votca {
namespace xtp {

void Calculatorfactory::RegisterAll(void) {

  Calculators().Register<Neighborlist>("neighborlist");
  Calculators().Register<MapChecker>("mapchecker");
  Calculators().Register<IAnalyze>("ianalyze");
  Calculators().Register<EAnalyze>("eanalyze");
  Calculators().Register<EInternal>("einternal");
  Calculators().Register<KMCLifetime>("kmclifetime");
  Calculators().Register<KMCMultiple>("kmcmultiple");
  Calculators().Register<VAverage>("vaverage");
}
}  // namespace xtp
}  // namespace votca
