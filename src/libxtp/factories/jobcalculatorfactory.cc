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
#include "votca/xtp/jobcalculatorfactory.h"

// Local private VOTCA includes
#include "libxtp/jobcalculators/eqm.h"
#include "libxtp/jobcalculators/iexcitoncl.h"
#include "libxtp/jobcalculators/iqm.h"
#include "libxtp/jobcalculators/qmmm.h"

namespace votca {
namespace xtp {

void JobCalculatorfactory::RegisterAll(void) {
  JobCalculators().Register<IQM>("iqm");
  JobCalculators().Register<EQM>("eqm");
  JobCalculators().Register<IEXCITON>("iexcitoncl");
  JobCalculators().Register<QMMM>("qmmm");
}

}  // namespace xtp
}  // namespace votca
