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
#include "votca/xtp/qmpackagefactory.h"

// Local private VOTCA includes
#include "qmpackages/orca.h"
#include "qmpackages/xtpdft.h"

namespace votca {
namespace xtp {

void QMPackageFactory::RegisterAll(void) {
  QMPackages().Register<Orca>("orca");
  QMPackages().Register<XTPDFT>("xtp");
}
}  // namespace xtp
}  // namespace votca
