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
#include "votca/xtp/toolfactory.h"

// Local private VOTCA includes
#include "libxtp/tools/apdft.h"
#include "libxtp/tools/coupling.h"
#include "libxtp/tools/densityanalysis.h"
#include "libxtp/tools/dftgwbse.h"
#include "libxtp/tools/excitoncoupling.h"
#include "libxtp/tools/gencube.h"
#include "libxtp/tools/log2mps.h"
#include "libxtp/tools/molpol.h"
#include "libxtp/tools/partialcharges.h"
#include "libxtp/tools/qmsandbox.h"
#include "libxtp/tools/spectrum.h"

namespace votca {
namespace xtp {

void QMToolFactory::RegisterAll(void) {

  QMTools().Register<Log2Mps>("log2mps");
  QMTools().Register<DftGwBse>("dftgwbse");
  QMTools().Register<QMSandbox>("qmsandbox");
  QMTools().Register<Spectrum>("spectrum");
  QMTools().Register<ExcitonCoupling>("excitoncoupling");
  QMTools().Register<GenCube>("gencube");
  QMTools().Register<Partialcharges>("partialcharges");
  QMTools().Register<DensityAnalysis>("densityanalysis");
  QMTools().Register<Coupling>("coupling");
  QMTools().Register<MolPol>("molpol");
  QMTools().Register<APDFT>("apdft");
}

}  // namespace xtp
}  // namespace votca
