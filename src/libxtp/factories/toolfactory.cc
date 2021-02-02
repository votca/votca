/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#include "tools/apdft.h"
#include "tools/coupling.h"
#include "tools/densityanalysis.h"
#include "tools/dftgwbse.h"
#include "tools/excitoncoupling.h"
#include "tools/gencube.h"
#include "tools/log2mps.h"
#include "tools/mol2orb.h"
#include "tools/molpol.h"
#include "tools/orb2fchk.h"
#include "tools/orb2mol.h"
#include "tools/partialcharges.h"
#include "tools/qmsandbox.h"
#include "tools/spectrum.h"

namespace votca {
namespace xtp {

void QMToolFactory::RegisterAll(void) {

  this->Register<Log2Mps>("log2mps");
  this->Register<DftGwBse>("dftgwbse");
  this->Register<QMSandbox>("qmsandbox");
  this->Register<Spectrum>("spectrum");
  this->Register<ExcitonCoupling>("excitoncoupling");
  this->Register<GenCube>("gencube");
  this->Register<Partialcharges>("partialcharges");
  this->Register<DensityAnalysis>("densityanalysis");
  this->Register<Coupling>("coupling");
  this->Register<MolPol>("molpol");
  this->Register<APDFT>("apdft");
  this->Register<Mol2Orb>("mol2orb");
  this->Register<Orb2Mol>("orb2mol");
  this->Register<Orb2Fchk>("orb2fchk");
}

}  // namespace xtp
}  // namespace votca
