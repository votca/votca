/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#include "tools/diabatization.h"
#include "tools/excitoncoupling.h"
#include "tools/gencube.h"
#include "tools/gpu_benchmark.h"
#include "tools/log2mps.h"
#include "tools/mol2orb.h"
#include "tools/molpol.h"
#include "tools/orb2fchk.h"
#include "tools/orb2mol.h"
#include "tools/partialcharges.h"
#include "tools/spectrum.h"

namespace votca {
namespace xtp {

void QMToolFactory::RegisterAll(void) {

  QMTools().Register<Log2Mps>("log2mps");
  QMTools().Register<DftGwBse>("dftgwbse");
  QMTools().Register<Spectrum>("spectrum");
  QMTools().Register<ExcitonCoupling>("excitoncoupling");
  QMTools().Register<GenCube>("gencube");
  QMTools().Register<Partialcharges>("partialcharges");
  QMTools().Register<DensityAnalysis>("densityanalysis");
  QMTools().Register<Coupling>("coupling");
  QMTools().Register<MolPol>("molpol");
  QMTools().Register<APDFT>("apdft");
  QMTools().Register<Mol2Orb>("mol2orb");
  QMTools().Register<Orb2Mol>("orb2mol");
  QMTools().Register<Diabatization>("diabatization");
  QMTools().Register<Orb2Fchk>("orb2fchk");
  QMTools().Register<GPUBenchmark>("gpu_benchmark");
}

}  // namespace xtp
}  // namespace votca
