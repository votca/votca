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

#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/eigen.h"
#include <libint2/initialize.h>
namespace votca {
namespace xtp {
bool QMCalculator::EvaluateFrame(Topology& top) {
  libint2::initialize();
  OPENMP::setMaxThreads(_nThreads);
  std::cout << " Using " << OPENMP::getMaxThreads() << " threads" << std::flush;
  bool success = Evaluate(top);
  libint2::finalize();
  return success;
}

void QMCalculator::Initialize(const tools::Property& opt) {
  tools::Property options = LoadDefaultsAndUpdateWithUserOptions("xtp", opt);
  ParseOptions(options);
}

}  // namespace xtp
}  // namespace votca
