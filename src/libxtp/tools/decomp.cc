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

#include "decomp.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>
#include "votca/xtp/pmdecomposition.h"

namespace votca {
namespace xtp {

void Decomp::ParseOptions(const tools::Property&) {
  orbitals.ReadFromCpt(_job_name + ".orb");
}

bool Decomp::Run() {
  log.setReportLevel(Log::current_level);
  log.setMultithreading(true);
  log.setCommonPreface("\n... ...");
  XTP_LOG(Log::error,log) << "Starting decomp tool" << std::endl;
  PMDecomposition pmd(orbitals, log);
  pmd.compute();
  XTP_LOG(Log::error, log) << "There you go!!" << std::endl;
  return true;
}

}  // namespace xtp
}  // namespace votca
