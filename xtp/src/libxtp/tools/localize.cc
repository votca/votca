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

#include "localize.h"
#include <sstream>

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/pmlocalization.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Localize::ParseOptions(const tools::Property& options) {
  options_ = options;
  orbitals.ReadFromCpt(job_name_ + ".orb");
  std::string temp = options.get("activeatoms").as<std::string>();
  std::stringstream ss(temp);
  Index tmp;
  while (ss >> tmp) {
    activeatoms.push_back(tmp);
  }
  std::cout << "Atoms in active region: ";
  for (const auto& atom : activeatoms) {
    std::cout << atom << " ";
  }
  std::cout << std::endl;
}

bool Localize::Run() {
  log.setReportLevel(Log::current_level);
  log.setMultithreading(true);
  log.setCommonPreface("\n... ...");
  XTP_LOG(Log::error, log) << "Starting localization tool" << std::endl;
  PMLocalization pml(log, options_);
  pml.computePML(orbitals);
  XTP_LOG(Log::error, log) << "Computing Dmat_A now" << std::endl;
  ActiveDensityMatrix Dmat_A(orbitals, activeatoms);
  Dmat_A.compute_Dmat_A();
  return true;
}

}  // namespace xtp
}  // namespace votca
