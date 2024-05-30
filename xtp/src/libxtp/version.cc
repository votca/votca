/*
 *            Copyright 2009-2024 The VOTCA Development Team
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

// Standard includes
#include <iostream>

// VOTCA includes
#include <votca/tools/version.h>
#include <votca/tools/votca_tools_config.h>

// Local VOTCA includes
#include "votca/xtp/version.h"
#include "votca/xtp/votca_xtp_config.h"

namespace votca {
namespace xtp {

void HelpTextHeader(const std::string &tool_name) {
  std::cout << "==================================================\n"
            << "========   VOTCA (http://www.votca.org)   ========\n"
            << "==================================================\n\n"
            << "please read and cite: " PROJECT_CITATION "\n"
            << "and submit bugs to " PACKAGE_BUGREPORT "\n\n"
            << tool_name << ", version " << votca::tools::ToolsVersionStr()
            << "\n\n";
}

}  // namespace xtp
}  // namespace votca
