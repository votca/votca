/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <iostream>
#include <votca/csg/version.h>
#include <votca/tools/version.h>
#include <votca/xtp/version.h>
#include <votca/xtp/votca_config.h>

extern "C" {
void VotcaMd2QmFromC() {
  // do nothing - this is just that we have a c function for autotools
}
}

namespace votca {
namespace xtp {

// defines gitversion
#include "gitversion.h"
static const std::string version_str = std::string(VERSION) + " " + gitversion +
                                       " (compiled " __DATE__ ", " __TIME__ ")";

const std::string &XtpVersionStr() { return version_str; }

void HelpTextHeader(const std::string &tool_name) {
  std::cout << "==================================================\n"
            << "========   VOTCA (http://www.votca.org)   ========\n"
            << "==================================================\n\n"
            << "please submit bugs to " PACKAGE_BUGREPORT "\n\n"
            << tool_name << ", version " << XtpVersionStr()
            << "\nvotca_csg, version " << votca::csg::CsgVersionStr()
            << "\nvotca_tools, version " << votca::tools::ToolsVersionStr()
            << "\n\n";
}

}  // namespace xtp
}  // namespace votca
