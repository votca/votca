/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/tools/version.h"

// Local private VOTCA includes
#include "votca_tools_config.h"

namespace votca {
namespace tools {

// defines gitversion
#include "gitversion.h"
static const std::string version_str = std::string(PROJECT_VERSION) + " " +
                                       gitversion +
                                       " (compiled " __DATE__ ", " __TIME__ ")";

const std::string &ToolsVersionStr() { return version_str; }

}  // namespace tools
}  // namespace votca
