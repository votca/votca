/*
 * Copyright 2009-2024 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <iostream>

// Third party includes
#ifdef GMX_DOUBLE
#include <gromacs/utility/baseversion.h>
// this one is needed because of bool is defined in one of the headers included
// by gmx
#undef bool
#endif

// VOTCA includes
#include <votca/tools/version.h>
#include <votca/tools/votca_tools_config.h>

// Local VOTCA includes
#include "votca/csg/version.h"

// Local private VOTCA includes
#include "votca_csg_config.h"

namespace votca {
namespace csg {

void HelpTextHeader(const std::string &tool_name) {
  std::cout << "==================================================\n"
            << "========   VOTCA (http://www.votca.org)   ========\n"
            << "==================================================\n\n"
            << "please read and cite: " PROJECT_CITATION "\n"
            << "please submit bugs to " PROJECT_BUGREPORT "\n\n"
            << tool_name << ", version " << votca::tools::ToolsVersionStr()
#ifdef GMX_DOUBLE
            << "\ngromacs, " << gmx_version()
#if (GMX_DOUBLE == 1)
            << " (double precision)"
#else
            << " (single precision)"
#endif
#else
            << "\n"
#endif
            << "\n\n";
// sanity check if GMX is the write precision
#ifdef GMX_DOUBLE
#if (GMX_DOUBLE == 1)
  gmx_is_double_precision();
#else
  gmx_is_single_precision();
#endif
#endif
}

}  // namespace csg
}  // namespace votca
