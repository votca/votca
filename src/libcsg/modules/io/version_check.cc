/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include "votca_config.h"

#if GMX == 50
#include <gromacs/legacyheaders/copyrite.h>
#elif GMX == 45
#include <gromacs/copyrite.h>
#elif GMX == 40
    extern "C"
    {
        #include <copyrite.h>
    }
#endif

#ifdef GMX
// this one is needed because of bool is defined in one of the headers included by gmx
#undef bool
#endif

#include "version_check.h"
#include "gmx_version.h"
#include <iostream>

namespace votca { namespace csg { namespace gmx {
    using namespace votca::csg;
    
    void CheckVersion() {

        std::string GromacsVersionString =
#ifdef GMX
         GromacsVersion()
#ifdef GMX_DOUBLE
         + std::string(" (double precision)");
#else
         + std::string(" (single precision)");
#endif
#endif
        if(GmxVersionStr()==std::string("VERSION NOT SET")) {
            std::cout << "WARNING: GROMACS version string not set." << std::endl;
        }
        else if(GmxVersionStr() != GromacsVersionString) {
            std::cout << "WARNING: VOTCA was compiled using a different Gromacs library version\n"
                    << "compiled: " << GmxVersionStr() <<
                    "\nloaded:" << GromacsVersionString <<
                    "\nTry to source another GMXRC or be prepared for unexpectred behaviour." << std::endl;
        }
    }

}}}
