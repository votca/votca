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

#include <votca_config.h>
#include <votca/tools/version.h>
#include <votca/csg/version.h>
#include <iostream>
#include <votca/ctp/version.h>


namespace votca { namespace md2qm {

//defines hgversion
#include "hgversion.h"
static const std::string version_str = std::string(VERSION) + " " + hgversion + " (compiled " __DATE__ ", " __TIME__ ")";

const std::string &CtpVersionStr()
{
    return version_str;
}

void HelpTextHeader(const std::string &tool_name)
{
    std::cout
         << "==================================================\n"
         << "========   VOTCA (http://www.votca.org)   ========\n"
         << "==================================================\n\n"
	 << "please submit bugs to " PACKAGE_BUGREPORT "\n\n" 
	 << tool_name << ", version " << votca::md2qm::CtpVersionStr() 
         << "\nvotca_csg, version " << votca::csg::CsgVersionStr()
         << "\n\n";
}

}}


