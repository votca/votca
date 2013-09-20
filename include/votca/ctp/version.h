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

/**

\mainpage VOTCA-CTP

This is the code documentation of the VOTCA-CTP package (http://www.votca.org). 

\section installation_sec Installation 

Ideally the following bash script should do the job for you: 

\code
    prefix=~/votca
    mkdir -p ${prefix}/src
    cd ${prefix}/src
    wget http://votca.googlecode.com/hg/build.sh
    chmod +x build.sh
    ./build.sh --prefix ${prefix} --dev tools csg kmc moo ctp
\endcode

\section started_sec Getting started

To run the program, source the VOTCARC file, i.e. in bash

\code
    source ~/votca/bin/VOTCARC.bash
\endcode

 */

#ifndef __VOTCA_MD2QM_VERSION_H
#define	__VOTCA_MD2QM_VERSION_H

#include <string>

namespace votca { namespace ctp {
    const std::string & CtpVersionStr();
    void HelpTextHeader(const std::string &tool_name);
}}

#endif

