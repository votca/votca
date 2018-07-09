/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_QMSANDBOX_H
#define _VOTCA_XTP_QMSANDBOX_H

#include <stdio.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmtool.h>

#include <votca/xtp/qmpackagefactory.h>
#include<votca/xtp/aobasis.h>
#include<votca/xtp/aomatrix.h>

namespace votca { namespace xtp {

    
class QMSandbox : public xtp::QMTool
{
public:

    QMSandbox() { };
   ~QMSandbox() { };

    std::string Identify() { return "qmsandbox"; }

    void   Initialize(tools::Property *options);
    bool   Evaluate();


private:
    
    std::string      _orbfile;
    std::string      _output_file;
    
    xtp::Logger      _log;
 
    std::string      _logfile;

    std::string      _package;
    tools::Property    _package_options; 
    
    void CheckContent(  Orbitals& _orbitals );

};

void QMSandbox::Initialize(tools::Property* options) {

    // update options with the VOTCASHARE defaults   
    //UpdateWithDefaults( options, "xtp" );
 

    std::string key = "options." + Identify();

    // orbitals file or pure DFT output
    
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    
    QMPackageFactory::RegisterAll();
    
    
}

bool QMSandbox::Evaluate() {
    
 
    return true;
}




void QMSandbox::CheckContent( Orbitals& _orbitals ){


   
    
    return;
}

}}


#endif
