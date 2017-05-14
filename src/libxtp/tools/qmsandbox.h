/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <votca/ctp/logger.h>
#include <votca/ctp/qmtool.h>

#include <votca/xtp/qmpackagefactory.h>
#include<votca/xtp/aobasis.h>
#include<votca/xtp/aomatrix.h>
#include<boost/numeric/ublas/matrix.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class QMSandbox : public ctp::QMTool
{
public:

    QMSandbox() { };
   ~QMSandbox() { };

    string Identify() { return "qmsandbox"; }

    void   Initialize(Property *options);
    bool   Evaluate();


private:
    
    string      _orbfile;
    string      _output_file;
    
    ctp::Logger      _log;
 
    string      _logfile;

    string      _package;
    Property    _package_options; 
    
    void CheckContent(  Orbitals& _orbitals );

};

void QMSandbox::Initialize(Property* options) {

    // update options with the VOTCASHARE defaults   
    //UpdateWithDefaults( options, "xtp" );
 

    string key = "options." + Identify();

    // orbitals file or pure DFT output
    _logfile  = options->get(key + ".log").as<string> ();
    _package  = options->get(key + ".package").as<string> ();
    
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