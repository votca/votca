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

#ifndef _VOTCA_XTP_PARTIALCHARGES_H
#define _VOTCA_XTP_PARTIALCHARGES_H

#include <stdio.h>
#include <votca/xtp/esp2multipole.h>
#include <votca/xtp/logger.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    
class Partialcharges : public QMTool
{
public:

    Partialcharges () { };
   ~Partialcharges () { };

    std::string Identify() { return "partialcharges"; }

    void   Initialize(tools::Property *options);
    bool   Evaluate();
    // two access functions for egwbse interface
    

private:
    
    std::string      _orbfile;
    std::string      _output_file;
    tools::Property    _esp_options;
    
    Logger      _log;
    
    
};

void Partialcharges::Initialize(tools::Property* options) {
    
            // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();
 
    _orbfile      = options->get(key + ".input").as<std::string> ();
    _output_file  = options->get(key + ".output").as<std::string> ();
    std::string _esp2multipole_xml = options->get(key + ".esp_options").as<std::string> ();
    load_property_from_xml(_esp_options,_esp2multipole_xml.c_str());
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

}

bool Partialcharges::Evaluate() {
    
    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    XTP_LOG(logDEBUG, _log) << "Converting serialized QM data in " << _orbfile << flush;

    Orbitals orbitals;
    // load the QM data from serialized orbitals object

    XTP_LOG(logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    orbitals.ReadFromCpt(_orbfile);

    Esp2multipole esp2multipole=Esp2multipole(&_log);
    esp2multipole.Initialize(_esp_options);
    esp2multipole.Extractingcharges(orbitals);
    
    esp2multipole.WritetoFile(_output_file,orbitals);
    
    XTP_LOG(logDEBUG, _log) << "Written charges to " << _output_file << flush;
    
    return true;
}








}}


#endif
