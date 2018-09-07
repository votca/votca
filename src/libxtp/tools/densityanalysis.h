/* 
 *            Copyright 2016 The MUSCET Development Team
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

#ifndef _VOTCA_XTP_DENSITYANALYSIS_H
#define _VOTCA_XTP_DENSITYANALYSIS_H

#include <stdio.h>
#include <votca/xtp/gyration.h>
#include <votca/ctp/logger.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class DensityAnalysis : public ctp::QMTool
{
public:

    string Identify() { return "densityanalysis"; }

    void   Initialize(Property *options);
    bool   Evaluate();

private:
    
    string      _orbfile;
    string      _output_file;
    Property    _gyration_options;
    
    ctp::Logger      _log;
    
    
};

void DensityAnalysis::Initialize(Property* options) {
    
    string key = "options." + Identify(); 
    _orbfile      = options->get(key + ".input").as<string> ();

    string _gyration_xml = options->get(key + ".gyration_options").as<string> ();
    load_property_from_xml(_gyration_options,_gyration_xml.c_str());

    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

}

bool DensityAnalysis::Evaluate() {
    
    _log.setReportLevel( ctp::logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ..."); 

    Orbitals orbitals;
    CTP_LOG(ctp::logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    orbitals.ReadFromCpt(_orbfile);

    Density2Gyration density2gyration=Density2Gyration(&_log);
    density2gyration.Initialize(_gyration_options);
    density2gyration.AnalyzeDensity(orbitals);
     
    return true;
}








}}


#endif