/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include "exciton.h"
#include "votca/xtp/qminterface.h"
#include <votca/xtp/gwbseengine.h>


using namespace std;

namespace votca {
    namespace xtp {

        void Exciton::Initialize(tools::Property* options) {
         
            _do_optimize = false;
            // _do_guess=false; //Writing guess for dimer calculation

            std::string key = "options." + Identify();
            _archive_file = options->ifExistsReturnElseReturnDefault<string>(key + ".archive", "system.orb");
            _reporting    = options->ifExistsReturnElseReturnDefault<string>(key + ".reporting", "default");

            // job tasks
            string _tasks_string = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".tasks");
            if (_tasks_string.find("optimize") != std::string::npos) _do_optimize = true;

            // GWBSEENGINE options
            _gwbseengine_options = options->get(key + ".gwbse_engine"); 
            
            // options for dft package
            string _package_xml = options->get(key + ".dftpackage").as<string> ();
            load_property_from_xml(_package_options, _package_xml.c_str());
            _package = _package_options.get("package.name").as<string> ();

            // MOLECULE properties
            _xyzfile = options->ifExistsReturnElseThrowRuntimeError<string>(key + ".molecule");

            // XML OUTPUT
            _xml_output = options->ifExistsReturnElseReturnDefault<string>(key + ".output", "exciton.out.xml");
            
            // if optimization is chosen, get options for geometry_optimizer
            if (_do_optimize) _geoopt_options = options->get(key + ".geometry_optimization");

            // register all QM packages (Gaussian, TURBOMOLE, etc)
            QMPackageFactory::RegisterAll();

        }

        bool Exciton::Evaluate() {


            if (_reporting == "silent")  _log.setReportLevel(ctp::logERROR); // only output ERRORS, GEOOPT info, and excited state info for trial geometry
            if (_reporting == "noisy")   _log.setReportLevel(ctp::logDEBUG); // OUTPUT ALL THE THINGS
            if (_reporting == "default") _log.setReportLevel(ctp::logINFO); // 

            _log.setMultithreading(true);
            _log.setPreface(ctp::logINFO,    "\n... ...");
            _log.setPreface(ctp::logERROR,   "\n... ...");
            _log.setPreface(ctp::logWARNING, "\n... ...");
            _log.setPreface(ctp::logDEBUG,   "\n... ...");
            
            // Get orbitals object
            Orbitals _orbitals;

            // Read molecular geometry from xyz file and store in a segment (WHY SEGMENT?)
            std::vector <ctp::Segment* > _segments;
            ctp::Segment _segment(0, "mol");
            
            if(_orbitals.hasQMAtoms()){
              
            }
            CTP_LOG(ctp::logDEBUG, _log) << "Reading molecular coordinates from " << _xyzfile << flush;
            Orbitals temp;
            temp.LoadFromXYZ(_xyzfile);
            QMInterface qminterface;
            qminterface.Orbitals2Segment(&_segment,&temp);
            _segments.push_back(&_segment);

            // Get and initialize QMPackage for DFT ground state
            QMPackage *_qmpackage = QMPackages().Create(_package);
            _qmpackage->setLog(&_log);
            _qmpackage->Initialize(&_package_options);
            _qmpackage->setRunDir(".");

            // Get GWBSEENGINE Object and initialize
            GWBSEENGINE _gwbse_engine;
            _gwbse_engine.setLog(&_log);
            _gwbse_engine.Initialize(&_gwbseengine_options, _archive_file);

            if ( _do_optimize ) {
                // Run Geometry Optimization
                GeometryOptimization _geoopt(_gwbse_engine,_qmpackage, _segments, &_orbitals);
                _geoopt.setLog(&_log);
                _geoopt.Initialize(&_geoopt_options);
                _geoopt.Evaluate();
            } else {
                // Run GWBSE
                _gwbse_engine.ExcitationEnergies(_qmpackage, _segments, &_orbitals);
            }

            CTP_LOG(ctp::logDEBUG, _log) << "Saving data to " << _archive_file << flush;
            _orbitals.WriteToCpt(_archive_file);
            
            tools::Property _summary = _gwbse_engine.ReportSummary();
            if(_summary.exists("output")){  //only do gwbse summary output if we actually did gwbse
                tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
                CTP_LOG(ctp::logDEBUG, _log) << "Writing output to " << _xml_output << flush;
                std::ofstream ofout(_xml_output.c_str(), std::ofstream::out);
                ofout << (_summary.get("output"));
                ofout.close();
            }

            return true;
        }

       
    }
}
