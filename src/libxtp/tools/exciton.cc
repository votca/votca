/* 
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/xtp/gwbseengine.h>


using namespace std;

namespace votca {
    namespace xtp {

        void Exciton::Initialize(Property* options) {
         
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

            ctp::TLogLevel _ReportLevel = _log.getReportLevel(); // backup report level

            // Get orbitals object
            Orbitals _orbitals;

            // Read molecular geometry from xyz file and store in a segment (WHY SEGMENT?)
            std::vector <ctp::Segment* > _segments;
            ctp::Segment _segment(0, "mol");
            CTP_LOG(ctp::logDEBUG, _log) << "Reading molecular coordinates from " << _xyzfile << flush;
            ReadXYZ(&_segment, _xyzfile);
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






            /*
   
                   CTP_LOG(ctp::logINFO,_log) << " Summary of iteration " << _iteration +1 << flush;
                   CTP_LOG(ctp::logINFO,_log) << "    total energy of " << _spintype << " " << _opt_state << " is " << energy_new << " Hartree "  << flush;
                   CTP_LOG(ctp::logINFO,_log) << "    BFGS-TRM: trust radius was " << _old_TR << " Bohr" << flush;
                   CTP_LOG(ctp::logINFO,_log) << "    convergence: " << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       energy change: " << setprecision(12) << energy_delta << " a.u. " << Convergence( _energy_converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       RMS force:     " << setprecision(12) << _RMSForce << Convergence( _RMSForce_converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       Max force:     " << setprecision(12) << _MaxForce << Convergence( _MaxForce_converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       RMS step:      " << setprecision(12) << _RMSStep << Convergence( _RMSStep_converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       Max step:      " << setprecision(12) << _MaxStep << Convergence( _MaxStep_converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "       Geometry       " << Convergence( _converged ) << flush;
                   CTP_LOG(ctp::logINFO,_log) << "   " << flush;
                   _log.setReportLevel( _ReportLevel ); // go silent again

             */

            CTP_LOG(ctp::logDEBUG, _log) << "Saving data to " << _archive_file << flush;
            std::ofstream ofs((_archive_file).c_str());
            boost::archive::binary_oarchive oa(ofs);



            oa << _orbitals;
            ofs.close();

            return true;
        }



        /* Calculate forces on atoms numerically by central differences */
        /*void Exciton::NumForceCentral(double energy, std::vector<ctp::Atom*> _atoms, ub::matrix<double>& _force, 
                QMPackage* _qmpackage, std::vector<ctp::Segment*> _molecule, Orbitals* _orbitals){
    

            std::vector< ctp::Atom* > ::iterator ait;
            int _i_atom = 0;
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
               // get its current coordinates
               vec _current_pos = (*ait)->getQMPos(); // in nm
    
               // go through all cartesian components
               for ( int _i_cart = 0 ; _i_cart < 3; _i_cart++ ){
           
                   // get displacement vector in positive direction
                   vec _displaced(0, 0, 0);
                   if (_i_cart == 0) {
                      _displaced.setX(_displacement / tools::conv::nm2bohr ); // x, _displacement in Bohr
                      _current_xyz(_i_atom,_i_cart ) = _current_pos.getX() * tools::conv::nm2bohr; // in Bohr
                   }
                   if (_i_cart == 1) {
                      _current_xyz(_i_atom,_i_cart ) = _current_pos.getY() * tools::conv::nm2bohr; // in Bohr
                      _displaced.setY(_displacement / tools::conv::nm2bohr ); // y, _displacement in in Angstrom
                   }
                   if (_i_cart == 2) {
                      _current_xyz(_i_atom,_i_cart ) = _current_pos.getZ() * tools::conv::nm2bohr; // in Bohr
                      _displaced.setZ(_displacement / tools::conv::nm2bohr ); // z, _displacement in in Angstrom
                   }
                            
                   // update the coordinate
                   vec _pos_displaced = _current_pos + _displaced;
                   (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                   // run DFT and GW-BSE for this geometry
                   ExcitationEnergies(_qmpackage, _molecule, _orbitals);
                            
                   // get total energy for this excited state
                   double energy_displaced_plus = _orbitals->GetTotalEnergy(_spintype, _opt_state);
           
                   // get displacement vector in negative direction

                   // update the coordinate
                   _pos_displaced = _current_pos - 2.0 * _displaced;
                   (*ait)->setQMPos(_pos_displaced); // put updated coordinate into segment

                   // run DFT and GW-BSE for this geometry
                   ExcitationEnergies(_qmpackage, _molecule, _orbitals);
                            
                   // get total energy for this excited state
                   double energy_displaced_minus =  _orbitals->GetTotalEnergy(_spintype, _opt_state);
           
                   // calculate force and put into matrix
                   _force(_i_atom,_i_cart) = 0.5 * (energy_displaced_minus - energy_displaced_plus) / _displacement ; // force a.u./a.u.
                                                  
                   (*ait)->setQMPos(_current_pos); // restore original coordinate into segment
                }
               
              _i_atom++;
           }
    

    
            return;
        } */



       

       

        void Exciton::ReadXYZ(ctp::Segment* _segment, string filename) {

            string line;
            std::ifstream in;


            string label, type;
            vec pos;


            in.open(filename.c_str(), std::ios::in);
            if (!in) throw runtime_error(string("Error reading coordinates from: ")
                    + filename);


            int atomCount = 1;

            if (in.is_open()) {
                while (in.good()) {
                    std::getline(in, line);

                    std::vector< string > split;
                    Tokenizer toker(line, " \t");
                    toker.ToVector(split);
                    if (!split.size() ||
                            split.size() != 4 ||
                            split[0] == "#" ||
                            split[0].substr(0, 1) == "#") {
                        continue;
                    }

                    // Interesting information written here: e.g. 'C 0.000 0.000 0.000'
                    atomCount++;
                    string element = split[0];
                    double x = boost::lexical_cast<double>(split[1]) / 10.; //°A to NM
                    double y = boost::lexical_cast<double>(split[2]) / 10.;
                    double z = boost::lexical_cast<double>(split[3]) / 10.;
                    vec Pos = vec(x, y, z);
                    ctp::Atom *pAtom = new ctp::Atom(atomCount, element);
                    pAtom->setPos(Pos);
                    pAtom->setQMPart(atomCount, Pos);
                    pAtom->setElement(element);
                    _segment->AddAtom(pAtom);

                }
            } else {
                throw std::runtime_error("No such file: '" + filename + "'.");
            }

            return;
        }


    }
}
