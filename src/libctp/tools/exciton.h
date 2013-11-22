/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef _VOTCA_CTP_EXCITON_H
#define _VOTCA_CTP_EXCITON_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/ctp/mbpt.h>
#include <votca/ctp/qmpackagefactory.h>

namespace votca { namespace ctp {
    using namespace std;
    
class Exciton : public QMTool
{
public:

    Exciton() { };
   ~Exciton() { };

    string Identify() { return "exciton"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    MBPT _mbpt;
 

private:
    
    string      _orbfile;
    string      _logfile;

    string      _package;
    Property    _package_options; 
    
    string      _output_file;
    
    Logger      _log;

};

void Exciton::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options );
            // setting some defaults
            _mbpt.set_do_qp_diag(false);
            _mbpt.set_do_bse_singlets(false);
            _mbpt.set_do_bse_triplets(false);
            _mbpt.set_ranges("default");
            _mbpt.set_store_qp_pert(true);

            // _bse_nmax        = 100;

            string key = "options." + Identify();
            // _jobfile = options->get(key + ".file").as<string>();

            key = "options." + Identify();
            // getting level ranges 
            _mbpt.set_ranges(options->get(key + ".ranges").as<string> ());
            // now check validity, and get rpa, qp, and bse level ranges accordingly
            if (_mbpt.get_ranges() == "factor") {
                // get factors
                _mbpt.set_rpamaxfactor(options->get(key + ".rpamax").as<double> ());
                _mbpt.set_qpminfactor(options->get(key + ".qpmin").as<double> ());
                _mbpt.set_qpmaxfactor(options->get(key + ".qpmax").as<double> ());
                _mbpt.set_bseminfactor(options->get(key + ".bsemin").as<double> ());
                _mbpt.set_bsemaxfactor(options->get(key + ".bsemax").as<double> ());
            } else if (_mbpt.get_ranges() == "explicit") {
                //get explicit numbers
                _mbpt.set_rpamax(options->get(key + ".rpamax").as<unsigned int> ());
                _mbpt.set_qpmin(options->get(key + ".qpmin").as<unsigned int> ());
                _mbpt.set_qpmax(options->get(key + ".qpmax").as<unsigned int> ());
                _mbpt.set_bse_vmin(options->get(key + ".bsemin").as<unsigned int> ());
                _mbpt.set_bse_cmax(options->get(key + ".bsemax").as<unsigned int> ());
            } else if (_mbpt.get_ranges() == "") {
                _mbpt.set_ranges("default");
            } else {
                cerr << "\nSpecified range option " << _mbpt.get_ranges() << " invalid. ";
                throw std::runtime_error("\nValid options are: default,factor,explicit");
            }

            _mbpt.set_bse_nmax(options->get(key + ".exctotal").as<int> ());


            _mbpt.set_gwbasis_name(options->get(key + ".gwbasis").as<string> ());
            _mbpt.set_dftbasis_name(options->get(key + ".dftbasis").as<string> ());
            _mbpt.set_shift(options->get(key + ".shift").as<double> ());


            // possible tasks
            // diagQP, singlets, triplets, all
            string _tasks_string = options->get(key + ".tasks").as<string> ();
            if (_tasks_string.find("all") != std::string::npos) {
                _mbpt.set_do_qp_diag(true);
                _mbpt.set_do_bse_singlets(true);
                _mbpt.set_do_bse_triplets(true);
            }
            if (_tasks_string.find("qpdiag") != std::string::npos) _mbpt.set_do_qp_diag(true);
            if (_tasks_string.find("singlets") != std::string::npos) _mbpt.set_do_bse_singlets(true);
            if (_tasks_string.find("triplets") != std::string::npos) _mbpt.set_do_bse_triplets(true);

            // possible storage 
            // qpPert, qpdiag_energies, qp_diag_coefficients, bse_singlet_energies, bse_triplet_energies, bse_singlet_coefficients, bse_triplet_coefficients

            string _store_string = options->get(key + ".store").as<string> ();
            if ((_store_string.find("all") != std::string::npos) || (_store_string.find("") != std::string::npos)) {
                // store according to tasks choice
                if (_mbpt.get_do_qp_diag()) _mbpt.set_store_qp_diag(true);
                if (_mbpt.get_do_bse_singlets()) _mbpt.set_store_bse_singlets(true);
                if (_mbpt.get_do_bse_triplets()) _mbpt.set_store_bse_triplets(true);
            }
            if (_store_string.find("qpdiag") != std::string::npos) _mbpt.set_store_qp_diag(true);
            if (_store_string.find("singlets") != std::string::npos) _mbpt.set_store_bse_singlets(true);
            if (_store_string.find("triplets") != std::string::npos) _mbpt.set_store_bse_triplets(true);


            
            
            // something unique
           _package = options->get(key + ".package").as<string> ();
       
           // orbitals file or pure DFT output
           _orbfile  = options->get(key + ".molecule.orbitals").as<string> ();
           _logfile  = options->get(key + ".molecule.log").as<string> ();

            
    
  
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/qmpackages/") + _package + string("_idft_pair.xml");
    load_property_from_xml( _package_options, xmlFile );    

    // register all QM packages (Gaussian, TURBOMOLE, etc)
    QMPackageFactory::RegisterAll();
    
}

bool Exciton::Evaluate() {

    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
   _qmpackage->setLog( &_log );       
   _qmpackage->Initialize( &_package_options );

     Orbitals _orbitals;
  
    _qmpackage->setOrbitalsFileName( _orbfile );
    int _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitals );

    _qmpackage->setLogFileName( _logfile );
    int _parse_log_status = _qmpackage->ParseLogFile( &_orbitals );
    
//     MBPT _overlap; 
    _mbpt.setLogger(&_log);
     bool _evaluate = _mbpt.Evaluate( &_orbitals );
     std::cout << _log;
 
     
     // output the results
/*    Property _summary; 
    Property *_job_output = &_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    int HOMO_A = _orbitalsA.getNumberOfElectrons();
    int HOMO_B = _orbitalsB.getNumberOfElectrons();
    int LUMO_A = HOMO_A + 1;
    int LUMO_B = HOMO_B + 1;
    _pair_summary->setAttribute("homoA", HOMO_A);
    _pair_summary->setAttribute("homoB", HOMO_B);
    for (int levelA = HOMO_A - _levA +1; levelA <= LUMO_A + _levA - 1; ++levelA ) {
        for (int levelB = HOMO_B - _levB + 1; levelB <= LUMO_B + _levB -1 ; ++levelB ) {        
                double JAB = _overlap.getCouplingElement( levelA , levelB, &_orbitalsA, &_orbitalsB, &_JAB, _degeneracy );
                Property *_overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<string>(JAB)); 
                double energyA = _orbitalsA.getEnergy( levelA );
                double energyB = _orbitalsB.getEnergy( levelB );
                _overlap_summary->setAttribute("orbA", levelA);
                _overlap_summary->setAttribute("orbB", levelB);
                //_overlap_summary->setAttribute("jAB", JAB);
                _overlap_summary->setAttribute("eA", energyA);
                _overlap_summary->setAttribute("eB", energyB);
        }
    }
 */
    Property _summary; 
    Property *_job_output = &_summary.add("output","");
    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
     
    std::ofstream ofs (_output_file.c_str(), std::ofstream::out);
    ofs << *_job_output;    
    ofs.close();
    
    return true;
}


}}


#endif