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

#ifndef _VOTCA_XTP_TOOLS_COUPLINGH_H
#define _VOTCA_XTP_TOOLS_COUPLINGH_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/xtp/dftcoupling.h>
#include <votca/xtp/qmpackagefactory.h>

namespace votca { namespace xtp {
    
class Coupling : public ctp::QMTool
{
public:

    Coupling() { };
   ~Coupling() { };

    std::string Identify() { return "coupling"; }

    void   Initialize(tools::Property *options);
    bool   Evaluate();

 

private:
    
    std::string      _orbA, _orbB, _orbAB;
    std::string      _logA, _logB, _logAB;
    int         _levA, _levB;
    int         _trimA, _trimB;
    double      _degeneracy;

    std::string      _package;
    tools::Property    _package_options; 
    
    std::string      _output_file;
    
    ctp::Logger      _log;

};

void Coupling::Initialize(tools::Property* options) 
{

   // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "ctp" );
    std::string key = "options." + Identify();    
    
    _degeneracy = options->get(key + ".degeneracy").as<double> ();
    
            
    _orbA  = options->get(key + ".moleculeA.orbitals").as<std::string> ();
    _orbB  = options->get(key + ".moleculeB.orbitals").as<std::string> ();
    _orbAB = options->get(key + ".dimerAB.orbitals").as<std::string> ();

    _logA  = options->get(key + ".moleculeA.log").as<std::string> ();
    _logB  = options->get(key + ".moleculeB.log").as<std::string> ();
    _logAB = options->get(key + ".dimerAB.log").as<std::string> ();

    _levA  = options->get(key + ".moleculeA.levels").as<int> ();
    _levB  = options->get(key + ".moleculeB.levels").as<int> ();
 
    _trimA  = options->get(key + ".moleculeA.trim").as<int> ();
    _trimB  = options->get(key + ".moleculeB.trim").as<int> ();
   
    _output_file = options->get(key + ".output").as<std::string> ();

     string _package_xml = options->get(key + ".dftpackage").as<string> ();
    load_property_from_xml(_package_options, _package_xml.c_str());
     _package = _package_options.get("package.name").as<string> ();
    

    // register all QM packages (Gaussian, TURBOMOLE, etc)
    xtp::QMPackageFactory::RegisterAll();
    
}

bool Coupling::Evaluate() {

    _log.setReportLevel( ctp::logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ..."); 

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
   _qmpackage->setLog( &_log );       
   _qmpackage->Initialize( &_package_options );
    _qmpackage->setRunDir(".");
     Orbitals _orbitalsA, _orbitalsB, _orbitalsAB;
     

   
    _qmpackage->setLogFileName( _logA );
    bool _parse_logA_status = _qmpackage->ParseLogFile( _orbitalsA );
    if ( !_parse_logA_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read log of molecule A" << std::flush; }
    
    _qmpackage->setLogFileName( _logB );
    bool _parse_logB_status = _qmpackage->ParseLogFile( _orbitalsB );
    if ( !_parse_logB_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read log of molecule B" << std::flush; }
    
    _qmpackage->setLogFileName( _logAB );
    bool _parse_logAB_status = _qmpackage->ParseLogFile( _orbitalsAB );
    if ( !_parse_logAB_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read log of molecule AB" << std::flush; }
    
        _qmpackage->setOrbitalsFileName( _orbA );
    bool _parse_orbitalsA_status = _qmpackage->ParseOrbitalsFile( _orbitalsA );
    if ( !_parse_orbitalsA_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read orbitals of molecule A" << std::flush; }

    _qmpackage->setOrbitalsFileName( _orbB );   
    bool _parse_orbitalsB_status = _qmpackage->ParseOrbitalsFile( _orbitalsB );
    if ( !_parse_orbitalsB_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read orbitals of molecule B" << std::flush; }
    
    _qmpackage->setOrbitalsFileName( _orbAB );   
    bool _parse_orbitalsAB_status = _qmpackage->ParseOrbitalsFile( _orbitalsAB );
     if ( !_parse_orbitalsAB_status ) { CTP_LOG(ctp::logERROR,_log) << "Failed to read orbitals of dimer AB" << std::flush; }

    int _degAH = 1;
    int _degAL = 1;
    int _degBH = 1;
    int _degBL = 1;
    
    // trim monomers A and B to one level 
    if ((_trimA < 0) || (_trimB <0) ) { // any -1 overrides the specification of the other 
        
        // find degeneracy of HOMOs and LUMOs
        std::vector<int> list_levelsAH  = (_orbitalsA.getDegeneracy( _orbitalsA.getNumberOfElectrons()-1, _degeneracy ));
        _degAH = list_levelsAH.size();
        std::vector<int> list_levelsAL  = (_orbitalsA.getDegeneracy( _orbitalsA.getNumberOfElectrons(), _degeneracy ));
        _degAL = list_levelsAL.size();  
        
        std::vector<int> list_levelsBH  = (_orbitalsB.getDegeneracy( _orbitalsB.getNumberOfElectrons()-1, _degeneracy ));
        _degBH = list_levelsBH.size();
        std::vector<int> list_levelsBL  = (_orbitalsB.getDegeneracy( _orbitalsB.getNumberOfElectrons(), _degeneracy ));
        _degBL = list_levelsBL.size();  
        
        _orbitalsA.Trim(_degAH,_degAL);
        _orbitalsB.Trim(_degBH,_degBL);
    
    // trim by the factors   (_trimA-1)  and  (_trimB-1) 
    } else {
 
        if ( _orbitalsA.getNumberOfElectrons()*(_trimA-1) <  int(_orbitalsA.getNumberOfLevels()) - _orbitalsA.getNumberOfElectrons() ) {
            CTP_LOG(ctp::logDEBUG,_log) << "Trimming virtual orbitals A:" 
                    << _orbitalsA.getNumberOfLevels() - _orbitalsA.getNumberOfElectrons() << "->" 
                    << _orbitalsA.getNumberOfElectrons()*(_trimA-1) << std::flush;  
            _orbitalsA.Trim(_trimA);
        }
    
        if ( _orbitalsB.getNumberOfElectrons()*(_trimB-1) <   int(_orbitalsB.getNumberOfLevels()) - _orbitalsB.getNumberOfElectrons() ) {
            CTP_LOG(ctp::logDEBUG,_log) << "Trimming virtual orbitals B:" 
                    << _orbitalsB.getNumberOfLevels() - _orbitalsB.getNumberOfElectrons() << "->" 
                    << _orbitalsB.getNumberOfElectrons()*(_trimB-1) << std::flush;      
            _orbitalsB.Trim(_trimB);
        }
    
    }
    
     DFTcoupling dftcoupling; 
    dftcoupling.setLogger(&_log);
          
    Eigen::MatrixXd _JAB = dftcoupling.CalculateIntegrals( _orbitalsA, _orbitalsB, _orbitalsAB);  
    std::cout << _log;
 
     
     // output the results
    tools::Property _summary; 
    tools::Property *_job_output = &_summary.add("output","");
    tools::Property *_pair_summary = &_job_output->add("pair","");
    int HOMO_A = _orbitalsA.getNumberOfElectrons();
    int HOMO_B = _orbitalsB.getNumberOfElectrons();
    int LUMO_A = HOMO_A + 1;
    int LUMO_B = HOMO_B + 1;
    _pair_summary->setAttribute("homoA", HOMO_A);
    _pair_summary->setAttribute("homoB", HOMO_B);

    if ( (_trimA <0) || (_trimB <0) ) {

        // HOMO-HOMO coupling
        double JAB = dftcoupling.getCouplingElement(_degAH, _degBH , _orbitalsA, _orbitalsB, &_JAB, _degeneracy);
        tools::Property *_overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<std::string>(JAB));
        double energyA = _orbitalsA.getEnergy(_degAH);
        double energyB = _orbitalsB.getEnergy(_degBH);
        _overlap_summary->setAttribute("orbA", HOMO_A);
        _overlap_summary->setAttribute("orbB", HOMO_B);
        //_overlap_summary->setAttribute("jAB", JAB);
        _overlap_summary->setAttribute("eA", energyA);
        _overlap_summary->setAttribute("eB", energyB);
                
        // LUMO-LUMO coupling
        JAB = dftcoupling.getCouplingElement(_degAH+1, _degBH+1 , _orbitalsA, _orbitalsB, &_JAB, _degeneracy);
        _overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<std::string>(JAB));
        energyA = _orbitalsA.getEnergy(_degAH +1);
        energyB = _orbitalsB.getEnergy(_degBH +1);
        _overlap_summary->setAttribute("orbA", LUMO_A);
        _overlap_summary->setAttribute("orbB", LUMO_B);
        //_overlap_summary->setAttribute("jAB", JAB);
        _overlap_summary->setAttribute("eA", energyA);
        _overlap_summary->setAttribute("eB", energyB);                

    } else {
    
        for (int levelA = HOMO_A - _levA +1; levelA <= LUMO_A + _levA - 1; ++levelA ) {
            for (int levelB = HOMO_B - _levB + 1; levelB <= LUMO_B + _levB -1 ; ++levelB ) {        
                double JAB = dftcoupling.getCouplingElement( levelA , levelB, _orbitalsA, _orbitalsB, &_JAB, _degeneracy );
                tools::Property *_overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<std::string>(JAB)); 
                double energyA = _orbitalsA.getEnergy( levelA );
                double energyB = _orbitalsB.getEnergy( levelB );
                _overlap_summary->setAttribute("orbA", levelA);
                _overlap_summary->setAttribute("orbB", levelB);
                //_overlap_summary->setAttribute("jAB", JAB);
                _overlap_summary->setAttribute("eA", energyA);
                _overlap_summary->setAttribute("eB", energyB);
            }
        }
    }
    
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
     
    std::ofstream ofs (_output_file.c_str(), std::ofstream::out);
    ofs << *_job_output;    
    ofs.close();
    
    return true;
}


}}


#endif
