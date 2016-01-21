/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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
#include <votca/xtp/qmpackagefactory.h>
#include<votca/xtp/aobasis.h>
#include<votca/xtp/aomatrix.h>
#include<boost/numeric/ublas/matrix.hpp>

namespace votca { namespace xtp {
    using namespace std;
    
class QMSandbox : public QMTool
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
    
    Logger      _log;
 
    string      _logfile;

    string      _package;
    Property    _package_options; 
    
    void CheckContent(  Orbitals& _orbitals );

};

void QMSandbox::Initialize(Property* options) {

    // update options with the VOTCASHARE defaults   
    //UpdateWithDefaults( options );
 

    string key = "options." + Identify();

    // orbitals file or pure DFT output
    _logfile  = options->get(key + ".log").as<string> ();
    _package  = options->get(key + ".package").as<string> ();
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    // string xmlFile = string(getenv("VOTCASHARE")) + string("/xtp/qmpackages/") + _package + string("_idft_pair.xml");
    // load_property_from_xml( _package_options, xmlFile );    

    // register all QM packages (Gaussian, TURBOMOLE, etc)
    // QMPackageFactory::RegisterAll();
    
    // register all QM packages (Gaussian, TURBOMOLE, etc)
    QMPackageFactory::RegisterAll();
    
    
}

bool QMSandbox::Evaluate() {

    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    LOG(logDEBUG, _log) << "Analyzing serialized QM data in " << _orbfile << flush;

    Orbitals _orbitals;
    // load the QM data from serialized orbitals object

    //std::ifstream ifs( (_orbfile).c_str());
    //LOG(logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    //boost::archive::binary_iarchive ia(ifs);
    //ia >> _orbitals;
    //ifs.close();
    
    // output info about contents of serialized data
    //CheckContent( _orbitals );

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
   _qmpackage->setLog( &_log );     
    _qmpackage->setLogFileName( _logfile );
    //int _parse_log_status = _qmpackage->ParseLogFile( &_orbitals );
   _qmpackage->ParseLogFile( &_orbitals );
   //_qmpackage->Initialize( &_package_options );

    //_qmpackage->setOrbitalsFileName( _orbfile );
    //int _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitals );

                std::vector<QMAtom*> _atoms = _orbitals.QMAtoms();
            
            // load DFT basis set (element-wise information) from xml file
            BasisSet dftbs;
            dftbs.LoadBasisSet("ubecppol");
            _orbitals.setDFTbasis( "ubecppol" );
            LOG(logDEBUG, _log) << TimeStamp() << " Loaded DFT Basis Set "  << flush;

            // fill DFT AO basis by going through all atoms 
            AOBasis dftbasis;
            dftbasis.AOBasisFill(&dftbs, _atoms);
            LOG(logDEBUG, _log) << TimeStamp() << " Filled DFT Basis of size " << dftbasis._AOBasisSize << flush;

            // fill auxiliary DFT AO basis by going through all atoms
            AOOverlap _dftoverlap;
            // initialize overlap matrix
            _dftoverlap.Initialize(dftbasis._AOBasisSize);
            // Fill overlap
            _dftoverlap.Fill(&dftbasis);
            LOG(logDEBUG, _log) << TimeStamp() << " Filled AO Overlap matrix of dimension: " << _dftoverlap._aomatrix.size1() << flush;
            _dftoverlap.Print("Denis does not like this ");    
    
    
            boost::numeric::ublas::matrix<AOBasis*> n;
    // LOG(logDEBUG, _log) << "Written text data to " << _output_file << flush;
    
    
    return true;
}




void QMSandbox::CheckContent( Orbitals& _orbitals ){


   
    LOG(logDEBUG, _log) << "===== Summary of serialized content ===== " << flush;
    LOG(logDEBUG, _log) << "   Information about DFT:" << flush;
    
          
    // DFT atoms
    if ( _orbitals.hasQMAtoms() ) {
        LOG(logDEBUG, _log) << "      atoms:                  " << _orbitals.QMAtoms().size() << flush;
    } else {
        LOG(logDEBUG, _log) << "      atoms:                  not stored "<< flush;
    } 
    
    // QM package
    if ( _orbitals.hasQMpackage() ){
        LOG(logDEBUG, _log) << "      QM package:             " << _orbitals.getQMpackage() << flush;
    
    } else {
        LOG(logDEBUG, _log) << "      QM package:             not stored " << flush;
    }
    
    
    // DFT basis set
    if ( _orbitals.hasDFTbasis() ) {
        LOG(logDEBUG, _log) << "      basis set:              " << _orbitals.getDFTbasis() << flush;
    } else {
        LOG(logDEBUG, _log) << "      basis set:              not stored "<< flush;
    }

    // DFT basis set size
    if ( _orbitals.hasBasisSetSize() ) {
        LOG(logDEBUG, _log) << "      basis set size:         " << _orbitals.getBasisSetSize() << flush;
   
    } else {
        LOG(logDEBUG, _log) << "      basis set size:          not stored "<< flush;
    }

    // DFT number of electrons
    if ( _orbitals.hasNumberOfElectrons() ) {
        LOG(logDEBUG, _log) << "      number of electrons:    " << _orbitals.getNumberOfElectrons() << flush;
    } else {
         LOG(logDEBUG, _log) << "      number of electrons:    not stored "<< flush;
    }    
    
    // DFT number of levels
    if ( _orbitals.hasNumberOfLevels() ) {
        LOG(logDEBUG, _log) << "      number of levels:       " << _orbitals.getNumberOfLevels() << flush;
    } else {
         LOG(logDEBUG, _log) << "      number of levels:       not stored "<< flush;
    }    

    // DFT orbital energies
    if ( _orbitals.hasMOEnergies() ) {
        LOG(logDEBUG, _log) << "      MO energies:            " << _orbitals.getEnergies()->size() << flush;
    } else {
         LOG(logDEBUG, _log) << "      MO energies:            not stored "<< flush;
    }    

    // DFT orbital coefficients
    if ( _orbitals.hasMOCoefficients() ) {
        LOG(logDEBUG, _log) << "      MO coefficients:        " << _orbitals.MOCoefficients().size1() << " x " << _orbitals.MOCoefficients().size2() << flush;
    } else {
        LOG(logDEBUG, _log) << "      MO coefficients:        not stored "<< flush;
    }  
    
    // DFT AO overlap matrix
    if ( _orbitals.hasAOOverlap() ) {
        LOG(logDEBUG, _log) << "      AO overlap matrix:      " << _orbitals.getOverlap()->size1()  << " x " << _orbitals.getOverlap()->size2() << flush;
    } else {
        LOG(logDEBUG, _log) << "      AO overlap matrix:      not stored "<< flush;
    }    
    
    // DFT AO XC matrix
    if ( _orbitals.hasAOVxc() ) {
        LOG(logDEBUG, _log) << "      AO XC matrix:           " << _orbitals.AOVxc().size1()  << " x " << _orbitals.AOVxc().size2() << flush;
    } else {
        LOG(logDEBUG, _log) << "      AO XC matrix:           not stored "<< flush;
    }    

    // QM total energy
    if ( _orbitals.hasQMEnergy() ){
        LOG(logDEBUG, _log) << "      QM energy:              " << _orbitals.getQMEnergy() << flush;
    } else{
        LOG(logDEBUG, _log) << "      QM energy:              not stored " << flush;
    }
    
    // MM self-energy 
    if ( _orbitals.hasSelfEnergy() ){
        LOG(logDEBUG, _log) << "      MM self energy:         " << _orbitals.getSelfEnergy() << flush;
    } else{
        LOG(logDEBUG, _log) << "      MM self energy:         not stored " << flush;
    }    
    
    // DFT transfer integrals
    if ( _orbitals.hasMOCouplings() ) {
        LOG(logDEBUG, _log) << "      DFT transfer integrals: " << _orbitals.MOCouplings().size1() << " x " << _orbitals.MOCouplings().size2() << flush;

    } else {
         LOG(logDEBUG, _log) << "      DFT transfer integrals: not stored "<< flush;
    }    
    
    
    

    LOG(logDEBUG, _log) << "   Information about GWA:" << flush;
    
    // GWA basis set
    if ( _orbitals.hasGWbasis() ) {
        LOG(logDEBUG, _log) << "      basis set:              " << _orbitals.getGWbasis() << flush;
    } else {
        LOG(logDEBUG, _log) << "      basis set:              not stored "<< flush;
    }
    
    // RPA index range
    if ( _orbitals.hasRPAindices() ){
        LOG(logDEBUG, _log) << "      RPA level range:        " << _orbitals.getRPAmin() << " : " << _orbitals.getRPAmax() << flush;        
    } else {
        LOG(logDEBUG, _log) << "      RPA level range:        not stored" << flush;
    }

    // GWA index range
    if ( _orbitals.hasGWAindices() ){
        LOG(logDEBUG, _log) << "      GWA level range:        " << _orbitals.getGWAmin() << " : " << _orbitals.getGWAmax() << flush;        
    } else {
        LOG(logDEBUG, _log) << "      GWA level range:        not stored" << flush;
    }
    
    // perturbative QP energies
    if ( _orbitals.hasQPpert()){
        LOG(logDEBUG, _log) << "      number of QP levels:    " << _orbitals.QPpertEnergies().size1() << flush;
    } else {
        LOG(logDEBUG, _log) << "      number of QP levels:    not stored" << flush;
    }
    
    // diagonalized QP energies
    if ( _orbitals.hasQPdiag() ){
        LOG(logDEBUG, _log) << "      diagonalized QP levels: " << _orbitals.QPdiagEnergies().size() << flush;
    } else {
        LOG(logDEBUG, _log) << "      diagonalized QP levels: not stored" << flush;
    }


    LOG(logDEBUG, _log) << "   Information about BSE:" << flush;

    // BSE index range
    if ( _orbitals.hasBSEindices() ){
        LOG(logDEBUG, _log) << "      BSE level range:        [" << _orbitals.getBSEvmin() << " : " << _orbitals.getBSEvmax()  << "] x [ " << _orbitals.getBSEcmin() << " : " << _orbitals.getBSEcmax() << "]" << flush;        
    } else {
        LOG(logDEBUG, _log) << "      BSE level range:        not stored" << flush;
    }
    
    // BSE EH interaction
    if ( _orbitals.hasEHinteraction() ){
        LOG(logDEBUG, _log) << "      direct interaction:     " << _orbitals.eh_d().size1() << flush;
        LOG(logDEBUG, _log) << "      exchange interaction:   " << _orbitals.eh_x().size1() << flush;
    } else {
        LOG(logDEBUG, _log) << "      e-h interactions:       not stored" << flush;
    }
    
    // BSE singlet excitons
    if ( _orbitals.hasBSESinglets()){
        LOG(logDEBUG, _log) << "      BSE singlet excitons:   " << _orbitals.getBSESingletEnergies()->size() << flush;
    } else {
        LOG(logDEBUG, _log) << "      BSE singlet excitons:   not stored" << flush;
    }  
    
    // BSE triplet excitons
    if ( _orbitals.hasBSETriplets()){
        LOG(logDEBUG, _log) << "      BSE triplet excitons:   " << _orbitals.getBSETripletEnergies()->size() << flush;
    } else {
        LOG(logDEBUG, _log) << "      BSE triplet excitons:   not stored" << flush;
    }  

    
    
    
}

}}


#endif