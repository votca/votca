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

#ifndef _VOTCA_CTP_COUPLINGH_H
#define _VOTCA_CTP_COUPLINGH_H

#include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>

namespace votca { namespace ctp {
    using namespace std;
    
class Coupling : public QMTool
{
public:

    Coupling() { };
   ~Coupling() { };

    string Identify() { return "coupling"; }

    void   Initialize(Property *options);
    bool   Evaluate();

 

private:
    
    string      _orbA, _orbB, _orbAB;
    string      _logA, _logB, _logAB;
    int         _levA, _levB;
    int         _trimA, _trimB;
    string      _qmpackageXML;
    double      _degeneracy;

    string              _package;
    Property            _package_options; 
    
    Logger      _log;

};

void Coupling::Initialize(Property* options) 
{
    
   // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options );
    std::string key = "options." + Identify();    

    _degeneracy = options->get(key + ".degeneracy").as<double> ();
    _qmpackageXML = options->get(key + ".package").as<string> ();
            
    _orbA  = options->get(key + ".moleculeA.orbitals").as<string> ();
    _orbB  = options->get(key + ".moleculeB.orbitals").as<string> ();
    _orbAB = options->get(key + ".dimerAB.orbitals").as<string> ();

    _logA  = options->get(key + ".moleculeA.log").as<string> ();
    _logB  = options->get(key + ".moleculeB.log").as<string> ();
    _logAB = options->get(key + ".dimerAB.log").as<string> ();

    _levA  = options->get(key + ".moleculeA.levels").as<int> ();
    _levB  = options->get(key + ".moleculeB.levels").as<int> ();
 
    _trimA  = options->get(key + ".moleculeA.trim").as<int> ();
    _trimB  = options->get(key + ".moleculeB.trim").as<int> ();
   
     //key = "package";
    //_package = _package_options.get(key+".name").as<string> ();

    _package = "gaussian";
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/qmpackages/") + _package + string("/idft.pair.xml ");
    //load_property_from_xml( _package_options, xmlFile );    
    
}

bool Coupling::Evaluate() {

    LOG( logINFO, _log ) << "Loading orbitals " << _orbA << flush;
    LOG( logINFO, _log ) << "Loading orbitals " << _orbB << flush;
    LOG( logINFO, _log ) << "Loading orbitals " << _orbAB << flush;
   
    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
   _qmpackage->setLog( &_log );       
   //_qmpackage->Initialize( &_package_options );
   
    Orbitals _orbitalsA, _orbitalsB, _orbitalsAB;
    
    int _parse_logA_status = _qmpackage->ParseLogFile( &_orbitalsA );
    int _parse_logB_status = _qmpackage->ParseLogFile( &_orbitalsB );
    int _parse_logAB_status = _qmpackage->ParseLogFile( &_orbitalsAB );

    int _parse_orbitalsA_status = _qmpackage->ParseOrbitalsFile( &_orbitalsA );
    int _parse_orbitalsB_status = _qmpackage->ParseOrbitalsFile( &_orbitalsB );
    int _parse_orbitalsAB_status = _qmpackage->ParseOrbitalsFile( &_orbitalsAB );

    
    std::cout << _log;
    
    return true;
}


}}


#endif