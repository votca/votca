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

#ifndef _VOTCA_CTP_ESPFIT_TOOL_H
#define _VOTCA_CTP_ESPFIT_TOOL_H

#include <stdio.h>
#include <votca/ctp/espfit.h>
#include <votca/ctp/logger.h>
#include <boost/filesystem.hpp>
#include <votca/ctp/qmmachine.h>

namespace votca { namespace ctp {
    using namespace std;
    
class ESPFit_Tool : public QMTool
{
public:

    ESPFit_Tool () { };
   ~ESPFit_Tool () { };

    string Identify() { return "ESPFit"; }

    void   Initialize(Property *options);
    bool   Evaluate();


private:
    
    string      _orbfile;
    string      _output_file;
    int         _state_no;   
    string      _state;
    bool        _use_mps;
    bool        _use_pdb;
    
    Logger      _log;
    
    void FitESP( Orbitals& _orbitals );

};

void ESPFit_Tool::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options );
 

            string key = "options." + Identify();
            // _jobfile = options->get(key + ".file").as<string>();

            // key = "options." + Identify();
 
       
           // orbitals file or pure DFT output
           _orbfile      = options->get(key + ".input").as<string> ();
           _state     = options->get(key + ".state").as<string> ();
           _output_file  = options->get(key + ".output").as<string> ();
           _state_no     = options->get(key + ".state number").as<int> ();
           data_format  = boost::filesystem::extension(_output_file);
           
    if (data_format=="mps"){
        _use_mps=true; 
        _use_pdb=false;
    }
    else if(data_format=="pdb"){
       _use_mps=false; 
        _use_pdb=true;
    }
    else  throw std::runtime_error("Outputfile format not recognized. Export only to .pdb and .mps");
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    // string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/qmpackages/") + _package + string("_idft_pair.xml");
    // load_property_from_xml( _package_options, xmlFile );    

    // register all QM packages (Gaussian, TURBOMOLE, etc)
    // QMPackageFactory::RegisterAll();
   
}

bool ESPFit_Tool::Evaluate() {

    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    LOG(logDEBUG, _log) << "Converting serialized QM data in " << _orbfile << flush;

    Orbitals _orbitals;
    // load the QM data from serialized orbitals object

    std::ifstream ifs( (_orbfile).c_str());
    LOG(logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
    boost::archive::binary_iarchive ia(ifs);
    ia >> _orbitals;
    ifs.close();

    FitESP(_orbitals);

    if(_use_mps){
        QMMInterface Converter;
        PolarSeg* result=Converter.Convert(_orbitals.QMAtoms(););
        result->WriteMPS(_output_file);
        }
    else if(use_pdb){
        FILE *out;
        out = fopen(_output_file.c_str(), "w");
       _orbitals.WritePDB(out); 
    }
    
    
    LOG(logDEBUG, _log) << "Written charges to " << _output_file << flush;
    
    return true;
}




void ESPFit_Tool::FitESP( Orbitals& _orbitals ){


   
    LOG(logDEBUG, _log) << "===== Running ESPFIT ===== " << flush;

        vector< QMAtoms* >& Atomlist =_orbitals.QMAtoms();
        
        BasisSet dftbs;
        dftbs.LoadBasisSet(_orbitals.getDFTbasis());
        AOBasis basis;
        basis.AOBasisFill(&dftbs, Atomlist );
        basis.ReorderMOs(_orbitals.MOCoefficients(), _orbitals.getQMpackage(), "votca" );
        ub::matrix<double> &DMATGS=_orbitals.DensityMatrixGroundState();
        ub::matrix<double> DMAT_tot=DMATGS;
        if ( _state > 0 ){ 
	  std::vector<ub::matrix<double> > &DMAT = _orbitals.DensityMatrixExcitedState( _orbitals.MOCoefficients() , _orbitals.BSESingletCoefficients(), _state-1);
	  DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1]; // Ground state + hole_contribution + electron contribution
	}
        Espfit esp;
        esp.setLog(_log);
        esp.Fit2Density(Atomlist, DMAT_tot, basis);   
}

}}


#endif