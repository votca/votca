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


#ifndef _CALC_DFT_ENERGIES_H
#define	_CALC_DFT_ENERGIES_H

#include <votca/ctp/segment.h>
#include <votca/ctp/gaussian.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/parallelsitecalc.h>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace votca { namespace ctp {

/**
* \brief Site energies and orbitals for QM pairs
*
* QM orbitals and energies for all molecules
* Requires first-principles package, i.e. GAUSSIAN installation
*
* Callname: idft
*/

class EDFT : public ParallelSiteCalculator
{
public:

    EDFT() {};
   ~EDFT() {};

    string  Identify() { return "EDFT"; }
    void    Initialize(Topology *top, Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    void    EvalSite(Topology *top, Segment *seg, int slot);

    void    CleanUp();


private:

    bool   _maverick;
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;    

};

void EDFT::CleanUp() {

}

void EDFT::Initialize(Topology *top, Property *options) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;
    //cout << endl <<  "... ... Reading the input ";
    
    /* ---- OPTIONS.XML Structure -----
     *
     * <edft>
     *
     *      <package>gaussian.xml</package>
     *
     * </edft>
     *
     */

    this->ParseOrbitalsXML(top, options);

}


void EDFT::ParseOrbitalsXML(Topology *top, Property *opt) {

    string key = "options.edft";
    string _package_xml = opt->get(key+".package").as<string> ();
    cout << endl << "... ... Parsing " << _package_xml << endl ;

    load_property_from_xml( _package_options, _package_xml.c_str() );    
    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();
    
    //cout << endl << "... ... Using " << _package << " package" << endl ;    

}


void EDFT::EvalSite(Topology *top, Segment *seg, int slot) {


    FILE *out;
    Orbitals _orbitals;
    vector < Segment* > segments;
    
    segments.push_back( seg );
    
    _outParent = "frame" + boost::lexical_cast<string>(top->getDatabaseId());
    mkdir(_outParent.c_str(), 0755);

    string ID   = boost::lexical_cast<string>( seg->getId() );
    //string fileName = "mol_" + ID;
    string fileName = "monomer";
    string DIR  = _outParent + "/" + "mol_" + ID;
    string XYZ_FILE = fileName + ".xyz";
    string ORB_FILE = fileName + ".orb";
    string COM_FILE = fileName + ".com";
    string LOG_FILE = fileName + ".log"; 
    string SHL_FILE = fileName + ".sh";
    string GAUSSIAN_ORB_FILE = "fort.7" ;
        
    mkdir(DIR.c_str(), 0755);        
    out = fopen((DIR + "/" + XYZ_FILE).c_str(),"w");
    seg->WriteXYZ(out);
    fclose(out);
 
   if ( _package == "gaussian" ) { 
        
        Gaussian _gaussian( &_package_options );
               
        _gaussian.setRunDir( DIR );
        _gaussian.setInputFile( COM_FILE );

        // provide a separate scratch dir for every thread
        if ( ( _gaussian.getScratchDir() ).size() != 0 ) {
          _gaussian.setShellFile( SHL_FILE );
           string SCR_DIR  = _gaussian.getScratchDir() + "/mol_" + ID;
          _gaussian.setScratchDir( SCR_DIR );
          _gaussian.WriteShellScript ();
        } 
        
        _gaussian.WriteInputFile( segments );
        
        // Run the executable
        //_gaussian.Run( );
        
        // Collect information     
        _gaussian.setLogFile( DIR + "/" + LOG_FILE );
        _gaussian.ParseLogFile( &_orbitals );
        
        _gaussian.setOrbitalsFile( DIR + "/" + GAUSSIAN_ORB_FILE );
        _gaussian.ParseOrbitalsFile( &_orbitals );
 
        // save orbitals 
        std::ofstream ofs( (DIR + "/" + ORB_FILE).c_str() );
        boost::archive::binary_oarchive oa( ofs );
        oa << _orbitals;
        ofs.close();
       
        Orbitals _orbitals_new;
        std::ifstream ifs( (DIR +"/" + ORB_FILE).c_str() );
        boost::archive::binary_iarchive ia( ifs );
        ia >> _orbitals_new;
        ifs.close();
        exit(0);
        std::vector< QMAtom* >* _atoms = _orbitals_new.getAtoms();
        cout << (_atoms->front())->type;
        exit(0);
        
        //cout << "ELECTRONS " << _orbitals_new.getNumberOfElectrons() << endl ;
        //cout << "BASIS SIZE " << _orbitals_new.getBasisSetSize() << endl ;
 
        _gaussian.CleanUp( ID );

        //exit(0);
        
   }    
   
}

}}

#endif	/* _CALC_DFT_ENERGIES_H */
