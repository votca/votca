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

#include "votca/ctp/gaussian.h"
#include "votca/ctp/segment.h"
#include <stdio.h>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <sys/stat.h>

using namespace std;

namespace votca { namespace ctp {

Gaussian::Gaussian( tools::Property *opt ) { 
    
    string key = "package";

    string _name = opt->get(key+".name").as<string> ();
    
    if ( _name != "gaussian" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Package is not supported.");
    }
    
    _executable =       opt->get(key + ".executable").as<string> ();
    _functional =       opt->get(key + ".functional").as<string> ();
    _basis_set =        opt->get(key + ".basisset").as<string> ();
    _charge =           opt->get(key + ".charge").as<int> ();
    _spin =             opt->get(key + ".spin").as<int> ();
    _options =          opt->get(key + ".options").as<string> ();
    _memory =           opt->get(key + ".memory").as<string> ();
    _threads =          opt->get(key + ".threads").as<int> ();
    _checkpointfile =   opt->get(key + ".checkpoint").as<string> ();
        
};   
    
Gaussian::~Gaussian() { 
}  

/**
 * Prepares the com file from a vector of segments
 */
bool Gaussian::WriteInputFile( Segment *seg, string filename ) {

    vector< Atom* > ::iterator ait;

    int qmatoms = 0;

    vector< Atom* > _atoms;
    _atoms  = seg-> Atoms();

    ofstream _com_file;
    _com_file.open ( filename.c_str() );
    // header 
    if ( _checkpointfile.size() != 0 ) {
        _com_file << "%chk = " << _checkpointfile << endl;
    }
    
    if ( _memory.size() != 0 ) {
        _com_file << "%mem = " << _memory << endl ;
    }
    
    if ( _threads != 0 ) {
         _com_file << "%nprocshared = "  << _threads << endl;
    }
    
    if ( _functional.size() != 0 and _basis_set.size() != 0 ) {
        _com_file <<  "# " <<  _functional << "/" <<  _basis_set;
    }
    
    if ( _options.size() != 0 ) {
        _com_file <<  "  " << _options << endl ;
    }
    
    _com_file << endl << seg->getName() << endl << endl;
    _com_file << setw(2) << _charge << setw(2) << _spin << endl;
    
    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

        if ((*ait)->HasQMPart() == false) { continue; }

        vec     pos = (*ait)->getQMPos();
        string  name = (*ait)->getElement();

        //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
        _com_file << setw(3) << name.c_str() 
                  << setw(12) << setprecision(5) << pos.getX()*10
                  << setw(12) << setprecision(5) << pos.getY()*10
                  << setw(12) << setprecision(5) << pos.getZ()*10 
                  << endl;
    }
    
    _com_file << endl;
    _com_file.close();
}

/**
 * Runs the Gaussian job
 */
bool Gaussian::Run( string filename )
{
    //boost::filesystem::path _current_path = boost::filesystem::current_path();
    
    //boost::filesystem::path p( filename );
    //boost::filesystem::path dir = p.parent_path();
    //boost::filesystem::path file = p.filename();
            
    //chdir( (dir.string()).c_str() );
    //mkdir ( "temp", S_IRWXU|S_IRGRP|S_IXGRP ) ;

    //cout << endl << "... ... Running " << _executable << " package." << endl ;  
    //execlp( _executable.c_str(), _executable.c_str(), (file.string()).c_str(), NULL ); 
    //cout << filename << endl ;
    
    execlp( _executable.c_str(), _executable.c_str(), filename.c_str(), NULL );
    
    //execlp( _executable.c_str(), _executable.c_str(), filename.c_str(), NULL ); 
    cout << filename << endl ;
    //rmdir ( "temp" );
    
    //chdir( (_current_path.string()).c_str() );
}


}}
