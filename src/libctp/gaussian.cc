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

using namespace std;

namespace votca { namespace ctp {

Gaussian::Gaussian( tools::Property *opt ) { 
    
    string key = "package";

    string _name = opt->get(key+".name").as<string> ();
    
    if ( _name != "gaussian" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Package is not supported.");
    }
    
    _functional =       opt->get(key + ".functional").as<string> ();
    _basis_set =        opt->get(key + ".basisset").as<string> ();
    _options =          opt->get(key + ".options").as<string> ();
    _memory =           opt->get(key + ".memory").as<string> ();
    _threads =          opt->get(key + ".threads").as<string> ();
    _checkpointfile =   opt->get(key + ".checkpoint").as<string> ();
        
};   
    
Gaussian::~Gaussian() { 
}  

/**
 * Prepares the com file from a vector of segments
 */
bool Gaussian::WriteInputFile( Segment *seg, FILE *out ) {

    vector< Atom* > ::iterator ait;

    int qmatoms = 0;

    vector< Atom* > _atoms;
    _atoms  = seg-> Atoms();

    // header 
    out << "%chk=" << _checkpointfile <<
           "%mem=" << _memory <<
           "%nprocshared="  << _threads << 
           "# " <<  _functional << "/" <<  _basis_set << 
            << _options ;
    
    
    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

        if ((*ait)->HasQMPart() == false) { continue; }

        vec     pos = (*ait)->getQMPos();
        string  name = (*ait)->getElement();

        fprintf(out, "%2s %4.7f %4.7f %4.7f \n",
                        name.c_str(),
                        pos.getX()*10,
                        pos.getY()*10,
                        pos.getZ()*10);
    }
}

/**
 * Runs the Gaussian job
 */
bool Gaussian::Run( )
{
}


}}
