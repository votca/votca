/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#include <string>
#include <xtp/moo_orbitals.h>
#include <moo/mol_and_orb.h>

using namespace std;

int main(int argc, char **argv){
    mol_and_orb mol;
    basis_set bs;
    orb orbital;
    if ( string(argv[1]).compare (  string("--help") ) == 0 ) {
	cout << "Usage: namebasissetXML file, name xyz coordinate, name orbital file (fort.7) nameoutput " <<endl;	    
	return 1;
    } 
    string fbasis =string(argv[1]);
    string fxyz   =string(argv[2]);
    string forb   =string(argv[3]);

    bs.set_basis_set(fbasis);
    mol.define_bs(bs);
    mol.init(fxyz.c_str());
    mol.init_orbitals(orbital, forb.c_str());


    ofstream out(argv[4]);
    for (int i =0; i< orbital.getNBasis(); i++){
    	for (int j=0; j< orbital.getNBasis(); j++) {
	    out << scientific << setprecision(8) << orbital.getorb(i)[j] << "  ";
	}
	out << '\n';
    }

    return 0;
}
