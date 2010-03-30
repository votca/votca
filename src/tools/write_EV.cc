#include <string>
#include <moo/orbitals.h>
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
	out << scientific << setprecision(8) << orbital.getevl()[i] << "  ";
	out << '\n';
    }

    return 0;
}
