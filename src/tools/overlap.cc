//
// File:   overlap_integral.cc
// Author: James L
// calculates overlap integrals for molecular orbitals using libmoo
// usage:
//       overlap_integral --coord1 <coordinate file for mol1>
//			  --coord2 <coordinate file for mol2>
//                        --listcrg <xml file containing the definition for the charge unit type>
//                        
//

#include <boost/program_options.hpp>
#include <string>
#include "jcalc.h"
#include "crgunit.h"

using namespace std;

int main(int argc, char **argv) {
    // catches the exceptions if the program does not behave 
    // feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);

    // parser for the program arguments
    namespace po = boost::program_options;
    // declare the supported options
    po::options_description desc("Allowed options");
    string listcrg, pos1, pos2;
    desc.add_options()
            ("help,h", "produce help message")
            ("listcharges,l", po::value<string >(&listcrg)-> default_value("list_charges.xml"), "xml filename containing the list of the charge unit type")
            ("posor1,1", po::value<string > (&pos1) -> default_value("posmol1"), "list of charge unit type, position and orientation for mol 1")
	    ("posor2,2", po::value<string > (&pos2) -> default_value("posmol2"), " list of charge unit type, position and orientation for mol 2")
            ;
    

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } 
    catch (po::error &err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }

    try {
    cout << "Reading Crg unit Types"<<endl;
    JCalc jcalc(listcrg);
    cout << "Finished reading Crg unit Types" <<endl;

    ifstream in1(pos1.c_str());
    ifstream in2(pos2.c_str());
    while ( in1 && in2){
    	string name1, name2;
	vec com1;
	vec com2;
	matrix or1;
	matrix or2;
	in1 >> name1 >> com1.x() >> com1.y() >> com1.z() >> 
		or1[0][0] >> or1[0][1] >> or1[0][2] >>
		or1[1][0] >> or1[1][1] >> or1[1][2] >>
		or1[2][0] >> or1[2][1] >> or1[2][2] ;
	in2 >> name2 >> com2.x() >> com2.y() >> com2.z() >> 
		or2[0][0] >> or2[0][1] >> or2[0][2] >>
		or2[1][0] >> or2[1][1] >> or2[1][2] >>
		or2[2][0] >> or2[2][1] >> or2[2][2] ;

	// cout << "mol1: " << name1 << com1 << or1<<endl;
	// cout << "mol2: " << name2 << com2 << or2<<endl;


	CrgUnit * A = jcalc.CreateCrgUnit(0, name1);
	CrgUnit * B = jcalc.CreateCrgUnit(1, name2);
	A->SetPos(0, com1);
	B->SetPos(0, com2);
	A->SetNorm(0, or1[2]);
	A->SetPlane(0, or1[1]);
	B->SetNorm(0, or2[2]);
	B->SetPlane(0, or2[1]);

//	cout << "Compute J" <<endl;
	vector <double> Js = jcalc.CalcJ(*A,*B);
//	cout << "Finished computing J" <<endl;
	vector <double>::iterator itJ= Js.begin();
	for (; itJ!=Js.end(); ++itJ) cout << '\t'<< *itJ <<endl;

    }
    }
catch (std::exception &err) {
        cout << "error : " << err.what() << endl;
        return -1;
    }


}

