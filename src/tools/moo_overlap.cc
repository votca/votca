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
#include <votca/moo/jcalc.h>
#include <votca/moo/crgunit.h>
#include <votca/tools/application.h>

using namespace std;

class MOOApplication
	: public Application
{
public:
	    virtual string ProgramName() { return "moo_overlap"; }
	    virtual void HelpText(std::ostream &out) { };
	    void Initialize();
	    bool EvaluateOptions();
	    void Run(void);
protected: 		
	    string listcrg, pos1, pos2, pdbfile;
};

void MOOApplication::Initialize()
{
	AddProgramOptions("Options")
		("listcharges", boost::program_options::value<string>(), "xml filename with charge unit types")
		("mol1", boost::program_options::value<string>() -> default_value("mol1"), "position and orientation for mol 1")
		("mol2", boost::program_options::value<string>() -> default_value("mol2"), "position and orientation for mol 2")		
		("pdb", boost::program_options::value<string>() -> default_value("geometry.pdb"), "pdb file of two molecules")		
		;
}

bool MOOApplication::EvaluateOptions()
{
	CheckRequired("listcharges");
	//_runtime = OptionsMap()["time"].as<double>();
	return true;
}

void MOOApplication::Run(void)
{
        cout << "Reading Crg unit Types" << endl;
        JCalc jcalc(listcrg);
        cout << "Finished reading Crg unit Types" << endl;

        ifstream in1(pos1.c_str());
        ifstream in2(pos2.c_str());
        int written=0;
        while (in1 && in2) {
            string name1, name2;
            vec com1;
            vec com2;
            matrix or1;
            matrix or2;
            in1 >> name1 >> com1.x() >> com1.y() >> com1.z() >>
                    or1[0][0] >> or1[0][1] >> or1[0][2] >>
                    or1[1][0] >> or1[1][1] >> or1[1][2] >>
                    or1[2][0] >> or1[2][1] >> or1[2][2];
            in2 >> name2 >> com2.x() >> com2.y() >> com2.z() >>
                    or2[0][0] >> or2[0][1] >> or2[0][2] >>
                    or2[1][0] >> or2[1][1] >> or2[1][2] >>
                    or2[2][0] >> or2[2][1] >> or2[2][2];

            if (!in1 || !in2) break;
            //	cout << "mol1: " << name1 << com1 << or1<<endl;
            //	cout << "mol2: " << name2 << com2 << or2<<endl;


            CrgUnit * A = jcalc.CreateCrgUnit(0, name1);
            CrgUnit * B = jcalc.CreateCrgUnit(1, name2);
            A->SetPos(0, com1);
            B->SetPos(0, com2);
            A->SetNorm(0, or1[2]);
            A->SetPlane(0, or1[1]);
            B->SetNorm(0, or2[2]);
            B->SetPlane(0, or2[1]);


           //write pdb file
           mol_and_orb *molecule = ( A -> rotate_translate_beads() );
           (*molecule).write_pdb(pdbfile, "m1", written);
           written += (*molecule).getN();
           delete molecule;
           molecule = ( B -> rotate_translate_beads() );
           (*molecule).write_pdb(pdbfile, "m2", written);
           written += (*molecule).getN();
           delete molecule;
           ofstream fl;
           fl.open(pdbfile.c_str(), ios::app);
           fl.setf(ios::fixed);
           fl << "END" <<endl;

           
            //	cout << "Compute J" <<endl;
            vector <double> Js = jcalc.CalcJ(*A, *B);
            //	cout << "Finished computing J" <<endl;
            vector <double>::iterator itJ = Js.begin();
            for (; itJ != Js.end(); ++itJ) cout << '\t' << *itJ << endl;
	};
}

int main(int argc, char **argv)
{
	MOOApplication app;
	return app.Exec(argc, argv);
}
 

/*
int main(int argc, char **argv) {
    
    // catches the exceptions if the program does not behave
    // feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);

    // parser for the program arguments
    namespace po = boost::program_options;
    // declare the supported options
    po::options_description desc("Allowed options");
    string listcrg, pos1, pos2, pdbfile;
    desc.add_options()
            ("help,h", "produce help message")
            ("listcharges,l", po::value<string > (&listcrg)-> default_value("list_charges.xml"), "xml filename containing the list of the charge unit type")
            ("posor1,1", po::value<string > (&pos1) -> default_value("posmol1"), "list of charge unit type, position and orientation for mol 1")
            ("posor2,2", po::value<string > (&pos2) -> default_value("posmol2"), "list of charge unit type, position and orientation for mol 2")
            ("pdb,p", po::value<string > (&pdbfile) -> default_value("posmol.pdb"), "write pdb file with used molecule positions and orientations")
    ;


    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (po::error &err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }

    if (vm.count("help")) {
        //help_text();
        cout << desc << endl;
        return 0;
    }
    try {
        cout << "Reading Crg unit Types" << endl;
        JCalc jcalc(listcrg);
        cout << "Finished reading Crg unit Types" << endl;

        ifstream in1(pos1.c_str());
        ifstream in2(pos2.c_str());
        int written=0;
        while (in1 && in2) {
            string name1, name2;
            vec com1;
            vec com2;
            matrix or1;
            matrix or2;
            in1 >> name1 >> com1.x() >> com1.y() >> com1.z() >>
                    or1[0][0] >> or1[0][1] >> or1[0][2] >>
                    or1[1][0] >> or1[1][1] >> or1[1][2] >>
                    or1[2][0] >> or1[2][1] >> or1[2][2];
            in2 >> name2 >> com2.x() >> com2.y() >> com2.z() >>
                    or2[0][0] >> or2[0][1] >> or2[0][2] >>
                    or2[1][0] >> or2[1][1] >> or2[1][2] >>
                    or2[2][0] >> or2[2][1] >> or2[2][2];

            if (!in1 || !in2) break;
            //	cout << "mol1: " << name1 << com1 << or1<<endl;
            //	cout << "mol2: " << name2 << com2 << or2<<endl;


            CrgUnit * A = jcalc.CreateCrgUnit(0, name1);
            CrgUnit * B = jcalc.CreateCrgUnit(1, name2);
            A->SetPos(0, com1);
            B->SetPos(0, com2);
            A->SetNorm(0, or1[2]);
            A->SetPlane(0, or1[1]);
            B->SetNorm(0, or2[2]);
            B->SetPlane(0, or2[1]);


           //write pdb file
           mol_and_orb *molecule = ( A -> rotate_translate_beads() );
           (*molecule).write_pdb(pdbfile, "m1", written);
           written += (*molecule).getN();
           delete molecule;
           molecule = ( B -> rotate_translate_beads() );
           (*molecule).write_pdb(pdbfile, "m2", written);
           written += (*molecule).getN();
           delete molecule;
           ofstream fl;
           fl.open(pdbfile.c_str(), ios::app);
           fl.setf(ios::fixed);
           fl << "END" <<endl;

           
            //	cout << "Compute J" <<endl;
            vector <double> Js = jcalc.CalcJ(*A, *B);
            //	cout << "Finished computing J" <<endl;
            vector <double>::iterator itJ = Js.begin();
            for (; itJ != Js.end(); ++itJ) cout << '\t' << *itJ << endl;

        }
    } catch (std::exception &err) {
        cout << "error : " << err.what() << endl;
        return -1;
    }

}

*/