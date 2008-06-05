// 
// File:   csg_nemat.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:35 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "cgengine.h"
#include "libversion.h"
#include "nematicorder.h"

using namespace std;

class CGNematicOrder
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {};
    void EndCG() {
        cout << "nematic order of u: "<< endl
            << nemat.NematicU().eigenvalues[0] << ": " << nemat.NematicU().eigenvecs[0] << "\n"
            << nemat.NematicU().eigenvalues[1] << ": " << nemat.NematicU().eigenvecs[1] << "\n"
            << nemat.NematicU().eigenvalues[2] << ": " << nemat.NematicU().eigenvecs[2] << "\n";
        cout << "nematic order of v: "<< endl
            << nemat.NematicV().eigenvalues[0] << ": " << nemat.NematicV().eigenvecs[0] << "\n"
            << nemat.NematicV().eigenvalues[1] << ": " << nemat.NematicV().eigenvecs[1] << "\n"
            << nemat.NematicV().eigenvalues[2] << ": " << nemat.NematicV().eigenvecs[2] << "\n";
/*        cout << "nematic order of w: "<< endl
            << nemat.NematicW().eigenvalues[0] << ": " << nemat.NematicW().eigenvecs[0] << "\n"
            << nemat.NematicW().eigenvalues[1] << ": " << nemat.NematicW().eigenvecs[1] << "\n"
            << nemat.NematicW().eigenvalues[2] << ": " << nemat.NematicW().eigenvecs[2] << "\n";
*/
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
        nemat.Process(*conf);
        cout << conf->getTime() << " " <<  nemat.NematicU().eigenvalues[2] << " " 
            << nemat.NematicV().eigenvalues[2]<< endl;
    }
    
protected:
    NematicOrder nemat;
};


int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGNematicOrder no;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&no);


    // initialize the readers/writers,
    // this will be combined in an initialize function later
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    
    // lets read in some program options
    namespace po = boost::program_options;
        
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    
    
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }
    // did an error occour?
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    return 0;
}

