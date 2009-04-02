// 
// File:   template.cc
// Author: ruehle
//
// Created on June 8, 2008, 10:41 PM
//

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>
#include "imc.h"

using namespace std;

int main(int argc, char** argv)
{    
    // we have one observer
    Imc imc;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&imc);


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
   
    desc.add_options()
      ("options", po::value<string>(), "options file for coarse graining");
     
    // now read in the command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }
    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_imc, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_imc, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    if(!vm.count("options")) {
        cout << "need to specify options file\n";
        cout << desc << endl;
        return -1;
    }

    imc.LoadOptions(vm["options"].as<string>());
            
    // try to run the cg process, go through the frames, etc...
    //try {
        cg_engine.Run(desc, vm);
    //}
    // did an error occour?
    //catch(string error) {
     //   cerr << "An error occoured!" << endl << error << endl;
    //}
    return 0;
}


