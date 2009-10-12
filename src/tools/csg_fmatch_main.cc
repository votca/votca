/* 
 * File:   csg_fmatch_main.cc
 * Author: lukyanov
 *
 * Created on June 10, 2009, 5:00 PM
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>
#include "csg_fmatch.h"
/*
 * 
 */
int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGForceMatching fmatch;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&fmatch);


    // initialize the readers/writers,
    // this will be combined in an initialize function later
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    
    // lets read in some program options
    namespace po = boost::program_options;
        
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    

    desc.add_options()
      ("options", po::value<string>(), "options file for coarse graining");    
    
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    if(!vm.count("options")) {
        cout << "need to specify options file\n";
        cout << desc << endl;
        return -1;
    }
    
    fmatch.LoadOptions(vm["options"].as<string>());
    
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }

    // did an error occour?
    catch(exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
        return -1;
    }
    return 0;
}
