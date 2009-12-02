// 
// File:   template.cc
// Author: ruehle
//
// Created on June 8, 2008, 10:41 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <topologyreader.h>
#include "version.h"

using namespace std;

void help_text(void)
{
    votca::csg::HelpTextHeader("csg_dump");
    cout << "Print atoms that are read from topology file to help"
        " debugging atom naming.\n\n";  
}

int main(int argc, char** argv)
{    
    // initialize the readers/writers,
    // this will be combined in an initialize function later
//    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    // lets read in some program options
    namespace po = boost::program_options;
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    
    
    // let cg_engine add some program options
    desc.add_options()
    ("help", "produce this help message")
    //("version", "show version info")
    ("top", boost::program_options::value<string>(), "atomistic topology file");
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        help_text();
        cout << desc << endl;
        return 0;
    }

    if(!vm.count("top")) {
        cout << "no topology specified\n";
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        Topology top;
        TopologyReader *reader;
        reader = TopReaderFactory().Create(vm["top"].as<string>());
        if(reader == NULL) 
            throw std::runtime_error("input format not supported: " + vm["top"].as<string>());
        
        reader->ReadTopology(vm["top"].as<string>(), top);
        cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
        
        MoleculeContainer::iterator mol;
        for(mol=top.Molecules().begin(); mol!=top.Molecules().end();++mol) {
            cout << "molecule: " << (*mol)->getId() + 1 << " " << (*mol)->getName() 
              << " beads: " << (*mol)->BeadCount() << endl;
            for(int i=0; i<(*mol)->BeadCount(); ++i) {
                cout << (*mol)->getBeadId(i) << " " << 
                    (*mol)->getBeadName(i) << " " << (*mol)->getBead(i)->getType()->getName() << endl;
            }
        }
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

