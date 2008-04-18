// 
// File:   csg_nemat.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:35 PM
//

// 
// File:   main.cc
// Author: ruehle
//
// Created on April 5, 2007, 12:29 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "vec.h"
#include "cgmoleculedef.h"
#include "cgengine.h"
#include "molecule.h"
#include "topologyreader.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include "tools.h"
#include "matrix.h"
#include "libversion.h"
#include "nematicorder.h"

using namespace std;

int main(int argc, char** argv)
{    
    TopologyReader *reader;
    TrajectoryReader *traj_reader;
    Topology top;
    Topology top_cg;
    Configuration conf(&top);
    Configuration conf_cg(&top_cg);
    CGEngine cg_engine;
    TrajectoryWriter *writer;
    TopologyMap *map;
    bool bWrite = false;
    namespace po = boost::program_options;
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help", "produce this help message")
    ("version", "show version info")
    ("top", po::value<string>(), "atomistic topology file")
    ("trj", po::value<string>(), "atomistic trajectory file")
    ("cg", po::value<string>(), "coarse graining definitions (xml-file)")
    ;

    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    if (vm.count("version")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    if (!vm.count("top")) {
        cout << desc << endl;
        cout << "no topology file specified" << endl;
        return 1;
    }
    if (!vm.count("cg")) {
        cout << desc << endl;
        cout << "no coarse graining definition specified" << endl;
        return 1;
    }
    if (vm.count("out")) {
        if (!vm.count("trj")) {
            cout << desc << endl;
            cout << "no trajectory file specified" << endl;
            return 1;
        }

        writer = TrjWriterFactory().Create(vm["out"].as<string>());
        if(writer == NULL) {
            cerr << "output format not supported:" << vm["out"].as<string>() << endl;
            return 1;
        }
        bWrite = true;
        writer->Open(vm["out"].as<string>());
    }
        
    try {
        reader = TopReaderFactory().Create(vm["top"].as<string>());
        if(reader == NULL) {
            cerr << "input format not supported:" << vm["out"].as<string>() << endl;
            return 1;
        }
        reader->ReadTopology(vm["top"].as<string>(), top);
        conf.Initialize();
        cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
        //top.CreateMoleculesByResidue();    
        //top.CreateOneBigMolecule("PS1");    
        
        cg_engine.LoadMoleculeType(vm["cg"].as<string>());
        map = cg_engine.CreateCGTopology(top, top_cg);
        //cg_def.CreateMolecule(top_cg);
        conf_cg.Initialize();
              
        cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
        
        NematicOrder nemat;
        if (vm.count("trj")) {
            traj_reader = TrjReaderFactory().Create(vm["trj"].as<string>());
            if(traj_reader == NULL) {
                cerr << "input format not supported:" << vm["trj"].as<string>() << endl;
                return 1;
            }
            traj_reader->Open(vm["trj"].as<string>());
            traj_reader->FirstFrame(conf);    
        
            cg_engine.BeginCG(top_cg);
            bool bok=true;
            while(bok) {
                map->Apply(conf, conf_cg);
                cg_engine.EvalConfiguration(conf_cg);
                
                nemat.Process(conf_cg);
                cout << conf.getTime() << " " <<  nemat.NematicU().eigenvalues[2] << " " 
                   << nemat.NematicV().eigenvalues[2]<< endl;
                if(bWrite) writer->Write(&conf_cg);
                bok = traj_reader->NextFrame(conf);
            }
            cg_engine.EndCG();
            traj_reader->Close();
        }
        delete traj_reader;
        delete reader;
        delete map;            
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
    }
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    
    return 0;
}

