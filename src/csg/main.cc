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
#include <tools/vec.h>
#include "cgmoleculedef.h"
#include "cgengine.h"
#include "molecule.h"
#include "topologyreader.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include <tools/tokenizer.h>
#include <tools/matrix.h>
#include "analysistool.h"
#include "version.h"
#include <tools/rangeparser.h>
#include "bondedstatistics.h"
#include "libversion.h"

using namespace std;

int main(int argc, char** argv)
{    
    BondedStatistics bs;
    TopologyReader *reader;
    TrajectoryReader *traj_reader;
    Topology top;
    Topology top_cg;
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
    ("out", po::value<string>(), "write pdb cg trajectory")
    ("excl", po::value<string>(), "write exclusion list to file")
    ;
    
    cg_engine.AddObserver(&bs);
    
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();
    
    po::variables_map vm;    
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }

    if (vm.count("help")) {
        cout << "csg version " << VERSION_STR << "\n";                
        cout << "libcsg version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    if (vm.count("version")) {
        cout << "csg version " << VERSION_STR  << "\n";                        
        cout << "libcsg version " << LIB_VERSION_STR << "\n\n";                
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
            cerr << "input format not supported:" << vm["top"].as<string>() << endl;
            return 1;
        }
        reader->ReadTopology(vm["top"].as<string>(), top);
        cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
        //top.CreateMoleculesByResidue();    
        //top.CreateOneBigMolecule("PS1");    
        
        cg_engine.LoadMoleculeType(vm["cg"].as<string>());
        //cg_engine.LoadMoleculeType("Cl.xml");
        map = cg_engine.CreateCGTopology(top, top_cg);
        //cg_def.CreateMolecule(top_cg);
              
        cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
        
        if (vm.count("excl")) {
            ExclusionList *ex;
            CGMoleculeDef *cg_def;
            cg_def = cg_engine.getMoleculeDef(top.MoleculeByIndex(0)->getName());
            if(!cg_def)
                throw string("error, no molecule definition found to create exclusion list");
                            
            ex = cg_def->CreateExclusionList(*top.MoleculeByIndex(0));
            ofstream fl;
            fl.open(vm["excl"].as<string>().c_str());
            fl << *ex;        
            fl.close();
            delete ex;
            return 0;
        }

        if (vm.count("trj")) {
            traj_reader = TrjReaderFactory().Create(vm["trj"].as<string>());
            if(traj_reader == NULL) {
                cerr << "input format not supported:" << vm["trj"].as<string>() << endl;
                return 1;
            }
            traj_reader->Open(vm["trj"].as<string>());
            traj_reader->FirstFrame(top);    
        
            cg_engine.BeginCG(top_cg);
            bool bok=true;
            while(bok) {
                map->Apply();
                cg_engine.EvalConfiguration(top_cg);
                if(bWrite) writer->Write(&top_cg);
                bok = traj_reader->NextFrame(top);
            }
            cg_engine.EndCG();
            traj_reader->Close();
        }
        delete traj_reader;
        delete reader;
        delete map;
    
    }
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    if (vm.count("out")) {
        writer->Close();
        delete writer;
    }
    
    cout << "Interactive mode, expecting commands:" << endl;
    cout << "q: quit" << endl;
    cout << "list: list all available bonds" << endl;
    cout << "vals <file> <selection>: write values to file" << endl;
    cout << "hist <file> <selection>: create histogram" << endl;
    cout << "tab <file> <selection>: create tabulated potential" << endl;
    cout << "autocor <file> <selection>: calculate autocorrelation, only one row allowed in selection!" << endl;
    cout << "cor <file> <selection>: calculate correlations, first row is correlated with all other rows" << endl;
    
    while(1) {
        string line;
        cout << "> ";
        getline(cin, line);
        size_t start;
        size_t end;
        
        start = line.find_first_not_of(' ');
        if(start == string::npos) continue;
        end = line.find(' ', start);
        
        string cmd = line.substr(start, end-start);
        start = end+1;        
        try {
    
            if(cmd == "q") break;
            string arg_str(line.substr(start));
            vector<string> args;
            Tokenizer tok(arg_str, " \t");
            tok.ToVector(args);
            AnalysisTool *tool = AnalysisFactory().get(cmd);
            if(!tool) {
                cout << "error, command not found" << endl;
                continue;
            }
                tool->Command(bs, cmd, args);
        }
        catch(string error) {
            cerr << "An error occoured:" << endl << error << endl;
        }
    
    }

    return 0;
}

