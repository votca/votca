/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
// 
// File:   main.cc
// Author: ruehle
//
// Created on April 5, 2007, 12:29 PM
//
// TODO: This code need lots of cleaning up! please do not look at anything in here!
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <votca/tools/vec.h>
#include "cgmoleculedef.h"
#include "cgengine.h"
#include "molecule.h"
#include "topologyreader.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include <votca/tools/tokenizer.h>
#include <votca/tools/matrix.h>
#include "analysistool.h"
#include "version.h"
#include <votca/tools/rangeparser.h>
#include "bondedstatistics.h"
#include "version.h"
#include <map>
#include <string>
#include "tabulatedpotential.h"
#include "stdanalysis.h"

using namespace std;

void help_text()
{
    votca::csg::HelpTextHeader("csg_boltzmann");
    cout << "Perform tasks that are needed for simple boltzmann\n"
            "inversion in an interactive environment.\n\n";
}

ExclusionList *CreateExclusionList(Molecule &atomistic, Molecule &cg)
{
    list<int> exclude;
   
    ExclusionList *ex = new ExclusionList();
    ex->ExcludeAll(atomistic.BeadCount());

    // reintroduce bead internal nonbonded interaction
    for(int i=0; i<cg.BeadCount(); ++i) {
        exclude.clear();
        
        vector<int> &v = cg.getBead(i)->ParentBeads();
        exclude.insert(exclude.begin(), v.begin(), v.end());
        ex->Remove(exclude);
    }

    Topology *top_cg = cg.getParent();
    InteractionContainer::iterator iter;
    // reintroduce nonbonded interactions for bonded beads
    for(iter = top_cg->BondedInteractions().begin();
            iter!=top_cg->BondedInteractions().end(); ++iter) {
        Interaction *ic = *iter;
        exclude.clear();
        for(size_t i=0; i<ic->BeadCount(); i++) {
            vector<int> &v = top_cg->getBead(ic->getBeadId(i))->ParentBeads();
            exclude.insert(exclude.end(), v.begin(), v.end());
        }
        ex->Remove(exclude);
    }
    return ex;
}


int main(int argc, char** argv)
{    
    BondedStatistics bs;
    TopologyReader *reader;
    TrajectoryReader *traj_reader;
    Topology top;
    Topology top_cg;
    CGEngine cg_engine;
    TopologyMap *map;

    namespace po = boost::program_options;
    std::map<std::string, AnalysisTool *> cmds;
    TabulatedPotential tab;
    StdAnalysis std;
    tab.Register(cmds);
    std.Register(cmds);

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help", "produce this help message")
    ("top", po::value<string>(), "atomistic topology file")
    ("trj", po::value<string>(), "atomistic trajectory file")
    ("cg", po::value<string>(), "coarse graining definitions (xml-file)")
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
        help_text();
        cout << desc << endl;
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

    try {
        reader = TopReaderFactory().Create(vm["top"].as<string>());
        if(reader == NULL) {
            cerr << "input format not supported:" << vm["top"].as<string>() << endl;
            return 1;
        }
        reader->ReadTopology(vm["top"].as<string>(), top);
        cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;

        top.CheckMoleculeNaming();
        //top.CreateMoleculesByResidue();    
        //top.CreateOneBigMolecule("PS1");    
        
        cg_engine.LoadMoleculeType(vm["cg"].as<string>());
        //cg_engine.LoadMoleculeType("Cl.xml");
        map = cg_engine.CreateCGTopology(top, top_cg);
        //cg_def.CreateMolecule(top_cg);
              
        cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
        
        if (vm.count("excl")) {
            ExclusionList *ex;
            if(top.MoleculeCount() > 1)
                cout << "WARNING: cannot create exclusion list for topology with"
                "multiple molecules, using only first molecule\n";
            
            map->Apply();
            cout << "Writing exclusion list for atomistic molecule "
                    << top.MoleculeByIndex(0)->getName()
                    << " in coarse grained representation "
                    << top.MoleculeByIndex(0)->getName() << endl;
            ex = CreateExclusionList(*top.MoleculeByIndex(0), *top_cg.MoleculeByIndex(0));
            ofstream fl;
            fl.open(vm["excl"].as<string>().c_str());
            fl << "# atomistic: " << top.MoleculeByIndex(0)->getName()
               << " cg: " << top.MoleculeByIndex(0)->getName()
               << " cgmap: " << vm["cg"].as<string>() << endl;
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
                bok = traj_reader->NextFrame(top);
            }
            cg_engine.EndCG();
            traj_reader->Close();
        }
        delete traj_reader;
        delete reader;
        delete map;
    
    }
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
        exit(1);
    }
    
    string help_text = 
        "Interactive mode, expecting commands:\n"
        "help: show this help\n"
        "q: quit\n"
        "list: list all available bonds\n"
    	"vals <file> <selection>: write values to file\n"
    	"hist <file> <selection>: create histogram\n"
    	"tab <file> <selection>: create tabulated potential\n"
    	"autocor <file> <selection>: calculate autocorrelation, only one row allowed in selection!\n" 
    	"cor <file> <selection>: calculate correlations, first row is correlated with all other rows";

    cout << help_text << endl;
    
    while(1) {
        string line;
        cout << "> ";
        getline(cin, line);
        size_t start;
        size_t end;

        boost::trim(line);
        vector<string> args;
        Tokenizer tok(line, " \t");
        tok.ToVector(args);

        if(args.size() == 0) continue;

        string cmd = args.front();
        args.erase(args.begin());

        try {
    
            if(cmd == "q") break;            

            std::map<string, AnalysisTool *>::iterator tool;
            if(cmd == "help") {
                if(args.size() == 0) {
                    cout << help_text << endl;
                    continue;
                }
                cmd = args.front();
                args.erase(args.begin());
                tool = cmds.find(cmd);
                if(tool == cmds.end()) {
                    cout << "error, no help item found" << endl;
                    continue;
                }
                tool->second->Help(cmd, args);
                cout << endl;
                continue;
            }

            tool = cmds.find(cmd);
            if(tool == cmds.end()) {
                cout << "error, command not found" << endl;
                continue;
            }
            
            tool->second->Command(bs, cmd, args);
        }
        catch(std::exception &error) {
            cerr << "An error occoured:" << endl << error.what() << endl;
        }
    
    }

    return 0;
}

