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
#include <boost/format.hpp>
#include <topologyreader.h>
#include <libversion.h>
#include "cgengine.h"

void help_text(void)
{
    cout << "csg_gmxtopol, lib version " << LIB_VERSION_STR << "\n\n";
    cout << "Create skeleton for gromacs topology based on atomistic topology\n"
            "and a mapping file. Files still needs to be modified by the user.\n\n";     
}

using namespace std;

using boost::format;

void WriteAtoms(ostream &out, Molecule &cg)
{
    out << "[atoms]\n";
    out << "; nr type resnr residue atom cgnr charge mass\n";
    for(int i=0; i<cg.BeadCount(); ++i) {
        Bead *b=cg.getBead(i);
       
        out << format("%d %s 1 RES %s %d %f %f\n")
            % (i+1) % b->getType()->getName() % b->getName() % (i+1) % b->getQ() % b->getM();
    }
    out << endl;
}

void WriteInteractions(ostream &out, Molecule &cg)
{
    int nb=-1;
    
    Interaction *ic;
    vector<Interaction *>::iterator iter;
  
    InteractionContainer &ics=cg.getParent()->BondedInteractions();

    for(iter=ics.begin(); iter!=ics.end(); ++iter) {
        ic = *iter;
        if(ic->getMolecule() != cg.getId()) continue;
        if(nb != ic->BeadCount()) {
            nb=ic->BeadCount();
            switch(nb) {
                case 2:
                    out << "\n[ bonds ]\n";
                    break;
                case 3:
                    out << "\n[ angles ]\n";
                    break;
                case 4:
                    out << "\n[ dihedrals ]\n";
                    break;
                default:
                    throw runtime_error(string("cannot handle number of beads in interaction:") +
                            ic->getName());
            }
        }
        for(int i=0; i<nb; ++i)
            out << ic->getBeadId(i)+1 << " ";
        out << "  1  ; " << ic->getName() << endl;
    }
}

void WriteMolecule(ostream &out, Molecule &cg)
{
    out << "[ moleculetype ]\n";
    out << cg.getName() << " 3\n\n";

    WriteAtoms(out, cg);
    WriteInteractions(out, cg);
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
    ("top", boost::program_options::value<string>(), "atomistic topology file")
    ("cg", boost::program_options::value<string>(), "coarse grained mapping")
    ("out", boost::program_options::value<string>(), "output topology (will create .top and in future also .itp)");
    
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

    if(!vm.count("cg")) {
        cout << "no mapping specified\n";
        return 0;
    }

    if(!vm.count("out")) {
        cout << "no output topology specified\n";
        return 0;
    }

    if(!vm.count("top")) {
        cout << "no topology specified\n";
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        Topology top;
        Topology top_cg;
        CGEngine cg_engine;
        
        TopologyReader *reader;
        TopologyMap *map;
        
        reader = TopReaderFactory().Create(vm["top"].as<string>());
        if(reader == NULL) 
            throw string("input format not supported: ") + vm["top"].as<string>();
        
        reader->ReadTopology(vm["top"].as<string>(), top);
        cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;

        cg_engine.LoadMoleculeType(vm["cg"].as<string>());
        map = cg_engine.CreateCGTopology(top, top_cg);

        cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;

        if(top.MoleculeCount() > 1)
                cout << "WARNING: cannot create topology for topology with"
                "multiple molecules, using only first molecule\n";
        map->Apply();
	ofstream fl;
        fl.open((vm["out"].as<string>() + ".top").c_str());
                    
	WriteMolecule(fl, *top_cg.MoleculeByIndex(0));
        fl.close();
    }
    // did an error occour?
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    return 0;
}

