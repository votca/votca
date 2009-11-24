// 
// File:   cgengine.cc
// Author: ruehle
//
// Created on April 26, 2007, 4:20 PM
//

#include <fstream>
#include "cgengine.h"
#include <votca/tools/tokenizer.h>

CGEngine::~CGEngine()
{
    map<string, CGMoleculeDef *>::iterator i;
    for(i=_molecule_defs.begin(); i!=_molecule_defs.end();++i)
        delete (*i).second;
    _molecule_defs.clear();
}

/**
    \todo melts with different molecules
*/
TopologyMap *CGEngine::CreateCGTopology(Topology &in, Topology &out)
{
    MoleculeContainer &mols = in.Molecules();
    MoleculeContainer::iterator iter;
    TopologyMap *m = new TopologyMap(&in, &out);
    for(iter=mols.begin(); iter!=mols.end(); ++iter) {
        Molecule *mol = *iter;
        CGMoleculeDef *def = getMoleculeDef(mol->getName());
        if(!def) {
            cout << "unknown molecule " << mol->getName() << " in topology" << endl;
            continue;
        }
        Molecule *mcg = def->CreateMolecule(out);
        Map *map = def->CreateMap(*mol, *mcg);
        m->AddMoleculeMap(map);
    }
    out.RebuildExclusions();    
    return m;
}

void CGEngine::LoadMoleculeType(string filename)
{
    Tokenizer tok(filename, ";");
    Tokenizer::iterator file; 
   
    for(file=tok.begin(); file!=tok.end(); ++file) {
        CGMoleculeDef *def = new CGMoleculeDef();        
        def->Load(*file);
        _molecule_defs[def->getIdent()] = def;
    }
}
 
void CGEngine::BeginCG(Topology &top)
{
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->BeginCG(&top);
}

void CGEngine::EndCG()
{
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->EndCG();
}
    
void CGEngine::EvalConfiguration(Topology &top_cg)
{
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->EvalConfiguration(&top_cg);
}

void CGEngine::AddProgramOptions(boost::program_options::options_description &desc)
{
    desc.add_options()
    ("help", "produce this help message")
    ("version", "show version info")
    ("top", boost::program_options::value<string>(), "atomistic topology file")
    ("trj", boost::program_options::value<string>(), "atomistic trajectory file")
    ("cg", boost::program_options::value<string>(), "coarse graining definitions (xml-file)")
    ("begin", boost::program_options::value<double>(), "skip frames before this time")
    ("first-frame", boost::program_options::value<int>()->default_value(0), "start with this frame")
    ("nframes", boost::program_options::value<int>(), "process so many frames")
    ;
}

void CGEngine::Run(boost::program_options::options_description &desc, boost::program_options::variables_map &vm)
{
    TopologyReader *reader;
    TrajectoryReader *traj_reader;
    Topology top;
    Topology top_cg;
    TopologyMap *map;
    double begin;
    int first_frame;
    
    bool has_begin=false;
    if (!vm.count("top")) {
        cout << desc << endl;
        throw runtime_error("no topology file specified");
    }
    if (!vm.count("cg")) {
        cout << desc << endl;
        throw runtime_error("no coarse graining definition specified");
    }
    
    if(vm.count("begin")) {
        has_begin = true;
        begin = vm["begin"].as<double>();
    }
    
    
    int nframes = -1;
    if(vm.count("nframes")) {
        nframes = vm["nframes"].as<int>();
    }
    
    first_frame = vm["first-frame"].as<int>();
    
    // create reader for atomistic topology
    reader = TopReaderFactory().Create(vm["top"].as<string>());
    if(reader == NULL) 
        throw runtime_error(string("input format not supported: ") + vm["top"].as<string>());
        
    // read in the atomistic topology
    reader->ReadTopology(vm["top"].as<string>(), top);
    cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
    
    //top.CreateMoleculesByResidue();    
    //top.CreateOneBigMolecule("PS1");    
        
    // read in the coarse graining definitions (xml files)
    LoadMoleculeType(vm["cg"].as<string>());
    // create the mapping + cg topology
    map = CreateCGTopology(top, top_cg);
    //cg_def.CreateMolecule(top_cg);
              
    cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
   

    // if trj is given, process trajectory file
    if (vm.count("trj")) {
        // create reader for trajectory
        traj_reader = TrjReaderFactory().Create(vm["trj"].as<string>());
        if(traj_reader == NULL)
            throw runtime_error(string("input format not supported: ") + vm["trj"].as<string>());
        // open the trajectory
        traj_reader->Open(vm["trj"].as<string>());
        // read in first frame
        traj_reader->FirstFrame(top);    
        
        // notify all observer that coarse graining has begun
             
        map->Apply();
        BeginCG(top_cg);
        
        for(bool bok=true; bok==true; bok = traj_reader->NextFrame(top)) {
            if(top.getTime() < begin || first_frame > 0) {
                first_frame--;
                continue;                
            }
            if(nframes == 0 ) break;
            map->Apply();
            EvalConfiguration(top_cg);            
            nframes--;
        }
        EndCG();
        traj_reader->Close();
    }
    delete traj_reader;
    delete reader;
    delete map;                   
}
