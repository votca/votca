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

#include <fstream>
#include "cgengine.h"
#include <votca/tools/tokenizer.h>
#include "version.h"

namespace po = boost::program_options;

CGEngine::CGEngine()
    : _op_desc("Allowed options"), _op_desc_specific("Specific options")
{
    InitializeStandardOptions();
}

CGEngine::~CGEngine()
{
    map<string, CGMoleculeDef *>::iterator i;
    for(i=_molecule_defs.begin(); i!=_molecule_defs.end();++i)
        delete (*i).second;
    _molecule_defs.clear();
}

void CGEngine::Initialize()
{
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();
}

void CGEngine::ParseCommandLine(int argc, char **argv)
{
    _op_desc.add(_op_desc_specific);
    namespace po=boost::program_options;
    try {
        po::store(po::parse_command_line(argc, argv, _op_desc), _op_vm);
        po::notify(_op_vm);
    }
    catch(boost::program_options::error err) {
        throw runtime_error(string("error parsing command line: ") + err.what());
    }
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
            cout << "--------------------------------------\n"
                 << "WARNING: unknown molecule \"" << mol->getName() << "\" with id "
                 << mol->getId() << " in topology" << endl
                 << "molecule will not be mapped to CG representation\n"
                 << "Check weather a mapping file for all molecule exists, was specified in --cg "
                 << "separated by ; and the ident tag in xml-file matches the molecule name\n"
                 << "--------------------------------------\n";
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
    Tokenizer::iterator iter;
   
    for(iter=tok.begin(); iter!=tok.end(); ++iter) {
        CGMoleculeDef *def = new CGMoleculeDef();
        string  file = *iter;
        boost::trim(file);
        def->Load(file);
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

void CGEngine::InitializeStandardOptions()
{
    po::options_description trj("Trajectory options");
    _op_desc.add_options()
    ("help", "  produce this help message")
    ("top", boost::program_options::value<string>(), "  atomistic topology file")
    ("cg", boost::program_options::value<string>(), "  coarse graining definitions (xml-file)")
    ;
    
    trj.add_options()
    ("trj", boost::program_options::value<string>(), "  atomistic trajectory file")
    ("begin", boost::program_options::value<double>(), "  skip frames before this time")
    ("first-frame", boost::program_options::value<int>()->default_value(0), "  start with this frame")
    ("nframes", boost::program_options::value<int>(), "  process so many frames")
    ;

    _op_desc.add(trj);
}

void CGEngine::Run()
{
    TopologyReader *reader;
    TrajectoryReader *traj_reader;
    Topology top;
    Topology top_cg;
    TopologyMap *map;
    double begin;
    int first_frame;
    
    bool has_begin=false;

//    if (vm.count("version")) {
//        votca::csg::HelpTextHeader("csg");
//        exit(0);
//    }


    if (!_op_vm.count("top")) {
        cout << _op_desc << endl;
        throw runtime_error("no topology file specified");
    }
    if (!_op_vm.count("cg")) {
        cout << _op_desc << endl;
        throw runtime_error("no coarse graining definition specified");
    }
    
    if(_op_vm.count("begin")) {
        has_begin = true;
        begin = _op_vm["begin"].as<double>();
    }    
    
    int nframes = -1;
    if(_op_vm.count("nframes")) {
        nframes = _op_vm["nframes"].as<int>();
    }
    
    first_frame = _op_vm["first-frame"].as<int>();
    
    // create reader for atomistic topology
    reader = TopReaderFactory().Create(_op_vm["top"].as<string>());
    if(reader == NULL) 
        throw runtime_error(string("input format not supported: ") + _op_vm["top"].as<string>());
        
    // read in the atomistic topology
    reader->ReadTopology(_op_vm["top"].as<string>(), top);
    cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
    top.CheckMoleculeNaming();
    
    //top.CreateMoleculesByResidue();    
    //top.CreateOneBigMolecule("PS1");    
        
    // read in the coarse graining definitions (xml files)
    LoadMoleculeType(_op_vm["cg"].as<string>());
    // create the mapping + cg topology
    map = CreateCGTopology(top, top_cg);
    //cg_def.CreateMolecule(top_cg);
              
    cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
   

    // if trj is given, process trajectory file
    if (_op_vm.count("trj")) {
        // create reader for trajectory
        traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<string>());
        if(traj_reader == NULL)
            throw runtime_error(string("input format not supported: ") + _op_vm["trj"].as<string>());
        // open the trajectory
        traj_reader->Open(_op_vm["trj"].as<string>());
        // read in first frame
        traj_reader->FirstFrame(top);    
        
        // notify all observer that coarse graining has begun
             
        map->Apply();
        BeginCG(top_cg);
        
        for(bool bok=true; bok==true; bok = traj_reader->NextFrame(top)) {
            if(((top.getTime() < begin) && has_begin)|| first_frame > 0) {
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
