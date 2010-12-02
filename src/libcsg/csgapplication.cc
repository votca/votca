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

#include "csgapplication.h"
#include "trajectorywriter.h"
#include "trajectoryreader.h"
#include "topologyreader.h"
#include "topologymap.h"
#include "cgengine.h"
#include "version.h"

namespace votca { namespace csg {

CsgApplication::CsgApplication(void)
{
}

CsgApplication::~CsgApplication(void)
{
}

void CsgApplication::Initialize(void)
{
    // register all io plugins
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    AddProgramOptions()
        ("top", boost::program_options::value<string>(), "  atomistic topology file");
    if(DoMapping()) {
        if(DoMappingDefault()) {
            AddProgramOptions("Mapping options")
               ("cg", boost::program_options::value<string>(), "  coarse graining mapping definitions (xml-file)")
               ("no-map", "  disable mapping and act on original trajectory");
        } else {
            AddProgramOptions("Mapping options")
                ("cg", boost::program_options::value<string>(), "  [OPTIONAL] coarse graining mapping definitions\n"
                                                                "  (xml-file). If no file is given, program acts on original trajectory");
        }
    }
    

    if(DoTrajectory())
        AddProgramOptions("Trajectory options")
            ("trj", boost::program_options::value<string>(), "  atomistic trajectory file")
            ("begin", boost::program_options::value<double>()->default_value(0.0), "  skip frames before this time")
            ("first-frame", boost::program_options::value<int>()->default_value(0), "  start with this frame")
            ("nframes", boost::program_options::value<int>(), "  process so many frames")
            ;
}

bool CsgApplication::EvaluateOptions(void)
{
    _do_mapping = false;
    CheckRequired("top", "no topology file specified");
    
    // check for mapping options
    if (DoMapping()) {
        // default mapping is on
        if(DoMappingDefault()) {
            // if the user does not explicitly ask to turn it off, cg is needed
            if(OptionsMap().count("no-map")==0) {
                CheckRequired("cg", "no coarse graining definition specified");
                _do_mapping = true;
            }
        }
        // default mapping is off, if user gives cg, then do mapping
        else if(OptionsMap().count("cg")) {
            _do_mapping = true;
        }
    }
    return true;
}

void CsgApplication::ShowHelpText(std::ostream &out)
{
    string name =  ProgramName();
    if(VersionString() != "")
         name = name + ", version " + VersionString();
    
    HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}

void CsgApplication::Run(void)
{
    // first read in the topology
    TopologyReader *reader;
    Topology top;

    Topology top_cg;
    TopologyMap *map=0;
    CGEngine cg;

    // create reader for atomistic topology
    reader = TopReaderFactory().Create(_op_vm["top"].as<string>());
    if(reader == NULL)
        throw runtime_error(string("input format not supported: ") + _op_vm["top"].as<string>());

    // read in the topology
    reader->ReadTopology(_op_vm["top"].as<string>(), top);
    cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
    top.CheckMoleculeNaming();

    if(_do_mapping) {
        // read in the coarse graining definitions (xml files)
        cg.LoadMoleculeType(_op_vm["cg"].as<string>());
        // create the mapping + cg topology
        map = cg.CreateCGTopology(top, top_cg);

        cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
        map->Apply();
        if(!EvaluateTopology(&top_cg, &top))
            return;
    }
    else
        if(!EvaluateTopology(&top))
            return;

    // do we need to read a trajectory?
    if(DoTrajectory() && _op_vm.count("trj")) {
        TrajectoryReader *traj_reader;

        double begin;
        int first_frame;
        bool has_begin=false;

        if(_op_vm.count("begin")) {
            has_begin = true;
            begin = _op_vm["begin"].as<double>();
        }

        int nframes = -1;
        if(_op_vm.count("nframes")) {
            nframes = _op_vm["nframes"].as<int>();
        }

        first_frame = _op_vm["first-frame"].as<int>();        

        // create reader for trajectory
        traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<string>());
        if(traj_reader == NULL)
            throw runtime_error(string("input format not supported: ") + _op_vm["trj"].as<string>());
        // open the trajectory
        traj_reader->Open(_op_vm["trj"].as<string>());
        // read in first frame
        traj_reader->FirstFrame(top);

        // notify all observer that coarse graining has begun

        if(_do_mapping) {
            map->Apply();
            BeginEvaluate(&top_cg, &top);
        } else
            BeginEvaluate(&top);

        for(bool bok=true; bok==true; bok = traj_reader->NextFrame(top)) {
            if(((top.getTime() < begin) && has_begin)|| first_frame > 1) {
                first_frame--;
                continue;
            }
            if(nframes == 0 ) break;
            if(_do_mapping) {
                map->Apply();
                EvalConfiguration(&top_cg, &top);

            } else
                EvalConfiguration(&top);
            nframes--;
        }
        EndEvaluate();
        traj_reader->Close();

        delete traj_reader;
    }
    if(map)
        delete map;
    delete reader;
    
}

void CsgApplication::BeginEvaluate(Topology *top, Topology *top_ref) {
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->BeginCG(top, top_ref);
}

void CsgApplication::EndEvaluate()
{
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->EndCG();
}

void CsgApplication::EvalConfiguration(Topology *top, Topology *top_ref)
{
    list<CGObserver *>::iterator iter;
    for(iter=_observers.begin(); iter!=_observers.end(); ++iter)
        (*iter)->EvalConfiguration(top, top_ref);
}


}}
