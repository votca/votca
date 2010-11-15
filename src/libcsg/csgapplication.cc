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
#include "votca/tools/easylock.h"

namespace votca {
    namespace csg {

        CsgApplication::CsgApplication(void) {
        }

        CsgApplication::~CsgApplication(void) {
        }

        void CsgApplication::Initialize(void) {
            // register all io plugins
            TrajectoryWriter::RegisterPlugins();
            TrajectoryReader::RegisterPlugins();
            TopologyReader::RegisterPlugins();

            AddProgramOptions()
                    ("top", boost::program_options::value<string > (), "  atomistic topology file");
            if (DoMapping()) {
                if (DoMappingDefault()) {
                    AddProgramOptions("Mapping options")
                            ("cg", boost::program_options::value<string > (), "  coarse graining mapping definitions (xml-file)")
                            ("no-map", "  disable mapping and act on original trajectory");
                } else {
                    AddProgramOptions("Mapping options")
                            ("cg", boost::program_options::value<string > (), "  [OPTIONAL] coarse graining mapping definitions\n"
                            "  (xml-file). If no file is given, program acts on original trajectory");
                }
            }


            if (DoTrajectory())
                AddProgramOptions("Trajectory options")
                ("trj", boost::program_options::value<string > (), "  atomistic trajectory file")
                ("begin", boost::program_options::value<double>()->default_value(0.0), "  skip frames before this time")
                ("first-frame", boost::program_options::value<int>()->default_value(0), "  start with this frame")
                ("nframes", boost::program_options::value<int>(), "  process so many frames")
                ;

            if (DoThreaded())
                AddProgramOptions("Threading options")
                ("nt", boost::program_options::value<int>(), "  number of threads")
                ;

        }

        bool CsgApplication::EvaluateOptions(void) {
            _do_mapping = false;
            CheckRequired("top", "no topology file specified");

            // check for mapping options
            if (DoMapping()) {
                // default mapping is on
                if (DoMappingDefault()) {
                    // if the user does not explicitly ask to turn it off, cg is needed
                    if (OptionsMap().count("no-map") == 0) {
                        CheckRequired("cg", "no coarse graining definition specified");
                        _do_mapping = true;
                    }
                }// default mapping is off, if user gives cg, then do mapping
                else if (OptionsMap().count("cg")) {
                    _do_mapping = true;
                }
            }

            /* check threading options */
            if (DoThreaded()) {
                /* TODO
                 * does the number of threads make sense?
                 * which criteria should be used? smaller than system's cores?
                 */
            }

            return true;
        }

        void CsgApplication::ShowHelpText(std::ostream &out) {
            string name = ProgramName();
            if (VersionString() != "")
                name = name + ", version " + VersionString();

            HelpTextHeader(name);
            HelpText(out);
            out << "\n\n" << OptionsDesc() << endl;
        }

        void CsgApplication::Run(void) {
            if (DoThreaded()) {
                RunThreaded();
                return;
            }
            // first read in the topology
            TopologyReader *reader;
            Topology top;

            Topology top_cg;
            TopologyMap *map = 0;
            CGEngine cg;

            // create reader for atomistic topology
            reader = TopReaderFactory().Create(_op_vm["top"].as<string > ());
            if (reader == NULL)
                throw runtime_error(string("input format not supported: ") + _op_vm["top"].as<string > ());

            // read in the topology
            reader->ReadTopology(_op_vm["top"].as<string > (), top);
            cout << "I have " << top.BeadCount() << " beads in " << top.MoleculeCount() << " molecules" << endl;
            top.CheckMoleculeNaming();

            if (_do_mapping) {
                // read in the coarse graining definitions (xml files)
                cg.LoadMoleculeType(_op_vm["cg"].as<string > ());
                // create the mapping + cg topology
                map = cg.CreateCGTopology(top, top_cg);

                cout << "I have " << top_cg.BeadCount() << " beads in " << top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
                map->Apply();
                if (!EvaluateTopology(&top_cg, &top))
                    return;
            } else
                if (!EvaluateTopology(&top))
                return;

            // do we need to read a trajectory?
            if (DoTrajectory()) {
                TrajectoryReader *traj_reader;

                double begin;
                int first_frame;
                bool has_begin = false;

                if (_op_vm.count("begin")) {
                    has_begin = true;
                    begin = _op_vm["begin"].as<double>();
                }

                int nframes = -1;
                if (_op_vm.count("nframes")) {
                    nframes = _op_vm["nframes"].as<int>();
                }

                first_frame = _op_vm["first-frame"].as<int>();

                // create reader for trajectory
                traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<string > ());
                if (traj_reader == NULL)
                    throw runtime_error(string("input format not supported: ") + _op_vm["trj"].as<string > ());
                // open the trajectory
                traj_reader->Open(_op_vm["trj"].as<string > ());

                // read in first frame
                traj_reader->FirstFrame(top);

                // notify all observer that coarse graining has begun

                if (_do_mapping) {
                    map->Apply();
                    BeginEvaluate(&top_cg, &top);
                } else
                    BeginEvaluate(&top);

                for (bool bok = true; bok == true; bok = traj_reader->NextFrame(top)) {
                    if (((top.getTime() < begin) && has_begin) || first_frame > 1) {
                        first_frame--;
                        continue;
                    }
                    if (nframes == 0) break;
                    if (_do_mapping) {
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
            if (map)
                delete map;
            delete reader;

        }

        void CsgApplication::Worker::Initialize(int id, Topology * top, Topology * top_cg, TopologyMap * map,                \
                TrajectoryReader *traj_reader, bool do_mapping, int * nframes,                \
                Mutex * nframesMutex, Mutex * traj_readerMutex) {
            _myId = id;
            _top = top;
            _top_cg = top_cg;
            _map = map;
            _traj_reader = traj_reader;
            _do_mapping = do_mapping;
            _nframes = nframes;
            _nframesMutex = nframesMutex;
            _traj_readerMutex = traj_readerMutex;

        }

#define NUMBERTHREADS 8

        void CsgApplication::Worker::Run(void) {
            //move the whole nframes issue into getData?
            while (GetData(_traj_reader, _top)) {
                _nframesMutex->Lock();
                if (*_nframes > 0 || *_nframes <= -99) { // in getData
                    Decr(_nframes);
                    _nframesMutex->Unlock();
                    if (_do_mapping) {
                        _map->Apply();
                        EvalConfiguration(_top_cg, _top);

                    } else
                        EvalConfiguration(_top);
                } else {
                    _nframesMutex->Unlock();
                    return;
                }
            }
            return;
        }

        bool CsgApplication::Worker::GetData(TrajectoryReader *traj_reader, Topology * top) {
            //easyLock
            //EasyLock locker(_traj_readerMutex);
            //return traj_reader->NextFrame(*top);
            bool tmpRes;
            _traj_readerMutex->Lock();
            //std::cout << "lock: _traj_reader...";
            tmpRes = _traj_reader->NextFrame(*top);
            _traj_readerMutex->Unlock();
            //std::cout << "...unlock!";
            //std::cout << "         getData: " << tmpRes << std::endl;
            return tmpRes;
        }

        void CsgApplication::RunThreaded(void) {
            Worker * myWorker;
            // first read in the topology
            TopologyReader *reader;

            // create reader for atomistic topology
            reader = TopReaderFactory().Create(_op_vm["top"].as<string > ());
            if (reader == NULL)
                throw runtime_error(string("input format not supported: ") + _op_vm["top"].as<string > ());

            if (DoTrajectory()) {
                TrajectoryReader *traj_reader; // nur einmal fuer alle

                double begin;
                int first_frame;
                bool has_begin = false;

                if (_op_vm.count("begin")) {
                    has_begin = true;
                    begin = _op_vm["begin"].as<double>();
                }

                int nframes = -99;
                if (_op_vm.count("nframes")) {
                    nframes = _op_vm["nframes"].as<int>();
                }

                first_frame = _op_vm["first-frame"].as<int>();

                // create reader for trajectory
                traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<string > ());
                if (traj_reader == NULL)
                    throw runtime_error(string("input format not supported: ") + _op_vm["trj"].as<string > ());
                // open the trajectory
                traj_reader->Open(_op_vm["trj"].as<string > ());


                Mutex * nframesMutex = new Mutex; // als member in application ohne *
                Mutex * traj_readerMutex = new Mutex; // als member in application ohne *
                _myWorkers = new Worker *[NUMBERTHREADS]; // als std::vector workers, push_back
                //for every thread
                for (int thread = 0; thread < NUMBERTHREADS; thread++) {
                    myWorker = ForkWorker();
                    _myWorkers[thread] = myWorker;
                    Topology *top = new Topology;
                    Topology *top_cg = new Topology;
                    TopologyMap *map = 0; // in constructor
                    CGEngine cg; // nur einmal nicht fuer alle threads


                    // read in the topology
                    reader->ReadTopology(_op_vm["top"].as<string > (), *top);
                    //dirty check
                    if (thread == 0)
                        cout << "I have " << top->BeadCount() << " beads in " << top->MoleculeCount() << " molecules" << endl;
                    top->CheckMoleculeNaming();

                    if (_do_mapping) {
                        // read in the coarse graining definitions (xml files)
                        cg.LoadMoleculeType(_op_vm["cg"].as<string > ());
                        // create the mapping + cg topology
                        map = cg.CreateCGTopology(*top, *top_cg);
                        //dirty check
                        if (thread == 0)
                            cout << "I have " << top_cg->BeadCount() << " beads in " << top_cg->MoleculeCount() << " molecules for the coarsegraining" << endl;
                        map->Apply();

                        if (!EvaluateTopology(top_cg, top))
                            return;
                    } else

                        if (!EvaluateTopology(top))
                        return;

                    if (thread == 0) {
                        // read in first frame
                        traj_reader->FirstFrame(*top);

                        // notify all observer that coarse graining has begun

                        if (_do_mapping) {
                            map->Apply();
                            BeginEvaluate(top_cg, top); // einmal fuer application
                        } else
                            BeginEvaluate(top); // einmal fuer application
                    }

                    // ab hier wirds interessant fuer threads

                    myWorker->Initialize(thread, top, top_cg, map, traj_reader, _do_mapping, &nframes,              \
                            nframesMutex, traj_readerMutex);

                    //seek to first frame, let thread0 do that
                    if (thread == 0) {
                        for (bool bok = true; bok == true; bok = traj_reader->NextFrame(*top)) {
                            if (((top->getTime() < begin) && has_begin) || first_frame > 1) {
                                first_frame--;
                                continue;
                            }
                            break;

                        }
                    }
                    //start threads
                    myWorker->Start();
                }

                // hier noch nen block rein
                for (int thread = 0; thread < NUMBERTHREADS; thread++) {
                    //break;
                    _myWorkers[thread]->WaitDone();
                    MergeWorker(_myWorkers[thread]);
                    delete _myWorkers[thread];
                }
                //                    if (thread == 0) {
                EndEvaluate(); // einmal fuer application
                //                    }

                traj_reader->Close();
                delete traj_reader;
            }

            delete reader;

        }

        CsgApplication::Worker::~Worker() {
            if (_map)
                delete _map;
            delete _top;
            delete _top_cg;

        }

        void CsgApplication::BeginEvaluate(Topology *top, Topology * top_ref) {
            list<CGObserver *>::iterator iter;
            for (iter = _observers.begin(); iter != _observers.end(); ++iter)
                (*iter)->BeginCG(top, top_ref);
        }

        void CsgApplication::EndEvaluate() {
            list<CGObserver *>::iterator iter;
            for (iter = _observers.begin(); iter != _observers.end(); ++iter)
                (*iter)->EndCG();
        }

        void CsgApplication::EvalConfiguration(Topology *top, Topology * top_ref) {
            list<CGObserver *>::iterator iter;
            for (iter = _observers.begin(); iter != _observers.end(); ++iter)
                (*iter)->EvalConfiguration(top, top_ref);
        }


    }
}
