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
                ("nt", boost::program_options::value<int>()->default_value(1), "  number of threads")
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

        //void CsgApplication::Worker::Initialize(int id, Topology top, Topology top_cg, TopologyMap * map) {

        void CsgApplication::Worker::Initialize(int id) {
            _myId = id;
        }

        void CsgApplication::Worker::Run(void) {
            while (_app->ProcessData(this));
        }

        bool CsgApplication::ProcessData(Worker * worker) {
            //easyLock
            //EasyLock locker(_traj_readerMutex);
            //return traj_reader->NextFrame(*top);

            //not so easy lock
            _nframesMutex.Lock();
            if (_nframes > 0 || _nframes <= -99) {
                _nframes--;
                _nframesMutex.Unlock();
                bool tmpRes;
                _traj_readerMutex.Lock();

                tmpRes = _traj_reader->NextFrame(worker->_top);

                _traj_readerMutex.Unlock();
                if (_do_mapping) {
                    worker->_map->Apply();
                    worker->EvalConfiguration(&worker->_top_cg, &worker->_top);
                } else
                    worker->EvalConfiguration(&worker->_top);
                return tmpRes;
            } else {
                _nframesMutex.Unlock();
                return false;
            }
        }

        void CsgApplication::RunThreaded(void) {
            Worker * myWorker;
            // first read in the topology
            CGEngine cg;

            TopologyReader *reader;

            // create reader for atomistic topology
            reader = TopReaderFactory().Create(_op_vm["top"].as<string > ());
            if (reader == NULL)
                throw runtime_error(string("input format not supported: ") + _op_vm["top"].as<string > ());

            if (DoTrajectory()) {
                double begin;
                int first_frame;
                bool has_begin = false;

                if (_op_vm.count("begin")) {
                    has_begin = true;
                    begin = _op_vm["begin"].as<double>();
                }

                _nframes = -99;
                if (_op_vm.count("nframes")) {
                    _nframes = _op_vm["nframes"].as<int>();
                }

                first_frame = _op_vm["first-frame"].as<int>();

                // create reader for trajectory
                _traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<string > ());
                if (_traj_reader == NULL)
                    throw runtime_error(string("input format not supported: ") + _op_vm["trj"].as<string > ());
                // open the trajectory
                _traj_reader->Open(_op_vm["trj"].as<string > ());


                //for first thread(=0): initialize
                myWorker = ForkWorker();
                myWorker->setApplication(this);
                _myWorkers.push_back(myWorker);

                // read in the topology
                reader->ReadTopology(_op_vm["top"].as<string > (), myWorker->_top);

                cout << "I have " << myWorker->_top.BeadCount() << " beads in " << myWorker->_top.MoleculeCount() << " molecules" << endl;
                myWorker->_top.CheckMoleculeNaming();

                if (_do_mapping) {
                    // read in the coarse graining definitions (xml files)
                    cg.LoadMoleculeType(_op_vm["cg"].as<string > ());
                    // create the mapping + cg topology
                    myWorker->_map = cg.CreateCGTopology(myWorker->_top, myWorker->_top_cg);
                    cout << "I have " << myWorker->_top_cg.BeadCount() << " beads in " << myWorker->_top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
                    myWorker->_map->Apply();

                    if (!EvaluateTopology(&myWorker->_top_cg, &myWorker->_top))
                        return;
                } else
                    if (!EvaluateTopology(&myWorker->_top))
                    return;

                _traj_reader->FirstFrame(myWorker->_top);

                // notify all observer that coarse graining has begun
                if (_do_mapping) {
                    myWorker->_map->Apply();
                    BeginEvaluate(&myWorker->_top_cg, &myWorker->_top);
                } else
                    BeginEvaluate(&myWorker->_top);


                //myWorker->Initialize(thread, top, top_cg, map);
                long thread = 0;
                myWorker->Initialize(thread);

                //seek to first frame, let thread0 do that
                for (bool bok = true; bok == true; bok = _traj_reader->NextFrame(myWorker->_top)) {
                    if (((myWorker->_top.getTime() < begin) && has_begin) || first_frame > 1) {
                        first_frame--;
                        continue;
                    }
                    break;
                }

                //for other threads: initialize
                for (int thread = 1; thread < _op_vm["nt"].as<int > (); thread++) {
                    myWorker = ForkWorker();
                    myWorker->setApplication(this);
                    _myWorkers.push_back(myWorker);


                    // read in the topology
                    reader->ReadTopology(_op_vm["top"].as<string > (), myWorker->_top);

                    //dirty check
                    //if (thread == 0)
                    //    cout << "I have " << myWorker->_top.BeadCount() << " beads in " << myWorker->_top.MoleculeCount() << " molecules" << endl;
                    myWorker->_top.CheckMoleculeNaming();


                    if (_do_mapping) {
                        // read in the coarse graining definitions (xml files)
                        //cg.LoadMoleculeType(_op_vm["cg"].as<string > ());
                        // create the mapping + cg topology
                        myWorker->_map = cg.CreateCGTopology(myWorker->_top, myWorker->_top_cg);
                        //dirty check
                        //if (thread == 0)
                        //    cout << "I have " << myWorker->_top_cg.BeadCount() << " beads in " << myWorker->_top_cg.MoleculeCount() << " molecules for the coarsegraining" << endl;
                        myWorker->_map->Apply();

                        if (!EvaluateTopology(&myWorker->_top_cg, &myWorker->_top))
                            return;
                    } else
                        if (!EvaluateTopology(&myWorker->_top))
                        return;

                    //if (thread == 0) {
                    // read in first frame
                    //_traj_reader->FirstFrame(myWorker->_top);

                    // notify all observer that coarse graining has begun

                    //if (_do_mapping) {
                    //    myWorker->_map->Apply();
                    //    BeginEvaluate(&myWorker->_top_cg, &myWorker->_top); // einmal fuer application
                    //} else
                    //    BeginEvaluate(&myWorker->_top); // einmal fuer application
                    //}

                    //myWorker->Initialize(thread, top, top_cg, map);
                    myWorker->Initialize(thread);

                    //seek to first frame, let thread0 do that
                    //if (thread == 0) {
                    //    for (bool bok = true; bok == true; bok = _traj_reader->NextFrame(myWorker->_top)) {
                    //        if (((myWorker->_top.getTime() < begin) && has_begin) || first_frame > 1) {
                    //            first_frame--;
                    //            continue;
                    //        }
                    //        break;
                    //    }
                    //}
                }

                //start threads
                for (long thread = 0; thread < _myWorkers.size(); thread++) {
                    _myWorkers[thread]->Start();
                }
                // hier noch nen block rein
                for (long thread = 0; thread < _myWorkers.size(); thread++) {
                    _myWorkers[thread]->WaitDone();
                    MergeWorker(_myWorkers[thread]);
                    delete _myWorkers[thread];
                }
                EndEvaluate();

                _traj_reader->Close();
                delete _traj_reader;
            }

            delete reader;

        }

        CsgApplication::Worker::~Worker() {
            if (_map)
                delete _map;

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
