/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <boost/algorithm/string/trim.hpp>
#include <votca/csg/cgengine.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/topologymap.h>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/version.h>

namespace votca {
namespace csg {

void CsgApplication::Initialize() {
  // register all io plugins
  TrajectoryWriter::RegisterPlugins();
  TrajectoryReader::RegisterPlugins();
  TopologyReader::RegisterPlugins();

  if (NeedsTopology()) {
    AddProgramOptions()("top", boost::program_options::value<std::string>(),
                        "  atomistic topology file");
  }
  if (DoMapping()) {
    if (DoMappingDefault()) {
      AddProgramOptions("Mapping options")(
          "cg", boost::program_options::value<std::string>(),
          "  coarse graining mapping and bond definitions (xml-file)")(
          "map-ignore", boost::program_options::value<std::string>(),
          "  list of molecules to ignore separated by ;")(
          "no-map", "  disable mapping and act on original trajectory");
    } else {
      AddProgramOptions("Mapping options")(
          "cg", boost::program_options::value<std::string>(),
          "  [OPTIONAL] coarse graining mapping and bond definitions\n"
          "  (xml-file). If no file is given, program acts on original "
          "trajectory")(
          "map-ignore", boost::program_options::value<std::string>(),
          "  list of molecules to ignore if mapping is done separated by ;");
    }
  }

  if (DoTrajectory())
    AddProgramOptions("Trajectory options")(
        "trj", boost::program_options::value<std::string>(),
        "  atomistic trajectory file")(
        "begin", boost::program_options::value<double>()->default_value(0.0),
        "  skip frames before this time (only works for Gromacs files)")(
        "first-frame", boost::program_options::value<int>()->default_value(0),
        "  start with this frame")("nframes",
                                   boost::program_options::value<int>(),
                                   "  process the given number of frames");

  if (DoThreaded())
    /*
     * TODO default value of 1 for nt is not smart
     */
    AddProgramOptions("Threading options")(
        "nt", boost::program_options::value<int>()->default_value(1),
        "  number of threads");
}

bool CsgApplication::EvaluateOptions() {
  _do_mapping = false;
  if (NeedsTopology()) {
    CheckRequired("top", "no topology file specified");
  }

  // check for mapping options
  if (DoMapping()) {
    // default mapping is on
    if (DoMappingDefault()) {
      // if the user does not explicitly ask to turn it off, cg is needed
      if (OptionsMap().count("no-map") == 0) {
        CheckRequired("cg", "no coarse graining definition specified");
        _do_mapping = true;
      }
      if (OptionsMap().count("no-map") && OptionsMap().count("cg")) {
        ShowHelpText(std::cout);
        throw std::runtime_error(
            "no-map and cg options are mutually exclusive!");
      }
    }  // default mapping is off, if user gives cg, then do mapping
    else if (OptionsMap().count("cg")) {
      _do_mapping = true;
    }
  }

  /* check threading options */
  if (DoThreaded()) {
    _nthreads = _op_vm["nt"].as<int>();
    /* TODO
     * does the number of threads make sense?
     * which criteria should be used? smaller than system's cores?
     */
  }

  return true;
}

void CsgApplication::ShowHelpText(std::ostream &out) {
  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();

  HelpTextHeader(name);
  HelpText(out);

  out << "\n\n" << VisibleOptions() << std::endl;
}

void CsgApplication::Worker::Run() {
  while (_app->ProcessData(this)) {
    if (_app->SynchronizeThreads()) {
      int id = getId();
      _app->_threadsMutexesOut[id]->Lock();
      _app->MergeWorker(this);
      _app->_threadsMutexesOut[(id + 1) % _app->_nthreads]->Unlock();
    }
  }
}

bool CsgApplication::ProcessData(Worker *worker) {

  int id;
  id = worker->getId();

  if (SynchronizeThreads()) {
    // wait til its your turn
    _threadsMutexesIn[id]->Lock();
  }
  _traj_readerMutex.Lock();
  if (_nframes == 0) {
    _traj_readerMutex.Unlock();

    if (SynchronizeThreads()) {
      // done processing? don't forget to unlock next worker anyway
      _threadsMutexesIn[(id + 1) % _nthreads]->Unlock();
    }

    return false;
  }
  _nframes--;
  if (!_is_first_frame || worker->getId() != 0) {
    // get frame
    bool tmpRes = _traj_reader->NextFrame(worker->_top);
    if (!tmpRes) {
      _traj_readerMutex.Unlock();
      if (SynchronizeThreads())
        _threadsMutexesIn[(id + 1) % _nthreads]->Unlock();
      return false;
    }
  }
  if (worker->getId() == 0) _is_first_frame = false;

  _traj_readerMutex.Unlock();
  if (SynchronizeThreads()) {
    // unlock next frame for input
    _threadsMutexesIn[(id + 1) % _nthreads]->Unlock();
  }
  // evaluate
  if (_do_mapping) {
    worker->_map->Apply();
    worker->EvalConfiguration(&worker->_top_cg, &worker->_top);
  } else
    worker->EvalConfiguration(&worker->_top);

  return true;
}

void CsgApplication::Run(void) {
  TopologyReader *reader;
  // create reader for atomistic topology
  reader = TopReaderFactory().Create(_op_vm["top"].as<std::string>());
  if (reader == nullptr)
    throw std::runtime_error(std::string("input format not supported: ") +
                             _op_vm["top"].as<std::string>());

  class DummyWorker : public Worker {
   public:
    void EvalConfiguration(Topology *top, Topology *top_ref) override {
      _app->EvalConfiguration(top, top_ref);
    }
  };

  // create the master worker
  Worker *master = nullptr;
  if (DoThreaded())
    master = ForkWorker();
  else
    master = new DummyWorker();

  master->setApplication(this);
  master->setId(0);
  _myWorkers.push_back(master);

  CGEngine cg;

  //////////////////////////////////////////////////
  // read in the topology for master
  //////////////////////////////////////////////////
  reader->ReadTopology(_op_vm["top"].as<std::string>(), master->_top);
  std::cout << "I have " << master->_top.BeadCount() << " beads in "
            << master->_top.MoleculeCount() << " molecules" << std::endl;
  master->_top.CheckMoleculeNaming();

  if (_do_mapping) {
    // read in the coarse graining definitions (xml files)
    cg.LoadMoleculeType(_op_vm["cg"].as<std::string>());
    // create the mapping + cg topology

    if (_op_vm.count("map-ignore") != 0) {
      tools::Tokenizer tok(_op_vm["map-ignore"].as<std::string>(), ";");
      for (std::string str : tok) {
        boost::trim(str);
        if (str.length() > 0) cg.AddIgnore(str);
      }
    }

    master->_map = cg.CreateCGTopology(master->_top, master->_top_cg);

    std::cout << "I have " << master->_top_cg.BeadCount() << " beads in "
              << master->_top_cg.MoleculeCount()
              << " molecules for the coarsegraining" << std::endl;
    master->_map->Apply();
    if (!EvaluateTopology(&master->_top_cg, &master->_top)) return;
  } else if (!EvaluateTopology(&master->_top))
    return;

  //////////////////////////////////////////////////
  // Here trajectory parsing starts
  //////////////////////////////////////////////////
  if (DoTrajectory() && _op_vm.count("trj")) {
    double begin = 0;
    int first_frame;
    bool has_begin = false;

    if (_op_vm.count("begin")) {
      has_begin = true;
      begin = _op_vm["begin"].as<double>();
    }

    _nframes = -1;
    if (_op_vm.count("nframes")) {
      _nframes = _op_vm["nframes"].as<int>();
    }

    first_frame = _op_vm["first-frame"].as<int>();

    // create reader for trajectory
    _traj_reader = TrjReaderFactory().Create(_op_vm["trj"].as<std::string>());
    if (_traj_reader == nullptr)
      throw std::runtime_error(std::string("input format not supported: ") +
                               _op_vm["trj"].as<std::string>());
    // open the trajectory
    _traj_reader->Open(_op_vm["trj"].as<std::string>());

    //////////////////////////////////////////////////
    // Create all the workers
    /////////////////verbose/////////////////////////////////
    for (int thread = 1; thread < _nthreads && DoThreaded(); thread++) {
      Worker *myWorker = ForkWorker();
      myWorker->setApplication(this);
      myWorker->setId(thread);
      _myWorkers.push_back(myWorker);

      // this will be changed to CopyTopologyData
      // read in the topology

      reader->ReadTopology(_op_vm["top"].as<std::string>(), myWorker->_top);
      myWorker->_top.CheckMoleculeNaming();

      if (_do_mapping) {
        // create the mapping + cg topology
        myWorker->_map = cg.CreateCGTopology(myWorker->_top, myWorker->_top_cg);
      }
    }

    //////////////////////////////////////////////////
    // Proceed to first frame of interest
    //////////////////////////////////////////////////

    _traj_reader->FirstFrame(master->_top);
    if (master->_top.getBoxType() == BoundaryCondition::typeOpen) {
      std::cout
          << "NOTE: You are using OpenBox boundary conditions. Check if this "
             "is intended.\n"
          << std::endl;
    }
    // seek first frame, let thread0 do that
    bool bok;
    for (bok = true; bok == true; bok = _traj_reader->NextFrame(master->_top)) {
      if ((has_begin && (master->_top.getTime() < begin)) || first_frame > 1) {
        first_frame--;
        continue;
      }
      break;
    }
    if (!bok) {  // trajectory was too short and we did not proceed to first
                 // frame
      _traj_reader->Close();
      delete _traj_reader;

      throw std::runtime_error(
          "trajectory was too short, did not process a single frame");
    }

    // notify all observers that coarse graining has begun
    if (_do_mapping) {
      master->_map->Apply();
      BeginEvaluate(&master->_top_cg, &master->_top);
    } else
      BeginEvaluate(&master->_top);

    _is_first_frame = true;
    /////////////////////////////////////////////////////////////////////////
    // start threads
    if (DoThreaded()) {
      for (size_t thread = 0; thread < _myWorkers.size(); thread++) {

        if (SynchronizeThreads()) {
          tools::Mutex *myMutexIn = new tools::Mutex;
          _threadsMutexesIn.push_back(myMutexIn);
          // lock each worker for input
          myMutexIn->Lock();

          tools::Mutex *myMutexOut = new tools::Mutex;
          _threadsMutexesOut.push_back(myMutexOut);
          // lock each worker for output
          myMutexOut->Lock();
        }
      }
      for (size_t thread = 0; thread < _myWorkers.size(); thread++)
        _myWorkers[thread]->Start();

      if (SynchronizeThreads()) {
        // unlock first thread and start ordered input/output
        _threadsMutexesIn[0]->Unlock();
        _threadsMutexesOut[0]->Unlock();
      }
      // mutex needed for merging if SynchronizeThreads()==False
      tools::Mutex mergeMutex;
      for (size_t thread = 0; thread < _myWorkers.size(); thread++) {
        _myWorkers[thread]->WaitDone();
        if (!SynchronizeThreads()) {
          mergeMutex.Lock();
          MergeWorker(_myWorkers[thread]);
          mergeMutex.Unlock();
        }
        delete _myWorkers[thread];
      }
      for (size_t thread = 0; thread < _threadsMutexesIn.size(); ++thread) {
        delete _threadsMutexesIn[thread];
        delete _threadsMutexesOut[thread];
      }

    } else {
      master->Start();
      master->WaitDone();
      delete master;
    }

    EndEvaluate();

    _myWorkers.clear();
    _threadsMutexesIn.clear();
    _threadsMutexesOut.clear();
    _traj_reader->Close();

    delete _traj_reader;
  }

  delete reader;
}

CsgApplication::Worker::~Worker() {
  if (_map) delete _map;
}

void CsgApplication::BeginEvaluate(Topology *top, Topology *top_ref) {
  for (CGObserver *ob : _observers) {
    ob->BeginCG(top, top_ref);
  }
}

void CsgApplication::EndEvaluate() {
  for (CGObserver *ob : _observers) {
    ob->EndCG();
  }
}

void CsgApplication::EvalConfiguration(Topology *top, Topology *top_ref) {
  for (CGObserver *ob : _observers) {
    ob->EvalConfiguration(top, top_ref);
  }
}

CsgApplication::Worker *CsgApplication::ForkWorker(void) {
  throw std::runtime_error("ForkWorker not implemented in application");
  return nullptr;
}

void CsgApplication::MergeWorker(CsgApplication::Worker *worker) {
  throw std::runtime_error("MergeWorker not implemented in application");
}

}  // namespace csg
}  // namespace votca
