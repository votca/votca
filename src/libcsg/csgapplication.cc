/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// Third party includes
#include <boost/algorithm/string/trim.hpp>
#include <memory>

// Local VOTCA includes
#include "votca/csg/cgengine.h"
#include "votca/csg/csgapplication.h"
#include "votca/csg/topologymap.h"
#include "votca/csg/topologyreader.h"
#include "votca/csg/trajectoryreader.h"
#include "votca/csg/trajectorywriter.h"
#include "votca/csg/version.h"

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
          "  list of molecules to ignore separated by ;")("no-map",
                                                          "  disable mapping "
                                                          "and act on original "
                                                          "trajectory");
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

  if (DoTrajectory()) {
    AddProgramOptions("Trajectory options")(
        "trj", boost::program_options::value<std::string>(),
        "  atomistic trajectory file")(
        "begin", boost::program_options::value<double>()->default_value(0.0),
        "  skip frames before this time (only works for Gromacs files)")(
        "first-frame", boost::program_options::value<Index>()->default_value(0),
        "  start with this frame")("nframes",
                                   boost::program_options::value<Index>(),
                                   "  process the given number of frames");
  }

  if (DoThreaded()) {
    /*
     * TODO default value of 1 for nt is not smart
     */
    AddProgramOptions("Threading options")(
        "nt", boost::program_options::value<Index>()->default_value(1),
        "  number of threads");
  }
}

bool CsgApplication::EvaluateOptions() {
  do_mapping_ = false;
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
        do_mapping_ = true;
      }
      if (OptionsMap().count("no-map") && OptionsMap().count("cg")) {
        ShowHelpText(std::cout);
        throw std::runtime_error(
            "no-map and cg options are mutually exclusive!");
      }
    }  // default mapping is off, if user gives cg, then do mapping
    else if (OptionsMap().count("cg")) {
      do_mapping_ = true;
    }
  }

  /* check threading options */
  if (DoThreaded()) {
    nthreads_ = OptionsMap()["nt"].as<Index>();
    /* TODO
     * does the number of threads make sense?
     * which criteria should be used? smaller than system's cores?
     */
  }

  return true;
}

void CsgApplication::ShowHelpText(std::ostream &out) {
  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }

  HelpTextHeader(name);
  HelpText(out);

  out << "\n\n" << VisibleOptions() << std::endl;
}

void CsgApplication::Worker::Run() {
  while (app_->ProcessData(this)) {
    if (app_->SynchronizeThreads()) {
      Index id = getId();
      app_->threadsMutexesOut_[id]->Lock();
      app_->MergeWorker(this);
      app_->threadsMutexesOut_[(id + 1) % app_->nthreads_]->Unlock();
    }
  }
}

bool CsgApplication::ProcessData(Worker *worker) {

  Index id;
  id = worker->getId();

  if (SynchronizeThreads()) {
    // wait til its your turn
    threadsMutexesIn_[id]->Lock();
  }
  traj_readerMutex_.Lock();
  if (nframes_ == 0) {
    traj_readerMutex_.Unlock();

    if (SynchronizeThreads()) {
      // done processing? don't forget to unlock next worker anyway
      threadsMutexesIn_[(id + 1) % nthreads_]->Unlock();
    }

    return false;
  }
  nframes_--;
  if (!is_first_frame_ || worker->getId() != 0) {
    // get frame
    bool tmpRes = traj_reader_->NextFrame(worker->top_);
    if (!tmpRes) {
      traj_readerMutex_.Unlock();
      if (SynchronizeThreads()) {
        threadsMutexesIn_[(id + 1) % nthreads_]->Unlock();
      }
      return false;
    }
  }
  if (worker->getId() == 0) {
    is_first_frame_ = false;
  }

  traj_readerMutex_.Unlock();
  if (SynchronizeThreads()) {
    // unlock next frame for input
    threadsMutexesIn_[(id + 1) % nthreads_]->Unlock();
  }
  // evaluate
  if (do_mapping_) {
    worker->map_->Apply();
    worker->EvalConfiguration(&worker->top_cg_, &worker->top_);
  } else {
    worker->EvalConfiguration(&worker->top_);
  }

  return true;
}

void CsgApplication::Run(void) {
  // create reader for atomistic topology
  std::unique_ptr<TopologyReader> reader =
      TopReaderFactory().Create(OptionsMap()["top"].as<std::string>());
  if (reader == nullptr) {
    throw std::runtime_error(std::string("input format not supported: ") +
                             OptionsMap()["top"].as<std::string>());
  }

  class DummyWorker : public Worker {
   public:
    void EvalConfiguration(Topology *top, Topology *top_ref) override {
      app_->EvalConfiguration(top, top_ref);
    }
  };

  // create the master worker
  if (DoThreaded()) {
    myWorkers_.push_back(ForkWorker());
  } else {
    myWorkers_.emplace_back(std::make_unique<DummyWorker>());
  }
  myWorkers_.back()->setApplication(this);
  myWorkers_.back()->setId(0);
  Worker *master = (myWorkers_.back().get());

  CGEngine cg;

  //////////////////////////////////////////////////
  // read in the topology for master
  //////////////////////////////////////////////////
  reader->ReadTopology(OptionsMap()["top"].as<std::string>(), master->top_);
  // Ensure that the coarse grained topology will have the same boundaries
  master->top_cg_.setBox(master->top_.getBox());

  std::cout << "I have " << master->top_.BeadCount() << " beads in "
            << master->top_.MoleculeCount() << " molecules" << std::endl;
  master->top_.CheckMoleculeNaming();

  if (do_mapping_) {
    // read in the coarse graining definitions (xml files)
    cg.LoadMoleculeType(OptionsMap()["cg"].as<std::string>());
    // create the mapping + cg topology

    if (OptionsMap().count("map-ignore") != 0) {
      tools::Tokenizer tok(OptionsMap()["map-ignore"].as<std::string>(), ";");
      for (std::string str : tok) {
        boost::trim(str);
        if (str.length() > 0) {
          cg.AddIgnore(str);
        }
      }
    }

    master->map_ = cg.CreateCGTopology(master->top_, master->top_cg_);

    std::cout << "I have " << master->top_cg_.BeadCount() << " beads in "
              << master->top_cg_.MoleculeCount()
              << " molecules for the coarsegraining" << std::endl;

    // If the trajectory reader is off but mapping flag is specified do apply
    // the mapping, this switch is necessary in cases where xml files are
    // specified, which do not contain positional information. In such cases
    // it is not possible to apply the positional mapping, a trajectory file
    // must be read in.
    if (DoTrajectory() == false) {
      master->map_->Apply();
      if (!EvaluateTopology(&master->top_cg_, &master->top_)) {
        return;
      }
    }
  } else if (!EvaluateTopology(&master->top_)) {
    return;
  }

  //////////////////////////////////////////////////
  // Here trajectory parsing starts
  //////////////////////////////////////////////////
  if (DoTrajectory() && OptionsMap().count("trj")) {
    double begin = 0;
    Index first_frame;
    bool has_begin = false;

    if (OptionsMap().count("begin")) {
      has_begin = true;
      begin = OptionsMap()["begin"].as<double>();
    }

    nframes_ = -1;
    if (OptionsMap().count("nframes")) {
      nframes_ = OptionsMap()["nframes"].as<Index>();
    }

    first_frame = OptionsMap()["first-frame"].as<Index>();

    // create reader for trajectory
    traj_reader_ =
        TrjReaderFactory().Create(OptionsMap()["trj"].as<std::string>());
    if (traj_reader_ == nullptr) {
      throw std::runtime_error(std::string("input format not supported: ") +
                               OptionsMap()["trj"].as<std::string>());
    }
    // open the trajectory
    traj_reader_->Open(OptionsMap()["trj"].as<std::string>());

    //////////////////////////////////////////////////
    // Create all the workers
    /////////////////verbose/////////////////////////////////
    for (Index thread = 1; thread < nthreads_ && DoThreaded(); thread++) {
      myWorkers_.push_back(ForkWorker());
      myWorkers_.back()->setApplication(this);
      myWorkers_.back()->setId(thread);

      // this will be changed to CopyTopologyData
      // read in the topology

      reader->ReadTopology(OptionsMap()["top"].as<std::string>(),
                           myWorkers_.back()->top_);
      myWorkers_.back()->top_.CheckMoleculeNaming();

      if (do_mapping_) {
        // create the mapping + cg topology
        myWorkers_.back()->map_ = cg.CreateCGTopology(
            myWorkers_.back()->top_, myWorkers_.back()->top_cg_);
      }
    }

    //////////////////////////////////////////////////
    // Proceed to first frame of interest
    //////////////////////////////////////////////////

    traj_reader_->FirstFrame(master->top_);
    if (master->top_.getBoxType() == BoundaryCondition::typeOpen) {
      std::cout
          << "NOTE: You are using OpenBox boundary conditions. Check if this "
             "is intended.\n"
          << std::endl;
    }
    // seek first frame, let thread0 do that
    bool bok;
    for (bok = true; bok == true; bok = traj_reader_->NextFrame(master->top_)) {
      if ((has_begin && (master->top_.getTime() < begin)) || first_frame > 1) {
        first_frame--;
        continue;
      }
      break;
    }
    if (!bok) {  // trajectory was too short and we did not proceed to first
                 // frame
      traj_reader_->Close();

      throw std::runtime_error(
          "trajectory was too short, did not process a single frame");
    }

    // notify all observers that coarse graining has begun
    if (do_mapping_) {
      master->map_->Apply();
      BeginEvaluate(&master->top_cg_, &master->top_);
    } else {
      BeginEvaluate(&master->top_);
    }

    is_first_frame_ = true;
    /////////////////////////////////////////////////////////////////////////
    // start threads
    if (DoThreaded()) {
      for (size_t thread = 0; thread < myWorkers_.size(); thread++) {

        if (SynchronizeThreads()) {
          threadsMutexesIn_.push_back(std::make_unique<tools::Mutex>());
          // lock each worker for input
          threadsMutexesIn_.back()->Lock();

          threadsMutexesOut_.push_back(std::make_unique<tools::Mutex>());
          // lock each worker for output
          threadsMutexesOut_.back()->Lock();
        }
      }
      for (auto &myWorker_ : myWorkers_) {
        myWorker_->Start();
      }

      if (SynchronizeThreads()) {
        // unlock first thread and start ordered input/output
        threadsMutexesIn_[0]->Unlock();
        threadsMutexesOut_[0]->Unlock();
      }
      // mutex needed for merging if SynchronizeThreads()==False
      tools::Mutex mergeMutex;
      for (auto &myWorker : myWorkers_) {
        myWorker->WaitDone();
        if (!SynchronizeThreads()) {
          mergeMutex.Lock();
          MergeWorker(myWorker.get());
          mergeMutex.Unlock();
        }
      }

    } else {
      master->Start();
      master->WaitDone();
    }

    EndEvaluate();

    myWorkers_.clear();
    threadsMutexesIn_.clear();
    threadsMutexesOut_.clear();
    traj_reader_->Close();
  }
}

void CsgApplication::BeginEvaluate(Topology *top, Topology *top_ref) {
  for (CGObserver *ob : observers_) {
    ob->BeginCG(top, top_ref);
  }
}

void CsgApplication::EndEvaluate() {
  for (CGObserver *ob : observers_) {
    ob->EndCG();
  }
}

void CsgApplication::EvalConfiguration(Topology *top, Topology *top_ref) {
  for (CGObserver *ob : observers_) {
    ob->EvalConfiguration(top, top_ref);
  }
}

std::unique_ptr<CsgApplication::Worker> CsgApplication::ForkWorker(void) {
  throw std::runtime_error("ForkWorker not implemented in application");
  return nullptr;
}

void CsgApplication::MergeWorker(CsgApplication::Worker *) {
  throw std::runtime_error("MergeWorker not implemented in application");
}

}  // namespace csg
}  // namespace votca
