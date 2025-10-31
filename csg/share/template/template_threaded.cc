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

#include <cstdlib>
#include <memory>

#include <votca/csg/beadlist.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/histogram.h>

using namespace std;
using namespace votca::csg;

// comments were mainly added to explain the "overhead" needed for threaded
// calculations/analyzations

// to sum it up: instead of having one "thread" doing all your work (the whole
// tracetory), you may split it into single frames and distribute it among many
// "workers". a solid choice is: number of cores = number of workers. you, as
// the user, are required to define how to initialize and merge your workers.
// the main part of the program, EvalConfiguration, is shifted to the Worker
// class but other than that stays untouched compared to a non-threaded version

class CsgParallelTestApp : public CsgApplication {

  string ProgramName() override { return "template_threaded_rdf"; }

  void HelpText(ostream &out) override {
    out << "template for threaded rdf calculations";
  }

  void Initialize() override;

  bool DoTrajectory() override { return true; }

  // explicitly turn on threaded mode by overriding DoThreaded() and returning
  // true note that threads will be started and merged in an ordered way by
  // default this has the disadvantage of slowing everything down a bit (you
  // will likely not notice a decrease of performance), but the advantage of
  // processing frames in their original order in most cases, you want that in
  // some cases, where reading and writing/merging does not have to occur in
  // order, you may consider switching SynchronizeThreads() off in this example,
  // where an rdf-like value is calculated, ordered reading/writing is not
  // neccessary. however, leave it untouched to prevent future mistakes

  bool DoThreaded() override { return true; }
  // are you sure? really?
  // bool SynchronizeThreads() {
  //      return false;
  //  }

  void BeginEvaluate(Topology *top, Topology *top_ref) override;
  void EndEvaluate() override;

  // ForkWorker is the function you need to override and initialize your workers
  std::unique_ptr<CsgApplication::Worker> ForkWorker(void) override;

  // MergeWorker needs you to define how to merge different workers and their
  // data
  void MergeWorker(Worker *worker) override;

 protected:
  // data belonging to the main class CsgParallelTestApp
  votca::tools::HistogramNew rdf_;
  double cut_off_;
};

// derive from CsgApplication::Worker and define your worker
class RDFWorker : public CsgApplication::Worker {
 public:
  ~RDFWorker() override = default;
  // override EvalConfiguration with your analysis routine
  void EvalConfiguration(Topology *, Topology *) override;
  // data belonging to this particular worker
  votca::tools::HistogramNew rdf_;
  double cut_off_;
};

int main(int argc, char **argv) {
  CsgParallelTestApp app;

  return app.Exec(argc, argv);
}

void CsgParallelTestApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("RDF options")(
      "c", boost::program_options::value<double>()->default_value(1.0),
      "the cutoff");
}

void CsgParallelTestApp::BeginEvaluate(Topology *, Topology *) {
  cut_off_ = OptionsMap()["c"].as<double>();
  rdf_.Initialize(0, cut_off_, 50);
}

// create and initialize single workers
// ForkWorker() will be called as often as the parameter '--nt NTHREADS'
// it creates a new worker and the user is required to initialize variables etc.
// (if needed)
std::unique_ptr<CsgApplication::Worker> CsgParallelTestApp::ForkWorker() {
  auto worker = std::make_unique<RDFWorker>();
  // initialize
  worker->cut_off_ = OptionsMap()["c"].as<double>();
  worker->rdf_.Initialize(0, worker->cut_off_, 50);
  return worker;
}

// EvalConfiguration does the actual calculation
// you won't see any explicit threaded stuff here
void RDFWorker::EvalConfiguration(Topology *top, Topology *) {
  BeadList b;
  b.Generate(*top, "*");
  NBListGrid nb;
  nb.setCutoff(cut_off_);
  nb.Generate(b);
  for (auto &pair : nb) {
    rdf_.Process(pair->dist());
  }
}

// the user is required to define how to merge the single data
// belonging to each thread into the main data belonging to CsgParallelTestApp
void CsgParallelTestApp::MergeWorker(Worker *worker) {
  RDFWorker *myRDFWorker;
  // cast generel Worker into your derived worker class(here RDFWorker)
  myRDFWorker = dynamic_cast<RDFWorker *>(worker);

  // the next comment block explains how mutexes are used internally for this
  // function: mutexes are used to exclusively work on data e.g., if you read or
  // write global data, make sure that nobody else (i.e. no other worker) works
  // on that very same piece of data at the same time; otherwise, you will end
  // up with wrong results that you struggle to understand the parent class
  // handles a "merging mutex" for you internally; this is what happens: first,
  // a mutex is created, e.g.
  //        Mutex rdfMutex;
  // then, for each worker, the mutex is first locked
  //        rdfMutex.Lock())
  // and MergeWorker(worker) is called (i.e. the code you define here is
  // executed) after MergeWorker exits, the mutex is unlocked
  //        rdfMutex.Unlock();
  // and allows other threads to get a lock and start merging

  // now follows your code

  // merging of data in this simple example is easy and does not have to follow
  // the original order of frames (since plain summing is commutative)
  rdf_.data().y() = rdf_.data().y() + myRDFWorker->rdf_.data().y();
}

void CsgParallelTestApp::EndEvaluate() {
  rdf_.data().y() = rdf_.data().y().cwiseQuotient(rdf_.data().x().cwiseAbs2());
  rdf_.data().Save("rdf.dat");
}
