/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include <stdlib.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/histogramnew.h>

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

class CsgTestApp : public CsgApplication {

  string ProgramName() { return "template_threaded_rdf"; }

  void HelpText(ostream &out) {
    out << "template for threaded rdf calculations";
  }

  void Initialize();

  bool DoTrajectory() { return true; }

  // explicitly turn on threaded mode by overriding DoThreaded() and returning
  // true note that threads will be started and merged in an ordered way by
  // default this has the disadvantage of slowing everything down a bit (you
  // will likely not notice a decrease of performance), but the advantage of
  // processing frames in their original order in most cases, you want that in
  // some cases, where reading and writing/merging does not have to occur in
  // order, you may consider switching SynchronizeThreads() off in this example,
  // where an rdf-like value is calculated, ordered reading/writing is not
  // neccessary. however, leave it untouched to prevent future mistakes

  bool DoThreaded() { return true; }
  // are you sure? really?
  // bool SynchronizeThreads() {
  //      return false;
  //  }

  void BeginEvaluate(Topology *top, Topology *top_ref);
  void EndEvaluate();

  // ForkWorker is the function you need to override and initialize your workers
  CsgApplication::Worker *ForkWorker(void);

  // MergeWorker needs you to define how to merge different workers and their
  // data
  void MergeWorker(Worker *worker);

 protected:
  // data belonging to the main class CsgTestApp
  HistogramNew _rdf;
  double _cut_off;
};

// derive from CsgApplication::Worker and define your worker
class RDFWorker : public CsgApplication::Worker {
 public:
  ~RDFWorker(){};
  // override EvalConfiguration with your analysis routine
  void EvalConfiguration(Topology *top, Topology *top_ref);
  // data belonging to this particular worker
  HistogramNew _rdf;
  double _cut_off;
};

int main(int argc, char **argv) {
  CsgTestApp app;

  return app.Exec(argc, argv);
}

void CsgTestApp::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("RDF options")(
      "c", boost::program_options::value<double>()->default_value(1.0),
      "the cutoff");
}

void CsgTestApp::BeginEvaluate(Topology *top, Topology *top_ref) {
  _cut_off = OptionsMap()["c"].as<double>();
  _rdf.Initialize(0, _cut_off, 50);
}

// create and initialize single workers
// ForkWorker() will be called as often as the parameter '--nt NTHREADS'
// it creates a new worker and the user is required to initialize variables etc.
// (if needed)
CsgApplication::Worker *CsgTestApp::ForkWorker() {
  RDFWorker *worker;
  worker = new RDFWorker();
  // initialize
  worker->_cut_off = OptionsMap()["c"].as<double>();
  worker->_rdf.Initialize(0, worker->_cut_off, 50);
  return worker;
}

// EvalConfiguration does the actual calculation
// you won't see any explicit threaded stuff here
void RDFWorker::EvalConfiguration(Topology *top, Topology *top_ref) {
  BeadList b;
  b.Generate(*top, "*");
  NBListGrid nb;
  nb.setCutoff(_cut_off);
  nb.Generate(b);
  NBList::iterator i;
  for (i = nb.begin(); i != nb.end(); ++i) {
    _rdf.Process((*i)->dist());
  }
}

// the user is required to define how to merge the single data
// belonging to each thread into the main data belonging to CsgTestApp
void CsgTestApp::MergeWorker(Worker *worker) {
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
  _rdf.data().y() = _rdf.data().y() + myRDFWorker->_rdf.data().y();
}

void CsgTestApp::EndEvaluate() {
  _rdf.data().y() = _rdf.data().y().cwiseQuotient(_rdf.data().x().cwiseAbs2());
  _rdf.data().Save("rdf.dat");
}
