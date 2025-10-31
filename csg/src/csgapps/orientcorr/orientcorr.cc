/*
 * Copyright 2009-2025 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <cstdlib>
#include <memory>

// VOTCA includes
#include <votca/tools/histogram.h>

// Local VOTCA includes
#include <votca/csg/beadlist.h>
#include <votca/csg/csgapplication.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>

using namespace std;
using namespace votca::csg;

class OrientCorrApp : public CsgApplication {

  string ProgramName() override { return "orientcorr"; }

  void HelpText(ostream &out) override {
    out << "Calculates the orientational correlation function\n"
           "<3/2*u(0)*u(r) - 1/2>\n"
           "for a polymer melt, where u is the vector pointing along a bond "
           "and \n"
           "r the distance between bond segments (centered on middle of "
           "bond).\n\n"
           "The output is correlation.dat (with intra-molecular contributions) "
           "and\n"
           "correlation_excl.dat, where inter-molecular contributions are "
           "excluded.";
  }

  void Initialize() override;

  bool DoTrajectory() override { return true; }
  // do a threaded analyzis, splitting in time domain
  bool DoThreaded() override { return true; }
  // we don't care about the order of the analyzed frames
  bool SynchronizeThreads() override { return false; }

  // called before groing through all frames
  void BeginEvaluate(Topology *top, Topology *top_ref) override;
  // called after all frames were parsed
  void EndEvaluate() override;

  // creates a worker for a thread
  std::unique_ptr<CsgApplication::Worker> ForkWorker(void) override;
  // merge data of worker into main
  void MergeWorker(Worker *worker) override;

 public:
  // helper class to choose nbsearch algorithm
  static std::unique_ptr<NBList> CreateNBSearch();

 protected:
  votca::tools::Histogram cor_;
  votca::tools::Histogram count_;
  votca::tools::Histogram cor_excl_;
  votca::tools::Histogram count_excl_;
  static string nbmethod_;
  double cut_off_;
  votca::Index nbins_;
};

string OrientCorrApp::nbmethod_;

// Earch thread has a worker and analysis data
class MyWorker : public CsgApplication::Worker {
 public:
  ~MyWorker() override = default;

  // evaluate the current frame
  void EvalConfiguration(Topology *top, Topology *top_ref) override;

  // callback if neighborsearch finds a pair
  bool FoundPair(Bead *b1, Bead *b2, const Eigen::Vector3d &r,
                 const double dist);

  // accumulator of the 3/2*u(0)u(r) - 1/2
  votca::tools::Histogram cor_;
  // number of hits for each bin
  votca::tools::Histogram count_;
  // accumulator of the 3/2*u(0)u(r) - 1/2, only inter-molecular
  votca::tools::Histogram cor_excl_;
  // number of hits for each bin, only inter-molecular
  votca::tools::Histogram count_excl_;
  double cut_off_;
};

int main(int argc, char **argv) {
  OrientCorrApp app;
  return app.Exec(argc, argv);
}

// add some program options
void OrientCorrApp::Initialize() {
  // add all standard application options
  CsgApplication::Initialize();
  // some application specific options
  AddProgramOptions("Neighbor search options")(
      "cutoff,c",
      boost::program_options::value<double>(&cut_off_)->default_value(1.0),
      "cutoff for the neighbor search")(
      "nbins",
      boost::program_options::value<votca::Index>(&nbins_)->default_value(40),
      "number of bins for the grid")(
      "nbmethod",
      boost::program_options::value<string>(&nbmethod_)->default_value("grid"),
      "neighbor search algorithm (simple or "
      "grid)");
}

std::unique_ptr<NBList> OrientCorrApp::CreateNBSearch() {
  if (nbmethod_ == "simple") {
    return std::make_unique<NBList>();
  }
  if (nbmethod_ == "grid") {
    return std::make_unique<NBListGrid>();
  }

  throw std::runtime_error(
      "unknown neighbor search method, use simple or grid");
  return nullptr;
}

// initialize the histograms
void OrientCorrApp::BeginEvaluate(Topology *, Topology *) {
  cor_.Initialize(0, cut_off_, nbins_);
  count_.Initialize(0, cut_off_, nbins_);
  cor_excl_.Initialize(0, cut_off_, nbins_);
  count_excl_.Initialize(0, cut_off_, nbins_);
}

// creates worker for each thread
std::unique_ptr<CsgApplication::Worker> OrientCorrApp::ForkWorker() {
  auto worker = std::make_unique<MyWorker>();
  worker->cut_off_ = cut_off_;
  worker->cor_.Initialize(0, worker->cut_off_, nbins_);
  worker->count_.Initialize(0, worker->cut_off_, nbins_);
  worker->cor_excl_.Initialize(0, worker->cut_off_, nbins_);
  worker->count_excl_.Initialize(0, worker->cut_off_, nbins_);
  return worker;
}

// evaluates a frame
void MyWorker::EvalConfiguration(Topology *top, Topology *) {

  // first genearate a mapped topology
  // the beads are sitting on the bonds and have an orientation which
  // is pointing along bond direction
  Topology mapped;
  cout << "generating mapped topology...";

  // copy box size
  mapped.setBox(top->getBox());

  // loop over all molecules
  for (const auto &mol_src : top->Molecules()) {
    // create a molecule in mapped topology
    Molecule *mol = mapped.CreateMolecule(mol_src.getName());
    // loop over beads in molecule
    for (votca::Index i = 0; i < mol_src.BeadCount() - 1; ++i) {
      // create a bead in mapped topology
      string bead_type = "A";
      if (mapped.BeadTypeExist(bead_type) == false) {
        mapped.RegisterBeadType(bead_type);
      }
      Bead *b =
          mapped.CreateBead(Bead::ellipsoidal, "A", bead_type, 1, 0.0, 0.0);
      Eigen::Vector3d p1 = mol_src.getBead(i)->getPos();
      Eigen::Vector3d p2 = mol_src.getBead(i + 1)->getPos();
      // position is in middle of bond
      Eigen::Vector3d pos = 0.5 * (p1 + p2);
      // orientation pointing along bond
      Eigen::Vector3d v = p2 - p1;
      v.normalize();
      b->setPos(pos);
      b->setV(v);
      mol->AddBead(b, "A");
    }
  }
  cout << "done\n";

  // the neighbor search only finds pairs, add self-self correlation parts here
  cor_.Process(0.0f, (double)mapped.BeadCount());
  count_.Process(0.0f, (double)mapped.BeadCount());

  // search for all beads
  BeadList b;
  b.Generate(mapped, "*");

  // create/initialize neighborsearch
  std::unique_ptr<NBList> nb = OrientCorrApp::CreateNBSearch();
  nb->setCutoff(cut_off_);

  // set callback for each pair found
  nb->SetMatchFunction(this, &MyWorker::FoundPair);

  // execute the search
  nb->Generate(b);
}

// process a pair, since return value is falsed, pairs are not cached which
// saves a lot of memory for the big systems
bool MyWorker::FoundPair(Bead *b1, Bead *b2, const Eigen::Vector3d &,
                         const double dist) {
  double tmp = b1->getV().dot(b2->getV());
  double P2 = 3. / 2. * tmp * tmp - 0.5;

  // calculate average without exclusions
  cor_.Process(dist, P2);
  count_.Process(dist);

  if (b1->getMoleculeId() == b2->getMoleculeId()) {
    return false;
  }

  // calculate average with excluding intramolecular contributions
  cor_excl_.Process(dist, P2);
  count_excl_.Process(dist);

  return false;
}

// merge analysed data of a worker into main applications
void OrientCorrApp::MergeWorker(Worker *worker) {
  MyWorker *myWorker = dynamic_cast<MyWorker *>(worker);
  cor_.data().y() += myWorker->cor_.data().y();
  count_.data().y() += myWorker->count_.data().y();
  cor_excl_.data().y() += myWorker->cor_excl_.data().y();
  count_excl_.data().y() += myWorker->count_excl_.data().y();
}

// write out the data
void OrientCorrApp::EndEvaluate() {
  cor_.data().y() = cor_.data().y().cwiseQuotient(count_.data().y());
  cor_.data().Save("correlation.dat");

  cor_excl_.data().y() =
      cor_excl_.data().y().cwiseQuotient(count_excl_.data().y());
  cor_excl_.data().Save("correlation_excl.dat");
}
