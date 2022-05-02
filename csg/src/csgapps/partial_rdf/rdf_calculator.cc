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

// Standard includes
#include <fstream>
#include <iomanip>
#include <iostream>

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/rangeparser.h>

// Local VOTCA includes
#include <votca/csg/beadlist.h>
#include <votca/csg/imcio.h>
#include <votca/csg/nblistgrid.h>

// Local private includes
#include "rdf_calculator.h"

namespace votca {
namespace csg {

RDFCalculator::RDFCalculator()
    : write_every_(0),
      do_blocks_(false),
      nblock_(0),
      subvol_rad_(0),
      do_vol_corr_(false),
      processed_some_frames_(false) {}

RDFCalculator::~RDFCalculator() = default;

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void RDFCalculator::Initialize() {
  // do some output
  std::cout << "begin to calculate distribution functions\n";
  std::cout << "# of bonded interactions: " << bonded_.size() << std::endl;
  std::cout << "# of non-bonded interactions: " << nonbonded_.size()
            << std::endl;

  if (bonded_.size() + nonbonded_.size() == 0) {
    throw std::runtime_error(
        "No interactions defined in options xml-file - nothing to be done");
  }

  // initialize non-bonded structures
  for (Property *prop : nonbonded_) {
    interaction_t *i = AddInteraction(prop);
    i->is_bonded_ = false;
  }
}

void RDFCalculator::BeginEvaluate(Topology *top, Topology *) {
  Eigen::Matrix3d box = top->getBox();
  boxc_ = box.rowwise().sum() / 2.0;

  std::cout << "Using center of box: " << boxc_ << std::endl;
  // we didn't process any frames so far
  nframes_ = 0;
  nblock_ = 0;
  processed_some_frames_ = false;

  // initialize non-bonded structures
  for (Property *prop : nonbonded_) {
    std::string name = prop->get("name").value();

    interaction_t &i = *interactions_[name];

    // count total species for ideal densities

    BeadList allbeads1, allbeads2;
    allbeads1.Generate(*top, prop->get("type1").value());
    allbeads2.Generate(*top, prop->get("type2").value());

    if (allbeads1.size() == 0) {
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type1").value() +
                               "\"\n"
                               "This was specified in type1 of interaction \"" +
                               name + "\"");
    }
    if (allbeads2.size() == 0) {
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type2").value() +
                               "\"\n"
                               "This was specified in type2 of interaction \"" +
                               name + "\"");
    }
    // calculate normalization factor for rdf

    if (do_vol_corr_) {
      std::cout << "Volume correction on" << std::endl;
      i.norm_ = 1. / (4.0 * votca::tools::conv::Pi * i.step_);
    } else {
      std::cout << "Volume correction off" << std::endl;
      i.norm_ = 1. / (4. * votca::tools::conv::Pi * i.step_);
    }
  }
}
// create an entry for interactions
RDFCalculator::interaction_t *RDFCalculator::AddInteraction(Property *p) {
  std::string name = p->get("name").value();
  std::string group;

  group = "none";

  auto inter = std::make_unique<interaction_t>();
  interaction_t *i = inter.get();
  inter->index_ = interactions_.size();
  interactions_[name] = std::move(inter);
  getGroup(group)->interactions_.push_back(i);

  i->step_ = p->get("step").as<double>();
  i->min_ = p->get("min").as<double>();
  i->max_ = p->get("max").as<double>();
  i->norm_ = 1.0;
  i->p_ = p;

  // initialize the current and average histogram
  Index n = (Index)((i->max_ - i->min_) / i->step_ + 1.000000001);

  i->average_.Initialize(i->min_, i->max_ + i->step_, n);

  return i;
}

// end of trajectory, post processing data
void RDFCalculator::EndEvaluate() {
  if (nframes_ > 0) {
    if (!do_blocks_) {
      WriteDist();
    }
  }
  // clear interactions and groups
  interactions_.clear();
  groups_.clear();
  if (!processed_some_frames_) {
    throw std::runtime_error(
        "no frames were processed. Please check your input");
  }
}

// load options from xml file
void RDFCalculator::LoadOptions(const std::string &file) {
  options_.LoadFromXML(file);
  bonded_ = options_.Select("cg.bonded");
  nonbonded_ = options_.Select("cg.non-bonded");
}

// evaluate current conformation
void RDFCalculator::Worker::EvalConfiguration(Topology *top, Topology *) {
  cur_vol_ = 4.0 / 3.0 * votca::tools::conv::Pi * rdfcalculator_->subvol_rad_ *
             rdfcalculator_->subvol_rad_ * rdfcalculator_->subvol_rad_;
  // process non-bonded interactions
  DoNonbonded(top);
  // process bonded interactions
  DoBonded(top);
}

void RDFCalculator::ClearAverages() {

  nframes_ = 0;
  for (auto &interaction_ : interactions_) {
    interaction_.second->average_.Clear();
  }

  for (auto &group_ : groups_) {
    group_.second->corr_.setZero();
  }
}

class IMCNBSearchHandler {
 public:
  IMCNBSearchHandler(HistogramNew *hist, double subvol_rad,
                     Eigen::Vector3d boxc, bool do_vol_corr)
      : hist_(hist),
        subvol_rad_(subvol_rad),
        boxc_(boxc),
        do_vol_corr_(do_vol_corr) {}

  HistogramNew *hist_;
  double subvol_rad_;
  Eigen::Vector3d boxc_;  // center of box
  bool do_vol_corr_;

  bool FoundPair(Bead *b1, Bead *, const Eigen::Vector3d &, const double dist) {

    if (do_vol_corr_) {
      double dr = (b1->Pos() - boxc_).norm();
      if (dist + dr > subvol_rad_) {
        // 2.0 is because everything is normalized to 4 PI
        hist_->Process(dist, 2.0 / SurfaceRatio(dist, dr));
      } else {
        hist_->Process(dist);
      }

    } else {
      hist_->Process(dist);
    }
    return false;
  }

  double SurfaceRatio(double dist, double r) {
    // r: distance of particle from ex center
    // dist: distance between particles
    return (1.0 + (subvol_rad_ * subvol_rad_ - r * r - dist * dist) /
                      (2.0 * dist * r));
  }
};

// process non-bonded interactions for current frame
void RDFCalculator::Worker::DoNonbonded(Topology *top) {
  for (Property *prop : rdfcalculator_->nonbonded_) {
    std::string name = prop->get("name").value();

    interaction_t &i = *rdfcalculator_->interactions_[name];

    // generate the bead lists
    BeadList beads1, beads2;

    beads1.GenerateInSphericalSubvolume(*top, prop->get("type1").value(),
                                        rdfcalculator_->boxc_,
                                        rdfcalculator_->subvol_rad_);
    beads2.GenerateInSphericalSubvolume(*top, prop->get("type2").value(),
                                        rdfcalculator_->boxc_,
                                        rdfcalculator_->subvol_rad_);

    cur_beadlist_1_count_ = (double)beads1.size();
    cur_beadlist_2_count_ = (double)beads2.size();

    // same types, so put factor 1/2 because of already counted interactions
    if (prop->get("type1").value() == prop->get("type2").value()) {
      cur_beadlist_2_count_ /= 2.0;
    }

    // generate the neighbour list
    std::unique_ptr<NBList> nb;

    bool gridsearch = true;

    if (rdfcalculator_->options_.exists("cg.nbsearch")) {
      if (rdfcalculator_->options_.get("cg.nbsearch").as<std::string>() ==
          "grid") {
        gridsearch = true;
      } else if (rdfcalculator_->options_.get("cg.nbsearch")
                     .as<std::string>() == "simple") {
        gridsearch = false;
      } else {
        throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
      }
    }
    if (gridsearch) {
      nb = std::make_unique<NBListGrid>();
    } else {
      nb = std::make_unique<NBList>();
    }

    nb->setCutoff(i.max_ + i.step_);

    // clear the current histogram
    current_hists_[i.index_].Clear();

    IMCNBSearchHandler h(&(current_hists_[i.index_]),
                         rdfcalculator_->subvol_rad_, rdfcalculator_->boxc_,
                         rdfcalculator_->do_vol_corr_);
    nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

    // is it same types or different types?
    if (prop->get("type1").value() == prop->get("type2").value()) {
      nb->Generate(beads1);
    } else {
      nb->Generate(beads1, beads2);
    }

    // store particle number in subvolume for each interaction
    i.avg_beadlist_1_count_.Process(cur_beadlist_1_count_);
    i.avg_beadlist_2_count_.Process(cur_beadlist_2_count_);
  }
}

// process non-bonded interactions for current frame
void RDFCalculator::Worker::DoBonded(Topology *top) {
  for (Property *prop : rdfcalculator_->bonded_) {
    std::string name = prop->get("name").value();

    interaction_t &i = *rdfcalculator_->interactions_[name];

    // clear the current histogram
    current_hists_[i.index_].Clear();

    // now fill with new data
    std::vector<Interaction *> vec = top->InteractionsInGroup(name);

    for (auto ic : vec) {
      double v = ic->EvaluateVar(*top);
      current_hists_[i.index_].Process(v);
    }
  }
}

// returns a group, creates it if doesn't exist
RDFCalculator::group_t *RDFCalculator::getGroup(const std::string &name) {
  std::map<std::string, std::unique_ptr<group_t>>::iterator iter;
  iter = groups_.find(name);
  if (iter == groups_.end()) {
    return (groups_[name] = std::make_unique<group_t>()).get();
  }
  return (*iter).second.get();
}

// write the distribution function
void RDFCalculator::WriteDist(const std::string &suffix) {

  // for all interactions
  for (auto &interaction_ : interactions_) {
    // calculate the rdf
    Table &t = interaction_.second->average_.data();
    Table dist(t);

    interaction_.second->norm_ /=
        (interaction_.second->avg_beadlist_1_count_.getAvg() *
         interaction_.second->avg_beadlist_2_count_.getAvg());
    dist.y() = avg_vol_.getAvg() * interaction_.second->norm_ *
               dist.y().cwiseQuotient(dist.x().cwiseAbs2());

    dist.Save((interaction_.first) + suffix + ".dist.new");
    std::cout << "written " << (interaction_.first) + suffix + ".dist.new\n";

    std::cout << "Avg. number of particles in subvol for "
              << (interaction_.first) << std::endl;
    std::cout << "beadlist 1: "
              << interaction_.second->avg_beadlist_1_count_.getAvg()
              << std::endl;
    std::cout << "beadlist 2: "
              << interaction_.second->avg_beadlist_2_count_.getAvg()
              << std::endl;
  }

  std::cout << "Volume used for normalization: " << avg_vol_.getAvg()
            << std::endl;
}

std::unique_ptr<CsgApplication::Worker> RDFCalculator::ForkWorker() {
  auto worker = std::make_unique<RDFCalculator::Worker>();

  worker->current_hists_.resize(interactions_.size());
  worker->rdfcalculator_ = this;

  for (auto &interaction_ : interactions_) {
    interaction_t *i = interaction_.second.get();
    worker->current_hists_[i->index_].Initialize(
        i->average_.getMin(), i->average_.getMax(), i->average_.getNBins());
  }
  return worker;
}

void RDFCalculator::MergeWorker(CsgApplication::Worker *worker_) {
  processed_some_frames_ = true;
  RDFCalculator::Worker *worker =
      dynamic_cast<RDFCalculator::Worker *>(worker_);
  // update the average

  ++nframes_;

  avg_vol_.Process(worker->cur_vol_);

  for (auto &interaction_ : interactions_) {
    interaction_t *i = interaction_.second.get();
    i->average_.data().y() =
        (((double)nframes_ - 1.0) * i->average_.data().y() +
         worker->current_hists_[i->index_].data().y()) /
        (double)nframes_;
  }

  if (write_every_ != 0) {
    if ((nframes_ % write_every_) == 0) {
      nblock_++;
      std::string suffix =
          std::string("_") + boost::lexical_cast<std::string>(nblock_);
      WriteDist(suffix);
      if (do_blocks_) {
        ClearAverages();
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
