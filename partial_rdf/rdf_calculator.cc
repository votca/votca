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

#include "rdf_calculator.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <votca/csg/beadlist.h>
#include <votca/csg/imcio.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/constants.h>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

RDFCalculator::RDFCalculator()
    : _write_every(0),
      _do_blocks(false),
      _nblock(0),
      _subvol_rad(0),
      _do_vol_corr(false),
      _processed_some_frames(false) {}

RDFCalculator::~RDFCalculator() = default;

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void RDFCalculator::Initialize() {
  // do some output
  std::cout << "begin to calculate distribution functions\n";
  std::cout << "# of bonded interactions: " << _bonded.size() << std::endl;
  std::cout << "# of non-bonded interactions: " << _nonbonded.size()
            << std::endl;

  if (_bonded.size() + _nonbonded.size() == 0) {
    throw std::runtime_error(
        "No interactions defined in options xml-file - nothing to be done");
  }

  // initialize non-bonded structures
  for (Property *prop : _nonbonded) {
    interaction_t *i = AddInteraction(prop);
    i->_is_bonded = false;
  }
}

void RDFCalculator::BeginEvaluate(Topology *top, Topology *) {
  Eigen::Matrix3d box = top->getBox();
  Eigen::Vector3d a = box.col(0);
  Eigen::Vector3d b = box.col(1);
  Eigen::Vector3d c = box.col(2);
  _boxc = box.rowwise().sum() / 2.0;

  std::cout << "Using center of box: " << _boxc << std::endl;
  // we didn't process any frames so far
  _nframes = 0;
  _nblock = 0;
  _processed_some_frames = false;

  // initialize non-bonded structures
  for (Property *prop : _nonbonded) {
    std::string name = prop->get("name").value();

    interaction_t &i = *_interactions[name];

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

    if (_do_vol_corr) {
      std::cout << "Volume correction on" << std::endl;
      i._norm = 1. / (4.0 * votca::tools::conv::Pi * i._step);
    } else {
      std::cout << "Volume correction off" << std::endl;
      i._norm = 1. / (4. * votca::tools::conv::Pi * i._step);
    }
  }
}
// create an entry for interactions
RDFCalculator::interaction_t *RDFCalculator::AddInteraction(Property *p) {
  std::string name = p->get("name").value();
  std::string group;

  group = "none";

  interaction_t *i = new interaction_t;
  i->_index = _interactions.size();
  _interactions[name] = i;
  getGroup(group)->_interactions.push_back(i);

  i->_step = p->get("step").as<double>();
  i->_min = p->get("min").as<double>();
  i->_max = p->get("max").as<double>();
  i->_norm = 1.0;
  i->_p = p;

  // initialize the current and average histogram
  int n = (int)((i->_max - i->_min) / i->_step + 1.000000001);

  i->_average.Initialize(i->_min, i->_max + i->_step, n);

  return i;
}

// end of trajectory, post processing data
void RDFCalculator::EndEvaluate() {
  if (_nframes > 0) {
    if (!_do_blocks) {
      WriteDist();
    }
  }
  // clear interactions and groups
  _interactions.clear();
  _groups.clear();
  if (!_processed_some_frames) {
    throw std::runtime_error(
        "no frames were processed. Please check your input");
  }
}

// load options from xml file
void RDFCalculator::LoadOptions(const std::string &file) {
  _options.LoadFromXML(file);
  _bonded = _options.Select("cg.bonded");
  _nonbonded = _options.Select("cg.non-bonded");
}

// evaluate current conformation
void RDFCalculator::Worker::EvalConfiguration(Topology *top, Topology *) {
  _cur_vol = 4.0 / 3.0 * votca::tools::conv::Pi * _rdfcalculator->_subvol_rad *
             _rdfcalculator->_subvol_rad * _rdfcalculator->_subvol_rad;
  // process non-bonded interactions
  DoNonbonded(top);
  // process bonded interactions
  DoBonded(top);
}

void RDFCalculator::ClearAverages() {

  _nframes = 0;
  for (auto &_interaction : _interactions) {
    _interaction.second->_average.Clear();
  }

  for (auto &_group : _groups) {
    _group.second->_corr.setZero();
  }
}

class IMCNBSearchHandler {
 public:
  IMCNBSearchHandler(HistogramNew *hist, double subvol_rad,
                     Eigen::Vector3d boxc, bool do_vol_corr)
      : _hist(hist),
        _subvol_rad(subvol_rad),
        _boxc(boxc),
        _do_vol_corr(do_vol_corr) {}

  HistogramNew *_hist;
  double _subvol_rad;
  Eigen::Vector3d _boxc;  // center of box
  bool _do_vol_corr;

  bool FoundPair(Bead *b1, Bead *, const Eigen::Vector3d &, const double dist) {

    if (_do_vol_corr) {
      double dr = (b1->Pos() - _boxc).norm();
      if (dist + dr > _subvol_rad) {
        // 2.0 is because everything is normalized to 4 PI
        _hist->Process(dist, 2.0 / SurfaceRatio(dist, dr));
      } else {
        _hist->Process(dist);
      }

    } else {
      _hist->Process(dist);
    }
    return false;
  }

  double SurfaceRatio(double dist, double r) {
    // r: distance of particle from ex center
    // dist: distance between particles
    return (1.0 + (_subvol_rad * _subvol_rad - r * r - dist * dist) /
                      (2.0 * dist * r));
  }
};

// process non-bonded interactions for current frame
void RDFCalculator::Worker::DoNonbonded(Topology *top) {
  for (Property *prop : _rdfcalculator->_nonbonded) {
    std::string name = prop->get("name").value();

    interaction_t &i = *_rdfcalculator->_interactions[name];

    // generate the bead lists
    BeadList beads1, beads2;

    beads1.GenerateInSphericalSubvolume(*top, prop->get("type1").value(),
                                        _rdfcalculator->_boxc,
                                        _rdfcalculator->_subvol_rad);
    beads2.GenerateInSphericalSubvolume(*top, prop->get("type2").value(),
                                        _rdfcalculator->_boxc,
                                        _rdfcalculator->_subvol_rad);

    _cur_beadlist_1_count = (double)beads1.size();
    _cur_beadlist_2_count = (double)beads2.size();

    // same types, so put factor 1/2 because of already counted interactions
    if (prop->get("type1").value() == prop->get("type2").value()) {
      _cur_beadlist_2_count /= 2.0;
    }

    // generate the neighbour list
    std::unique_ptr<NBList> nb;

    bool gridsearch = true;

    if (_rdfcalculator->_options.exists("cg.nbsearch")) {
      if (_rdfcalculator->_options.get("cg.nbsearch").as<std::string>() ==
          "grid") {
        gridsearch = true;
      } else if (_rdfcalculator->_options.get("cg.nbsearch")
                     .as<std::string>() == "simple") {
        gridsearch = false;
      } else {
        throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
      }
    }
    if (gridsearch) {
      nb = std::make_unique<NBList>(NBListGrid());
    } else {
      nb = std::make_unique<NBList>(NBList());
    }

    nb->setCutoff(i._max + i._step);

    // clear the current histogram
    _current_hists[i._index].Clear();

    IMCNBSearchHandler h(&(_current_hists[i._index]),
                         _rdfcalculator->_subvol_rad, _rdfcalculator->_boxc,
                         _rdfcalculator->_do_vol_corr);
    nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

    // is it same types or different types?
    if (prop->get("type1").value() == prop->get("type2").value()) {
      nb->Generate(beads1);
    } else {
      nb->Generate(beads1, beads2);
    }

    // store particle number in subvolume for each interaction
    i._avg_beadlist_1_count.Process(_cur_beadlist_1_count);
    i._avg_beadlist_2_count.Process(_cur_beadlist_2_count);
  }
}

// process non-bonded interactions for current frame
void RDFCalculator::Worker::DoBonded(Topology *top) {
  for (Property *prop : _rdfcalculator->_bonded) {
    std::string name = prop->get("name").value();

    interaction_t &i = *_rdfcalculator->_interactions[name];

    // clear the current histogram
    _current_hists[i._index].Clear();

    // now fill with new data
    std::list<Interaction *> list = top->InteractionsInGroup(name);

    for (auto ic : list) {
      double v = ic->EvaluateVar(*top);
      _current_hists[i._index].Process(v);
    }
  }
}

// returns a group, creates it if doesn't exist
RDFCalculator::group_t *RDFCalculator::getGroup(const std::string &name) {
  std::map<std::string, group_t *>::iterator iter;
  iter = _groups.find(name);
  if (iter == _groups.end()) {
    return _groups[name] = new group_t;
  }
  return (*iter).second;
}

// write the distribution function
void RDFCalculator::WriteDist(const std::string &suffix) {

  // for all interactions
  for (auto &_interaction : _interactions) {
    // calculate the rdf
    Table &t = _interaction.second->_average.data();
    Table dist(t);

    _interaction.second->_norm /=
        (_interaction.second->_avg_beadlist_1_count.getAvg() *
         _interaction.second->_avg_beadlist_2_count.getAvg());
    dist.y() = _avg_vol.getAvg() * _interaction.second->_norm *
               dist.y().cwiseQuotient(dist.x().cwiseAbs2());

    dist.Save((_interaction.first) + suffix + ".dist.new");
    std::cout << "written " << (_interaction.first) + suffix + ".dist.new\n";

    std::cout << "Avg. number of particles in subvol for "
              << (_interaction.first) << std::endl;
    std::cout << "beadlist 1: "
              << _interaction.second->_avg_beadlist_1_count.getAvg()
              << std::endl;
    std::cout << "beadlist 2: "
              << _interaction.second->_avg_beadlist_2_count.getAvg()
              << std::endl;
  }

  std::cout << "Volume used for normalization: " << _avg_vol.getAvg()
            << std::endl;
}

CsgApplication::Worker *RDFCalculator::ForkWorker() {
  RDFCalculator::Worker *worker;
  worker = new RDFCalculator::Worker;

  worker->_current_hists.resize(_interactions.size());
  worker->_rdfcalculator = this;

  for (auto &_interaction : _interactions) {
    interaction_t *i = _interaction.second;
    worker->_current_hists[i->_index].Initialize(
        i->_average.getMin(), i->_average.getMax(), i->_average.getNBins());
  }
  return worker;
}

void RDFCalculator::MergeWorker(CsgApplication::Worker *worker_) {
  _processed_some_frames = true;
  RDFCalculator::Worker *worker =
      dynamic_cast<RDFCalculator::Worker *>(worker_);
  // update the average

  ++_nframes;

  _avg_vol.Process(worker->_cur_vol);

  for (auto &_interaction : _interactions) {
    interaction_t *i = _interaction.second;
    i->_average.data().y() =
        (((double)_nframes - 1.0) * i->_average.data().y() +
         worker->_current_hists[i->_index].data().y()) /
        (double)_nframes;
  }

  if (_write_every != 0) {
    if ((_nframes % _write_every) == 0) {
      _nblock++;
      std::string suffix =
          std::string("_") + boost::lexical_cast<std::string>(_nblock);
      WriteDist(suffix);
      if (_do_blocks) {
        ClearAverages();
      }
    }
  }
}

}  // namespace csg
}  // namespace votca
