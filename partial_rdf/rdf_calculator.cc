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

#include "rdf_calculator.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <votca/csg/beadlist.h>
#include <votca/csg/imcio.h>
#include <votca/csg/nblistgrid.h>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

RDFCalculator::RDFCalculator()
    : _write_every(0),
      _do_blocks(false),
      _do_vol_corr(false),
      _processed_some_frames(false) {}

RDFCalculator::~RDFCalculator() {}

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void RDFCalculator::Initialize() {
  // do some output
  cout << "begin to calculate distribution functions\n";
  cout << "# of bonded interactions: " << _bonded.size() << endl;
  cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;

  if (_bonded.size() + _nonbonded.size() == 0)
    throw std::runtime_error(
        "No interactions defined in options xml-file - nothing to be done");

  // initialize non-bonded structures
  for (Property *prop:_nonbonded) {
    interaction_t *i = AddInteraction(prop);
    i->_is_bonded = false;
  }
};

void RDFCalculator::BeginEvaluate(Topology *top, Topology *top_atom) {
  Eigen::Matrix3d box = top->getBox();
  Eigen::Vector3d a = box.col(0);
  Eigen::Vector3d b = box.col(1);
  Eigen::Vector3d c = box.col(2);
  _boxc = a / 2 + b / 2 + c / 2;

  cout << "Using center of box: " << _boxc << endl;
  // we didn't process any frames so far
  _nframes = 0;
  _nblock = 0;
  _processed_some_frames = false;

  // initialize non-bonded structures
   for (Property *prop:_nonbonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_interactions[name];

    // count total species for ideal densities

    BeadList allbeads1, allbeads2;
    allbeads1.Generate(*top, prop->get("type1").value());
    allbeads2.Generate(*top, prop->get("type2").value());

    if (allbeads1.size() == 0)
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type1").value() +
                               "\"\n"
                               "This was specified in type1 of interaction \"" +
                               name + "\"");
    if (allbeads2.size() == 0)
      throw std::runtime_error("Topology does not have beads of type \"" +
                               prop->get("type2").value() +
                               "\"\n"
                               "This was specified in type2 of interaction \"" +
                               name + "\"");
    // calculate normalization factor for rdf

    /*if ((*iter)->get("type1").value() == (*iter)->get("type2").value())
        i._norm = 1. / (4. * M_PI * i._step * beads1.size()*(beads2.size() - 1.)
    / 2.); else i._norm = 1. / (4. * M_PI * i._step * beads1.size() *
    beads2.size());*/
    /*if ((*iter)->get("type1").value() == (*iter)->get("type2").value())
        i._norm = 1. / ( (beads2.size() - 1.) / 2.);
    else
        i._norm = 1. / (  beads2.size());*/

    if (_do_vol_corr) {
      cout << "Volume correction on" << endl;
      i._norm = 1. / (4.0 * M_PI * i._step);
    } else {
      cout << "Volume correction off" << endl;
      i._norm = 1. / (4. * M_PI * i._step);
    }
  }
}
// create an entry for interactions
RDFCalculator::interaction_t *RDFCalculator::AddInteraction(Property *p) {
  string name = p->get("name").value();
  string group;

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
  if (!_processed_some_frames)
    throw std::runtime_error(
        "no frames were processed. Please check your input");
}

// load options from xml file
void RDFCalculator::LoadOptions(const string &file) {
  load_property_from_xml(_options, file);
  _bonded = _options.Select("cg.bonded");
  _nonbonded = _options.Select("cg.non-bonded");
}

// evaluate current conformation
void RDFCalculator::Worker::EvalConfiguration(Topology *top,
                                              Topology *top_atom) {
  _cur_vol = 4.0 / 3.0 * M_PI * _rdfcalculator->_subvol_rad *
             _rdfcalculator->_subvol_rad * _rdfcalculator->_subvol_rad;
  // process non-bonded interactions
  DoNonbonded(top);
  // process bonded interactions
  DoBonded(top);
}

void RDFCalculator::ClearAverages() {
  map<string, interaction_t *>::iterator ic_iter;
  map<string, group_t *>::iterator group_iter;

  _nframes = 0;
  for (ic_iter = _interactions.begin(); ic_iter != _interactions.end();
       ++ic_iter)
    ic_iter->second->_average.Clear();

  for (group_iter = _groups.begin(); group_iter != _groups.end(); ++group_iter)
    group_iter->second->_corr.setZero();
}

class IMCNBSearchHandler {
 public:
  IMCNBSearchHandler(HistogramNew *hist, double subvol_rad, Eigen::Vector3d boxc,
                     bool do_vol_corr)
      : _hist(hist),
        _subvol_rad(subvol_rad),
        _boxc(boxc),
        _do_vol_corr(do_vol_corr) {}

  HistogramNew *_hist;
  double _subvol_rad;
  Eigen::Vector3d _boxc;  // center of box
  bool _do_vol_corr;

  bool FoundPair(Bead *b1, Bead *b2, const Eigen::Vector3d &r, const double dist) {

    if (_do_vol_corr) {
      double dr = (b1->Pos() - _boxc).norm();
      if (dist + dr > _subvol_rad)
        // 2.0 is because everything is normalized to 4 PI
        _hist->Process(dist, 2.0 / SurfaceRatio(dist, dr));
      else
        _hist->Process(dist);

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
  for (Property *prop: _rdfcalculator->_nonbonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_rdfcalculator->_interactions[name];

    // generate the bead lists
    BeadList beads1, beads2;

    beads1.GenerateInSphericalSubvolume(*top, prop->get("type1").value(),
                                        _rdfcalculator->_boxc,
                                        _rdfcalculator->_subvol_rad);
    beads2.GenerateInSphericalSubvolume(*top, prop->get("type2").value(),
                                        _rdfcalculator->_boxc,
                                        _rdfcalculator->_subvol_rad);

    _cur_beadlist_1_count = beads1.size();
    _cur_beadlist_2_count = beads2.size();

    // same types, so put factor 1/2 because of already counted interactions
    if (prop->get("type1").value() == prop->get("type2").value()) {
      _cur_beadlist_2_count /= 2.0;
    }

    // generate the neighbour list
    NBList *nb;

    bool gridsearch = true;

    if (_rdfcalculator->_options.exists("cg.nbsearch")) {
      if (_rdfcalculator->_options.get("cg.nbsearch").as<string>() == "grid")
        gridsearch = true;
      else if (_rdfcalculator->_options.get("cg.nbsearch").as<string>() ==
               "simple")
        gridsearch = false;
      else
        throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }
    if (gridsearch)
      nb = new NBListGrid();
    else
      nb = new NBList();

    nb->setCutoff(i._max + i._step);

    // clear the current histogram
    _current_hists[i._index].Clear();

    IMCNBSearchHandler h(&(_current_hists[i._index]),
                         _rdfcalculator->_subvol_rad, _rdfcalculator->_boxc,
                         _rdfcalculator->_do_vol_corr);
    nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

    // is it same types or different types?
    if (prop->get("type1").value() == prop->get("type2").value())
      nb->Generate(beads1);
    else
      nb->Generate(beads1, beads2);

    // store particle number in subvolume for each interaction
    i._avg_beadlist_1_count.Process(_cur_beadlist_1_count);
    i._avg_beadlist_2_count.Process(_cur_beadlist_2_count);

    // process all pairs
    /*NBList::iterator pair_iter;
    for(pair_iter = nb->begin(); pair_iter!=nb->end();++pair_iter) {
            _current_hists[i._index].Process((*pair_iter)->dist());
    }*/

    delete nb;
  }
}

// process non-bonded interactions for current frame
void RDFCalculator::Worker::DoBonded(Topology *top) {
  for (Property *prop:_rdfcalculator->_bonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_rdfcalculator->_interactions[name];

    // clear the current histogram
    _current_hists[i._index].Clear();

    // now fill with new data
    std::list<Interaction *> list = top->InteractionsInGroup(name);

    std::list<Interaction *>::iterator ic_iter;
    for (ic_iter = list.begin(); ic_iter != list.end(); ++ic_iter) {
      Interaction *ic = *ic_iter;
      double v = ic->EvaluateVar(*top);
      _current_hists[i._index].Process(v);
    }
  }
}

// returns a group, creates it if doesn't exist
RDFCalculator::group_t *RDFCalculator::getGroup(const string &name) {
  map<string, group_t *>::iterator iter;
  iter = _groups.find(name);
  if (iter == _groups.end()) {
    return _groups[name] = new group_t;
  }
  return (*iter).second;
}

// write the distribution function
void RDFCalculator::WriteDist(const string &suffix) {
  map<string, interaction_t *>::iterator iter;

  // for all interactions
  for (iter = _interactions.begin(); iter != _interactions.end(); ++iter) {
    // calculate the rdf
    Table &t = iter->second->_average.data();
    Table dist(t);

    iter->second->_norm /= (iter->second->_avg_beadlist_1_count.getAvg() *
                            iter->second->_avg_beadlist_2_count.getAvg());
    dist.y() = _avg_vol.getAvg() * iter->second->_norm *
               dist.y().cwiseQuotient(dist.x().cwiseAbs2());

    dist.Save((iter->first) + suffix + ".dist.new");
    cout << "written " << (iter->first) + suffix + ".dist.new\n";

    cout << "Avg. number of particles in subvol for " << (iter->first) << endl;
    cout << "beadlist 1: " << iter->second->_avg_beadlist_1_count.getAvg()
         << endl;
    cout << "beadlist 2: " << iter->second->_avg_beadlist_2_count.getAvg()
         << endl;
  }

  cout << "Volume used for normalization: " << _avg_vol.getAvg() << endl;
}

CsgApplication::Worker *RDFCalculator::ForkWorker() {
  RDFCalculator::Worker *worker;
  worker = new RDFCalculator::Worker;
  map<string, interaction_t *>::iterator ic_iter;

  worker->_current_hists.resize(_interactions.size());
  worker->_rdfcalculator = this;

  for (ic_iter = _interactions.begin(); ic_iter != _interactions.end();
       ++ic_iter) {
    interaction_t *i = ic_iter->second;
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
  map<string, interaction_t *>::iterator ic_iter;
  // map<string, group_t *>::iterator group_iter;

  ++_nframes;

  _avg_vol.Process(worker->_cur_vol);

  for (ic_iter = _interactions.begin(); ic_iter != _interactions.end();
       ++ic_iter) {
    interaction_t *i = ic_iter->second;
    i->_average.data().y() =
        (((double)_nframes - 1.0) * i->_average.data().y() +
         worker->_current_hists[i->_index].data().y()) /
        (double)_nframes;
  }

  if (_write_every != 0) {
    if ((_nframes % _write_every) == 0) {
      _nblock++;
      string suffix = string("_") + boost::lexical_cast<string>(_nblock);
      WriteDist(suffix);
      if (_do_blocks) ClearAverages();
    }
  }
}

}  // namespace csg
}  // namespace votca
