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

#include "csg_stat_imc.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <votca/csg/beadlist.h>
#include <votca/csg/imcio.h>
#include <votca/csg/nblistgrid.h>
#include <votca/csg/nblistgrid_3body.h>
#include <votca/tools/rangeparser.h>

namespace votca {
namespace csg {

using namespace std;

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void Imc::Initialize() {
  // do some output
  if (_do_imc)
    cout << "begin to calculate inverse monte carlo parameters\n";
  else
    cout << "begin to calculate distribution functions\n";
  cout << "# of bonded interactions: " << _bonded.size() << endl;
  cout << "# of non-bonded interactions: " << _nonbonded.size() << endl;

  if (_bonded.size() + _nonbonded.size() == 0)
    throw std::runtime_error(
        "No interactions defined in options xml-file - nothing to be done");

  // initialize non-bonded structures
  for (tools::Property *prop : _nonbonded) {
    interaction_t *i = AddInteraction(prop);
    i->_is_bonded = false;
  }

  // initialize bonded structures
  for (tools::Property *prop : _bonded) {
    interaction_t *i = AddInteraction(prop);
    i->_is_bonded = true;
  }

  // initialize the group structures
  if (_do_imc) InitializeGroups();
};

void Imc::BeginEvaluate(Topology *top, Topology *top_atom) {
  // we didn't process any frames so far
  _nframes = 0;
  _nblock = 0;
  _processed_some_frames = false;

  // initialize non-bonded structures
  for (tools::Property *prop : _nonbonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_interactions[name];

    // Preliminary: Quickest way to incorporate 3 body correlations
    if (i._threebody) {

      // generate the bead lists
      BeadList beads1, beads2, beads3;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());
      beads3.Generate(*top, prop->get("type3").value());

      if (beads1.size() == 0)
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type1").value() +
            "\"\n"
            "This was specified in type1 of interaction \"" +
            name + "\"");
      if (beads2.size() == 0)
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type2").value() +
            "\"\n"
            "This was specified in type2 of interaction \"" +
            name + "\"");
      if (beads3.size() == 0)
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type3").value() +
            "\"\n"
            "This was specified in type3 of interaction \"" +
            name + "\"");
    }
    // 2body
    if (!i._threebody) {

      // generate the bead lists
      BeadList beads1, beads2;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());

      if (beads1.size() == 0)
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type1").value() +
            "\"\n"
            "This was specified in type1 of interaction \"" +
            name + "\"");
      if (beads2.size() == 0)
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type2").value() +
            "\"\n"
            "This was specified in type2 of interaction \"" +
            name + "\"");

      // calculate normalization factor for rdf
      if (prop->get("type1").value() == prop->get("type2").value())
        i._norm = 1. / (beads1.size() * (beads2.size()) / 2.);
      else
        i._norm = 1. / (beads1.size() * beads2.size());
    }
  }

  for (tools::Property *prop : _bonded) {
    string name = prop->get("name").value();

    std::list<Interaction *> list = top->InteractionsInGroup(name);
    if (list.empty())
      throw std::runtime_error(
          "Bonded interaction '" + name +
          "' defined in options xml-file, but not in topology - check name "
          "definition in the mapping file again");
  }
}

// create an entry for interactions
Imc::interaction_t *Imc::AddInteraction(tools::Property *p) {
  string name = p->get("name").value();
  string group;
  if (_do_imc)
    group = p->get("inverse.imc.group").value();
  else
    group = "none";

  int index = _interactions.size();
  auto success = _interactions.insert(std::make_pair(
      name, std::unique_ptr<interaction_t>(new interaction_t())));
  interaction_t *i = success.first->second.get();
  i->_index = index;
  getGroup(group)->_interactions.push_back(i);

  i->_step = p->get("step").as<double>();
  i->_min = p->get("min").as<double>();
  i->_max = p->get("max").as<double>();
  i->_norm = 1.0;
  i->_p = p;

  // if option threebody does not exist, replace it by default of 0
  i->_threebody = p->ifExistsReturnElseReturnDefault<bool>("threebody", 0);

  // if option force does not exist, replace it by default of 0
  i->_force = p->ifExistsReturnElseReturnDefault<bool>("force", 0);

  // if option cut does not exist, replace it by default of 0.37 nm
  i->_cut = p->ifExistsReturnElseReturnDefault<double>("cut", 0.37);

  // initialize the current and average histogram
  int n = static_cast<int>((i->_max - i->_min) / i->_step + 1.000000001);

  i->_average.Initialize(i->_min, i->_max, n);
  if (i->_force) {
    i->_average_force.Initialize(i->_min, i->_max, n);
  }

  return i;
}

// end of trajectory, post processing data
void Imc::EndEvaluate() {
  if (_nframes > 0) {
    if (_block_length == 0) {
      string suffix = string(".") + _extension;
      WriteDist(suffix);
      if (_do_imc) WriteIMCData();
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
void Imc::LoadOptions(const string &file) {
  _options.LoadFromXML(file);
  _bonded = _options.Select("cg.bonded");
  _nonbonded = _options.Select("cg.non-bonded");
}

// evaluate current conformation
void Imc::Worker::EvalConfiguration(Topology *top, Topology *top_atom) {

  _cur_vol = top->BoxVolume();
  // process non-bonded interactions
  DoNonbonded(top);
  // process bonded interactions
  DoBonded(top);
}

void Imc::ClearAverages() {
  _nframes = 0;

  for (auto &inter : _interactions) {
    inter.second->_average.Clear();
    if (inter.second->_force) {
      inter.second->_average_force.Clear();
    }
  }
  for (auto &group : _groups) {
    group.second->_corr.setZero();
  }
}

class IMCNBSearchHandler {
 public:
  explicit IMCNBSearchHandler(votca::tools::HistogramNew *hist)
      : _hist(*hist) {}

  votca::tools::HistogramNew &_hist;

  bool FoundPair(Bead *b1, Bead *b2, const Eigen::Vector3d &r,
                 const double dist) {
    _hist.Process(dist);
    return false;
  }
};

// process non-bonded interactions for current frame
void Imc::Worker::DoNonbonded(Topology *top) {
  for (tools::Property *prop : _imc->_nonbonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_imc->_interactions[name];

    // clear the current histogram
    _current_hists[i._index].Clear();
    _current_hists_force[i._index].Clear();

    bool gridsearch = true;

    if (_imc->_options.exists("cg.nbsearch")) {
      if (_imc->_options.get("cg.nbsearch").as<string>() == "grid")
        gridsearch = true;
      else if (_imc->_options.get("cg.nbsearch").as<string>() == "simple")
        gridsearch = false;
      else
        throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
    }

    // Preleminary: Quickest way to incorporate 3 body correlations
    if (i._threebody) {

      // generate the bead lists
      BeadList beads1, beads2, beads3;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());
      beads3.Generate(*top, prop->get("type3").value());

      // generate the neighbour list
      std::unique_ptr<NBList_3Body> nb;

      if (gridsearch)
        nb = std::unique_ptr<NBList_3Body>(new NBListGrid_3Body());
      else
        nb = std::unique_ptr<NBList_3Body>(new NBList_3Body());

      nb->setCutoff(i._cut);  // implement different cutoffs for different
                              // interactions!
      // Here, a is the distance between two beads of a triple, where the 3-body
      // interaction is zero

      // check if type1 and type2 are the same
      if (prop->get("type1").value() == prop->get("type2").value()) {
        // if all three types are the same
        if (prop->get("type2").value() == prop->get("type3").value()) {
          nb->Generate(beads1, true);
        }
        // if type2 and type3 are different, use the Generate function for 2
        // bead types
        if (prop->get("type2").value() != prop->get("type3").value()) {
          nb->Generate(beads1, beads3, true);
        }
      }
      // if type1 != type2
      if (prop->get("type1").value() != prop->get("type2").value()) {
        // if the last two types are the same, use Generate function with them
        // as the first two bead types Neighborlist_3body is constructed in a
        // way that the two equal bead types have two be the first 2 types
        if (prop->get("type2").value() == prop->get("type3").value()) {
          nb->Generate(beads1, beads2, true);
        }
        if (prop->get("type2").value() != prop->get("type3").value()) {
          // type1 = type3 !=type2
          if (prop->get("type1").value() == prop->get("type3").value()) {
            nb->Generate(beads2, beads1, true);
          }
          // type1 != type2 != type3
          if (prop->get("type1").value() != prop->get("type3").value()) {
            nb->Generate(beads1, beads2, beads3, true);
          }
        }
      }

      for (auto &triple : *nb) {
        Eigen::Vector3d rij = triple->r12();
        Eigen::Vector3d rik = triple->r13();
        double var = std::acos(rij.dot(rik) /
                               sqrt(rij.squaredNorm() * rik.squaredNorm()));
        _current_hists[i._index].Process(var);
      }
    }
    // 2body interaction
    if (!i._threebody) {

      // generate the bead lists
      BeadList beads1, beads2;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());

      {
        // generate the neighbour list
        std::unique_ptr<NBList> nb;
        if (gridsearch)
          nb = std::unique_ptr<NBList>(new NBListGrid());
        else
          nb = std::unique_ptr<NBList>(new NBListGrid());

        nb->setCutoff(i._max + i._step);

        IMCNBSearchHandler h(&(_current_hists[i._index]));

        nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

        // is it same types or different types?
        if (prop->get("type1").value() == prop->get("type2").value())
          nb->Generate(beads1);
        else
          nb->Generate(beads1, beads2);
      }

      // if one wants to calculate the mean force
      if (i._force) {
        std::unique_ptr<NBList> nb_force;
        if (gridsearch)
          nb_force = std::unique_ptr<NBList>(new NBListGrid());
        else
          nb_force = std::unique_ptr<NBList>(new NBListGrid());

        nb_force->setCutoff(i._max + i._step);

        // is it same types or different types?
        if (prop->get("type1").value() == prop->get("type2").value())
          nb_force->Generate(beads1);
        else
          nb_force->Generate(beads1, beads2);

        // process all pairs to calculate the projection of the
        // mean force on bead 1 on the pair distance: F1 * r12
        for (auto &pair : *nb_force) {
          Eigen::Vector3d F2 = pair->second()->getF();
          Eigen::Vector3d F1 = pair->first()->getF();
          Eigen::Vector3d r12 = pair->r();
          r12.normalize();
          double var = pair->dist();
          double scale = 0.5 * (F2 - F1).dot(r12);
          _current_hists_force[i._index].Process(var, scale);
        }
      }
    }
  }
}

// process non-bonded interactions for current frame
void Imc::Worker::DoBonded(Topology *top) {
  for (tools::Property *prop : _imc->_bonded) {
    string name = prop->get("name").value();

    interaction_t &i = *_imc->_interactions[name];

    // clear the current histogram
    _current_hists[i._index].Clear();

    // now fill with new data
    for (Interaction *ic : top->InteractionsInGroup(name)) {
      double v = ic->EvaluateVar(*top);
      _current_hists[i._index].Process(v);
    }
  }
}

// returns a group, creates it if doesn't exist
Imc::group_t *Imc::getGroup(const string &name) {
  map<string, std::unique_ptr<group_t> >::iterator iter;
  iter = _groups.find(name);
  if (iter == _groups.end()) {
    auto success = _groups.insert(
        std::make_pair(name, std::unique_ptr<group_t>(new group_t())));
    return success.first->second.get();
  }
  return (*iter).second.get();
}

// initialize the groups after interactions are added
void Imc::InitializeGroups() {
  if (!_do_imc) return;
  map<string, group_t *>::iterator group_iter;

  // clear all the pairs

  // iterator over all groups
  for (auto &group : _groups) {
    auto &grp = group.second;
    grp->_pairs.clear();

    auto &interactions = grp->_interactions;
    // count number of bins needed in matrix
    int n = std::accumulate(
        interactions.begin(), interactions.end(), 0,
        [](int j, interaction_t *i) { return j + i->_average.getNBins(); });

    // handy access to matrix
    group_matrix &M = grp->_corr;

    // initialize matrix with zeroes
    M = Eigen::MatrixXd::Zero(n, n);

    // now create references to the sub matrices
    // iterate over all possible cominations of pairs
    for (int i = 0; i < int(interactions.size()); i++) {
      int n1 = interactions[i]->_average.getNBins();
      for (int j = i; j < int(interactions.size()); j++) {
        int n2 = interactions[j]->_average.getNBins();
        // create matrix proxy with sub-matrix
        pair_matrix corr = M.block(i, j, n1, n2);
        // add the pair
        grp->_pairs.push_back(
            pair_t(interactions[i], interactions[j], i, j, corr));
      }
    }
  }
}

// update the correlation matrix
void Imc::DoCorrelations(Imc::Worker *worker) {
  if (!_do_imc) return;

  for (auto &group : _groups) {
    auto &grp = group.second;
    // update correlation for all pairs
    for (auto &pair : grp->_pairs) {
      Eigen::VectorXd &a = worker->_current_hists[pair._i1->_index].data().y();
      Eigen::VectorXd &b = worker->_current_hists[pair._i2->_index].data().y();
      pair_matrix &M = pair._corr;

      M = ((((double)_nframes - 1.0) * M) + a * b.transpose()) /
          (double)_nframes;
    }
  }
}

// write the distribution function
void Imc::WriteDist(const string &suffix) {
  map<string, interaction_t *>::iterator iter;

  cout << std::endl;  // Cosmetic, put \n before printing names of distribution
                      // files.
  // for all interactions
  for (auto &pair : _interactions) {
    // calculate the rdf
    auto &interaction = pair.second;
    votca::tools::Table &t = interaction->_average.data();

    // if no average force calculation, dummy table
    votca::tools::Table force;
    // if average force calculation, table force contains force data
    if (interaction->_force) {
      force = interaction->_average_force.data();
    }

    votca::tools::Table dist(t);
    if (!interaction->_is_bonded) {
      // Quickest way to incorporate 3 body correlations
      if (interaction->_threebody) {
        // \TODO normalize bond and angle differently....
        double norm = dist.y().cwiseAbs().sum();
        if (norm > 0) {
          dist.y() =
              interaction->_norm * dist.y() / (norm * interaction->_step);
        }
      }

      // 2body
      if (!interaction->_threebody) {
        // force normalization
        // normalize by number of pairs found at a specific distance
        for (unsigned int i = 0; i < force.y().size(); ++i) {
          // check if any number of pairs has been found at this distance, then
          // normalize
          if (dist.y()[i] != 0) {
            force.y()[i] /= dist.y()[i];
          }
          // else set to zero
          else {
            force.y()[i] = 0;
          }
        }

        // normalization is calculated using exact shell volume (difference of
        // spheres)
        for (unsigned int i = 0; i < dist.y().size(); ++i) {
          double x1 = dist.x()[i] - 0.5 * interaction->_step;
          double x2 = x1 + interaction->_step;
          if (x1 < 0) {
            dist.y()[i] = 0;
          } else {
            dist.y()[i] = _avg_vol.getAvg() * interaction->_norm * dist.y()[i] /
                          (4. / 3. * M_PI * (x2 * x2 * x2 - x1 * x1 * x1));
          }
        }
      }

    } else {
      // \TODO normalize bond and angle differently....
      double norm = dist.y().cwiseAbs().sum();
      if (norm > 0) {
        dist.y() = interaction->_norm * dist.y() / (norm * interaction->_step);
      }
    }

    dist.Save((pair.first) + suffix);
    cout << "written " << (pair.first) + suffix << "\n";

    // preliminary
    if (interaction->_force) {
      force.Save((pair.first) + ".force.new");
      cout << "written " << (pair.first) + ".force.new"
           << "\n";
    }
  }
}

/**
 *  Here the inverse monte carlo matrix is calculated and written out
 *
 *  steps:
 *      - calculate th
 */
void Imc::WriteIMCData(const string &suffix) {
  if (!_do_imc) return;
  // map<string, interaction_t *>::iterator ic_iter;
  map<string, group_t *>::iterator group_iter;

  // iterate over all groups
  for (auto &group : _groups) {
    auto &grp = group.second;
    string grp_name = group.first;

    // number of total bins for all interactions in group is matrix dimension
    int n = grp->_corr.rows();

    // build full set of equations + copy some data to make
    // code better to read
    group_matrix gmc(grp->_corr);
    Eigen::VectorXd dS(n);
    Eigen::VectorXd r(n);
    // the next two variables are to later extract the individual parts
    // from the whole data after solving equations
    vector<votca::tools::RangeParser> ranges;  // sizes of the individual
                                               // interactions
    vector<string> names;                      // names of the interactions

    // copy all averages+r of group to one vector
    n = 0;
    int begin = 1;
    for (interaction_t *ic : grp->_interactions) {

      // sub vector for dS
      Eigen::VectorBlock<Eigen::VectorXd> sub_dS =
          dS.segment(n, ic->_average.getNBins());

      // sub vector for r
      Eigen::VectorBlock<Eigen::VectorXd> sub_r =
          r.segment(n, ic->_average.getNBins());

      // read in target and calculate dS
      CalcDeltaS(ic, sub_dS);

      // copy r
      sub_r = ic->_average.data().x();

      // save size
      votca::tools::RangeParser rp;
      int end = begin + ic->_average.getNBins() - 1;
      rp.Add(begin, end);
      ranges.push_back(rp);
      begin = end + 1;
      // save name
      names.push_back(ic->_p->get("name").as<string>());

      // shift subrange by size of current
      n += ic->_average.getNBins();
    }

    // now we need to calculate the
    // A_ij = <S_i*S_j> - <S_i>*<S_j>
    for (pair_t &pair : grp->_pairs) {
      interaction_t *i1 = pair._i1;
      interaction_t *i2 = pair._i2;

      // make reference to <S_i>
      Eigen::VectorXd &a = i1->_average.data().y();
      // make reference to <S_j>
      Eigen::VectorXd &b = i2->_average.data().y();

      int i = pair._offset_i;
      int j = pair._offset_j;
      int n1 = i1->_average.getNBins();
      int n2 = i2->_average.getNBins();

      pair_matrix M = gmc.block(i, j, n1, n2);
      M = -(M - a * b.transpose());
      // matrix is symmetric
      gmc.block(j, i, n2, n1) = M.transpose().eval();
    }

    imcio_write_dS(grp_name + suffix + ".imc", r, dS);
    imcio_write_matrix(grp_name + suffix + ".gmc", gmc);
    imcio_write_index(grp_name + suffix + ".idx", names, ranges);
  }
}

// calculate deviation from target vectors
void Imc::CalcDeltaS(interaction_t *interaction,
                     Eigen::VectorBlock<Eigen::VectorXd> &dS) {
  const string &name = interaction->_p->get("name").as<string>();

  tools::Table target;
  target.Load(name + ".dist.tgt");

  if (!interaction->_is_bonded) {
    for (unsigned int i = 0; i < target.y().size(); ++i) {
      double x1 = target.x()[i] - 0.5 * interaction->_step;
      double x2 = x1 + interaction->_step;
      if (x1 < 0) x1 = x2 = 0;
      target.y()[i] = 1. / (_avg_vol.getAvg() * interaction->_norm) *
                      target.y()[i] *
                      (4. / 3. * M_PI * (x2 * x2 * x2 - x1 * x1 * x1));
    }
  } else {
    target.y() = (1.0 / interaction->_norm) * target.y();
  }
  if (target.y().size() != interaction->_average.data().y().size())
    throw std::runtime_error(
        "number of grid points in target does not match the grid");

  dS = interaction->_average.data().y() - target.y();
}

void Imc::WriteIMCBlock(const string &suffix) {

  if (!_do_imc) return;

  // iterate over all groups
  for (auto &group : _groups) {
    auto &grp = group.second;
    string grp_name = group.first;
    list<interaction_t *>::iterator iter;

    // number of total bins for all interactions in group is matrix dimension
    int n = grp->_corr.rows();

    // build full set of equations + copy some data to make code better to read
    group_matrix gmc(grp->_corr);
    Eigen::VectorXd dS(n);
    Eigen::VectorXd r(n);
    // the next two variables are to later extract the individual parts
    // from the whole data after solving equations
    vector<int> sizes;     // sizes of the individual interactions
    vector<string> names;  // names of the interactions

    // copy all averages+r of group to one vector
    n = 0;
    for (interaction_t *ic : grp->_interactions) {
      // sub vector for dS
      Eigen::VectorBlock<Eigen::VectorXd> sub_dS =
          dS.segment(n, ic->_average.getNBins());

      // sub vector for r
      Eigen::VectorBlock<Eigen::VectorXd> sub_r =
          r.segment(n, ic->_average.getNBins());

      // read in target and calculate dS
      sub_dS = ic->_average.data().y();
      // copy r
      sub_r = ic->_average.data().x();
      // save size
      sizes.push_back(ic->_average.getNBins());
      // save name
      names.push_back(ic->_p->get("name").as<string>());

      // shift subrange by size of current
      n += ic->_average.getNBins();
    }

    // write the dS
    ofstream out_dS;
    string name_dS = grp_name + suffix + ".S";
    out_dS.open(name_dS);
    out_dS << setprecision(8);
    if (!out_dS)
      throw runtime_error(string("error, cannot open file ") + name_dS);

    for (int i = 0; i < dS.size(); ++i) {
      out_dS << r[i] << " " << dS[i] << endl;
    }

    out_dS.close();
    cout << "written " << name_dS << endl;

    // write the correlations
    ofstream out_cor;
    string name_cor = grp_name + suffix + ".cor";
    out_cor.open(name_cor);
    out_cor << setprecision(8);

    if (!out_cor)
      throw runtime_error(string("error, cannot open file ") + name_cor);

    for (int i = 0; i < grp->_corr.rows(); ++i) {
      for (int j = 0; j < grp->_corr.cols(); ++j) {
        out_cor << grp->_corr(i, j) << " ";
      }
      out_cor << endl;
    }
    out_cor.close();
    cout << "written " << name_cor << endl;
  }
}

CsgApplication::Worker *Imc::ForkWorker() {

  Imc::Worker *worker = new Imc::Worker;
  worker->_current_hists.resize(_interactions.size());
  worker->_current_hists_force.resize(_interactions.size());
  worker->_imc = this;

  for (auto &interaction : _interactions) {
    auto &i = interaction.second;
    worker->_current_hists[i->_index].Initialize(
        i->_average.getMin(), i->_average.getMax(), i->_average.getNBins());
    // preliminary
    if (interaction.second->_force) {
      worker->_current_hists_force[i->_index].Initialize(
          i->_average_force.getMin(), i->_average_force.getMax(),
          i->_average_force.getNBins());
    }
  }
  return worker;
}

void Imc::MergeWorker(CsgApplication::Worker *worker_) {
  _processed_some_frames = true;
  Imc::Worker *worker = dynamic_cast<Imc::Worker *>(worker_);
  // update the average
  map<string, interaction_t *>::iterator ic_iter;
  // map<string, group_t *>::iterator group_iter;

  ++_nframes;
  _avg_vol.Process(worker->_cur_vol);
  for (auto &interaction : _interactions) {
    auto &i = interaction.second;
    i->_average.data().y() =
        (((double)_nframes - 1.0) * i->_average.data().y() +
         worker->_current_hists[i->_index].data().y()) /
        (double)_nframes;
    // preliminary
    if (i->_force) {
      i->_average_force.data().y() =
          (((double)_nframes - 1.0) * i->_average_force.data().y() +
           worker->_current_hists_force[i->_index].data().y()) /
          (double)_nframes;
    }
  }

  // update correlation matrices
  if (_do_imc) DoCorrelations(worker);

  if (_block_length != 0) {
    if ((_nframes % _block_length) == 0) {
      _nblock++;
      string suffix = string("_") + boost::lexical_cast<string>(_nblock) +
                      string(".") + _extension;
      WriteDist(suffix);
      WriteIMCData(suffix);
      WriteIMCBlock(suffix);
      ClearAverages();
    }
  }
}

}  // namespace csg
}  // namespace votca
