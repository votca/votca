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
#include <memory>
#include <numeric>

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/rangeparser.h>

// Local VOTCA includes
#include "votca/csg/beadlist.h"
#include "votca/csg/imcio.h"
#include "votca/csg/nblistgrid.h"
#include "votca/csg/nblistgrid_3body.h"

// Local private VOTCA includes
#include "csg_stat_imc.h"

namespace votca {
namespace csg {

using namespace std;

// begin the coarse graining process
// here the data structures are prepared to handle all the data
void Imc::Initialize() {
  // do some output
  if (do_imc_) {
    cout << "begin to calculate inverse monte carlo parameters\n";
    if (include_intra_) {
      throw runtime_error("error, can not have --do-imc and --include-intra");
    }
  } else {
    cout << "begin to calculate distribution functions\n";
  }
  cout << "# of bonded interactions: " << bonded_.size() << endl;
  cout << "# of non-bonded interactions: " << nonbonded_.size() << endl;

  if (bonded_.size() + nonbonded_.size() == 0) {
    throw std::runtime_error(
        "No interactions defined in options xml-file - nothing to be done");
  }

  // initialize non-bonded structures
  for (tools::Property *prop : nonbonded_) {
    bool bonded = false;
    AddInteraction(prop, bonded);
  }

  // initialize bonded structures
  for (tools::Property *prop : bonded_) {
    bool bonded = true;
    AddInteraction(prop, bonded);
  }

  // initialize the group structures
  if (do_imc_) {
    InitializeGroups();
  }
}

void Imc::BeginEvaluate(Topology *top, Topology *) {
  // we didn't process any frames so far
  nframes_ = 0;
  nblock_ = 0;
  processed_some_frames_ = false;

  // initialize non-bonded structures
  for (tools::Property *prop : nonbonded_) {
    string name = prop->get("name").value();

    interaction_t &i = *interactions_[name];

    // Preliminary: Quickest way to incorporate 3 body correlations
    if (i.threebody_) {

      // generate the bead lists
      BeadList beads1, beads2, beads3;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());
      beads3.Generate(*top, prop->get("type3").value());

      if (beads1.size() == 0) {
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type1").value() +
            "\"\n"
            "This was specified in type1 of interaction \"" +
            name + "\"");
      }
      if (beads2.size() == 0) {
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type2").value() +
            "\"\n"
            "This was specified in type2 of interaction \"" +
            name + "\"");
      }
      if (beads3.size() == 0) {
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type3").value() +
            "\"\n"
            "This was specified in type3 of interaction \"" +
            name + "\"");
      }
    }
    // 2body
    if (!i.threebody_) {

      double max_dist = 0.5 * top->ShortestBoxSize();
      double max = i.average_.getMax();
      if (max > max_dist) {
        throw std::runtime_error("The max of interaction \"" + name +
                                 "\" bigger is than half the box.");
      }

      // generate the bead lists
      BeadList beads1, beads2;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());

      if (beads1.size() == 0) {
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type1").value() +
            "\"\n"
            "This was specified in type1 of interaction \"" +
            name + "\"");
      }
      if (beads2.size() == 0) {
        throw std::runtime_error(
            "Topology does not have beads of type \"" +
            prop->get("type2").value() +
            "\"\n"
            "This was specified in type2 of interaction \"" +
            name + "\"");
      }

      // calculate normalization factor for rdf
      if (prop->get("type1").value() == prop->get("type2").value()) {
        i.norm_ = 2. / (double)(beads1.size() * beads2.size());
      } else {
        i.norm_ = 1. / (double)(beads1.size() * beads2.size());
      }
    }
  }

  for (tools::Property *prop : bonded_) {
    string name = prop->get("name").value();
    std::vector<Interaction *> vec = top->InteractionsInGroup(name);
    if (vec.empty()) {
      throw std::runtime_error(
          "Bonded interaction '" + name +
          "' defined in options xml-file, but not in topology - check name "
          "definition in the mapping file again");
    }
  }
}

// create an entry for interactions
Imc::interaction_t *Imc::AddInteraction(tools::Property *p, bool is_bonded) {
  string name = p->get("name").value();
  string group;
  if (do_imc_) {
    group = p->get("inverse.imc.group").value();
  } else {
    group = "none";
  }

  votca::Index index = Index(interactions_.size());
  auto success = interactions_.insert(
      std::make_pair(name, std::make_unique<interaction_t>()));
  interaction_t *i = success.first->second.get();
  i->index_ = index;
  if (group != "none") {
    getGroup(group)->interactions_.push_back(i);
  }

  i->is_bonded_ = is_bonded;
  i->step_ = p->get("step").as<double>();
  i->min_ = p->get("min").as<double>();
  i->max_ = p->get("max").as<double>();
  if (include_intra_ && (!i->is_bonded_)) {
    i->max_ = p->get("max_intra").as<double>();
  } else {
    i->max_ = p->get("max").as<double>();
  }

  i->norm_ = 1.0;
  i->p_ = p;

  // if option threebody does not exist, replace it by default of 0
  i->threebody_ = p->ifExistsReturnElseReturnDefault<bool>("threebody", 0);

  // if option force does not exist, replace it by default of 0
  i->force_ = p->ifExistsReturnElseReturnDefault<bool>("force", 0);

  // if option cut does not exist, replace it by default of 0.37 nm
  i->cut_ = p->ifExistsReturnElseReturnDefault<double>("cut", 0.37);

  // initialize the current and average histogram
  votca::Index n =
      static_cast<votca::Index>((i->max_ - i->min_) / i->step_ + 1.000000001);

  i->average_.Initialize(i->min_, i->max_, n);
  if (i->force_) {
    i->average_force_.Initialize(i->min_, i->max_, n);
  }

  return i;
}

// end of trajectory, post processing data
void Imc::EndEvaluate() {
  if (nframes_ > 0) {
    if (block_length_ == 0) {
      string suffix = string(".") + extension_;
      WriteDist(suffix);
      if (do_imc_) {
        WriteIMCData();
      }
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
void Imc::LoadOptions(const string &file) {
  options_.LoadFromXML(file);
  bonded_ = options_.Select("cg.bonded");
  nonbonded_ = options_.Select("cg.non-bonded");
}

// evaluate current conformation
void Imc::Worker::EvalConfiguration(Topology *top, Topology *) {

  cur_vol_ = top->BoxVolume();
  // process non-bonded interactions
  DoNonbonded(top);
  // process bonded interactions
  DoBonded(top);
}

void Imc::ClearAverages() {
  nframes_ = 0;

  for (auto &inter : interactions_) {
    inter.second->average_.Clear();
    if (inter.second->force_) {
      inter.second->average_force_.Clear();
    }
  }
  for (auto &group : groups_) {
    group.second->corr_.setZero();
  }
}

class IMCNBSearchHandler {
 public:
  explicit IMCNBSearchHandler(votca::tools::HistogramNew *hist)
      : hist_(*hist) {}

  votca::tools::HistogramNew &hist_;

  bool FoundPair(Bead *, Bead *, const Eigen::Vector3d &, const double dist) {
    hist_.Process(dist);
    return false;
  }
};

// process non-bonded interactions for current frame
void Imc::Worker::DoNonbonded(Topology *top) {
  for (tools::Property *prop : imc_->nonbonded_) {
    string name = prop->get("name").value();

    interaction_t &i = *imc_->interactions_[name];

    // clear the current histogram
    current_hists_[i.index_].Clear();
    current_hists_force_[i.index_].Clear();

    bool gridsearch = true;

    if (imc_->options_.exists("cg.nbsearch")) {
      if (imc_->options_.get("cg.nbsearch").as<string>() == "grid") {
        gridsearch = true;
      } else if (imc_->options_.get("cg.nbsearch").as<string>() == "simple") {
        gridsearch = false;
      } else {
        throw std::runtime_error("cg.nbsearch invalid, can be grid or simple");
      }
    }

    // Preleminary: Quickest way to incorporate 3 body correlations
    if (i.threebody_) {

      // generate the bead lists
      BeadList beads1, beads2, beads3;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());
      beads3.Generate(*top, prop->get("type3").value());

      // generate the neighbour list
      std::unique_ptr<NBList_3Body> nb;

      if (gridsearch) {
        nb = std::unique_ptr<NBList_3Body>(new NBListGrid_3Body());
      } else {
        nb = std::make_unique<NBList_3Body>();
      }

      nb->setCutoff(i.cut_);  // implement different cutoffs for different
                              // interactions!
      // Here, a is the distance between two beads of a triple, where the 3-body
      // interaction is zero

      // check if type2 and type3 are the same
      if (prop->get("type2").value() == prop->get("type3").value()) {
        // if then type2 and type1 are the same, all three types are the same
        // use the Generate function for this case
        if (prop->get("type1").value() == prop->get("type2").value()) {
          nb->Generate(beads1, true);
        }
        // else use the Generate function for type2 being equal to type3 (and type1 being different)
        if (prop->get("type1").value() != prop->get("type2").value()) {
          nb->Generate(beads1, beads2, true);
        }
      }
      // If type2 and type3 are not the same, use the Generate function for three different bead types
      // (Even if type1 and type2 or type1 and type3 are the same, the Generate function for two different beadtypes
      // is only applicable for the case that type2 is equal to type3
      if (prop->get("type2").value() != prop->get("type3").value()) {
        nb->Generate(beads1, beads2, beads3, true);
      }

      for (auto &triple : *nb) {
        Eigen::Vector3d rij = triple->r12();
        Eigen::Vector3d rik = triple->r13();
        double var = std::acos(rij.dot(rik) /
                               sqrt(rij.squaredNorm() * rik.squaredNorm()));
        current_hists_[i.index_].Process(var);
      }
    }
    // 2body interaction
    if (!i.threebody_) {

      // generate the bead lists
      BeadList beads1, beads2;

      beads1.Generate(*top, prop->get("type1").value());
      beads2.Generate(*top, prop->get("type2").value());

      {
        // generate the neighbour list
        std::unique_ptr<NBList> nb;
        if (gridsearch) {
          nb = std::unique_ptr<NBList>(new NBListGrid());
        } else {
          nb = std::unique_ptr<NBList>(new NBListGrid());
        }

        nb->setCutoff(i.max_ + i.step_);

        IMCNBSearchHandler h(&(current_hists_[i.index_]));

        nb->SetMatchFunction(&h, &IMCNBSearchHandler::FoundPair);

        // is it same types or different types?
        if (prop->get("type1").value() == prop->get("type2").value()) {
          nb->Generate(beads1, !(imc_->include_intra_));
        } else {
          nb->Generate(beads1, beads2, !(imc_->include_intra_));
        }
      }

      // if one wants to calculate the mean force
      if (i.force_) {
        std::unique_ptr<NBList> nb_force;
        if (gridsearch) {
          nb_force = std::unique_ptr<NBList>(new NBListGrid());
        } else {
          nb_force = std::unique_ptr<NBList>(new NBListGrid());
        }

        nb_force->setCutoff(i.max_ + i.step_);

        // is it same types or different types?
        if (prop->get("type1").value() == prop->get("type2").value()) {
          nb_force->Generate(beads1);
        } else {
          nb_force->Generate(beads1, beads2);
        }

        // process all pairs to calculate the projection of the
        // mean force on bead 1 on the pair distance: F1 * r12
        for (auto &pair : *nb_force) {
          Eigen::Vector3d F2 = pair->second()->getF();
          Eigen::Vector3d F1 = pair->first()->getF();
          Eigen::Vector3d r12 = pair->r();
          r12.normalize();
          double var = pair->dist();
          double scale = 0.5 * (F2 - F1).dot(r12);
          current_hists_force_[i.index_].Process(var, scale);
        }
      }
    }
  }
}

// process non-bonded interactions for current frame
void Imc::Worker::DoBonded(Topology *top) {
  for (tools::Property *prop : imc_->bonded_) {
    string name = prop->get("name").value();

    interaction_t &i = *imc_->interactions_[name];

    // clear the current histogram
    current_hists_[i.index_].Clear();

    // now fill with new data
    for (Interaction *ic : top->InteractionsInGroup(name)) {
      double v = ic->EvaluateVar(*top);
      current_hists_[i.index_].Process(v);
    }
  }
}

// returns a group, creates it if doesn't exist
Imc::group_t *Imc::getGroup(const string &name) {
  map<string, std::unique_ptr<group_t> >::iterator iter;
  iter = groups_.find(name);
  if (iter == groups_.end()) {
    auto success =
        groups_.insert(std::make_pair(name, std::make_unique<group_t>()));
    return success.first->second.get();
  }
  return (*iter).second.get();
}

// initialize the groups after interactions are added
void Imc::InitializeGroups() {
  if (!do_imc_) {
    return;
  }
  map<string, group_t *>::iterator group_iter;

  // clear all the pairs

  // iterator over all groups
  for (auto &group : groups_) {
    auto &grp = group.second;
    grp->pairs_.clear();

    auto &interactions = grp->interactions_;
    // count number of bins needed in matrix
    votca::Index n = std::accumulate(interactions.begin(), interactions.end(),
                                     0, [](votca::Index j, interaction_t *i) {
                                       return j + i->average_.getNBins();
                                     });

    // handy access to matrix
    group_matrix &M = grp->corr_;

    // initialize matrix with zeroes
    M = Eigen::MatrixXd::Zero(n, n);

    // now create references to the sub matrices and offsets
    votca::Index offset_i = 0;
    votca::Index offset_j = 0;
    // iterate over all possible cominations of pairs
    for (votca::Index i = 0; i < votca::Index(interactions.size()); i++) {
      votca::Index n1 = interactions[i]->average_.getNBins();
      offset_j = offset_i;
      for (votca::Index j = i; j < votca::Index(interactions.size()); j++) {
        votca::Index n2 = interactions[j]->average_.getNBins();
        // create matrix proxy with sub-matrix
        pair_matrix corr = M.block(offset_i, offset_j, n1, n2);
        // add the pair
        grp->pairs_.push_back(
            pair_t(interactions[i], interactions[j], offset_i, offset_j, corr));
        offset_j += n2;
      }
      offset_i += n1;
    }
  }
}

// update the correlation matrix
void Imc::DoCorrelations(Imc::Worker *worker) {
  if (!do_imc_) {
    return;
  }

  for (auto &group : groups_) {
    auto &grp = group.second;
    // update correlation for all pairs
    for (auto &pair : grp->pairs_) {
      Eigen::VectorXd &a = worker->current_hists_[pair.i1_->index_].data().y();
      Eigen::VectorXd &b = worker->current_hists_[pair.i2_->index_].data().y();
      pair_matrix &M = pair.corr_;

      M = ((((double)nframes_ - 1.0) * M) + a * b.transpose()) /
          (double)nframes_;
    }
  }
}

// write the distribution function
void Imc::WriteDist(const string &suffix) {
  map<string, interaction_t *>::iterator iter;

  cout << std::endl;  // Cosmetic, put \n before printing names of distribution
                      // files.
  // for all interactions
  for (auto &pair : interactions_) {
    // calculate the rdf
    auto &interaction = pair.second;
    votca::tools::Table &t = interaction->average_.data();

    // if no average force calculation, dummy table
    votca::tools::Table force;
    // if average force calculation, table force contains force data
    if (interaction->force_) {
      force = interaction->average_force_.data();
    }

    votca::tools::Table dist(t);
    if (!interaction->is_bonded_) {
      // Quickest way to incorporate 3 body correlations
      if (interaction->threebody_) {
        // \TODO normalize bond and angle differently....
        double norm = dist.y().cwiseAbs().sum();
        if (norm > 0) {
          dist.y() =
              interaction->norm_ * dist.y() / (norm * interaction->step_);
        }
      }

      // 2body
      if (!interaction->threebody_) {
        // force normalization
        // normalize by number of pairs found at a specific distance
        for (votca::Index i = 0; i < force.y().size(); ++i) {
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
        for (votca::Index i = 0; i < dist.y().size(); ++i) {
          double x1 = dist.x()[i] - 0.5 * interaction->step_;
          double x2 = x1 + interaction->step_;
          if (x1 < 0) {
            dist.y()[i] = 0;
          } else {
            dist.y()[i] = avg_vol_.getAvg() * interaction->norm_ * dist.y()[i] /
                          (4. / 3. * M_PI * (x2 * x2 * x2 - x1 * x1 * x1));
          }
        }
      }

    } else {
      // \TODO normalize bond and angle differently....
      double norm = dist.y().cwiseAbs().sum();
      if (norm > 0) {
        dist.y() = interaction->norm_ * dist.y() / (norm * interaction->step_);
      }
    }

    dist.Save((pair.first) + suffix);
    cout << "written " << (pair.first) + suffix << "\n";

    // preliminary
    if (interaction->force_) {
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
  if (!do_imc_) {
    return;
  }
  // map<string, interaction_t *>::iterator ic_iter;
  map<string, group_t *>::iterator group_iter;

  // iterate over all groups
  for (auto &group : groups_) {
    auto &grp = group.second;
    string grp_name = group.first;

    // number of total bins for all interactions in group is matrix dimension
    votca::Index n = grp->corr_.rows();

    // build full set of equations + copy some data to make
    // code better to read
    group_matrix gmc(grp->corr_);
    tools::Table dS;
    dS.resize(n);
    // the next two variables are to later extract the individual parts
    // from the whole data after solving equations
    vector<std::pair<std::string, votca::tools::RangeParser> >
        ranges;  // names of the interactions & sizes of the individual
                 // interactions

    // copy all averages+r of group to one vector
    n = 0;
    votca::Index begin = 1;
    for (interaction_t *ic : grp->interactions_) {

      // sub vector for dS
      Eigen::VectorBlock<Eigen::VectorXd> sub_dS =
          dS.y().segment(n, ic->average_.getNBins());

      // sub vector for r
      Eigen::VectorBlock<Eigen::VectorXd> sub_r =
          dS.x().segment(n, ic->average_.getNBins());

      // read in target and calculate dS
      CalcDeltaS(ic, sub_dS);

      // copy r
      sub_r = ic->average_.data().x();

      // save size
      votca::tools::RangeParser rp;
      votca::Index end = begin + ic->average_.getNBins() - 1;
      rp.Add(begin, end);
      ranges.push_back(std::pair<std::string, votca::tools::RangeParser>(
          ic->p_->get("name").as<string>(), rp));
      begin = end + 1;
      // save name

      // shift subrange by size of current
      n += ic->average_.getNBins();
    }

    // now we need to calculate the
    // A_ij = <S_i*S_j> - <S_i>*<S_j>
    for (pair_t &pair : grp->pairs_) {
      interaction_t *i1 = pair.i1_;
      interaction_t *i2 = pair.i2_;

      // make reference to <S_i>
      Eigen::VectorXd &a = i1->average_.data().y();
      // make reference to <S_j>
      Eigen::VectorXd &b = i2->average_.data().y();

      votca::Index i = pair.offset_i_;
      votca::Index j = pair.offset_j_;
      votca::Index n1 = i1->average_.getNBins();
      votca::Index n2 = i2->average_.getNBins();

      pair_matrix M = gmc.block(i, j, n1, n2);
      M = -(M - a * b.transpose());
      // matrix is symmetric
      gmc.block(j, i, n2, n1) = M.transpose().eval();
    }

    imcio_write_dS(grp_name + suffix + ".imc", dS);
    imcio_write_matrix(grp_name + suffix + ".gmc", gmc);
    imcio_write_index(grp_name + suffix + ".idx", ranges);
  }
}

// calculate deviation from target vectors
void Imc::CalcDeltaS(interaction_t *interaction,
                     Eigen::VectorBlock<Eigen::VectorXd> &dS) {
  const string &name = interaction->p_->get("name").as<string>();

  tools::Table target;
  target.Load(name + ".dist.tgt");

  if (!interaction->is_bonded_) {
    for (votca::Index i = 0; i < target.y().size(); ++i) {
      double x1 = target.x()[i] - 0.5 * interaction->step_;
      double x2 = x1 + interaction->step_;
      if (x1 < 0) {
        x1 = x2 = 0;
      }
      target.y()[i] = 1. / (avg_vol_.getAvg() * interaction->norm_) *
                      target.y()[i] *
                      (4. / 3. * M_PI * (x2 * x2 * x2 - x1 * x1 * x1));
    }
  } else {
    target.y() = (1.0 / interaction->norm_) * target.y();
  }
  if (target.y().size() != interaction->average_.data().y().size()) {
    throw std::runtime_error(
        "number of grid points in target does not match the grid");
  }

  dS = interaction->average_.data().y() - target.y();
}

void Imc::WriteIMCBlock(const string &suffix) {

  if (!do_imc_) {
    return;
  }

  // iterate over all groups
  for (auto &group : groups_) {
    auto &grp = group.second;
    string grp_name = group.first;
    list<interaction_t *>::iterator iter;

    // number of total bins for all interactions in group is matrix dimension
    votca::Index n = grp->corr_.rows();

    // build full set of equations + copy some data to make code better to read
    group_matrix gmc(grp->corr_);
    Eigen::VectorXd dS(n);
    Eigen::VectorXd r(n);
    // the next two variables are to later extract the individual parts
    // from the whole data after solving equations
    vector<votca::Index> sizes;  // sizes of the individual interactions
    vector<string> names;        // names of the interactions

    // copy all averages+r of group to one vector
    n = 0;
    for (interaction_t *ic : grp->interactions_) {
      // sub vector for dS
      Eigen::VectorBlock<Eigen::VectorXd> sub_dS =
          dS.segment(n, ic->average_.getNBins());

      // sub vector for r
      Eigen::VectorBlock<Eigen::VectorXd> sub_r =
          r.segment(n, ic->average_.getNBins());

      // read in target and calculate dS
      sub_dS = ic->average_.data().y();
      // copy r
      sub_r = ic->average_.data().x();
      // save size
      sizes.push_back(ic->average_.getNBins());
      // save name
      names.push_back(ic->p_->get("name").as<string>());

      // shift subrange by size of current
      n += ic->average_.getNBins();
    }

    // write the dS
    ofstream out_dS;
    string name_dS = grp_name + suffix + ".S";
    out_dS.open(name_dS);
    out_dS << setprecision(8);
    if (!out_dS) {
      throw runtime_error(string("error, cannot open file ") + name_dS);
    }

    for (votca::Index i = 0; i < dS.size(); ++i) {
      out_dS << r[i] << " " << dS[i] << endl;
    }

    out_dS.close();
    cout << "written " << name_dS << endl;

    // write the correlations
    ofstream out_cor;
    string name_cor = grp_name + suffix + ".cor";
    out_cor.open(name_cor);
    out_cor << setprecision(8);

    if (!out_cor) {
      throw runtime_error(string("error, cannot open file ") + name_cor);
    }

    for (votca::Index i = 0; i < grp->corr_.rows(); ++i) {
      for (votca::Index j = 0; j < grp->corr_.cols(); ++j) {
        out_cor << grp->corr_(i, j) << " ";
      }
      out_cor << endl;
    }
    out_cor.close();
    cout << "written " << name_cor << endl;
  }
}

std::unique_ptr<CsgApplication::Worker> Imc::ForkWorker() {

  auto worker = std::make_unique<Imc::Worker>();
  worker->current_hists_.resize(interactions_.size());
  worker->current_hists_force_.resize(interactions_.size());
  worker->imc_ = this;

  for (auto &interaction : interactions_) {
    auto &i = interaction.second;
    worker->current_hists_[i->index_].Initialize(
        i->average_.getMin(), i->average_.getMax(), i->average_.getNBins());
    // preliminary
    if (interaction.second->force_) {
      worker->current_hists_force_[i->index_].Initialize(
          i->average_force_.getMin(), i->average_force_.getMax(),
          i->average_force_.getNBins());
    }
  }
  return worker;
}

void Imc::MergeWorker(CsgApplication::Worker *worker_) {
  processed_some_frames_ = true;
  Imc::Worker *worker = dynamic_cast<Imc::Worker *>(worker_);
  // update the average
  map<string, interaction_t *>::iterator ic_iter;
  // map<string, group_t *>::iterator group_iter;

  ++nframes_;
  avg_vol_.Process(worker->cur_vol_);
  for (auto &interaction : interactions_) {
    auto &i = interaction.second;
    i->average_.data().y() =
        (((double)nframes_ - 1.0) * i->average_.data().y() +
         worker->current_hists_[i->index_].data().y()) /
        (double)nframes_;
    // preliminary
    if (i->force_) {
      i->average_force_.data().y() =
          (((double)nframes_ - 1.0) * i->average_force_.data().y() +
           worker->current_hists_force_[i->index_].data().y()) /
          (double)nframes_;
    }
  }

  // update correlation matrices
  if (do_imc_) {
    DoCorrelations(worker);
  }

  if (block_length_ != 0) {
    if ((nframes_ % block_length_) == 0) {
      nblock_++;
      string suffix = string("_") + boost::lexical_cast<string>(nblock_) +
                      string(".") + extension_;
      WriteDist(suffix);
      WriteIMCData(suffix);
      WriteIMCBlock(suffix);
      ClearAverages();
    }
  }
}

}  // namespace csg
}  // namespace votca
