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

#ifndef VOTCA_CSG_CSG_STAT_IMC_H
#define VOTCA_CSG_CSG_STAT_IMC_H

// VOTCA includes
#include <votca/tools/average.h>
#include <votca/tools/histogramnew.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"

namespace votca {
namespace csg {
/**
 * \brief class to calculate distribution functions and cross correlations for
 * inverse monte carlo
 *
 * This class calculates distribution functions as well as cross-correlations
 * for specific groups of interactions based on a given trajectory.
 *
 */
class Imc {
 public:
  void Initialize(void);

  /// load cg definitions file
  void LoadOptions(const std::string &file);

  /// begin coarse graining a trajectory
  void BeginEvaluate(Topology *top, Topology *top_atom);

  /// end coarse graining a trajectory
  void EndEvaluate();

  void BlockLength(votca::Index length) { block_length_ = length; }
  void DoImc(bool do_imc) { do_imc_ = do_imc; }
  void OnlyIntraNB(bool only_intra_nb) { only_intra_nb_ = only_intra_nb; }
  void Extension(std::string ext) { extension_ = ext; }

 protected:
  tools::Average<double> avg_vol_;

  using group_matrix = Eigen::MatrixXd;
  using pair_matrix = Eigen::Block<group_matrix>;

  /// struct to store collected information for interactions
  struct interaction_t {
    votca::Index index_;
    tools::Property *p_;
    tools::HistogramNew average_;
    tools::HistogramNew average_force_;
    double min_, max_, step_;
    double norm_;
    double cut_;
    bool is_bonded_;
    bool threebody_;
    bool force_;
  };

  // a pair of interactions which are correlated
  struct pair_t {
    interaction_t *i1_;
    interaction_t *i2_;
    votca::Index offset_i_, offset_j_;
    pair_matrix corr_;
    pair_t(interaction_t *i1, interaction_t *i2, votca::Index offset_i,
           votca::Index offset_j, const pair_matrix &corr);
  };

  /// struct to store collected information for groups (e.g. crosscorrelations)
  struct group_t {
    std::vector<interaction_t *> interactions_;
    group_matrix corr_;
    std::vector<pair_t> pairs_;
  };

  /// the options parsed from cg definition file
  tools::Property options_;
  // length of the block to write out and averages are clear after every write
  votca::Index block_length_ = 0;
  // calculate the inverse monte carlos parameters (cross correlations)
  bool do_imc_ = false;
  // only do the intramolecular non-bonded
  bool only_intra_nb_ = false;

  // file extension for the distributions
  std::string extension_;

  // number of frames we processed
  votca::Index nframes_;
  votca::Index nblock_;

  /// list of bonded interactions
  std::vector<tools::Property *> bonded_;
  /// list of non-bonded interactions
  std::vector<tools::Property *> nonbonded_;

  /// map interaction-name to interaction
  std::map<std::string, std::unique_ptr<interaction_t> > interactions_;
  /// map group-name to group
  std::map<std::string, std::unique_ptr<group_t> > groups_;

  /// create a new interaction entry based on given options
  interaction_t *AddInteraction(tools::Property *p, bool is_bonded);

  /// get group by name, creates one if it doesn't exist
  group_t *getGroup(const std::string &name);

  /// initializes the group structs after interactions were added
  void InitializeGroups();

  void WriteDist(const std::string &suffix = "");
  void WriteIMCData(const std::string &suffix = "");
  void WriteIMCBlock(const std::string &suffix);

  void CalcDeltaS(interaction_t *interaction,
                  Eigen::VectorBlock<Eigen::VectorXd> &dS);

  void ClearAverages();

  class Worker : public CsgApplication::Worker {
   public:
    std::vector<tools::HistogramNew> current_hists_;
    std::vector<tools::HistogramNew> current_hists_force_;
    Imc *imc_;
    double cur_vol_;

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom) override;
    /// process non-bonded interactions for given frame
    void DoNonbonded(Topology *top);
    /// process bonded interactions for given frame
    void DoBonded(Topology *top);
  };
  /// update the correlations after interations were processed
  void DoCorrelations(Imc::Worker *worker);

  bool processed_some_frames_ = false;

 public:
  std::unique_ptr<CsgApplication::Worker> ForkWorker();
  void MergeWorker(CsgApplication::Worker *worker_);
};

inline Imc::pair_t::pair_t(Imc::interaction_t *i1, Imc::interaction_t *i2,
                           votca::Index offset_i, votca::Index offset_j,
                           const pair_matrix &corr)
    : i1_(i1), i2_(i2), offset_i_(offset_i), offset_j_(offset_j), corr_(corr) {}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CSG_STAT_IMC_H
