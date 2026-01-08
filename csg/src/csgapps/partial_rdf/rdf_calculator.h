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

#ifndef VOTCA_CSG_RDF_CALCULATOR_H
#define VOTCA_CSG_RDF_CALCULATOR_H

// Standard includes
#include <cmath>
#include <memory>

// Third party includes
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

// VOTCA includes
#include <votca/tools/average.h>
#include <votca/tools/histogram.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include <votca/csg/csgapplication.h>

namespace votca {
namespace csg {
using namespace votca::tools;

/**
 * \brief class to calculate distribution functions and cross correlations for
 * inverse monte carlo
 *
 * This class calculates distribution functions as well as cross-correlations
 * for specific groups of interactions based on a given trajectory.
 *
 */
class RDFCalculator {
 public:
  RDFCalculator();
  ~RDFCalculator();

  void Initialize(void);

  /// load cg definitions file
  void LoadOptions(const std::string &file);

  /// begin coarse graining a trajectory
  void BeginEvaluate(Topology *top, Topology *top_atom);

  /// end coarse graining a trajectory
  void EndEvaluate();

  void WriteEvery(Index write_every) { write_every_ = write_every; }
  void DoBlocks(bool do_blocks) { do_blocks_ = do_blocks; }
  void DoVolumeCorrection(bool do_vol_corr) { do_vol_corr_ = do_vol_corr; }
  void SetSubvolRadius(double r) { subvol_rad_ = r; }
  double AnalyticVolumeCorrection(double t) {

    std::cout << "DBG " << t << " "
              << 1.0 / 24.0 *
                     (16.0 * t * t - 12.0 * t * t * t + t * t * t * t * t)
              << std::endl;
    return 1.0 / 24.0 * (16.0 * t * t - 12.0 * t * t * t + t * t * t * t * t);
  }

 protected:
  Average<double> avg_vol_;

  using group_matrix = Eigen::MatrixXd;
  using pair_matrix = Eigen::Block<group_matrix>;

  /// struct to store collected information for interactions
  struct interaction_t {
    Index index_;
    Property *p_;
    Histogram average_;
    double min_, max_, step_;
    double norm_;
    bool is_bonded_;
    Average<double> avg_beadlist_1_count_;
    Average<double> avg_beadlist_2_count_;
  };

  // a pair of interactions which are correlated
  struct pair_t {
    interaction_t *i1_;
    interaction_t *i2_;
    Index offset_i_, offset_j_;
    pair_matrix corr_;
    pair_t(interaction_t *i1, interaction_t *i2, Index offset_i, Index offset_j,
           const pair_matrix &corr);
  };

  /// struct to store collected information for groups (e.g. crosscorrelations)
  struct group_t {
    std::list<interaction_t *> interactions_;
    group_matrix corr_;
    std::vector<pair_t> pairs_;
  };

  /// the options parsed from cg definition file
  Property options_;
  // we want to write out every so many frames
  Index write_every_;
  // we want do do block averaging -> clear averagings every write out
  bool do_blocks_;

  // number of frames we processed
  Index nframes_;
  Index nblock_;
  double subvol_rad_;
  Eigen::Vector3d boxc_;  // center of box
  bool do_vol_corr_;

  /// list of bonded interactions
  std::vector<Property *> bonded_;
  /// list of non-bonded interactions
  std::vector<Property *> nonbonded_;

  /// std::map ineteractionm-name to interaction
  std::map<std::string, std::unique_ptr<interaction_t>> interactions_;
  /// std::map group-name to group
  std::map<std::string, std::unique_ptr<group_t>> groups_;

  /// create a new interaction entry based on given options
  interaction_t *AddInteraction(Property *p);

  /// get group by name, creates one if it doesn't exist
  group_t *getGroup(const std::string &name);

  void WriteDist(const std::string &suffix = "");

  void ClearAverages();

  class Worker : public CsgApplication::Worker {
   public:
    std::vector<Histogram> current_hists_;
    RDFCalculator *rdfcalculator_;
    double cur_vol_;
    double cur_beadlist_1_count_;  // need to normalize to avg density for
                                   // subvol
    double cur_beadlist_2_count_;

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom) override;
    /// process non-bonded interactions for given frame
    void DoNonbonded(Topology *top);
    /// process bonded interactions for given frame
    void DoBonded(Topology *top);
  };
  /// update the correlations after interations were processed
  void DoCorrelations(RDFCalculator::Worker *worker);

  bool processed_some_frames_;

 public:
  std::unique_ptr<CsgApplication::Worker> ForkWorker();
  void MergeWorker(CsgApplication::Worker *worker_);
};

inline RDFCalculator::pair_t::pair_t(RDFCalculator::interaction_t *i1,
                                     RDFCalculator::interaction_t *i2,
                                     Index offset_i, Index offset_j,
                                     const pair_matrix &corr)
    : i1_(i1), i2_(i2), offset_i_(offset_i), offset_j_(offset_j), corr_(corr) {}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_RDF_CALCULATOR_H
