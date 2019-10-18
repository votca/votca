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

#ifndef _RDFCALCULATOR_H
#define _RDFCALCULATOR_H

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <math.h>
#include <votca/csg/csgapplication.h>
#include <votca/tools/average.h>
#include <votca/tools/histogramnew.h>
#include <votca/tools/property.h>

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

  void WriteEvery(int write_every) { _write_every = write_every; }
  void DoBlocks(bool do_blocks) { _do_blocks = do_blocks; }
  void DoVolumeCorrection(bool do_vol_corr) { _do_vol_corr = do_vol_corr; }
  void SetSubvolRadius(double r) { _subvol_rad = r; }
  double AnalyticVolumeCorrection(double t) {

    std::cout << "DBG " << t << " "
              << 1.0 / 24.0 *
                     (16.0 * t * t - 12.0 * t * t * t + t * t * t * t * t)
              << std::endl;
    return 1.0 / 24.0 * (16.0 * t * t - 12.0 * t * t * t + t * t * t * t * t);
  }

 protected:
  Average<double> _avg_vol;

  using group_matrix = Eigen::MatrixXd;
  using pair_matrix = Eigen::Block<group_matrix>;

  /// struct to store collected information for interactions
  struct interaction_t {
    int _index;
    Property *_p;
    HistogramNew _average;
    double _min, _max, _step;
    double _norm;
    bool _is_bonded;
    Average<double> _avg_beadlist_1_count;
    Average<double> _avg_beadlist_2_count;
  };

  // a pair of interactions which are correlated
  struct pair_t {
    interaction_t *_i1;
    interaction_t *_i2;
    int _offset_i, _offset_j;
    pair_matrix _corr;
    pair_t(interaction_t *i1, interaction_t *i2, int offset_i, int offset_j,
           const pair_matrix &corr);
  };

  /// struct to store collected information for groups (e.g. crosscorrelations)
  struct group_t {
    std::list<interaction_t *> _interactions;
    group_matrix _corr;
    std::vector<pair_t> _pairs;
  };

  /// the options parsed from cg definition file
  Property _options;
  // we want to write out every so many frames
  int _write_every;
  // we want do do block averaging -> clear averagings every write out
  bool _do_blocks;

  // number of frames we processed
  int _nframes;
  int _nblock;
  double _subvol_rad;
  Eigen::Vector3d _boxc;  // center of box
  bool _do_vol_corr;

  /// list of bonded interactions
  std::vector<Property *> _bonded;
  /// list of non-bonded interactions
  std::vector<Property *> _nonbonded;

  /// std::map ineteractionm-name to interaction
  std::map<std::string, interaction_t *> _interactions;
  /// std::map group-name to group
  std::map<std::string, group_t *> _groups;

  /// create a new interaction entry based on given options
  interaction_t *AddInteraction(Property *p);

  /// get group by name, creates one if it doesn't exist
  group_t *getGroup(const std::string &name);

  void WriteDist(const std::string &suffix = "");

  void ClearAverages();

  class Worker : public CsgApplication::Worker {
   public:
    std::vector<HistogramNew> _current_hists;
    RDFCalculator *_rdfcalculator;
    double _cur_vol;
    double _cur_beadlist_1_count;  // need to normalize to avg density for
                                   // subvol
    double _cur_beadlist_2_count;

    /// evaluate current conformation
    void EvalConfiguration(Topology *top, Topology *top_atom) override;
    /// process non-bonded interactions for given frame
    void DoNonbonded(Topology *top);
    /// process bonded interactions for given frame
    void DoBonded(Topology *top);
  };
  /// update the correlations after interations were processed
  void DoCorrelations(RDFCalculator::Worker *worker);

  bool _processed_some_frames;

 public:
  CsgApplication::Worker *ForkWorker();
  void MergeWorker(CsgApplication::Worker *worker);
};

inline RDFCalculator::pair_t::pair_t(RDFCalculator::interaction_t *i1,
                                     RDFCalculator::interaction_t *i2,
                                     int offset_i, int offset_j,
                                     const pair_matrix &corr)
    : _i1(i1), _i2(i2), _offset_i(offset_i), _offset_j(offset_j), _corr(corr) {}

}  // namespace csg
}  // namespace votca

#endif /* _RDFCALCULATOR_H */
