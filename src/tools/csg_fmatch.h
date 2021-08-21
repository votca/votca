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

#ifndef VOTCA_CSG_CSG_FMATCH_H
#define VOTCA_CSG_CSG_FMATCH_H

// VOTCA includes
#include <votca/tools/cubicspline.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"
#include "votca/csg/trajectoryreader.h"

using namespace votca::csg;

using namespace std;

/**
    \brief Implements force matching algorithm using cubic spline basis set
 *
 *  Force matching method to obtain a coarse-grained force field is implemented
 *  using cubic spline basis set. Block averaging over trajectory blocks
 *  is used for calculating CG forces and their errors.
 *
 * \todo force matching needs a big cleanup!
 **/

class CGForceMatching : public CsgApplication {
 public:
  string ProgramName() override { return "csg_fmatch"; }
  void HelpText(ostream &out) override {
    out << "Perform force matching (also called multiscale coarse-graining)";
  }

  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }

  void Initialize(void) override;
  bool EvaluateOptions() override;

  /// \brief called before the first frame
  void BeginEvaluate(Topology *top, Topology *top_atom) override;
  /// \brief called after the last frame
  void EndEvaluate() override;
  /// \brief called for each frame which is mapped
  void EvalConfiguration(Topology *conf,
                         Topology *conf_atom = nullptr) override;
  /// \brief load options from the input file
  void LoadOptions(const string &file);

 protected:
  /// \brief structure, which contains CubicSpline object with related
  /// parameters
  struct SplineInfo {
    /// \brief constructor
    SplineInfo(votca::Index index, bool bonded_interaction,
               votca::Index matr_pos_col, votca::tools::Property *options);
    /// \brief number of spline functions
    votca::Index num_splinefun;
    /// \brief number of spline grid points
    votca::Index num_gridpoints;
    /// \brief number of grid points for output
    votca::Index num_outgrid;
    /// \brief interaction index
    votca::Index splineIndex;
    /// \brief true for bonded interactions, false for non-bonded
    bool bonded;
    /// \brief true if tabulated forces are periodic (e.g. for dihedral
    /// interactions)
    bool periodic;
    /// \brief true if non-bonded interaction is threebody interaction
    bool threebody;
    /// \brief additional variables for treating cutoff of
    double a;
    double sigma;
    double gamma;
    /// \brief CubicSpline object
    votca::tools::CubicSpline Spline;
    /// \brief position in the  A_ matrix (first coloumn which is occupied with
    /// this particular spline)
    votca::Index matr_pos;
    /// \brief dx for output. Calculated in the code
    double dx_out;
    /// \brief only for non-bonded interactions (seems like it is not used?)
    pair<votca::Index, votca::Index> beadTypes;

    /// \brief Result of 1 block calculation for f
    Eigen::VectorXd block_res_f;
    /// \brief Result of 1 block calculation for f''
    Eigen::VectorXd block_res_f2;
    /// \brief Final result: average over all blocks
    Eigen::VectorXd result;
    /// \brief accuracy of the final result
    Eigen::VectorXd error;
    /// \brief sum of all block_res (used to calculate error)
    Eigen::VectorXd resSum;
    /// \brief sum of all squares of block_res (used to calculate error)
    Eigen::VectorXd resSum2;

    // only needed for 3body nonbonded interactions as here force and potential
    // are calculated simultaneously
    /// \brief Final result of derivatives: average over all blocks
    Eigen::VectorXd resultDer;
    /// \brief accuracy of the final result
    Eigen::VectorXd errorDer;
    /// \brief sum of all block_res
    Eigen::VectorXd resSumDer;
    /// \brief sum of all squares of block_res (used to calculate error)
    Eigen::VectorXd resSumDer2;

    /// \brief Spline Name
    string splineName;
    /// \brief for non-bonded interactions: types of beads involved (type3 only
    /// used if threebody interaction)
    string type1, type2, type3;  //

    /// \brief pointer to Property object to handle input options
    votca::tools::Property *options_;
  };
  /// \brief Property object to handle input options
  votca::tools::Property options_;
  /// \brief list of bonded interactions
  std::vector<votca::tools::Property *> bonded_;
  /// \brief list of non-bonded interactions
  std::vector<votca::tools::Property *> nonbonded_;

  using SplineContainer = vector<SplineInfo>;
  /// \brief vector of SplineInfo * for all interactions
  SplineContainer splines_;

  /// \brief matrix used to store force matching equations
  Eigen::MatrixXd A_;
  /// \brief vector used to store reference forces on CG beads (from atomistic
  /// simulations)
  Eigen::VectorXd b_;
  /// \brief Solution of matrix equation  A_ *  x_ =  b_ : CG force-field
  /// parameters
  Eigen::VectorXd x_;  //
  /// \brief Additional matrix to handle constrained least squares fit
  /// contains constraints, which allow to get a real (smooth) spline (see VOTCA
  /// paper)
  Eigen::MatrixXd B_constr_;

  /// \brief Counter for trajectory frames
  votca::Index frame_counter_;
  /// \brief Number of CG beads
  votca::Index nbeads_;

  /// \brief accuracy for evaluating the difference in bead positions
  double dist_;

  /// \brief Flag: true for constrained least squares, false for simple least
  /// squares
  bool constr_least_sq_;
  /// \brief used in EvalConf to distinguish constrained and simple least
  /// squares
  votca::Index least_sq_offset_;
  /// \brief Number of frames used in one block for block averaging
  votca::Index nframes_;
  /// \brief Current number of blocks
  votca::Index nblocks_;

  /// \brief Counters for lines and columns in  B_constr_
  votca::Index line_cntr_, col_cntr_;

  bool has_existing_forces_;

  /// \brief Solves FM equations for one block and stores the results for
  /// further processing
  void FmatchAccumulateData();
  /// \brief Assigns smoothing conditions to matrices  A_ and  B_constr_
  void FmatchAssignSmoothCondsToMatrix(Eigen::MatrixXd &Matrix);
  /// \brief For each trajectory frame writes equations for bonded interactions
  /// to matrix  A_
  void EvalBonded(Topology *conf, SplineInfo *sinfo);
  /// \brief For each trajectory frame writes equations for non-bonded
  /// interactions to matrix  A_
  void EvalNonbonded(Topology *conf, SplineInfo *sinfo);
  /// \brief For each trajectory frame writes equations for non-bonded threebody
  /// interactions to matrix  A_
  void EvalNonbonded_Threebody(Topology *conf, SplineInfo *sinfo);
  /// \brief Write results to output files
  void WriteOutFiles();

  void OpenForcesTrajectory();

  Topology top_force_;
  std::unique_ptr<TrajectoryReader> trjreader_force_;
};

#endif  // VOTCA_CSG_CSG_FMATCH_H
