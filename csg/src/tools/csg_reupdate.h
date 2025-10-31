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

#ifndef VOTCA_CSG_CSG_REUPDATE_H
#define VOTCA_CSG_CSG_REUPDATE_H

// Third party includes
#include <boost/program_options.hpp>

// VOTCA includes
#include <votca/tools/histogram.h>
#include <votca/tools/property.h>
#include <votca/tools/table.h>

// Local VOTCA includes
#include "votca/csg/csgapplication.h"
#include "votca/csg/potentialfunctions/potentialfunction.h"
#include "votca/csg/potentialfunctions/potentialfunctioncbspl.h"
#include "votca/csg/potentialfunctions/potentialfunctionlj126.h"
#include "votca/csg/potentialfunctions/potentialfunctionljg.h"
#include "votca/csg/topologyreader.h"

using namespace votca::csg;
using namespace votca::tools;

struct PotentialInfo {

  PotentialInfo(votca::Index index, bool bonded_, votca::Index vec_pos_,
                std::string &param_in_ext_, Property *options,
                bool gentable = false);

  votca::Index potentialIndex;
  bool bonded;
  PotentialFunction *ucg;
  votca::Index vec_pos;
  std::pair<votca::Index, votca::Index> beadTypes;

  std::string potentialName;
  std::string potentialFunction;
  std::string type1, type2;

  double rmin, rcut;

  Property *options_;
};

class CsgREupdate : public CsgApplication {
 public:
  std::string ProgramName() override { return "csg_reupdate"; }
  void HelpText(std::ostream &out) override {
    out << "computes relative entropy update.";
  }

  bool DoTrajectory() override { return true; }

  bool DoMapping() override { return false; }

  bool DoThreaded() override { return true; }
  bool SynchronizeThreads() override { return false; }

  bool NeedsTopology() override { return false; }

  void Initialize() override;
  bool EvaluateOptions() override;
  void BeginEvaluate(Topology *top, Topology *top_atom = nullptr) override;
  void LoadOptions(const std::string &file);

  void Run() override;

  void EndEvaluate() override;
  std::unique_ptr<CsgApplication::Worker> ForkWorker(void) override;
  void MergeWorker(Worker *worker) override;

 private:
 protected:
  Property options_;
  std::vector<Property *> nonbonded_;

  using PotentialContainer = std::vector<PotentialInfo *>;
  PotentialContainer potentials_;

  votca::Index nlamda_;
  Eigen::VectorXd lamda_;
  //  HS_ is a symmetric matrix
  Eigen::MatrixXd HS_;
  Eigen::VectorXd DS_;
  Eigen::VectorXd dUFrame_;
  bool hessian_check_;

  double UavgAA_;
  double UavgCG_;
  double beta_;
  double relax_;
  votca::Index nframes_;

  bool gentable_;

  std::vector<Table *> aardfs_;
  std::vector<double *> aardfnorms_;

  // file extension for the inputs/outputs
  std::string param_in_ext_, param_out_ext_;
  std::string pot_out_ext_;
  std::string rdf_ext_;

  void WriteOutFiles();
  void EvalBonded(Topology *conf, PotentialInfo *potinfo);
  void EvalNonbonded(Topology *conf, PotentialInfo *potinfo);

  // Compute Avg U, dU, and d2U values in reference AA ensemble
  void AAavgBonded(PotentialInfo *potinfo);
  void AAavgNonbonded(PotentialInfo *potinfo);

  // Formulates  HS_ dlamda = -  DS_ system of Lin Eq.
  void REFormulateLinEq();

  // Solve  HS_ dlamda = -  DS_ and update  lamda_
  void REUpdateLamda();
};

class CsgREupdateWorker : public CsgApplication::Worker {
 public:
  ~CsgREupdateWorker() override = default;

  Property options_;
  std::vector<Property *> nonbonded_;

  using PotentialContainer = std::vector<PotentialInfo *>;
  PotentialContainer potentials_;

  votca::Index nlamda_;
  Eigen::VectorXd lamda_;
  Eigen::MatrixXd HS_;
  Eigen::VectorXd DS_;
  Eigen::VectorXd dUFrame_;

  double UavgCG_;
  double beta_;
  votca::Index nframes_;

  void EvalConfiguration(Topology *conf, Topology *conf_atom) override;
  void EvalBonded(Topology *conf, PotentialInfo *potinfo);
  void EvalNonbonded(Topology *conf, PotentialInfo *potinfo);
};

#endif  // VOTCA_CSG_CSG_REUPDATE_H
