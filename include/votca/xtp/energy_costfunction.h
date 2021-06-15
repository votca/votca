/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_ENERGY_COSTFUNCTION_H
#define VOTCA_XTP_ENERGY_COSTFUNCTION_H

// Local VOTCA includes
#include "forces.h"
#include "optimiser_costfunction.h"
#include "orbitals.h"

namespace votca {
namespace xtp {

class Energy_costfunction : public Optimiser_costfunction {
 public:
  struct conv_paras {
    double deltaE;
    double RMSForce;
    double MaxForce;
    double RMSStep;
    double MaxStep;
    Index maxforceindex = 0;
    Index maxstepindex = 0;
  };

  Energy_costfunction(GWBSEEngine& gwbse_engine, StateTracker& tracker,
                      Orbitals& orbitals, Forces& force_engine)
      : gwbse_engine_(gwbse_engine),
        tracker_(tracker),
        orbitals_(orbitals),
        force_engine_(force_engine){};

  double EvaluateCost(const Eigen::VectorXd& parameters) override;

  Eigen::VectorXd EvaluateGradient(const Eigen::VectorXd& parameters) override;

  Index NumParameters() const override {
    return Index(orbitals_.QMAtoms().size() * 3);
  };

  bool Converged(const Eigen::VectorXd& delta_parameters, double delta_cost,
                 const Eigen::VectorXd& gradient) override;

  void ForcesReport() const { return force_engine_.Report(); }

  const conv_paras& getConvParas() const { return convpara_; }

  void setConvergenceParameters(const conv_paras& convergence) {
    convpara_ = convergence;
  }

  void setLog(Logger* pLog) { pLog_ = pLog; }

  void Report(const conv_paras& val);
  static void Vector2QMAtoms(const Eigen::VectorXd& pos, QMMolecule& atoms);
  static Eigen::VectorXd QMAtoms2Vector(QMMolecule& atoms);
  static Eigen::VectorXd Write3XMatrixToVector(const Eigen::MatrixX3d& matrix);

 private:
  static std::string Converged(double val, double limit);
  GWBSEEngine& gwbse_engine_;
  StateTracker& tracker_;
  Orbitals& orbitals_;
  Forces& force_engine_;
  Index iteration_ = 0;
  double energy_;

  conv_paras convpara_;

  Logger* pLog_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_ENERGY_COSTFUNCTION_H
