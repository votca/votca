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
#ifndef VOTCA_XTP_POLARREGION_H
#define VOTCA_XTP_POLARREGION_H

// Local VOTCA includes
#include "eeinteractor.h"
#include "energy_terms.h"
#include "hist.h"
#include "mmregion.h"

/**
 * \brief defines a polar region and of interacting electrostatic and induction
 * segments
 *
 *
 *
 */

namespace votca {
namespace xtp {
class QMRegion;
class PolarRegion;
class StaticRegion;

class PolarRegion : public MMRegion<PolarSegment> {
 public:
  PolarRegion(Index id, Logger& log) : MMRegion<PolarSegment>(id, log) {}

  std::string identify() const override { return "polarregion"; }

  void Initialize(const tools::Property& prop) override;

  bool Converged() const override;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override;

  void Reset() override;

  double Etotal() const override { return E_hist_.back().Etotal(); }

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

 protected:
  void AppendResult(tools::Property& prop) const override;
  double InteractwithQMRegion(const QMRegion& region) override;
  double InteractwithPolarRegion(const PolarRegion& region) override;
  double InteractwithStaticRegion(const StaticRegion& region) override;

 private:
  void CalcInducedDipoles();
  double StaticInteraction();
  void PolarInteraction_scf();

  double PolarEnergy_extern() const;
  eeInteractor::E_terms PolarEnergy() const;
  Index CalcPolDoF() const;

  Eigen::VectorXd CalcInducedDipoleInsideSegments() const;
  Eigen::VectorXd ReadInducedDipolesFromLastIteration() const;

  Eigen::VectorXd CalcInducedDipolesViaPCG(
      const Eigen::VectorXd& initial_guess);
  void WriteInducedDipolesToSegments(const Eigen::VectorXd& x);

  hist<Energy_terms> E_hist_;
  double deltaE_ = 1e-5;
  double deltaD_ = 1e-5;
  Index max_iter_ = 100;
  double exp_damp_ = 0.39;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_POLARREGION_H
