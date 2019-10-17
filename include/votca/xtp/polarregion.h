/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/xtp/eeinteractor.h>
#include <votca/xtp/energy_terms.h>
#include <votca/xtp/hist.h>
#include <votca/xtp/mmregion.h>
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
  PolarRegion(int id, Logger& log) : MMRegion<PolarSegment>(id, log) {}

  std::string identify() const override { return "polar"; }

  void Initialize(const tools::Property& prop) override;

  bool Converged() const override;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override;

  void Reset() override;

  double Etotal() const override { return _E_hist.back().Etotal(); }

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

  hist<Energy_terms> _E_hist;
  double _deltaE = 1e-5;
  double _deltaD = 1e-5;
  int _max_iter = 100;
  double _exp_damp = 0.39;
};

}  // namespace xtp
}  // namespace votca

#endif /* VOTCA_XTP_MMREGION_H */
