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
#ifndef __VOTCA_XTP_OVERLAP_FILTER_H
#define __VOTCA_XTP_OVERLAP_FILTER_H

#include <votca/xtp/statefilter_base.h>

namespace votca {
namespace xtp {

/**
    \brief overlap_filter
    tracks states according to their overlap with a previous state
 */

class Overlap_filter : public StateFilter_base {
 public:
  std::string Identify() const final { return "overlap"; }

  void Initialize(const tools::Property& options) final;

  void Info(Logger& log) const final;

  void UpdateHist(const Orbitals& orb, QMState state) final;

  bool NeedsInitialState() const final { return true; }

  std::vector<Index> CalcIndeces(const Orbitals& orb,
                                 QMStateType type) const final;

  void WriteToCpt(CheckpointWriter& w) final;

  void ReadFromCpt(CheckpointReader& r) final;

 private:
  Eigen::VectorXd CalculateOverlap(const Orbitals& orb, QMStateType type) const;
  Eigen::MatrixXd CalcAOCoeffs(const Orbitals& orb, QMStateType type) const;

  Eigen::MatrixXd CalcExcitonAORepresentation(const Orbitals& orb,
                                              QMStateType type) const;
  double _threshold = 0.0;

  Eigen::VectorXd _laststatecoeff;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H */
