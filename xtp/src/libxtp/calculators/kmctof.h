/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_KMCTOF_H
#define VOTCA_XTP_KMCTOF_H

// Standard includes
#include <cstdint>
#include <string>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/kmcengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Time-of-flight kinetic Monte Carlo: carriers are injected in the
 * source region (one or more pulses), hop with site exclusion and leave at
 * the sink. Records per run the first arrival and the mean transit time per
 * carrier (the observables "first" and "transit" of grw) and the transient
 * current (Ramo-Shockley) and collected fraction on a logarithmic time axis.
 */
class KMCToF final : public QMCalculator {
 public:
  KMCToF() = default;
  ~KMCToF() override = default;

  std::string Identify() const override { return "kmctof"; }
  bool WriteToStateFile() const override { return false; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Evaluate(Topology& top) override;

 private:
  Logger log_;
  QMStateType carriertype_;
  double temperature_ = 0.0;                         // kT in Hartree
  Eigen::Vector3d field_ = Eigen::Vector3d::Zero();  // Hartree/bohr
  Eigen::Vector3d field_Vm_ = Eigen::Vector3d::Zero();
  ToFGraph::Options tof_;
  std::uint64_t seed_ = 1;
  Index runs_ = 1000;
  std::vector<KMCEngine::Pulse> pulses_;
  KMCEngine::Limits limits_;
  std::string waitingtime_ = "exponential";
  std::string runfile_;
  std::string currentfile_;
  std::string resultfile_;
  double tmin_ = 1e-15;
  double tmax_ = 1.0;
  Index bins_per_decade_ = 10;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCTOF_H
