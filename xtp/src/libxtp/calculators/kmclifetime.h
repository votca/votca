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
#ifndef VOTCA_XTP_KMCLIFETIME_H
#define VOTCA_XTP_KMCLIFETIME_H

// Standard includes
#include <cstdint>
#include <string>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/lifetimesimulation.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Kinetic Monte Carlo of excitons (or other carriers) with
 * site-dependent lifetimes in a periodic box: lifetime, distance travelled
 * and diffusion length, where the excitons decay, site occupations
 * (LifetimeSimulation on a bulk ToFGraph with KMCEngine). For one exciton
 * also the exact averages without sampling.
 */
class KMCLifetime final : public QMCalculator {
 public:
  KMCLifetime() = default;
  ~KMCLifetime() override = default;

  std::string Identify() const override { return "kmclifetime"; }
  bool WriteToStateFile() const override { return false; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Evaluate(Topology& top) override;

 private:
  Logger log_;
  QMStateType carriertype_;
  double temperature_ = 0.0;  // kT in Hartree
  std::string injection_pattern_ = "*";
  std::string ignore_types_;
  LifetimeSimulation::Settings settings_;
  Index insertions_ = 0;
  bool exact_ = true;
  bool carrierenergy_ = false;
  double alpha_ = 0.05;
  std::string waitingtime_ = "exponential";
  std::string lifetimefile_;
  std::string trajectoryfile_;
  std::string ratefile_;
  std::string occfile_;
  std::string probfile_;
  std::string energyfile_;
  std::string resultfile_;

  std::vector<double> ReadDecayRates(const ToFGraph& graph) const;
  void WriteRates(const Topology& top);
  void WriteDecayProbabilities(const ToFGraph& graph,
                               const std::vector<double>& decay,
                               const LifetimeSimulation::Result* kmc,
                               const LifetimeSimulation::ExactResult* exact);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCLIFETIME_H
