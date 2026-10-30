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
#ifndef VOTCA_XTP_KMCMULTIPLE_H
#define VOTCA_XTP_KMCMULTIPLE_H

// Standard includes
#include <cstdint>
#include <string>

// Local VOTCA includes
#include "votca/xtp/bulksimulation.h"
#include "votca/xtp/kmcengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Kinetic Monte Carlo of several electrons or holes with site
 * exclusion in a periodic box: drift velocity, mobility, diffusion tensor
 * and site occupations (BulkSimulation on a bulk ToFGraph with
 * KMCEngine).
 *
 * The run length (runtime) is a simulated time in seconds if it is at most
 * 100, otherwise a number of KMC events; outputtime is read in the same
 * unit.
 */
class KMCMultiple final : public QMCalculator {
 public:
  KMCMultiple() = default;
  ~KMCMultiple() override = default;

  std::string Identify() const override { return "kmcmultiple"; }
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
  std::string injection_pattern_ = "*";
  std::string ignore_types_;
  BulkSimulation::Settings settings_;
  bool steps_ = false;  // runtime given as number of events
  double runtime_ = 0.0;
  double outputtime_ = 0.0;
  std::string waitingtime_ = "exponential";
  std::string trajectoryfile_;
  std::string timefile_;
  std::string ratefile_;
  std::string occfile_;
  std::string resultfile_;
  bool exact_ = true;

  void WriteRates(const Topology& top);
  void WriteOccupation(const ToFGraph& graph,
                       const BulkSimulation::Result& result);
  void WriteTrajectory(const BulkSimulation::Result& result);
  void WriteTimeFile(const BulkSimulation::Result& result);
  /// sign of the carrier charge (0 for neutral excitations)
  double Charge() const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCMULTIPLE_H
