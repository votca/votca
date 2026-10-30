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
#ifndef VOTCA_XTP_GRW_H
#define VOTCA_XTP_GRW_H

// Local VOTCA includes
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/tofgraph.h"

namespace votca {
namespace xtp {

/**
 * \brief Time of flight and mobility from the graph random walk: expected
 * hitting times of the sink region, solved as a sparse linear system instead
 * of sampled by kinetic Monte Carlo.
 */
class GRW final : public QMCalculator {
 public:
  GRW() = default;
  ~GRW() override = default;

  std::string Identify() const override { return "grw"; }
  bool WriteToStateFile() const override { return false; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Evaluate(Topology& top) override;

 private:
  void WriteSources(const ToFGraph& graph,
                    const GraphRandomWalk::Result& result) const;
  void WriteOccupation(const ToFGraph& graph,
                       const GraphRandomWalk::Result& result) const;

  Logger log_;
  QMStateType carriertype_;
  double temperature_ = 0.0;                         // kT in Hartree
  Eigen::Vector3d field_ = Eigen::Vector3d::Zero();  // Hartree/bohr
  Eigen::Vector3d field_Vm_ = Eigen::Vector3d::Zero();
  ToFGraph::Options tof_;
  Index carriers_ = 1;
  GraphRandomWalk::Solver solver_ = GraphRandomWalk::Solver::Auto;
  GraphRandomWalk::Settings settings_;
  GraphRandomWalk::Observable observable_ = GraphRandomWalk::Observable::First;
  std::string sourcefile_;
  std::string occupationfile_;
  std::string resultfile_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GRW_H
