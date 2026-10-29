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
#ifndef VOTCA_XTP_NETWORKGEN_H
#define VOTCA_XTP_NETWORKGEN_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"
#include "votca/xtp/stochasticnetwork.h"

namespace votca {
namespace xtp {

/**
 * \brief Writes a state file with a stochastic network model (see
 * StochasticNetwork) of any size from a parameter file of networkfit.
 */
class NetworkGen final : public QMTool {
 public:
  NetworkGen() = default;
  ~NetworkGen() override = default;

  std::string Identify() const override { return "networkgen"; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Run() override;

 private:
  std::string parameterfile_;
  std::string statefile_;
  std::string energyfile_;
  Index sites_ = 0;
  std::vector<double> box_;  // nm, empty: from the density of the file
  Index seed_ = 0;           // 0: from the file
  double temperature_ = 300.0;
  double trap_factor_ = 1.0;
  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_NETWORKGEN_H
