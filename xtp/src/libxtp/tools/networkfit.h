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
#ifndef VOTCA_XTP_NETWORKFIT_H
#define VOTCA_XTP_NETWORKFIT_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"
#include "votca/xtp/stochasticnetworkfit.h"
#include "votca/xtp/trapfinder.h"

namespace votca {
namespace xtp {

/**
 * \brief Fits a stochastic network model (StochasticNetworkFit) to a state
 * file and writes the parameter file for networkgen.
 */
class NetworkFit final : public QMTool {
 public:
  NetworkFit() = default;
  ~NetworkFit() override = default;

  std::string Identify() const override { return "networkfit"; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Run() override;

 private:
  std::string statefile_;
  Index frame_ = -1;
  std::string parameterfile_;
  std::string statsfile_;
  std::string energyfile_;
  StochasticNetworkFit::Settings settings_;
  Index validate_ = 4;  // model networks for the validation
  bool find_traps_ = true;
  TrapFinder::Settings traps_;
  double temperature_ = 300.0;
  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_NETWORKFIT_H
