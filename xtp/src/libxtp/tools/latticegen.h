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
#ifndef VOTCA_XTP_LATTICEGEN_H
#define VOTCA_XTP_LATTICEGEN_H

// Local VOTCA includes
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

/**
 * \brief Writes a state file with a surrogate lattice hopping model
 * (see LatticeModel), e.g. for time-of-flight tests.
 */
class LatticeGen final : public QMTool {
 public:
  LatticeGen() = default;
  ~LatticeGen() override = default;

  std::string Identify() const override { return "latticegen"; }

 protected:
  void ParseOptions(const tools::Property& options) override;
  bool Run() override;

 private:
  LatticeModel::Options lattice_;
  std::string statefile_;
  std::string energyfile_;
  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LATTICEGEN_H
