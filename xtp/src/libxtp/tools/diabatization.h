/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#ifndef VOTCA_XTP_DIABATIZATION_H
#define VOTCA_XTP_DIABATIZATION_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmfragment.h"
#include "votca/xtp/qmtool.h"
#include <votca/tools/types.h>
#include <votca/xtp/erdiabatization.h>
#include <votca/xtp/fcddiabatization.h>
#include <votca/xtp/gmhdiabatization.h>

namespace votca {
namespace xtp {

class Diabatization final : public QMTool {
 public:
  Diabatization() = default;

  ~Diabatization() = default;

  std::string Identify() const { return "diabatization"; }

 protected:
  void ParseOptions(const tools::Property& user_options) final;
  bool Run() final;

 private:
  std::string orbfile1_;
  std::string orbfile2_;
  Logger log_;
  std::string method_;
  std::string qmtype_;

  Index state_idx_1_;
  Index state_idx_2_;

  bool isQMMM_;

  double E1_;
  double E2_;

  // for use in ER
  bool useRI_;
  // for use in FCD
  std::vector<QMFragment<BSE_Population> > fragments_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DIABATIZATION_H
