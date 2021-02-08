/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#ifndef VOTCA_XTP_ERDIABATIZATIONFRAME_H
#define VOTCA_XTP_ERDIABATIZATIONFRAME_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmtool.h"
#include <votca/tools/types.h>
#include <votca/xtp/erdiabatization.h>
namespace votca {
namespace xtp {

class ERDiabatizationFrame final : public QMTool {
 public:
  ERDiabatizationFrame() = default;

  ~ERDiabatizationFrame() = default;

  std::string Identify() { return "erdiabatization"; }

 protected:
  void ParseOptions(const tools::Property& user_options) final;
  bool Run() final;

 private:
  std::string _orbfile1;
  std::string _orbfile2;
  QMStateType _qmtype;
  std::string _xml_output;  // .xml output
  Logger _log;
  ERDiabatization::options_erdiabatization _options;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERDIABATIZATIONFRAME_H
