/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#ifndef VOTCA_XTP_DECOMP_H
#define VOTCA_XTP_DECOMP_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

class Decomp final : public QMTool {
 public:
  Decomp() = default;

  ~Decomp() = default;

  std::string Identify() { return "decomp"; }

  private:
    Orbitals orbitals;
    Logger log;


 protected:
  void ParseOptions(const tools::Property& options);
  bool Run();
};
 
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORB2FCHK_H
