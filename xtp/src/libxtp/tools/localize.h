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
#ifndef VOTCA_XTP_LOCALIZE_H
#define VOTCA_XTP_LOCALIZE_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class Localize final : public QMTool {
 public:
  Localize() = default;

  ~Localize() = default;

  std::string Identify() const { return "localize"; }

 private:
  Orbitals orbitals;
  Logger log;
  std::vector<Index> activeatoms;
  tools::Property options_;

 protected:
  void ParseOptions(const tools::Property& options);
  bool Run();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORB2FCHK_H
