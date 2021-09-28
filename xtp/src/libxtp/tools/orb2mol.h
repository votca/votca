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
#ifndef VOTCA_XTP_ORB2MOL_H
#define VOTCA_XTP_ORB2MOL_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class Orb2Mol final : public QMTool {
 public:
  Orb2Mol() = default;

  ~Orb2Mol() = default;

  std::string Identify() { return "orb2mol"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  // clang-format off
  std::string _moldenfile;
  std::string _orbfile;
  std::string _xyzfile;
  std::string _basisset_name;
  std::string _aux_basisset_name;
  Logger _log;
};

}  // namespace xtp
}  // namespace votca

#endif // VOTCA_XTP_ORB2MOL_H
