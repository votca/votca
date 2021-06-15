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
#ifndef VOTCA_XTP_GENCUBE_H
#define VOTCA_XTP_GENCUBE_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {
class AOBasis;

class GenCube final : public QMTool {
 public:
  GenCube() = default;

  ~GenCube() = default;

  std::string Identify() { return "gencube"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  void calculateCube();
  void subtractCubes();

  std::string orbfile_;
  std::string output_file_;
  std::string infile1_;
  std::string infile2_;

  bool dostateonly_;

  double padding_;
  Eigen::Array<Index, 3, 1> steps_;
  QMState state_;
  std::string mode_;
  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GENCUBE_H
