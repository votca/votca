/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#ifndef _VOTCA_XTP_GENCUBE_H
#define _VOTCA_XTP_GENCUBE_H

#include <votca/xtp/logger.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/qmtool.h>

namespace votca {
namespace xtp {
class AOBasis;

class GenCube : public QMTool {
 public:
  GenCube() = default;

  ~GenCube() override = default;

  std::string Identify() final { return "gencube"; }

  void Initialize(tools::Property& options) final;
  bool Evaluate() final;

 private:
  Eigen::VectorXd EvaluateBasisAtPosition(const AOBasis& dftbasis,
                                          const Eigen::Vector3d& pos);

  void calculateCube();
  void subtractCubes();

  std::string _orbfile;
  std::string _output_file;
  std::string _infile1;
  std::string _infile2;

  bool _dostateonly;

  double _padding;
  Index _xsteps;
  Index _ysteps;
  Index _zsteps;
  QMState _state;
  std::string _mode;
  Logger _log;
};

}  // namespace xtp
}  // namespace votca

#endif
