/*
 *            Copyright 2016 The MUSCET Development Team
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
#ifndef VOTCA_XTP_GYRATION_H
#define VOTCA_XTP_GYRATION_H

#include <stdio.h>
#include <votca/xtp/density_integration.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
namespace votca {
namespace xtp {

class Density2Gyration {
 public:
  Density2Gyration(Logger& log) : _log(log){};

  std::string Identify() { return "density2gyration"; }

  void Initialize(tools::Property& options);

  void AnalyzeDensity(const Orbitals& orbitals);

 private:
  void ReportAnalysis(std::string label, const Gyrationtensor& gyro,
                      const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& es);
  void AnalyzeGeometry(const QMMolecule& atoms);

  QMState _state;
  bool _dostateonly;
  std::string _integrationmethod;
  std::string _gridsize;
  Logger& _log;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GYRATION_H
