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
#ifndef _VOTCA_XTP_MOLPOL_H
#define _VOTCA_XTP_MOLPOL_H

#include <stdio.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmtool.h>

namespace votca {
namespace xtp {
class PolarRegion;
class MolPol : public QMTool {
 public:
  MolPol() : _input("", 0){};

  ~MolPol(){};

  std::string Identify() { return "molpol"; }

  void Initialize(tools::Property& options);
  bool Evaluate();

 private:
  void PrintPolarisation(const Eigen::Matrix3d& result) const;

  Eigen::Matrix3d CalcClassicalPol(const PolarSegment& input) const;
  Eigen::Vector3d Polarize(PolarRegion& pol,
                           const Eigen::Vector3d& ext_field) const;
  Logger _log;

  std::string _mps_output;
  PolarSegment _input;
  Eigen::Matrix3d _polarisation_target;

  Eigen::VectorXd _weights;

  tools::Property _polar_options;

  double _tolerance_pol = 1e-4;
  int _max_iter = 1000;
  double _alpha = 1.0;
};

}  // namespace xtp
}  // namespace votca

#endif
