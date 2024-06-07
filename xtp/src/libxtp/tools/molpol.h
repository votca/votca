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
#ifndef VOTCA_XTP_MOLPOL_H
#define VOTCA_XTP_MOLPOL_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {
class PolarRegion;
class MolPol final : public QMTool {
 public:
  MolPol() : input_("", 0) {};

  ~MolPol() = default;

  std::string Identify() const { return "molpol"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  void Printpolarization(const Eigen::Matrix3d& result) const;

  Eigen::Matrix3d CalcClassicalPol(const PolarSegment& input) const;
  Eigen::Vector3d CalcInducedDipole(const PolarSegment& input,
                                    const Eigen::Vector3d& ext_field) const;
  Logger log_;

  std::string mps_output_;
  PolarSegment input_;
  Eigen::Matrix3d polarization_target_;

  Eigen::VectorXd weights_;

  tools::Property polar_options_;

  double tolerance_pol_ = 1e-4;
  Index max_iter_ = 1000;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLPOL_H
