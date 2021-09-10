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
#ifndef VOTCA_XTP_BACKGROUND_H
#define VOTCA_XTP_BACKGROUND_H
#include <vector>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"

// Private VOTCA includes
#include "ewd_segment.h"
#include "unitcell.h"
#include "kspace.h"
#include "rspace.h"
#include "ewaldoptions.h"

namespace votca {
namespace xtp {

class Background {
 public:
  Background(Logger& log, const Eigen::Matrix3d& uc_matrix, const EwaldOptions options,
             std::vector<PolarSegment>& polar_background);

  Background(Logger& log, std::string& state_file) : log_(log) {
    readFromStateFile(state_file);
  }

  ~Background() = default;

  void Polarize();

  void writeToStateFile(std::string& state_file) { ; }

  void readFromStateFile(std::string& state_file) { ; }

 private:
  Index computeSystemSize(std::vector<EwdSegment>& ewaldSegments) const;
  Eigen::VectorXd solveLinearSystem(Eigen::MatrixXd A, Eigen::VectorXd b,
                                    Eigen::VectorXd guess);
  Logger& log_;
  UnitCell unit_cell_;
  EwaldOptions options_;
  std::vector<EwdSegment> ewald_background_;

};
}  // namespace xtp
}  // namespace votca

#endif