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
#include "votca/xtp/region.h"
#include "votca/xtp/segid.h"

// Private VOTCA includes
#include "votca/xtp/ewaldoptions.h"
#include "votca/xtp/ewd_segment.h"
#include "votca/xtp/kspace.h"
#include "votca/xtp/rspace.h"
#include "votca/xtp/unitcell.h"

namespace votca {
namespace xtp {

class Background {
 public:
  Background(Logger& log, const Eigen::Matrix3d& uc_matrix,
             const EwaldOptions options,
             std::vector<PolarSegment>& polar_background);

  Background(Logger& log) : log_(log), rspace(log), kspace(log) {}

  ~Background() = default;

  Index size() const { return ewald_background_.size(); }

  void Polarize();

  void ApplyBackgroundFields(
      std::vector<std::unique_ptr<votca::xtp::Region>>& regions,
      const std::vector<std::vector<SegId>>& region_seg_ids);

  double interactionEnergy(std::vector<std::unique_ptr<Region>>& regions,
                           std::vector<std::vector<SegId>>& region_seg_ids,
                           tools::Property& results);

  void writeToStateFile(std::string state_file);

  void readFromStateFile(const std::string state_file);

  bool operator==(const Background& other) {
    if (other.unit_cell_.getMatrix() != this->unit_cell_.getMatrix()) {
      std::cout << "UnitCellWrong" << std::endl;
      return false;
    }
    if ((this->options_ != other.options_)) {
      std::cout << "optionsWrong" << std::endl;
      return false;
    }
    if (other.ewald_background_.size() != this->ewald_background_.size()) {
      std::cout << "BackgroundSizesWrong" << std::endl;
      return false;

    } else {
      for (Index i = 0; i < Index(this->ewald_background_.size()); ++i) {
        if (this->ewald_background_[i] != other.ewald_background_[i]) {
          std::cout << "Site " << i << " is wrong!!!" << std::endl;
          std::cout << "THIS: " <<  ewald_background_[i] << std::endl;
          std::cout << "OTHER: " <<  other.ewald_background_[i] << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  bool empty() const { return ewald_background_.size() == 0; }

 private:
  Index computeSystemSize(std::vector<EwdSegment>& ewaldSegments) const;
  Eigen::VectorXd solveLinearSystem(Eigen::MatrixXd A, Eigen::VectorXd b,
                                    Eigen::VectorXd guess);
  void bgFieldAtSegment(PolarSegment& seg, std::vector<SegId> pCloud_indices);
  Logger& log_;
  UnitCell unit_cell_;
  EwaldOptions options_;
  std::vector<EwdSegment> ewald_background_;
  RSpace rspace;
  KSpace kspace;
};
}  // namespace xtp
}  // namespace votca

#endif