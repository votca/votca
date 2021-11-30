/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local private VOTCA includes
#include "votca/xtp/background.h"
#include "votca/tools/timer.h"

namespace votca {
namespace xtp {
Background::Background(Logger& log, const Eigen::Matrix3d& uc_matrix,
                       const EwaldOptions options,
                       std::vector<PolarSegment>& polar_background)
    : log_(log), unit_cell_(uc_matrix), options_(options) {
  for (const PolarSegment& pseg : polar_background) {
    BGSegment bgseg(pseg);
    background_segments_.emplace_back(bgseg);
  }
  // Place atoms in the simulation box
  for (BGSegment& seg : background_segments_) {
    for (BGSite& site : seg) {
      site.updatePos(unit_cell_.placeCoordInBox(site.getPos()));
    }
    // Should be removed before merge in master
    seg.calcPos();
  }

  // Now we print some basic info about the unit cell
  XTP_LOG(Log::error, log_) << unit_cell_ << std::endl;
  setupNeighbourList();
}

void Background::Polarize() {
  Timer timer;
  computeStaticField();

  XTP_LOG(Log::error, log_) << "Succesfully polarized the background in "
                            << timer.elapsedTimeAsString() << std::endl;
}

void Background::computeStaticField() { 
  Timer timer;
// #pragma omp parallel for
  for (Index segId = 0; segId < Index(background_segments_.size()); ++segId) {
    BGSegment& currentSeg = background_segments_[segId];
    for (const Neighbour& neighbour : nbList_.getNeighboursOf(segId)) {
      BGSegment& nbSeg = background_segments_[neighbour.getId()];
      for (BGSite& site : currentSeg) {
        for (BGSite& nbSite : nbSeg) {
          ;//site.addToStaticField(fieldAtBy(site.mp(), nbSite.mp(), neighbour.getDr()));
        }
      }
    }
  }

  XTP_LOG(Log::error, log_) << "Computed permanent fields in "
                            << timer.elapsedTimeAsString() << std::endl;
 }

void Background::setupNeighbourList() {
  Timer timer;
  std::array<Index, 3> maxCopies =
      unit_cell_.getNrOfRealSpaceCopiesForCutOff(options_.r_cutoff);
  XTP_LOG(Log::error, log_) << "Setting up neighbourlist with [" << maxCopies[0]
                            << "," << maxCopies[1] << "," << maxCopies[2]
                            << "] copies of the unit cell." << std::endl;
  nbList_.setSize(background_segments_.size());

#pragma omp parallel for
  for (Index segId = 0; segId < static_cast<Index>(background_segments_.size());
       ++segId) {
    BGSegment& currentSeg = background_segments_[segId];
    for (const BGSegment& seg : background_segments_) {
      Eigen::Vector3d minImage_dr = unit_cell_.minImage(seg, currentSeg);
      Eigen::Vector3d dr_dir = seg.getPos() - currentSeg.getPos();
      // triple for-loop is over all unitcell copies
      for (Index n1 = -maxCopies[0]; n1 < maxCopies[0]; ++n1) {
        for (Index n2 = -maxCopies[1]; n2 < maxCopies[1]; ++n2) {
          for (Index n3 = -maxCopies[2]; n3 < maxCopies[2]; ++n3) {
            if (n1 == 0 && n2 == 0 && n3 == 0 &&
                currentSeg.getId() == seg.getId()) {
              continue;
            }
            // LVector is the vector pointing to the n1,n2,n3th box
            Eigen::Vector3d lvector = unit_cell_.getLVector(n1, n2, n3);
            Eigen::Vector3d dr_l = minImage_dr + lvector;
            Eigen::Vector3d shift = dr_l - dr_dir;
            double dist = dr_l.norm();
            if (dist < options_.r_cutoff) {
              nbList_.addNeighbourTo(segId,
                                     Neighbour(seg.getId(), dr_l, shift, dist));
            }
          }
        }
      }
    }
    nbList_.sortOnDistance(segId);
  }
  XTP_LOG(Log::error, log_) << "Setup neighbourlist in "
                            << timer.elapsedTimeAsString() << std::endl;
}

}  // namespace xtp

}  // namespace votca
