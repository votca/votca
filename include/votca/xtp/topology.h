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
#ifndef VOTCA_XTP_TOPOLOGY_H
#define VOTCA_XTP_TOPOLOGY_H

// VOTCA includes
#include <votca/csg/boundarycondition.h>
#include <votca/csg/openbox.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/triclinicbox.h>

// Local VOTCA includes
#include "qmnblist.h"

namespace votca {
namespace xtp {

class Segment;
/**
 * \brief Container for segments and box
 * and atoms.
 */
class Topology {
 public:
  Topology() = default;

  Topology(const Topology &top);

  Topology &operator=(const Topology &top);

  // I do not have to manually make a move constructor or move assignment
  // operator or destructor because I only have to reassign pointers in qmnblist
  // object

  Segment &AddSegment(std::string segment_name);

  Segment &getSegment(Index id) { return segments_[id]; }
  const Segment &getSegment(Index id) const { return segments_[id]; }

  std::vector<Segment> &Segments() { return segments_; }
  const std::vector<Segment> &Segments() const { return segments_; }

  // Periodic boundary: Can be 'open', 'orthorhombic', 'triclinic'
  Eigen::Vector3d PbShortestConnect(const Eigen::Vector3d &r1,
                                    const Eigen::Vector3d &r2) const;
  const Eigen::Matrix3d &getBox() const { return bc_->getBox(); }
  double BoxVolume() const { return bc_->BoxVolume(); }
  void setBox(const Eigen::Matrix3d &box,
              csg::BoundaryCondition::eBoxtype boxtype =
                  csg::BoundaryCondition::typeAuto);

  QMNBList &NBList() { return nblist_; }
  const QMNBList &NBList() const { return nblist_; }

  // Trajectory meta data: step number, time, frame (= Db ID)

  Index getStep() const { return step_; }
  void setStep(Index step) { step_ = step; }
  double getTime() const { return time_; }
  void setTime(double time) { time_ = time; }

  void WriteToCpt(CheckpointWriter &w) const;

  void WriteToPdb(std::string filename) const;

  void ReadFromCpt(CheckpointReader &r);

  double GetShortestDist(const Segment &seg1, const Segment &seg2) const;

  std::vector<const Segment *> FindAllSegmentsOnMolecule(
      const Segment &seg1, const Segment &seg2) const;

 private:
  std::vector<Segment> segments_;

  std::unique_ptr<csg::BoundaryCondition> bc_ = nullptr;
  QMNBList nblist_;

  double time_;
  Index step_;

  csg::BoundaryCondition::eBoxtype AutoDetectBoxType(
      const Eigen::Matrix3d &box);

  static constexpr int topology_version() { return 1; }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOPOLOGY_H
