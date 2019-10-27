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
#ifndef VOTCA_XTP_TOPOLOGY_H
#define VOTCA_XTP_TOPOLOGY_H

#include <votca/csg/boundarycondition.h>
#include <votca/csg/openbox.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/triclinicbox.h>
#include <votca/xtp/qmnblist.h>

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

  Segment &getSegment(int id) { return _segments[id]; }
  const Segment &getSegment(int id) const { return _segments[id]; }

  std::vector<Segment> &Segments() { return _segments; }
  const std::vector<Segment> &Segments() const { return _segments; }

  // Periodic boundary: Can be 'open', 'orthorhombic', 'triclinic'
  Eigen::Vector3d PbShortestConnect(const Eigen::Vector3d &r1,
                                    const Eigen::Vector3d &r2) const;
  const Eigen::Matrix3d &getBox() const { return _bc->getBox(); }
  double BoxVolume() const { return _bc->BoxVolume(); }
  void setBox(const Eigen::Matrix3d &box,
              csg::BoundaryCondition::eBoxtype boxtype =
                  csg::BoundaryCondition::typeAuto);

  QMNBList &NBList() { return _nblist; }
  const QMNBList &NBList() const { return _nblist; }

  // Trajectory meta data: step number, time, frame (= Db ID)

  long getStep() const { return _step; }
  void setStep(long step) { _step = step; }
  double getTime() const { return _time; }
  void setTime(double time) { _time = time; }

  void WriteToCpt(CheckpointWriter &w) const;

  void WriteToPdb(std::string filename) const;

  void ReadFromCpt(CheckpointReader &r);

  double GetShortestDist(const Segment &seg1, const Segment &seg2) const;

  std::vector<const Segment *> FindAllSegmentsOnMolecule(
      const Segment &seg1, const Segment &seg2) const;

 protected:
  std::vector<Segment> _segments;

  std::unique_ptr<csg::BoundaryCondition> _bc = nullptr;
  QMNBList _nblist;

  double _time;
  long _step;

  csg::BoundaryCondition::eBoxtype AutoDetectBoxType(
      const Eigen::Matrix3d &box);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOPOLOGY_H
