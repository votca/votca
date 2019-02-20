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
 * \brief Container for molecules, conjugated segments, rigid fragments,
 * and atoms.
 */
class Topology {
 public:
  Segment &AddSegment(std::string segment_name);

  Segment &getSegment(int id) { return _segments[id]; }

  std::vector<Segment> &Segments() { return _segments; }

  // Periodic boundary: Can be 'open', 'orthorhombic', 'triclinic'
  Eigen::Vector3d PbShortestConnect(const Eigen::Vector3d &r1,
                                    const Eigen::Vector3d &r2) const;
  const Eigen::Matrix3d &getBox() { return _bc->getBox().ToEigenMatrix(); }
  double BoxVolume() { return _bc->BoxVolume(); }
  void setBox(const Eigen::Matrix3d box,
              csg::BoundaryCondition::eBoxtype boxtype =
                  csg::BoundaryCondition::typeAuto);

  QMNBList &NBList() { return _nblist; }

  // Trajectory meta data: step number, time, frame (= Db ID)

  const int &getStep() { return _step; }
  void setStep(int step) { _step = step; }
  const double &getTime() { return _time; }
  void setTime(double time) { _time = time; }

  void WriteToCpt(CheckpointWriter &w) const;

  void ReadFromCpt(CheckpointReader &r);

 protected:
  std::vector<Segment> _segments;

  std::unique_ptr<csg::BoundaryCondition> _bc = nullptr;
  QMNBList _nblist;

  double _time;
  int _step;

  csg::BoundaryCondition::eBoxtype AutoDetectBoxType(
      const Eigen::Matrix3d &box);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOPOLOGY_H
