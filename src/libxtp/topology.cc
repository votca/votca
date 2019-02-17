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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#include <votca/xtp/atom.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/topology.h>

#include <boost/lexical_cast.hpp>
#include <votca/tools/globals.h>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace xtp {

// +++++++++++++++++++++ //
// Clean-Up, Destruct    //
// +++++++++++++++++++++ //

void Topology::CleanUp() {

  _segments.clear();
  _atoms.clear();

  _bc.reset(nullptr);
  _nblist.Cleanup();
}

Segment &Topology::AddSegment(string segment_name) {
  int segment_id = _segments.size();
  Segment segment = Segment(segment_id, segment_name);
  segment.setTopology(this);
  _segments.push_back(segment);
  return _segments.back();
}

Atom &Topology::AddAtom(string atom_name) {
  int atom_id = _atoms.size();
  Atom atom = Atom(atom_id, atom_name);
  _atoms.push_back(atom);
  return _atoms.back();
}

// +++++++++++++++++ //
// Periodic Boundary //
// +++++++++++++++++ //

void Topology::setBox(const Eigen::Matrix3d &box,
                      csg::BoundaryCondition::eBoxtype boxtype) {

  // Determine box type automatically in case boxtype == typeAuto
  if (boxtype == csg::BoundaryCondition::typeAuto) {
    boxtype = AutoDetectBoxType(box);
  }

  if (_bc != nullptr) {
    if (votca::tools::globals::verbose) {
      cout << "Removing periodic box. Creating new... " << endl;
    }
  }
  _bc.reset(nullptr);
  switch (boxtype) {
    case csg::BoundaryCondition::typeTriclinic:
      _bc.reset(new csg::TriclinicBox());
      break;
    case csg::BoundaryCondition::typeOrthorhombic:
      _bc.reset(new csg::OrthorhombicBox());
      break;
    default:
      _bc.reset(new csg::OpenBox());
      break;
  }

  _bc->setBox(box);
}

csg::BoundaryCondition::eBoxtype Topology::AutoDetectBoxType(
    const Eigen::Matrix3d &box) {

  // Set box type to OpenBox in case "box" is the zero matrix,
  // to OrthorhombicBox in case "box" is a diagonal matrix,
  // or to TriclinicBox otherwise

  if (box.isApprox(Eigen::Matrix3d::Zero(), 1e-8)) {
    cout << "WARNING: No box vectors specified in trajectory."
            "Using open-box boundary conditions. "
         << endl;
    return csg::BoundaryCondition::typeOpen;
  }

  else if ((box - box.diagonal().asDiagonal())
               .isApprox(Eigen::Matrix3d::Zero(), 1e-8)) {
    return csg::BoundaryCondition::typeOrthorhombic;
  }

  else {
    return csg::BoundaryCondition::typeTriclinic;
  }

  return csg::BoundaryCondition::typeOpen;
}

Eigen::Vector3d Topology::PbShortestConnect(const Eigen::Vector3d &r1,
                                            const Eigen::Vector3d &r2) const {
  return _bc->BCShortestConnection(r1, r2).toEigen();
}

}  // namespace xtp
}  // namespace votca
