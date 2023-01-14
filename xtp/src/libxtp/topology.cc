/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

// Third party includes
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/csg/pdbwriter.h>
#include <votca/tools/globals.h>
#include <votca/tools/version.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/checkpointwriter.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {

Topology::Topology(const Topology &top) {
  segments_ = top.segments_;
  time_ = top.time_;
  step_ = top.step_;
  this->setBox(top.getBox());
  for (const QMPair *pair : top.NBList()) {
    const Segment &seg1 = segments_[pair->Seg1()->getId()];
    const Segment &seg2 = segments_[pair->Seg2()->getId()];
    nblist_.Add(seg1, seg2, pair->R());
  }
}

Topology &Topology::operator=(const Topology &top) {
  if (&top != this) {
    segments_ = top.segments_;
    time_ = top.time_;
    step_ = top.step_;
    this->setBox(top.getBox());
    nblist_.Cleanup();
    for (const QMPair *pair : top.NBList()) {
      const Segment &seg1 = segments_[pair->Seg1()->getId()];
      const Segment &seg2 = segments_[pair->Seg2()->getId()];
      nblist_.Add(seg1, seg2, pair->R());
    }
  }
  return *this;
}

Segment &Topology::AddSegment(std::string segment_name) {
  Index segment_id = Index(segments_.size());
  segments_.push_back(Segment(segment_name, segment_id));
  return segments_.back();
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

  if (bc_ != nullptr) {
    if (Log::verbose()) {
      std::cout << "Removing periodic box. Creating new... " << std::endl;
    }
  }
  bc_.reset(nullptr);
  switch (boxtype) {
    case csg::BoundaryCondition::typeTriclinic:
      bc_.reset(new csg::TriclinicBox());
      break;
    case csg::BoundaryCondition::typeOrthorhombic:
      bc_.reset(new csg::OrthorhombicBox());
      break;
    default:
      bc_.reset(new csg::OpenBox());
      break;
  }

  bc_->setBox(box);
}

csg::BoundaryCondition::eBoxtype Topology::AutoDetectBoxType(
    const Eigen::Matrix3d &box) {

  // Set box type to OpenBox in case "box" is the zero matrix,
  // to OrthorhombicBox in case "box" is a diagonal matrix,
  // or to TriclinicBox otherwise

  if (box.isApproxToConstant(0)) {
    std::cout << "WARNING: No box vectors specified in trajectory."
                 "Using open-box boundary conditions. "
              << std::endl;
    return csg::BoundaryCondition::typeOpen;
  }

  else if ((box - Eigen::Matrix3d(box.diagonal().asDiagonal()))
               .isApproxToConstant(0)) {
    return csg::BoundaryCondition::typeOrthorhombic;
  }

  else {
    return csg::BoundaryCondition::typeTriclinic;
  }

  return csg::BoundaryCondition::typeOpen;
}

Eigen::Vector3d Topology::PbShortestConnect(const Eigen::Vector3d &r1,
                                            const Eigen::Vector3d &r2) const {
  return bc_->BCShortestConnection(r1, r2);
}

double Topology::GetShortestDist(const Segment &seg1,
                                 const Segment &seg2) const {
  double R2 = std::numeric_limits<double>::max();
  for (const Atom &atom1 : seg1) {
    for (const Atom &atom2 : seg2) {
      double R2_test =
          PbShortestConnect(atom1.getPos(), atom2.getPos()).squaredNorm();
      if (R2_test < R2) {
        R2 = R2_test;
      }
    }
  }
  return std::sqrt(R2);
}

std::vector<const Segment *> Topology::FindAllSegmentsOnMolecule(
    const Segment &seg1, const Segment &seg2) const {
  const std::vector<Index> &ids1 = seg1.getMoleculeIds();
  const std::vector<Index> &ids2 = seg2.getMoleculeIds();
  std::vector<Index> common_elements;
  std::set_intersection(ids1.begin(), ids1.end(), ids2.begin(), ids2.end(),
                        std::back_inserter(common_elements));
  std::vector<const Segment *> results;
  if (common_elements.empty() || common_elements.size() > 1) {
    return results;
  }
  Index molid = common_elements[0];

  for (const Segment &seg : segments_) {
    if (std::find(seg.getMoleculeIds().begin(), seg.getMoleculeIds().end(),
                  molid) != seg.getMoleculeIds().end()) {
      results.push_back(&seg);
    }
  }
  return results;
}

void Topology::WriteToPdb(std::string filename) const {

  csg::PDBWriter writer;
  writer.Open(filename, false);
  writer.WriteHeader("Frame:" + std::to_string(this->getStep()));
  writer.WriteBox(this->getBox() * tools::conv::bohr2ang);
  for (const Segment &seg : segments_) {
    writer.WriteContainer(seg);
  }
  writer.Close();
}

void Topology::WriteToCpt(CheckpointWriter &w) const {
  w(votca::tools::ToolsVersionStr(), "XTPVersion");
  w(topology_version(), "version");
  w(time_, "time");
  w(step_, "step");
  w(this->getBox(), "box");
  CheckpointWriter ww = w.openChild("segments");
  for (const Segment &seg : segments_) {
    CheckpointWriter u = ww.openChild("segment" + std::to_string(seg.getId()));
    seg.WriteToCpt(u);
  }
  CheckpointWriter www = w.openChild("neighborlist");
  nblist_.WriteToCpt(www);
}

void Topology::ReadFromCpt(CheckpointReader &r) {
  r(time_, "time");
  r(step_, "step");
  Eigen::Matrix3d box;
  r(box, "box");
  setBox(box);
  CheckpointReader v = r.openChild("segments");
  segments_.clear();
  Index count = v.getNumDataSets();
  segments_.reserve(count);
  for (Index i = 0; i < count; i++) {
    CheckpointReader w = v.openChild("segment" + std::to_string(i));
    segments_.push_back(Segment(w));
  }
  CheckpointReader rr = r.openChild("neighborlist");
  nblist_.ReadFromCpt(rr, segments_);
}

}  // namespace xtp
}  // namespace votca
