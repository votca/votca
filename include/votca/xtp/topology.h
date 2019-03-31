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

#ifndef __VOTCA_XTP_TOPOLOGY_H
#define __VOTCA_XTP_TOPOLOGY_H

#include <string>
#include <vector>

#include <votca/tools/property.h>

#include <votca/csg/boundarycondition.h>
#include <votca/csg/openbox.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/triclinicbox.h>

#include <votca/xtp/qmnblist.h>

namespace CSG = votca::csg;

namespace votca {
namespace xtp {

class Molecule;
class Segment;
class SegmentType;
class Fragment;
class Atom;
class APolarSite;
/**
 * \brief Container for molecules, conjugated segments, rigid fragments,
 * and atoms.
 */
class Topology {
 public:
  Topology();
  ~Topology();

  // Population: Molecules, segments, fragments, atoms

  Molecule *AddMolecule(std::string molecule_name);
  Segment *AddSegment(std::string segment_name);
  Atom *AddAtom(std::string atom_name);
  Fragment *AddFragment(std::string fragment_name);
  APolarSite *AddAPolarSite(std::string siteName);
  SegmentType *AddSegmentType(std::string typeName);

  Molecule *getMolecule(int id) { return _molecules[id - 1]; }
  Segment *getSegment(int id) { return _segments[id - 1]; }
  Fragment *getFragment(int id) { return _fragments[id - 1]; }
  Atom *getAtom(int id) { return _atoms[id - 1]; }
  APolarSite *getAPolarSite(int id) { return _apolarSites[id - 1]; }
  SegmentType *getSegmentType(int id) { return _segmentTypes[id - 1]; }

  std::vector<Atom *> &Atoms() { return _atoms; }
  std::vector<Fragment *> &Fragments() { return _fragments; }
  std::vector<Segment *> &Segments() { return _segments; }
  std::vector<Molecule *> &Molecules() { return _molecules; }
  std::vector<APolarSite *> &APolarSites() { return _apolarSites; }
  std::vector<SegmentType *> &SegmentTypes() { return _segmentTypes; }

  bool Rigidify();
  void setCanRigidify(bool yesno) { _canRigidify = yesno; }
  const bool &canRigidify() { return _canRigidify; }
  const bool &isRigid() { return _isRigid; }
  void setIsEStatified(bool yesno) { _isEStatified = yesno; }
  const bool &isEStatified() { return _isEStatified; }

  // Periodic boundary: Can be 'open', 'orthorhombic', 'triclinic'

  votca::tools::vec PbShortestConnect(const votca::tools::vec &r1,
                                      const votca::tools::vec &r2) const;
  const votca::tools::matrix &getBox() { return _bc->getBox(); }
  double BoxVolume() { return _bc->BoxVolume(); }
  void setBox(const votca::tools::matrix &box,
              CSG::BoundaryCondition::eBoxtype boxtype =
                  CSG::BoundaryCondition::typeAuto);

  QMNBList &NBList() { return _nblist; }

  // Trajectory meta data: step number, time, frame (= Db ID)

  const int &getStep() { return _step; }
  void setStep(int step) { _step = step; }
  const double &getTime() { return _time; }
  void setTime(double time) { _time = time; }

  int getDatabaseId() { return _db_id; };
  void setDatabaseId(int id) { _db_id = id; }
  void CleanUp();

 protected:
  std::vector<Molecule *> _molecules;
  std::vector<Segment *> _segments;
  std::vector<Fragment *> _fragments;
  std::vector<Atom *> _atoms;
  std::vector<APolarSite *> _apolarSites;
  std::vector<SegmentType *> _segmentTypes;

  int _db_id;
  bool _hasPb;
  CSG::BoundaryCondition *_bc;
  QMNBList _nblist;

  bool _canRigidify;
  bool _isRigid;
  bool _isEStatified;

  double _time;
  int _step;

  CSG::BoundaryCondition::eBoxtype AutoDetectBoxType(
      const votca::tools::matrix &box);
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_TOPOLOGY_H */
