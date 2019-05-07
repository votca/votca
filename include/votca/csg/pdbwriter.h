/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _PDBWRITER_H
#define _PDBWRITER_H

#include <stdio.h>
#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/constants.h>

namespace votca {
namespace csg {

class PDBWriter : public TrajectoryWriter {
 public:
  void Open(std::string file, bool bAppend = false);
  void Close();

  void RegisteredAt(ObjectFactory<std::string, TrajectoryWriter> &factory) {}

  void Write(Topology *conf);

  template <class T>
  void WriteContainer(T &container);

  void WriteHeader(std::string header);

  void WriteBox(const Eigen::Matrix3d &box);

 private:
  template <class Atom>
  std::string getName(Atom &atom) {
    return atom.getElement();
  }

  std::string getName(Bead *bead) { return bead->getName(); }

  template <class T, class Atom>
  std::string getResname(T &container, Atom &atom) {
    return container.getName();
  }
  std::string getResname(Topology &conf, Bead *bead);

  template <class Atom>
  int getId(Atom &atom) {
    return atom.getId();
  }
  int getId(Bead *bead) { return bead->getId(); }

  template <class T, class Atom>
  int getResId(T &container, Atom &atom) {
    return container.getId();
  }
  int getResId(Topology &conf, Bead *bead) { return bead->getResnr() + 1; }

  template <class Atom>
  void writeSymmetry(Atom &atom) {
    return;
  }
  void writeSymmetry(Bead *bead);

  template <class Atom>
  Eigen::Vector3d getPos(Atom &atom) {
    return atom.getPos() * tools::conv::bohr2ang;
  }

  Eigen::Vector3d getPos(Bead *bead) {
    return bead->Pos() * tools::conv::nm2ang;
  }

  template <class T>
  T &getIterable(T &container) {
    return container;
  }

  BeadContainer &getIterable(Topology &top) { return top.Beads(); }

  std::ofstream _out;
};

template <class T>
inline void PDBWriter::WriteContainer(T &container) {
  boost::format atomfrmt(
      "ATOM  %1$5d %2$-4s %3$-3s %4$1s%5$4d    %6$8.3f%7$8.3f%8$8.3f\n");

  for (auto &atom : getIterable(container)) {
    Eigen::Vector3d r = getPos(atom);
    string resname = getResname(container, atom);
    string atomname = getName(atom);
    if (resname.size() > 3) {
      resname = resname.substr(0, 3);
    }
    if (atomname.size() > 4) {
      atomname = atomname.substr(0, 4);
    }

    _out << atomfrmt % (getId(atom) % 100000)  // atom serial number
                % atomname % resname % " "     // chain identifier 1 char
                % getResId(container, atom)    // residue sequence number
                % r.x() % r.y() % r.z();
    // we skip the charge
    writeSymmetry(atom);
  }
  _out << std::flush;
}
}  // namespace csg
}  // namespace votca

#endif /* _PDBWRITER_H */
