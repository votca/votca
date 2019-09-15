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
#ifndef VOTCA_XTP_QMATOM_H
#define VOTCA_XTP_QMATOM_H

#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>

namespace votca {
namespace xtp {

/**
 *    \brief container for QM atoms
 *
 *    Stores atom type, coordinates, charge
 */
class QMAtom {
  friend class ECPAOBasis;

 public:
  struct data {
    int index;
    char* element;
    double x;
    double y;
    double z;
    int nuccharge;
    int ecpcharge;
  };

  QMAtom(int index, std::string element, Eigen::Vector3d pos);

  QMAtom(data& d);

  const Eigen::Vector3d& getPos() const { return _pos; }

  void Translate(const Eigen::Vector3d& shift) { _pos += shift; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos);

  void setPos(const Eigen::Vector3d& position) { _pos = position; }

  const std::string& getElement() const { return _element; }

  int getId() const { return _index; }

  int getNuccharge() const { return _nuccharge - _ecpcharge; }

  std::string identify() const { return "qmatom"; }

  friend std::ostream& operator<<(std::ostream& out, const QMAtom& atom) {
    out << atom.getId() << " " << atom.getElement();
    out << " " << atom.getPos().x() << "," << atom.getPos().y() << ","
        << atom.getPos().z() << " " << atom.getNuccharge() << "\n";
    return out;
  }

 private:
  int _index;
  std::string _element;
  Eigen::Vector3d _pos;  // Bohr
  int _nuccharge = 0;    // nuc charge is set in aobasis fill and ecpfill
  int _ecpcharge = 0;

 public:
  void SetupCptTable(CptTable& table) const;

  void WriteData(data& d) const;

  void ReadData(data& d);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMATOM_H
