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
#ifndef VOTCA_XTP_QMATOM_H
#define VOTCA_XTP_QMATOM_H

// VOTCA includes
#include "eigen.h"
#include <votca/tools/elements.h>
#include <votca/tools/types.h>

namespace votca {
namespace xtp {
class CptTable;
/**
 *    \brief container for QM atoms
 *
 *    Stores atom type, coordinates, charge
 */
class QMAtom {
  friend class ECPAOBasis;

 public:
  struct data {
    Index index;
    char* element;
    double x;
    double y;
    double z;
    Index nuccharge;
    Index ecpcharge;
  };

  QMAtom(Index index, std::string element, Eigen::Vector3d pos);

  QMAtom(const data& d);

  const Eigen::Vector3d& getPos() const { return pos_; }

  void Translate(const Eigen::Vector3d& shift) { pos_ += shift; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos);

  void setPos(const Eigen::Vector3d& position) { pos_ = position; }

  const std::string& getElement() const { return element_; }

  Index getId() const { return index_; }

  Index getNuccharge() const { return nuccharge_ - ecpcharge_; }

  Index getElementNumber() const { return nuccharge_; }

  std::string identify() const { return "qmatom"; }

  friend std::ostream& operator<<(std::ostream& out, const QMAtom& atom) {
    out << atom.getId() << " " << atom.getElement();
    out << " " << atom.getPos().x() << "," << atom.getPos().y() << ","
        << atom.getPos().z() << " " << atom.getNuccharge() << "\n";
    return out;
  }

 private:
  Index index_;
  std::string element_;
  Eigen::Vector3d pos_;  // Bohr
  Index nuccharge_ = 0;
  Index ecpcharge_ = 0;  // ecp charge is set in ecpaobasis.fill

 public:
  static void SetupCptTable(CptTable& table);

  void WriteData(data& d) const;

  void ReadData(const data& d);
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMATOM_H
