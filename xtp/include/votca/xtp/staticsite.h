
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
nn * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_STATICSITE_H
#define VOTCA_XTP_STATICSITE_H

// Local VOTCA includes
#include "eigen.h"
#include "qmatom.h"

namespace votca {
namespace xtp {

/**
\brief Class to represent Atom/Site in electrostatic

 The units are atomic units, e.g. Bohr, Hartree.
*/
class StaticSite {

 public:
  struct data {
    Index id;
    char* element;
    double posX;
    double posY;
    double posZ;

    Index rank;

    double Q00;
    double Q11c;
    double Q11s;
    double Q10;
    double Q20;
    double Q21c;
    double Q21s;
    double Q22c;
    double Q22s;
  };
  StaticSite(Index id, std::string element, Eigen::Vector3d pos)
      : id_(id), element_(element), pos_(pos){};

  StaticSite(Index id, std::string element)
      : StaticSite(id, element, Eigen::Vector3d::Zero()){};

  StaticSite(const data& d) { ReadData(d); }

  StaticSite(const QMAtom& atom, double charge)
      : StaticSite(atom.getId(), atom.getElement(), atom.getPos()) {
    setCharge(charge);
  }
  virtual ~StaticSite() = default;

 protected:
  StaticSite() = default;

 public:
  Index getId() const { return id_; }
  Index getRank() const { return rank_; }
  const std::string& getElement() const { return element_; }
  const Eigen::Vector3d& getPos() const { return pos_; }

  void setMultipole(const Vector9d& multipole, Index rank) {
    Q_ = multipole;
    rank_ = rank;
  }
  // sets rank to 0 as well
  void setCharge(double q) {
    Q_(0) = q;
    rank_ = 0;
  }

  void setPos(const Eigen::Vector3d& position) { pos_ = position; }

  // COORDINATES TRANSFORMATION
  void Translate(const Eigen::VectorXd& shift);
  virtual void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& refPos);

  // MULTIPOLES DEFINITION
  double getCharge() const { return Q_(0); }
  const Vector9d& Q() const {
    return Q_;
  }  // Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s,...[NOT following Stone
     // order for dipoles]

  virtual Eigen::Vector3d getDipole() const { return Q_.segment<3>(1); }

  Eigen::Matrix3d CalculateCartesianMultipole() const;

  static Eigen::VectorXd CalculateSphericalMultipole(
      const Eigen::Matrix3d& quad_cart);

  std::string WriteMpsLine(std::string unit = "bohr") const;

  static void SetupCptTable(CptTable& table);

  void WriteData(data& d) const;
  void ReadData(const data& d);
  virtual void setpolarization(const Eigen::Matrix3d&) { return; }

  virtual std::string identify() const { return "staticsite"; }

  friend std::ostream& operator<<(std::ostream& out, const StaticSite& site) {
    out << site.getId() << " " << site.getElement() << " " << site.getRank();
    out << " " << site.getPos().transpose() << "\n";
    return out;
  }

 protected:
  virtual std::string writepolarization() const;

  Index id_ = -1;
  std::string element_ = "";
  Eigen::Vector3d pos_ = Eigen::Vector3d::Zero();
  Index rank_ = 0;

  Vector9d Q_ = Vector9d::Zero();  // Q00,Q11c,Q11s,Q10,Q20,Q21c,Q21s,Q22c,Q22s
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATICSITE_H
