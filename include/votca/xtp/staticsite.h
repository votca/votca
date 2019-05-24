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
nn * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef __VOTCA_XTP_STATICSITE_H
#define __VOTCA_XTP_STATICSITE_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/qmatom.h>

namespace votca {
namespace xtp {

/**
\brief Class to represent Atom/Site in electrostatic

 The units are atomic units, e.g. Bohr, Hartree.
*/
class StaticSite {

 public:
  struct data {
    int id;
    char* element;
    double posX;
    double posY;
    double posZ;

    int rank;

    double multipoleQ00;
    double multipoleQ11c;
    double multipoleQ11s;
    double multipoleQ10;
    double multipoleQ20;
    double multipoleQ21c;
    double multipoleQ21s;
    double multipoleQ22c;
    double multipoleQ22s;

    double fieldX;
    double fieldY;
    double fieldZ;

    double phi;
  };
  StaticSite(int id, std::string element, Eigen::Vector3d pos)
      : _id(id), _element(element), _pos(pos){};

  StaticSite(int id, std::string element)
      : StaticSite(id, element, Eigen::Vector3d::Zero()){};

  StaticSite(data& d) { ReadData(d); }

  StaticSite(const QMAtom& atom, double charge)
      : StaticSite(atom.getId(), atom.getElement(), atom.getPos()) {
    setCharge(charge);
  }
  virtual ~StaticSite(){};

 protected:
  StaticSite(){};

 public:
  int getId() const { return _id; }
  int getRank() const { return _rank; }
  const std::string& getElement() const { return _element; }
  const Eigen::Vector3d& getPos() const { return _pos; }

  void setMultipole(const Vector9d& multipole, int rank) {
    _multipole = multipole;
    _rank = rank;
  }

  void setCharge(double q) { _multipole(0) = q; }

  void setPos(const Eigen::Vector3d& position) { _pos = position; }

  // COORDINATES TRANSFORMATION
  void Translate(const Eigen::VectorXd& shift);
  virtual void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos);

  // MULTIPOLES DEFINITION

  double getCharge() const { return _multipole(0); }
  const Vector9d& getPermMultipole() const {
    return _multipole;
  }  // Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s,...[NOT following Stone
     // order for dipoles]

  virtual Eigen::Vector3d getDipole() const { return _multipole.segment<3>(1); }

  Eigen::Matrix3d CalculateCartesianMultipole() const;
  static Eigen::VectorXd CalculateSphericalMultipole(
      const Eigen::Matrix3d& quadrupole_cartesian);

  const Eigen::Vector3d& getField() const { return _localpermanentField; }
  Eigen::Vector3d& getField() { return _localpermanentField; }

  const double& getPotential() const { return _phi; }
  double& getPotential() { return _phi; }

  std::string WriteMpsLine(std::string unit = "bohr") const;

  virtual void SetupCptTable(CptTable& table) const;

  void WriteData(data& d) const;
  void ReadData(data& d);
  virtual void setPolarisation(const Eigen::Matrix3d pol) { return; }

  virtual std::string identify() const { return "staticsite"; }

  friend std::ostream& operator<<(std::ostream& out, const StaticSite& site) {
    out << site.getId() << " " << site.getElement() << " " << site.getRank();
    out << " " << site.getPos().x() << "," << site.getPos().y() << ","
        << site.getPos().z() << "\n";
    return out;
  }

 protected:
  virtual std::string writePolarisation() const { return "P 0.0 0.0 0.0"; };

  int _id = -1;
  std::string _element = "";
  Eigen::Vector3d _pos = Eigen::Vector3d::Zero();
  int _rank = 0;

  Vector9d _multipole =
      Vector9d::Zero();  // Q00,Q11c,Q11s,Q10,Q20, Q21c, Q21s, Q22c, Q22s

  Eigen::Vector3d _localpermanentField = Eigen::Vector3d::Zero();
  double _phi = 0.0;  // Electric potential (due to perm.)
};
}  // namespace xtp
}  // namespace votca

#endif
