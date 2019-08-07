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
#ifndef __VOTCA_XTP_POLARSITE_H
#define __VOTCA_XTP_POLARSITE_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/staticsite.h>

namespace votca {
namespace xtp {

/**
\brief Class to represent Atom/Site in electrostatic+polarisation

 The units are atomic units, e.g. Bohr, Hartree.
*/
class PolarSite : public StaticSite {

 public:
  // delete these two functions because we do not want to be able to read
  // StaticSite::data but PolarSite::data
  void WriteData(StaticSite::data& d) const = delete;
  void ReadData(StaticSite::data& d) = delete;

  PolarSite(int id, std::string element, Eigen::Vector3d pos);
  PolarSite(int id, std::string element)
      : PolarSite(id, element, Eigen::Vector3d::Zero()){};

  ~PolarSite(){};

  void setPolarisation(const Eigen::Matrix3d& pol) override;

  Eigen::Matrix3d getPolarisation() const { return _pinv.inverse(); }

  const Eigen::Matrix3d& getPInv() const { return _pinv; }

  // MULTIPOLES DEFINITION
  Eigen::Vector3d getDipole() const override;

  double getSqrtInvEigenDamp() const { return _eigendamp_invsqrt; }

  void Rotate(const Eigen::Matrix3d& R,
              const Eigen::Vector3d& ref_pos) override {
    StaticSite::Rotate(R, ref_pos);
    _pinv = R.transpose() * _pinv * R;
  }

  const Eigen::Vector3d& V() const { return _V; }

  Eigen::Vector3d& V() { return _V; }

  const Eigen::Vector3d& V_noE() const { return _V_noE; }

  Eigen::Vector3d& V_noE() { return _V_noE; }

  void Reset() {
    _V.setZero();
    _V_noE.setZero();
  }

  double deltaQ_V_ext() const { return _induced_dipole.dot(_V); }

  double InternalEnergy() const {
    return 0.5 * _induced_dipole.transpose() * _pinv * _induced_dipole;
  }

  const Eigen::Vector3d& Induced_Dipole() const { return _induced_dipole; }
  void setInduced_Dipole(const Eigen::Vector3d& induced_dipole) {
    _induced_dipole = induced_dipole;
  }

  struct data {
    int id;
    char* element;
    double posX;
    double posY;
    double posZ;

    int rank;

    double Q00;
    double Q11c;
    double Q11s;
    double Q10;
    double Q20;
    double Q21c;
    double Q21s;
    double Q22c;
    double Q22s;

    double Vx;
    double Vy;
    double Vz;

    double Vx_noE;
    double Vy_noE;
    double Vz_noE;

    double pxx;
    double pxy;
    double pxz;
    double pyy;
    double pyz;
    double pzz;

    double d_x_ind;
    double d_y_ind;
    double d_z_ind;
  };
  // do not move up has to be below data definition
  PolarSite(data& d);

  double DipoleChange() const;

  void SetupCptTable(CptTable& table) const override;
  void WriteData(data& d) const;
  void ReadData(data& d);

  std::string identify() const override { return "polarsite"; }

  friend std::ostream& operator<<(std::ostream& out, const PolarSite& site) {
    out << site.getId() << " " << site.getElement() << " " << site.getRank();
    out << " " << site.getPos().transpose() << " "
        << site.Induced_Dipole().transpose() << "\n";
    return out;
  }

 private:
  std::string writePolarisation() const override;

  // PolarSite has two external fields,
  // the first is used for interaction with regions, which are further out, i.e.
  // the interaction energy with it is included in the polar region energy
  Eigen::Vector3d _V = Eigen::Vector3d::Zero();
  // the second is used for interaction with regions, which are further inside,
  // i.e. the interaction energy with it is included in the other region's
  // energy
  Eigen::Vector3d _V_noE = Eigen::Vector3d::Zero();

  Eigen::Vector3d _induced_dipole = Eigen::Vector3d::Zero();
  Eigen::Matrix3d _pinv = Eigen::Matrix3d::Zero();
  double _eigendamp_invsqrt = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif
