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
#ifndef VOTCA_XTP_POLARSITE_H
#define VOTCA_XTP_POLARSITE_H

// Local VOTCA includes
#include "eigen.h"
#include "staticsite.h"

namespace votca {
namespace xtp {

/**
\brief Class to represent Atom/Site in electrostatic+polarization

 The units are atomic units, e.g. Bohr, Hartree.
*/
class PolarSite final : public StaticSite {

 public:
  // delete these two functions because we do not want to be able to read
  // StaticSite::data but PolarSite::data
  void WriteData(StaticSite::data& d) const = delete;
  void ReadData(StaticSite::data& d) = delete;

  PolarSite(Index id, std::string element, Eigen::Vector3d pos);
  PolarSite(Index id, std::string element)
      : PolarSite(id, element, Eigen::Vector3d::Zero()) {};

  ~PolarSite() final = default;

  void setpolarization(const Eigen::Matrix3d& pol) final;

  Eigen::Matrix3d getpolarization() const { return pinv_.inverse(); }

  const Eigen::Matrix3d& getPInv() const { return pinv_; }

  // MULTIPOLES DEFINITION
  Eigen::Vector3d getDipole() const final;

  double getSqrtInvEigenDamp() const { return eigendamp_invsqrt_; }

  void Rotate(const Eigen::Matrix3d& R, const Eigen::Vector3d& ref_pos) final {
    StaticSite::Rotate(R, ref_pos);
    pinv_ = R.transpose() * pinv_ * R;
  }

  const Eigen::Vector3d& V() const { return V_; }

  Eigen::Vector3d& V() { return V_; }

  const Eigen::Vector3d& V_noE() const { return V_noE_; }

  Eigen::Vector3d& V_noE() { return V_noE_; }

  void Reset() {
    V_.setZero();
    V_noE_.setZero();
  }

  double deltaQ_V_ext() const { return induced_dipole_.dot(V_); }

  double InternalEnergy() const {
    return 0.5 * induced_dipole_.transpose() * pinv_ * induced_dipole_;
  }

  const Eigen::Vector3d& Induced_Dipole() const { return induced_dipole_; }
  void setInduced_Dipole(const Eigen::Vector3d& induced_dipole) {
    induced_dipole_ = induced_dipole;
  }

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
  PolarSite(const data& d);

  double DipoleChange() const;

  static void SetupCptTable(CptTable& table);
  void WriteData(data& d) const;
  void ReadData(const data& d);

  std::string identify() const final { return "polarsite"; }

  friend std::ostream& operator<<(std::ostream& out, const PolarSite& site) {
    out << site.getId() << " " << site.getElement() << " " << site.getRank();
    out << " " << site.getPos().transpose() << " "
        << site.Induced_Dipole().transpose() << "\n";
    return out;
  }

 private:
  std::string writepolarization() const final;

  // PolarSite has two external fields,
  // the first is used for interaction with regions, which are further out, i.e.
  // the interaction energy with it is included in the polar region energy
  Eigen::Vector3d V_ = Eigen::Vector3d::Zero();
  // the second is used for interaction with regions, which are further inside,
  // i.e. the interaction energy with it is included in the other region's
  // energy
  Eigen::Vector3d V_noE_ = Eigen::Vector3d::Zero();

  Eigen::Vector3d induced_dipole_ = Eigen::Vector3d::Zero();
  Eigen::Matrix3d pinv_ = Eigen::Matrix3d::Zero();
  double eigendamp_invsqrt_ = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_POLARSITE_H
