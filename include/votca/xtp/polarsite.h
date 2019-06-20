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

  void setPolarisation(const Eigen::Matrix3d pol) override;
  void ResetInduction();

  const Eigen::Matrix3d& getPolarisation() const { return _Ps; }

  // MULTIPOLES DEFINITION
  Eigen::Vector3d getDipole() const override;

  const Vector9d& V_ind() const { return _V_ind; }
  Vector9d& V_ind() { return _V_ind; }

  double getEigenDamp() const { return _eigendamp; }

  void Rotate(const Eigen::Matrix3d& R,
              const Eigen::Vector3d& ref_pos) override {
    StaticSite::Rotate(R, ref_pos);
    _Ps = R * _Ps * R.transpose();
  }
  void calcDIIS_InducedDipole();

  const Eigen::Vector3d& getInduced_Dipole() const { return _induced_dipole; }

  double Energy() const { return FieldEnergy() + InternalEnergy(); }

  double FieldEnergy() const;

  double InternalEnergy() const;

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

    double V00;
    double V11c;
    double V11s;
    double V10;
    double V20;
    double V21c;
    double V21s;
    double V22c;
    double V22s;

    double pxx;
    double pxy;
    double pxz;
    double pyy;
    double pyz;
    double pzz;

    double V00_ind;
    double V11c_ind;
    double V11s_ind;
    double V10_ind;
    double V20_ind;
    double V21c_ind;
    double V21s_ind;
    double V22c_ind;
    double V22s_ind;
  };
  // do not move up has to be below data definition
  PolarSite(data& d);

  void Reset() override {
    StaticSite::Reset();
    ResetInduction();
    _dipole_hist.clear();
  }

  double DipoleChange() const;

  void SetupCptTable(CptTable& table) const override;
  void WriteData(data& d) const;
  void ReadData(data& d);

  std::string identify() const override { return "polarsite"; }

  friend std::ostream& operator<<(std::ostream& out, const PolarSite& site) {
    out << site.getId() << " " << site.getElement() << " " << site.getRank();
    out << " " << site.getPos().x() << "," << site.getPos().y() << ","
        << site.getPos().z() << "\n";
    return out;
  }

 private:
  std::string writePolarisation() const override;

  Eigen::Matrix3d _Ps = Eigen::Matrix3d::Zero();
  Vector9d _V_ind = Vector9d::Zero();

  // cached data
  Eigen::Vector3d _induced_dipole = Eigen::Vector3d::Zero();
  std::vector<Eigen::Vector3d> _dipole_hist;
  double _eigendamp = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif
