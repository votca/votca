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
#ifndef VOTCA_XTP_EWALDSITE_H
#define VOTCA_XTP_EWALDSITE_H

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/logger.h"

#include <iomanip>

namespace votca {
namespace xtp {

class EwdSite {
 public:
  EwdSite(const PolarSite& pol);

  ~EwdSite() = default;

  friend std::ostream& operator<<(std::ostream& out, const EwdSite& site) {
    out << site.getId() << std::fixed << std::setprecision(4) << " "
        << site.getPos().transpose() << " " << site.getCharge() << std::endl
        << " " << site.getPolarizationMatrix() << std::endl
        << std::scientific << std::setprecision(4) << " "
        << site.getStaticField().transpose();
    return out;
  }

  bool operator==(const EwdSite& other)  const;

  bool operator!=(const EwdSite& other) { return !operator==(other); }

  Index getId() const { return id_; }

  const Eigen::Vector3d& getPos() const { return position_; }

  void updatePos(Eigen::Vector3d pos) { position_ = pos; }

  const Eigen::Vector3d& getStaticDipole() const { return dipole_static_; }

  const Eigen::Vector3d getTotalDipole() const {
    return dipole_static_ + dipole_induced_;
  }

  const Eigen::Matrix3d& getPolarizationMatrix() const { return polarization_; }

  const Eigen::Vector3d getInducedDipole() const { return dipole_induced_; }

  void setInducedDipole(Eigen::Vector3d dip) { dipole_induced_ = dip; }

  double getCharge() const { return charge_; }

  const Eigen::Matrix3d& getQuadrupole() const { return quadrupole_; }

  Index getRank() const { return rank_; }

  std::string getElement() const { return element_; }

  const Eigen::Vector3d& getStaticField() const { return field_static_; }

  void addToStaticField(Eigen::Vector3d field) { field_static_ += field; }

  void addToInducedField(Eigen::Vector3d field) { field_induced_ += field; }

  void induceDirect();

  void resetInductionField() { field_induced_ = Eigen::Vector3d::Zero(); }

  struct data {
    Index id;
    double posX;
    double posY;
    double posZ;

    Index rank;

    double charge;
    double dipX;
    double dipY;
    double dipZ;
    double quadXX;
    double quadXY;
    double quadXZ;
    double quadYY;
    double quadYZ;
    double quadZZ;

    double polXX;
    double polXY;
    double polXZ;
    double polYY;
    double polYZ;
    double polZZ;

    double d_x_ind;
    double d_y_ind;
    double d_z_ind;

    double fieldX_stat;
    double fieldY_stat;
    double fieldZ_stat;
    double fieldX_ind;
    double fieldY_ind;
    double fieldZ_ind;
    char* element;
  };

  EwdSite(const data& d) { ReadData(d); }

  void WriteData(data& d);
  void ReadData(const data& d);
  static void SetupCptTable(CptTable& table);

 private:
  Index id_;
  Index rank_;
  Eigen::Vector3d position_;
  double charge_;
  Eigen::Vector3d dipole_static_;
  Eigen::Vector3d dipole_induced_ = Eigen::Vector3d::Zero();
  Eigen::Matrix3d quadrupole_;
  Eigen::Vector3d field_static_ = Eigen::Vector3d::Zero();
  Eigen::Vector3d field_induced_ = Eigen::Vector3d::Zero();
  Eigen::Matrix3d polarization_;
  std::string element_;
};
}  // namespace xtp
}  // namespace votca
#endif