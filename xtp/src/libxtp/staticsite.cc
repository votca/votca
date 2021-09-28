/*
 *           Copyright 2009-2020 The VOTCA Development Team
 *                      (http://www.votca.org)
 *
 *     Licensed under the Apache License,Version 2.0 (the "License")
 *
 *You may not use this file except in compliance with the License.
 *You may obtain a copy of the License at
 *
 *             http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing,software
 *distributed under the License is distributed on an "AS IS" BASIS,
 *WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,either express or implied.
 *See the License for the specific language governing permissions and
 *limitations under the License.
 *
 */

// Standard includes
#include <fstream>
#include <string>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/checkpointtable.h"
#include "votca/xtp/staticsite.h"

namespace votca {
namespace xtp {

Eigen::Matrix3d StaticSite::CalculateCartesianMultipole() const {
  // We are transforming here just quadrupoles
  // const  Eigen::VectorXd& MP =  multipole_;
  const Vector9d& MP = Q_;
  Eigen::Matrix3d theta = Eigen::Matrix3d::Zero();
  if (rank_ > 1) {
    double sqr3 = std::sqrt(3);
    theta(0, 0) = 0.5 * (-MP(4) + sqr3 * MP(7));     // theta_xx
    theta(1, 1) = 0.5 * (-MP(4) - sqr3 * MP(7));     // theta_yy
    theta(2, 2) = MP(4);                             // theta_zz
    theta(0, 1) = theta(1, 0) = 0.5 * sqr3 * MP(8);  // theta_xy = theta_yx
    theta(0, 2) = theta(2, 0) = 0.5 * sqr3 * MP(5);  // theta_xz = theta_zx
    theta(1, 2) = theta(2, 1) = 0.5 * sqr3 * MP(6);  // theta_yz = theta_zy
  }
  return theta;
}

Eigen::VectorXd StaticSite::CalculateSphericalMultipole(
    const Eigen::Matrix3d& quad_cart) {
  Eigen::VectorXd quadrupole_polar = Eigen::VectorXd::Zero(5);
  const double sqr3 = std::sqrt(3);
  quadrupole_polar(0) = quad_cart(2, 2);
  quadrupole_polar(1) = (2. / sqr3) * quad_cart(0, 2);
  quadrupole_polar(2) = (2. / sqr3) * quad_cart(1, 2);
  quadrupole_polar(3) = (1. / sqr3) * (quad_cart(0, 0) - quad_cart(1, 1));
  quadrupole_polar(4) = (2. / sqr3) * quad_cart(0, 1);
  return quadrupole_polar;
}

void StaticSite::Rotate(const Eigen::Matrix3d& R,
                        const Eigen::Vector3d& refPos) {
  Eigen::Vector3d dir = pos_ - refPos;
  dir = R * dir;
  pos_ = refPos + dir;  // Rotated Position
  if (rank_ > 0) {
    const Eigen::Vector3d temp = R * Q_.segment<3>(1);
    Q_.segment<3>(1) = temp;
  }
  if (rank_ > 1) {
    Eigen::Matrix3d cartesianquad = CalculateCartesianMultipole();
    Eigen::Matrix3d rotated = R * cartesianquad * R.transpose();
    Q_.segment<5>(4) = CalculateSphericalMultipole(rotated);
  }
  return;
}

void StaticSite::Translate(const Eigen::VectorXd& shift) {
  pos_ += shift;
  return;
}

std::string StaticSite::writepolarization() const {
  tools::Elements e;
  double default_pol = 1;  // default is alway 1A^3
  try {
    default_pol =
        e.getPolarizability(element_) * std::pow(tools::conv::nm2ang, 3);
  } catch (const std::runtime_error&) {
    ;
  }
  return (boost::format("     P %1$+1.7f\n") % default_pol).str();
}

std::string StaticSite::WriteMpsLine(std::string unit) const {
  double conv_pos = 1.;
  if (unit == "angstrom") {
    conv_pos = tools::conv::bohr2ang;
  } else if (unit == "bohr") {
    conv_pos = 1.;
  } else {
    throw std::runtime_error(
        " StaticSite::WriteMpsLine: Unit conversion not known");
  }
  std::string output = "";
  output += (boost::format(" %1$2s %2$+1.7f %3$+1.7f %4$+1.7f Rank %5$d\n") %
             element_ % (pos_(0) * conv_pos) % (pos_(1) * conv_pos) %
             (pos_(2) * conv_pos) % rank_)
                .str();
  output += (boost::format("    %1$+1.7f\n") % getCharge()).str();
  if (rank_ > 0) {
    // Dipole z x y
    output += (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f\n") % Q_(3) %
               Q_(1) % Q_(2))
                  .str();
    if (rank_ > 1) {
      // Quadrupole 20 21c 21s 22c 22s
      output +=
          (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f\n") %
           Q_(4) % Q_(5) % Q_(6) % Q_(7) % Q_(8))
              .str();
    }
  }
  // Polarizability
  output += writepolarization();
  return output;
}

void StaticSite::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, id));
  table.addCol<std::string>("type", HOFFSET(data, element));

  table.addCol<double>("posX", HOFFSET(data, posX));
  table.addCol<double>("posY", HOFFSET(data, posY));
  table.addCol<double>("posZ", HOFFSET(data, posZ));

  table.addCol<Index>("rank", HOFFSET(data, rank));

  table.addCol<double>("Q00", HOFFSET(data, Q00));
  table.addCol<double>("Q11c", HOFFSET(data, Q11c));
  table.addCol<double>("Q11s", HOFFSET(data, Q11s));
  table.addCol<double>("Q10", HOFFSET(data, Q10));
  table.addCol<double>("Q20", HOFFSET(data, Q20));
  table.addCol<double>("Q21c", HOFFSET(data, Q21c));
  table.addCol<double>("Q21s", HOFFSET(data, Q21s));
  table.addCol<double>("Q22c", HOFFSET(data, Q22c));
  table.addCol<double>("Q22s", HOFFSET(data, Q22s));
}

void StaticSite::WriteData(data& d) const {
  d.id = id_;
  d.element = const_cast<char*>(element_.c_str());
  d.posX = pos_[0];
  d.posY = pos_[1];
  d.posZ = pos_[2];

  d.rank = rank_;

  d.Q00 = Q_[0];
  d.Q11c = Q_[1];
  d.Q11s = Q_[2];
  d.Q10 = Q_[3];
  d.Q20 = Q_[4];
  d.Q21c = Q_[5];
  d.Q21s = Q_[6];
  d.Q22c = Q_[7];
  d.Q22s = Q_[8];
}

void StaticSite::ReadData(const data& d) {
  id_ = d.id;
  element_ = std::string(d.element);
  free(d.element);
  pos_[0] = d.posX;
  pos_[1] = d.posY;
  pos_[2] = d.posZ;

  rank_ = d.rank;

  Q_[0] = d.Q00;
  Q_[1] = d.Q11c;
  Q_[2] = d.Q11s;
  Q_[3] = d.Q10;
  Q_[4] = d.Q20;
  Q_[5] = d.Q21c;
  Q_[6] = d.Q21s;
  Q_[7] = d.Q22c;
  Q_[8] = d.Q22s;
}

}  // namespace xtp
}  // namespace votca
