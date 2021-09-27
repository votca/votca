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
#include "votca/xtp/polarsite.h"

using namespace std;

namespace votca {
namespace xtp {

PolarSite::PolarSite(Index id, std::string element, Eigen::Vector3d pos)
    : StaticSite(id, element, pos) {
  tools::Elements e;
  double default_pol = std::pow(tools::conv::ang2bohr, 3);
  try {
    default_pol =
        e.getPolarizability(element) * std::pow(tools::conv::nm2bohr, 3);
  } catch (const std::runtime_error&) {
    ;
  }
  setpolarization(default_pol * Eigen::Matrix3d::Identity());
}

PolarSite::PolarSite(const data& d) { ReadData(d); }

Eigen::Vector3d PolarSite::getDipole() const {
  return Q_.segment<3>(1) + induced_dipole_;
}

void PolarSite::setpolarization(const Eigen::Matrix3d& pol) {
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(pol);
  pinv_ = es.eigenvectors() * es.eigenvalues().cwiseInverse().asDiagonal() *
          es.eigenvectors().transpose();
  eigendamp_invsqrt_ = 1.0 / std::sqrt(es.eigenvalues().maxCoeff());
}

std::string PolarSite::writepolarization() const {
  double conv_pol = std::pow(tools::conv::bohr2ang, 3);
  Eigen::MatrixX3d pol = pinv_.inverse() * conv_pol;
  return (boost::format("     P %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f "
                        "%6$+1.7f\n") %
          pol(0, 0) % pol(1, 0) % pol(2, 0) % pol(1, 1) % pol(1, 2) % pol(2, 2))
      .str();
}

void PolarSite::SetupCptTable(CptTable& table) {
  table.addCol<Index>("index", HOFFSET(data, id));
  table.addCol<std::string>("type", HOFFSET(data, element));

  table.addCol<double>("posX", HOFFSET(data, posX));
  table.addCol<double>("posY", HOFFSET(data, posY));
  table.addCol<double>("posZ", HOFFSET(data, posZ));
  // same type as rank
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

  table.addCol<double>("Vx", HOFFSET(data, Vx));
  table.addCol<double>("Vy", HOFFSET(data, Vy));
  table.addCol<double>("Vz", HOFFSET(data, Vz));

  table.addCol<double>("Vx_noE", HOFFSET(data, Vx_noE));
  table.addCol<double>("Vy_NoE", HOFFSET(data, Vy_noE));
  table.addCol<double>("Vz_noE", HOFFSET(data, Vz_noE));
  // P_inv and P have same polarization so no problem, as only data type counts
  table.addCol<double>("pxx", HOFFSET(data, pxx));
  table.addCol<double>("pxy", HOFFSET(data, pxy));
  table.addCol<double>("pxz", HOFFSET(data, pxz));
  table.addCol<double>("pyy", HOFFSET(data, pyy));
  table.addCol<double>("pyz", HOFFSET(data, pyz));
  table.addCol<double>("pzz", HOFFSET(data, pzz));

  table.addCol<double>("d_ind_x", HOFFSET(data, d_x_ind));
  table.addCol<double>("d_ind_y", HOFFSET(data, d_y_ind));
  table.addCol<double>("d_ind_z", HOFFSET(data, d_z_ind));
}

void PolarSite::WriteData(data& d) const {
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

  d.Vx = V_[0];
  d.Vy = V_[1];
  d.Vz = V_[2];

  d.Vx_noE = V_noE_[0];
  d.Vy_noE = V_noE_[1];
  d.Vz_noE = V_noE_[2];

  Eigen::Matrix3d P = pinv_.inverse();
  d.pxx = P(0, 0);
  d.pxy = P(0, 1);
  d.pxz = P(0, 2);
  d.pyy = P(1, 1);
  d.pyz = P(1, 2);
  d.pzz = P(2, 2);

  d.d_x_ind = induced_dipole_.x();
  d.d_y_ind = induced_dipole_.y();
  d.d_z_ind = induced_dipole_.z();
}

void PolarSite::ReadData(const data& d) {
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

  V_[0] = d.Vx;
  V_[1] = d.Vy;
  V_[2] = d.Vz;

  V_noE_[0] = d.Vx_noE;
  V_noE_[1] = d.Vy_noE;
  V_noE_[2] = d.Vz_noE;

  Eigen::Matrix3d Ps;
  Ps(0, 0) = d.pxx;
  Ps(0, 1) = d.pxy;
  Ps(1, 0) = d.pxy;
  Ps(0, 2) = d.pxz;
  Ps(2, 0) = d.pxz;
  Ps(1, 1) = d.pyy;
  Ps(1, 2) = d.pyz;
  Ps(2, 1) = d.pyz;
  Ps(2, 2) = d.pzz;

  this->setpolarization(Ps);

  induced_dipole_.x() = d.d_x_ind;
  induced_dipole_.y() = d.d_y_ind;
  induced_dipole_.z() = d.d_z_ind;
}

}  // namespace xtp
}  // namespace votca
