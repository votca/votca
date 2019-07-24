/*
 *           Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/format.hpp>
#include <fstream>
#include <string>
#include <votca/tools/constants.h>
#include <votca/xtp/staticsite.h>

using namespace std;

namespace votca {
namespace xtp {

Eigen::Matrix3d StaticSite::CalculateCartesianMultipole() const {
  // We are transforming here just quadrupoles
  // const  Eigen::VectorXd& MP = _multipole;
  const Vector9d& MP = _Q;
  Eigen::Matrix3d theta = Eigen::Matrix3d::Zero();
  if (_rank > 1) {
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
  Eigen::Vector3d dir = _pos - refPos;
  dir = R * dir;
  _pos = refPos + dir;  // Rotated Position
  if (_rank > 0) {
    const Eigen::Vector3d temp = R * _Q.segment<3>(1);
    _Q.segment<3>(1) = temp;
  }
  if (_rank > 1) {
    Eigen::Matrix3d cartesianquad = CalculateCartesianMultipole();
    Eigen::Matrix3d rotated = R * cartesianquad * R.transpose();
    _Q.segment<5>(4) = CalculateSphericalMultipole(rotated);
  }
  return;
}

void StaticSite::Translate(const Eigen::VectorXd& shift) {
  _pos += shift;
  return;
}

std::string StaticSite::WriteMpsLine(string unit) const {
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
             _element % (_pos(0) * conv_pos) % (_pos(1) * conv_pos) %
             (_pos(2) * conv_pos) % _rank)
                .str();
  output += (boost::format("    %1$+1.7f\n") % getCharge()).str();
  if (_rank > 0) {
    // Dipole z x y
    output += (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f\n") % _Q(3) %
               _Q(1) % _Q(2))
                  .str();
    if (_rank > 1) {
      // Quadrupole 20 21c 21s 22c 22s
      output +=
          (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f\n") %
           _Q(4) % _Q(5) % _Q(6) % _Q(7) % _Q(8))
              .str();
    }
  }
  // Polarizability
  output += writePolarisation();
  return output;
}

void StaticSite::SetupCptTable(CptTable& table) const {
  table.addCol(_id, "index", HOFFSET(data, id));
  table.addCol(_element, "type", HOFFSET(data, element));

  table.addCol(_pos[0], "posX", HOFFSET(data, posX));
  table.addCol(_pos[1], "posY", HOFFSET(data, posY));
  table.addCol(_pos[2], "posZ", HOFFSET(data, posZ));

  table.addCol(_rank, "rank", HOFFSET(data, rank));

  table.addCol(_Q[0], "Q00", HOFFSET(data, Q00));
  table.addCol(_Q[1], "Q11c", HOFFSET(data, Q11c));
  table.addCol(_Q[2], "Q11s", HOFFSET(data, Q11s));
  table.addCol(_Q[3], "Q10", HOFFSET(data, Q10));
  table.addCol(_Q[4], "Q20", HOFFSET(data, Q20));
  table.addCol(_Q[5], "Q21c", HOFFSET(data, Q21c));
  table.addCol(_Q[6], "Q21s", HOFFSET(data, Q21s));
  table.addCol(_Q[7], "Q22c", HOFFSET(data, Q22c));
  table.addCol(_Q[8], "Q22s", HOFFSET(data, Q22s));
}

void StaticSite::WriteData(data& d) const {
  d.id = _id;
  d.element = const_cast<char*>(_element.c_str());
  d.posX = _pos[0];
  d.posY = _pos[1];
  d.posZ = _pos[2];

  d.rank = _rank;

  d.Q00 = _Q[0];
  d.Q11c = _Q[1];
  d.Q11s = _Q[2];
  d.Q10 = _Q[3];
  d.Q20 = _Q[4];
  d.Q21c = _Q[5];
  d.Q21s = _Q[6];
  d.Q22c = _Q[7];
  d.Q22s = _Q[8];
}

void StaticSite::ReadData(data& d) {
  _id = d.id;
  _element = std::string(d.element);
  free(d.element);
  _pos[0] = d.posX;
  _pos[1] = d.posY;
  _pos[2] = d.posZ;

  _rank = d.rank;

  _Q[0] = d.Q00;
  _Q[1] = d.Q11c;
  _Q[2] = d.Q11s;
  _Q[3] = d.Q10;
  _Q[4] = d.Q20;
  _Q[5] = d.Q21c;
  _Q[6] = d.Q21s;
  _Q[7] = d.Q22c;
  _Q[8] = d.Q22s;
}

}  // namespace xtp
}  // namespace votca
