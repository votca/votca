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
  // spherical_multipoles Q = ( Q00,Q10,Q11c,Q11s,Q20,Q21c,Q21s,Q22c,Q22s )
  // We are trasforming here just quadrupoles
  // const  Eigen::VectorXd& MP = _multipole;
  const Vector9d& MP = _multipole;
  Eigen::Matrix3d theta = Eigen::Matrix3d::Zero();
  if (_rank > 1) {
    double sqr3 = std::sqrt(3);
    theta(0, 0) = 0.5 * (-MP(4) + sqr3 * MP(7));     // theta_xx
    theta(1, 1) = 0.5 * (-MP(4) + sqr3 * (-MP(7)));  // theta_yy
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
    _multipole.segment<3>(1) = R * _multipole.segment<3>(1);
  }
  if (_rank > 1) {
    Eigen::Matrix3d cartesianquad = CalculateCartesianMultipole();
    Eigen::Matrix3d rotated = R * cartesianquad * R.transpose();
    _multipole.segment<5>(4) = CalculateSphericalMultipole(rotated);
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
    output += (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f\n") %
               _multipole(1) % _multipole(2) % _multipole(3))
                  .str();
    if (_rank > 1) {
      // Quadrupole 20 21c 21s 22c 22s
      output +=
          (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f\n") %
           _multipole(4) % _multipole(5) % _multipole(6) % _multipole(7) %
           _multipole(8))
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

  table.addCol(_multipole[0], "multipoleQ00", HOFFSET(data, multipoleQ00));
  table.addCol(_multipole[1], "multipoleQ11c", HOFFSET(data, multipoleQ11c));
  table.addCol(_multipole[2], "multipoleQ11s", HOFFSET(data, multipoleQ11s));
  table.addCol(_multipole[3], "multipoleQ10", HOFFSET(data, multipoleQ10));
  table.addCol(_multipole[4], "multipoleQ20", HOFFSET(data, multipoleQ20));
  table.addCol(_multipole[5], "multipoleQ21c", HOFFSET(data, multipoleQ21c));
  table.addCol(_multipole[6], "multipoleQ21s", HOFFSET(data, multipoleQ21s));
  table.addCol(_multipole[7], "multipoleQ22c", HOFFSET(data, multipoleQ22c));
  table.addCol(_multipole[8], "multipoleQ22s", HOFFSET(data, multipoleQ22s));
}

void StaticSite::WriteToCpt(CptTable& table, const std::size_t& idx) const {

  data d;

  d.id = _id;
  d.element = const_cast<char*>(_element.c_str());
  d.posX = _pos[0];
  d.posY = _pos[1];
  d.posZ = _pos[2];

  d.rank = _rank;

  d.multipoleQ00 = _multipole[0];
  d.multipoleQ11c = _multipole[1];
  d.multipoleQ11s = _multipole[2];
  d.multipoleQ10 = _multipole[3];
  d.multipoleQ20 = _multipole[4];
  d.multipoleQ21c = _multipole[5];
  d.multipoleQ21s = _multipole[6];
  d.multipoleQ22c = _multipole[7];
  d.multipoleQ22s = _multipole[8];

  table.writeToRow(&d, idx);
}

void StaticSite::WriteData(data& d) const {
  d.id = _id;
  d.element = const_cast<char*>(_element.c_str());
  d.posX = _pos[0];
  d.posY = _pos[1];
  d.posZ = _pos[2];

  d.rank = _rank;

  d.multipoleQ00 = _multipole[0];
  d.multipoleQ11c = _multipole[1];
  d.multipoleQ11s = _multipole[2];
  d.multipoleQ10 = _multipole[3];
  d.multipoleQ20 = _multipole[4];
  d.multipoleQ21c = _multipole[5];
  d.multipoleQ21s = _multipole[6];
  d.multipoleQ22c = _multipole[7];
  d.multipoleQ22s = _multipole[8];
}

void StaticSite::ReadFromCpt(CptTable& table, const std::size_t& idx) {
  data d;
  table.readFromRow(&d, idx);

  _id = d.id;
  _element = std::string(d.element);
  free(d.element);
  _pos[0] = d.posX;
  _pos[1] = d.posY;
  _pos[2] = d.posZ;

  _rank = d.rank;

  _multipole[0] = d.multipoleQ00;
  _multipole[1] = d.multipoleQ11c;
  _multipole[2] = d.multipoleQ11s;
  _multipole[3] = d.multipoleQ10;
  _multipole[4] = d.multipoleQ20;
  _multipole[5] = d.multipoleQ21c;
  _multipole[6] = d.multipoleQ21s;
  _multipole[7] = d.multipoleQ22c;
  _multipole[8] = d.multipoleQ22s;
}

void StaticSite::ReadData(data& d) {
  _id = d.id;
  _element = std::string(d.element);
  free(d.element);
  _pos[0] = d.posX;
  _pos[1] = d.posY;
  _pos[2] = d.posZ;

  _rank = d.rank;

  _multipole[0] = d.multipoleQ00;
  _multipole[1] = d.multipoleQ11c;
  _multipole[2] = d.multipoleQ11s;
  _multipole[3] = d.multipoleQ10;
  _multipole[4] = d.multipoleQ20;
  _multipole[5] = d.multipoleQ21c;
  _multipole[6] = d.multipoleQ21s;
  _multipole[7] = d.multipoleQ22c;
  _multipole[8] = d.multipoleQ22s;
}

}  // namespace xtp
}  // namespace votca
