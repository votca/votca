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
#include <votca/xtp/polarsite.h>

using namespace std;

namespace votca {
namespace xtp {

PolarSite::PolarSite(int id, std::string element, Eigen::Vector3d pos)
    : StaticSite(id, element, pos) {
  tools::Elements e;
  double default_pol = std::pow(tools::conv::ang2bohr, 3);
  try {
    default_pol =
        e.getPolarizability(element) * std::pow(tools::conv::nm2bohr, 3);
  } catch (const std::invalid_argument&) {
    ;
  }
  setPolarisation(default_pol * Eigen::Matrix3d::Identity());
};

PolarSite::PolarSite(data& d) { ReadData(d); };

double PolarSite::FieldEnergy() const {
  double e = _V.dot(_Q);
  e += _V.segment<3>(1).dot(_induced_dipole);
  return e;
}

double PolarSite::InternalEnergy() const {
  return -0.5 * _induced_dipole.transpose() * _pinv *
         _induced_dipole;  // internal
}

Eigen::Vector3d PolarSite::getDipole() const {
  Eigen::Vector3d dipole = _Q.segment<3>(1);
  dipole += Induced_Dipole();
  return dipole;
}

void PolarSite::setPolarisation(const Eigen::Matrix3d pol) {
  _Ps = pol;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(_Ps);
  _pinv = es.eigenvectors() * es.eigenvalues().cwiseInverse().asDiagonal() *
          es.eigenvectors().transpose();
  _eigendamp = es.eigenvalues().maxCoeff();
}

std::string PolarSite::writePolarisation() const {
  double conv_pol = std::pow(tools::conv::bohr2ang, 3);
  Eigen::MatrixX3d pol = _Ps * conv_pol;
  return (boost::format("     P %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f "
                        "%6$+1.7f\n") %
          pol(0, 0) % pol(1, 0) % pol(2, 0) % pol(1, 1) % pol(1, 2) % pol(2, 2))
      .str();
}

void PolarSite::SetupCptTable(CptTable& table) const {
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

  table.addCol(_V[0], "V00", HOFFSET(data, V00));
  table.addCol(_V[1], "V11c", HOFFSET(data, V11c));
  table.addCol(_V[2], "V11s", HOFFSET(data, V11s));
  table.addCol(_V[3], "V10", HOFFSET(data, V10));
  table.addCol(_V[4], "V20", HOFFSET(data, V20));
  table.addCol(_V[5], "V21c", HOFFSET(data, V21c));
  table.addCol(_V[6], "V21s", HOFFSET(data, V21s));
  table.addCol(_V[7], "V22c", HOFFSET(data, V22c));
  table.addCol(_V[8], "V22s", HOFFSET(data, V22s));

  table.addCol(_Ps(0, 0), "pxx", HOFFSET(data, pxx));
  table.addCol(_Ps(0, 1), "pxy", HOFFSET(data, pxy));
  table.addCol(_Ps(0, 2), "pxz", HOFFSET(data, pxz));
  table.addCol(_Ps(1, 1), "pyy", HOFFSET(data, pyy));
  table.addCol(_Ps(1, 2), "pyz", HOFFSET(data, pyz));
  table.addCol(_Ps(2, 2), "pzz", HOFFSET(data, pzz));
}

void PolarSite::WriteData(data& d) const {
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

  d.V11c = _V[1];
  d.V11s = _V[2];
  d.V10 = _V[3];
  d.V20 = _V[4];
  d.V21c = _V[5];
  d.V21s = _V[6];
  d.V22c = _V[7];
  d.V22s = _V[8];

  d.pxx = _Ps(0, 0);
  d.pxy = _Ps(0, 1);
  d.pxz = _Ps(0, 2);
  d.pyy = _Ps(1, 1);
  d.pyz = _Ps(1, 2);
  d.pzz = _Ps(2, 2);
}

void PolarSite::ReadData(data& d) {
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

  _V[0] = d.V00;
  _V[1] = d.V11c;
  _V[2] = d.V11s;
  _V[3] = d.V10;
  _V[4] = d.V20;
  _V[5] = d.V21c;
  _V[6] = d.V21s;
  _V[7] = d.V22c;
  _V[8] = d.V22s;

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

  this->setPolarisation(Ps);
}

}  // namespace xtp
}  // namespace votca
