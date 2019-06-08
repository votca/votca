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

void PolarSite::calcDIIS_InducedDipole() {
  Eigen::Vector3d induced_dipole =
      _Ps * (_localpermanentField + _localinducedField);
  _dipole_hist.push_back(induced_dipole);

  int hist_size = _dipole_hist.size();
  if (hist_size == 1) {
    _induced_dipole = induced_dipole;
    return;
  }

  Eigen::MatrixXd B = -1 * Eigen::MatrixXd::Ones(hist_size + 1, hist_size + 1);
  B(hist_size, hist_size) = 0;
  Eigen::VectorXd a = Eigen::VectorXd::Zero(hist_size);
  a(hist_size) = -1;
  for (int i = 1; i < hist_size; i++) {
    const Eigen::Vector3d dmui = _dipole_hist[i] - _dipole_hist[i - 1];
    for (int j = 1; j <= i; j++) {
      const Eigen::Vector3d dmuj = _dipole_hist[j] - _dipole_hist[j - 1];
      B(j, i) = dmui.dot(dmuj);
      if (i != j) {
        B(i, j) = B(j, i);
      }
    }
  }
  Eigen::VectorXd coeffs = B.colPivHouseholderQr().solve(a);
  _induced_dipole = Eigen::Vector3d::Zero();
  for (int i = 0; i < hist_size; i++) {
    _induced_dipole += coeffs[i] * _dipole_hist[i];
  }
}

double PolarSite::DipoleChange() const {
  if (_dipole_hist.empty()) {
    return _induced_dipole.norm();
  }
  return (_induced_dipole - _dipole_hist.back()).norm();
}

Eigen::Vector3d PolarSite::getDipole() const {
  Eigen::Vector3d dipole = _multipole.segment<3>(1);
  dipole += getInduced_Dipole();
  return dipole;
}

void PolarSite::ResetInduction() {
  _phi_induced = 0.0;
  _localinducedField.setZero();
}

void PolarSite::setPolarisation(const Eigen::Matrix3d pol) {
  _Ps = pol;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(_Ps, Eigen::EigenvaluesOnly);
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

  table.addCol(_multipole[0], "multipoleQ00", HOFFSET(data, multipoleQ00));
  table.addCol(_multipole[1], "multipoleQ11c", HOFFSET(data, multipoleQ11c));
  table.addCol(_multipole[2], "multipoleQ11s", HOFFSET(data, multipoleQ11s));
  table.addCol(_multipole[3], "multipoleQ10", HOFFSET(data, multipoleQ10));
  table.addCol(_multipole[4], "multipoleQ20", HOFFSET(data, multipoleQ20));
  table.addCol(_multipole[5], "multipoleQ21c", HOFFSET(data, multipoleQ21c));
  table.addCol(_multipole[6], "multipoleQ21s", HOFFSET(data, multipoleQ21s));
  table.addCol(_multipole[7], "multipoleQ22c", HOFFSET(data, multipoleQ22c));
  table.addCol(_multipole[8], "multipoleQ22s", HOFFSET(data, multipoleQ22s));

  table.addCol(_localpermanentField[0], "localPermFieldX",
               HOFFSET(data, fieldX));
  table.addCol(_localpermanentField[1], "localPermFieldY",
               HOFFSET(data, fieldY));
  table.addCol(_localpermanentField[2], "localPermFieldZ",
               HOFFSET(data, fieldZ));
  table.addCol(_phi, "phi", HOFFSET(data, phi));

  table.addCol(_Ps(0, 0), "pxx", HOFFSET(data, pxx));
  table.addCol(_Ps(0, 1), "pxy", HOFFSET(data, pxy));
  table.addCol(_Ps(0, 2), "pxz", HOFFSET(data, pxz));
  table.addCol(_Ps(1, 1), "pyy", HOFFSET(data, pyy));
  table.addCol(_Ps(1, 2), "pyz", HOFFSET(data, pyz));
  table.addCol(_Ps(2, 2), "pzz", HOFFSET(data, pzz));

  table.addCol(_localinducedField[0], "localInducedFieldX",
               HOFFSET(data, fieldX_induced));
  table.addCol(_localinducedField[1], "localInducedFieldY",
               HOFFSET(data, fieldY_induced));
  table.addCol(_localinducedField[2], "localInducedFieldZ",
               HOFFSET(data, fieldZ_induced));

  table.addCol(_phi_induced, "phiInduced", HOFFSET(data, phi_induced));
}

void PolarSite::WriteData(data& d) const {
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

  d.fieldX = _localpermanentField[0];
  d.fieldY = _localpermanentField[1];
  d.fieldZ = _localpermanentField[2];
  d.phi = _phi;

  d.pxx = _Ps(0, 0);
  d.pxy = _Ps(0, 1);
  d.pxz = _Ps(0, 2);
  d.pyy = _Ps(1, 1);
  d.pyz = _Ps(1, 2);
  d.pzz = _Ps(2, 2);

  d.fieldX_induced = _localinducedField[0];
  d.fieldY_induced = _localinducedField[1];
  d.fieldZ_induced = _localinducedField[2];

  d.phi_induced = _phi_induced;
}

void PolarSite::ReadData(data& d) {
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

  _localpermanentField[0] = d.fieldX;
  _localpermanentField[1] = d.fieldY;
  _localpermanentField[2] = d.fieldZ;
  _phi = d.phi;

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

  _localinducedField[0] = d.fieldX_induced;
  _localinducedField[1] = d.fieldY_induced;
  _localinducedField[2] = d.fieldZ_induced;
  _phi_induced = d.phi_induced;
}

}  // namespace xtp
}  // namespace votca
