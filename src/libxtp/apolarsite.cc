/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#include <boost/format.hpp>
#include <boost/math/special_functions/round.hpp>
#include <fstream>
#include <string>
#include <votca/xtp/apolarsite.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace xtp {

APolarSite::APolarSite(APolarSite *templ, bool do_depolarize)
    : _id(templ->_id),
      _name(templ->_name),
      _isVirtual(templ->_isVirtual),
      _pos(templ->_pos),

      _locX(templ->_locX),
      _locY(templ->_locY),
      _locZ(templ->_locZ),

      _top(templ->_top),
      _seg(templ->_seg),
      _frag(templ->_frag),

      _Qs(templ->_Qs),
      _rank(templ->_rank),
      _Ps(templ->_Ps),
      Pxx(templ->Pxx),
      Pxy(templ->Pxy),
      Pxz(templ->Pxz),
      Pyy(templ->Pyy),
      Pyz(templ->Pyz),
      Pzz(templ->Pzz),
      pax(templ->pax),
      pay(templ->pay),
      paz(templ->paz),
      eigenpxx(templ->eigenpxx),
      eigenpyy(templ->eigenpyy),
      eigenpzz(templ->eigenpzz),
      eigendamp(templ->eigendamp),

      Q00(templ->Q00),
      Q1x(templ->Q1x),
      Q1y(templ->Q1y),
      Q1z(templ->Q1z),
      Q20(templ->Q20),
      Q21c(templ->Q21c),
      Q21s(templ->Q21s),
      Q22c(templ->Q22c),
      Q22s(templ->Q22s),
      Qxx(templ->Qxx),
      Qxy(templ->Qxy),
      Qxz(templ->Qxz),
      Qyy(templ->Qyy),
      Qyz(templ->Qyz),
      Qzz(templ->Qzz),

      U1x(templ->U1x),
      U1y(templ->U1y),
      U1z(templ->U1z),
      FPx(templ->FPx),
      FPy(templ->FPy),
      FPz(templ->FPz),
      FUx(templ->FUx),
      FUy(templ->FUy),
      FUz(templ->FUz),

      _resolution(templ->_resolution),
      PhiP(templ->PhiP),
      PhiU(templ->PhiU) {

  if (do_depolarize) this->Depolarize();
}

void APolarSite::ImportFrom(APolarSite *templ, string tag) {

  _pos = templ->getPos();
  _resolution = templ->getResolution();

  if (tag == "basic") {
    _Qs[0] = templ->getQs(-1);
    _Qs[1] = templ->getQs(0);
    _Qs[2] = templ->getQs(1);
    _rank = templ->_rank;
    _Ps = templ->_Ps;
  }

  if (tag == "full") {
    _locX = templ->_locX;
    _locY = templ->_locY;
    _locZ = templ->_locZ;
    _top = templ->_top;
    _seg = templ->_seg;
    _frag = templ->_frag;
    _rank = templ->_rank;
    _Qs = templ->_Qs;
    _Ps = templ->_Ps;
    _id = templ->_id;
    _name = templ->_name;
  }
}

void APolarSite::Rotate(const matrix &rot, const vec &refPos) {

  vec dir = _pos - refPos;
  dir = rot * dir;
  _pos = refPos + dir;

  _locX = rot.getCol(0);
  _locY = rot.getCol(1);
  _locZ = rot.getCol(2);

  // Rotate multipoles into global frame >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  for (int state = -1; state < 2; state++) {

    //   0    1    2    3    4    5    6    7    8    9    10   ...
    //   Q00  Q10  Q1c  Q1s  Q20  Q21c Q21s Q22c Q22s Q30  Q31c ...

    matrix R = rot;
    matrix R_T = matrix(R.getRow(0), R.getRow(1), R.getRow(2));

    // Transform polarizability tensor into global frame
    matrix P_Global = R * _Ps[state + 1] * R_T;
    _Ps[state + 1] = P_Global;

    // Any multipoles for this charge state available?
    if (_Qs[state + 1].size() < 1) {
      continue;
    }

    // Transform dipole moment into global frame
    if (_Qs[state + 1].size() > 1) {

      double Qz = _Qs[state + 1][1];
      double Qx = _Qs[state + 1][2];
      double Qy = _Qs[state + 1][3];

      vec d = vec(Qx, Qy, Qz);
      d = R * d;

      _Qs[state + 1][1] = d.getZ();
      _Qs[state + 1][2] = d.getX();
      _Qs[state + 1][3] = d.getY();
    }

    // Transform quadrupole moment into global frame
    if (_Qs[state + 1].size() > 4) {

      double Qzz = _Qs[state + 1][4];
      double Qxx = -0.5 * _Qs[state + 1][4] + 0.5 * sqrt(3) * _Qs[state + 1][7];
      double Qyy = -0.5 * _Qs[state + 1][4] - 0.5 * sqrt(3) * _Qs[state + 1][7];

      double Qxy = 0.5 * sqrt(3) * _Qs[state + 1][8];
      double Qxz = 0.5 * sqrt(3) * _Qs[state + 1][5];
      double Qyz = 0.5 * sqrt(3) * _Qs[state + 1][6];

      matrix Q =
          matrix(vec(Qxx, Qxy, Qxz), vec(Qxy, Qyy, Qyz), vec(Qxz, Qyz, Qzz));

      matrix Q_Global = R * Q * R_T;

      /* if (this->getId() == 1) {
          cout << endl;
          cout << "  " << Q_Global.get(0,0);
          cout << "  " << Q_Global.get(0,1);
          cout << "  " << Q_Global.get(0,2);
          cout << endl;
          cout << "  " << Q_Global.get(1,0);
          cout << "  " << Q_Global.get(1,1);
          cout << "  " << Q_Global.get(1,2);
          cout << endl;
          cout << "  " << Q_Global.get(2,0);
          cout << "  " << Q_Global.get(2,1);
          cout << "  " << Q_Global.get(2,2);
          cout << endl;
      }                                   */

      _Qs[state + 1][4] = Q_Global.get(2, 2);                // Q20
      _Qs[state + 1][5] = 2 / sqrt(3) * Q_Global.get(0, 2);  // Q21c
      _Qs[state + 1][6] = 2 / sqrt(3) * Q_Global.get(1, 2);  // Q21s
      _Qs[state + 1][7] =
          1 / sqrt(3) * (Q_Global.get(0, 0) - Q_Global.get(1, 1));  // Q22c
      _Qs[state + 1][8] = 2 / sqrt(3) * Q_Global.get(0, 1);         // Q22s
    }
  }
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

bool APolarSite::getIsActive(bool estatics_only) {
  bool is_active = false;
  bool is_charged = IsCharged();
  bool is_polarizable = IsPolarizable();

  if (estatics_only)
    is_active = is_charged;
  else
    is_active = is_charged || is_polarizable;

  return is_active;
}

bool APolarSite::IsCharged() {
  bool is_charged = false;
  // Tolerances
  double q_tol = 1e-9;  // [e]
  double d_tol = 1e-9;  // [enm]
  double Q_tol = 1e-9;  // [enm^2]
  // Magnitudes
  double q_mag = sqrt(Q00 * Q00);
  double d_mag = sqrt(Q1x * Q1x + Q1y * Q1y + Q1z * Q1z);
  double Q_mag =
      sqrt(Q20 * Q20 + Q22c * Q22c + Q22s * Q22s + Q21c * Q21c + Q21s * Q21s);
  // Compare
  if (q_mag > q_tol) is_charged = true;
  if (_rank > 0 && d_mag > d_tol) is_charged = true;
  if (_rank > 1 && Q_mag > Q_tol) is_charged = true;
  return is_charged;
}

bool APolarSite::IsPolarizable() {
  bool is_polarizable = false;
  // Tolerances
  double p_tol = 1e-9;  // [nm^3]
  // Compare
  if (getIsoP() > p_tol) is_polarizable = true;
  return is_polarizable;
}

void APolarSite::Translate(const vec &shift) {
  _pos += shift;
  return;
}

void APolarSite::Charge(int state) {
  int idx = state + 1;
  // Adjust polarizability to charge state
  Pxx = _Ps[idx].get(0, 0);
  Pxy = _Ps[idx].get(0, 1);
  Pxz = _Ps[idx].get(0, 2);
  Pyy = _Ps[idx].get(1, 1);
  Pyz = _Ps[idx].get(1, 2);
  Pzz = _Ps[idx].get(2, 2);
  // Calculate principal axes
  matrix polarity =
      matrix(vec(Pxx, Pxy, Pxz), vec(Pxy, Pyy, Pyz), vec(Pxz, Pyz, Pzz));
  matrix::eigensystem_t eigensystem_polarity;
  polarity.SolveEigensystem(eigensystem_polarity);
  pax = eigensystem_polarity.eigenvecs[0];
  pay = eigensystem_polarity.eigenvecs[1];
  paz = eigensystem_polarity.eigenvecs[2];
  eigenpxx = eigensystem_polarity.eigenvalues[0];
  eigenpyy = eigensystem_polarity.eigenvalues[1];
  eigenpzz = eigensystem_polarity.eigenvalues[2];
  eigendamp = (eigenpxx > eigenpyy)
                  ? ((eigenpxx > eigenpzz) ? eigenpxx : eigenpzz)
                  : ((eigenpyy > eigenpzz) ? eigenpyy : eigenpzz);
  // cout << endl << "eigensystem ...";
  // cout << endl << pax << " --- " << eigenpxx;
  // cout << endl << pay << " --- " << eigenpyy;
  // cout << endl << paz << " --- " << eigenpzz;

  // Adjust multipole moments to charge state
  Q00 = _Qs[idx][0];

  if (_rank > 0) {
    Q1z = _Qs[idx][1];  // |
    Q1x = _Qs[idx][2];  // |-> NOTE: order z - x - y
    Q1y = _Qs[idx][3];  // |
  } else {
    Q1z = Q1x = Q1y = 0.0;
  }
  if (_rank > 1) {
    // Spherical tensor
    Q20 = _Qs[idx][4];
    Q21c = _Qs[idx][5];
    Q21s = _Qs[idx][6];
    Q22c = _Qs[idx][7];
    Q22s = _Qs[idx][8];

    // Cartesian tensor * 1/3
    Qzz = Q20;
    Qxx = -0.5 * Q20 + 0.5 * sqrt(3) * Q22c;
    Qyy = -0.5 * Q20 - 0.5 * sqrt(3) * Q22c;
    Qxy = +0.5 * sqrt(3) * Q22s;
    Qxz = +0.5 * sqrt(3) * Q21c;
    Qyz = +0.5 * sqrt(3) * Q21s;

    Qzz *= 1. / 3.;
    Qxx *= 1. / 3.;
    Qyy *= 1. / 3.;
    Qxy *= 1. / 3.;
    Qxz *= 1. / 3.;
    Qyz *= 1. / 3.;
  } else {
    Q20 = Q21c = Q21s = Q22c = Q22s = 0.0;
    Qzz = Qxx = Qyy = Qxy = Qxz = Qyz = 0.0;
  }
  return;
}

void APolarSite::ChargeDelta(int state1, int state2) {
  int idx1 = state1 + 1;
  int idx2 = state2 + 1;
  // Adjust polarizability to charge state
  Pxx = 0.0;
  Pxy = 0.0;
  Pxz = 0.0;
  Pyy = 0.0;
  Pyz = 0.0;
  Pzz = 0.0;
  // Adjust multipole moments to charge state
  Q00 = _Qs[idx2][0] - _Qs[idx1][0];
  if (_rank > 0) {
    Q1z = _Qs[idx2][1] - _Qs[idx1][1];  // |
    Q1x = _Qs[idx2][2] - _Qs[idx1][2];  // |-> NOTE: order z - x - y
    Q1y = _Qs[idx2][3] - _Qs[idx1][3];  // |
  }
  if (_rank > 1) {
    Q20 = _Qs[idx2][4] - _Qs[idx1][4];
    Q21c = _Qs[idx2][5] - _Qs[idx1][5];
    Q21s = _Qs[idx2][6] - _Qs[idx1][6];
    Q22c = _Qs[idx2][7] - _Qs[idx1][7];
    Q22s = _Qs[idx2][8] - _Qs[idx1][8];
  }
  return;
}

double APolarSite::getProjP(vec &dir) {
  double alpha_proj = fabs(dir * pax) * fabs(dir * pax) * eigenpxx +
                      fabs(dir * pay) * fabs(dir * pay) * eigenpyy +
                      fabs(dir * paz) * fabs(dir * paz) * eigenpzz;

  // Mix ?
  // alpha_proj = 0.5* (1./3. * (Pxx+Pyy+Pzz)) + 0.5*alpha_proj;

  // assert("APS::getProjP THIS FEATURE SHOULD NOT BE USED" == "" && false);
  return alpha_proj;
}

void APolarSite::Induce(double wSOR) {
  // SUCCESSIVE OVERRELAXATION
  U1_Hist.push_back(vec(U1x, U1y, U1z));  // Remember all previous moments
  // U1_Hist[0] = vec(U1x,U1y,U1z);           // Remember previous moment
  U1x = (1 - wSOR) * U1x +
        wSOR * (-Pxx * (FPx + FUx) - Pxy * (FPy + FUy) - Pxz * (FPz + FUz));
  U1y = (1 - wSOR) * U1y +
        wSOR * (-Pxy * (FPx + FUx) - Pyy * (FPy + FUy) - Pyz * (FPz + FUz));
  U1z = (1 - wSOR) * U1z +
        wSOR * (-Pxz * (FPx + FUx) - Pyz * (FPy + FUy) - Pzz * (FPz + FUz));
  return;

  /*
  // ANDERSON + SUCCESSIVE OVERRELAXATION
  U1_Hist.push_back( vec(U1x,U1y,U1z) );

  // Back up U(N,i)
  U1_i.push_back(vec(U1x, U1y, U1z));
  // Compute U(N,o)
  U1x = - Pxx * (FPx + FUx) - Pxy * (FPy + FUy) - Pxz * (FPz + FUz);
  U1y = - Pxy * (FPx + FUx) - Pyy * (FPy + FUy) - Pyz * (FPz + FUz);
  U1z = - Pxz * (FPx + FUx) - Pyz * (FPy + FUy) - Pzz * (FPz + FUz);
  // Back up U(N,o)
  U1_o.push_back(vec(U1x, U1y, U1z));

  // Compile Anderson in/out values from U1 record
  vec U1_i_mix = vec(0,0,0);
  vec U1_o_mix = vec(0,0,0);
  if (U1_i.size() > 2 && U1_o.size() > 2) {
      // Compute Anderson mixing factor from minimal rms deviation
      vec D_N   = U1_o[U1_o.size()-1] - U1_i[U1_i.size()-1];
      vec D_N_1 = U1_o[U1_o.size()-2] - U1_i[U1_i.size()-2];
      vec D_N_N_1 = D_N - D_N_1;
      double wAND = - D_N_1*D_N_N_1/(D_N_N_1*D_N_N_1);
      // Mingle new and old via wAND ...
      U1_i_mix = wAND*U1_i[U1_i.size()-1] + (1-wAND)*U1_i[U1_i.size()-2];
      U1_o_mix = wAND*U1_o[U1_o.size()-1] + (1-wAND)*U1_o[U1_o.size()-2];
  }
  else {
      U1_i_mix = U1_i.back();
      U1_o_mix = U1_o.back();
  }

  // Apply SOR
  vec U1 = (1-wSOR)*U1_i_mix + wSOR*U1_o_mix;
  U1x = U1.getX(); U1y = U1.getY(); U1z = U1.getZ();
  return;
  */
}

void APolarSite::InduceDirect() {
  U1_Hist.push_back(vec(0., 0., 0.));
  U1x = -Pxx * FPx - Pxy * FPy - Pxz * FPz;
  U1y = -Pxy * FPx - Pyy * FPy - Pyz * FPz;
  U1z = -Pxz * FPx - Pyz * FPy - Pzz * FPz;
  return;
}

double APolarSite::HistdU() {
  vec U0 = U1_Hist.back();
  vec U1 = vec(U1x, U1y, U1z);
  vec dU = U1 - U0;

  double abs_U0 = votca::tools::abs(U0);
  double abs_U1 = votca::tools::abs(U1);
  double abs_dU = votca::tools::abs(dU);
  double small = 1e-10;  // in e*nm
  double abs_U = (abs_U1 > abs_U0 && small > abs_U0) ? abs_U1 : abs_U0;

  double dU_U = 1.;
  if (small > abs_U)
    dU_U = abs_U;  // i.e.: U1, U0 << 1 enm
  else
    dU_U = abs_dU / abs_U;
  return dU_U;
}

double APolarSite::HistdU2() {
  vec U0 = U1_Hist.back();
  vec U1 = vec(U1x, U1y, U1z);
  double dU2 = (U1 - U0) * (U1 - U0);
  return dU2;
}

void APolarSite::Depolarize() {
  // Zero out induced moments
  U1x = U1y = U1z = 0.0;
  U1_Hist.clear();
  // Zero out fields
  FPx = FPy = FPz = 0.0;
  FUx = FUy = FUz = 0.0;
  return;
}

void APolarSite::PrintInfo(std::ostream &out) {

  printf(
      "\n%2s %2d POS %+2.3f %+2.3f %+2.3f "
      "POL %+2.3f %+2.3f %+2.3f %+2.3f %+2.3f %+2.3f ",
      _name.c_str(), _id, _pos.getX(), _pos.getY(), _pos.getZ(), Pxx * 1000,
      Pxy * 1000, Pxz * 1000, Pyy * 1000, Pyz * 1000, Pzz * 1000);
}

void APolarSite::PrintTensorPDB(FILE *out, int state) {

  this->Charge(state);

  double scale = 30.;

  vec rxyz = _pos * scale;

  double rx = rxyz.getX();
  double ry = rxyz.getY();
  double rz = rxyz.getZ();

  vec end_pax = rxyz + 1.0 * this->pax * this->eigenpxx / this->eigendamp;
  vec end_pay = rxyz + 1.0 * this->pay * this->eigenpyy / this->eigendamp;
  vec end_paz = rxyz + 1.0 * this->paz * this->eigenpzz / this->eigendamp;

  fprintf(out, "%-5s %+4.7f %+4.7f %+4.7f\n", this->_name.c_str(), rx, ry, rz);
  fprintf(out, "X     %+4.7f %+4.7f %+4.7f\n", end_pax.getX(), end_pax.getY(),
          end_pax.getZ());
  fprintf(out, "Y     %+4.7f %+4.7f %+4.7f\n", end_pay.getX(), end_pay.getY(),
          end_pay.getZ());
  fprintf(out, "Z     %+4.7f %+4.7f %+4.7f\n", end_paz.getX(), end_paz.getY(),
          end_paz.getZ());
}

void APolarSite::WritePdbLine(FILE *out, const string &tag) {

  fprintf(out,
          "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
          "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%4.7f\n",
          _id % 100000,      // Atom serial number           %5d
          _name.c_str(),     // Atom name                    %4s
          " ",               // alternate location indicator.%1s
          tag.c_str(),       // Residue name.                %3s
          "A",               // Chain identifier             %1s
          _id % 10000,       // Residue sequence number      %4d
          " ",               // Insertion of residues.       %1s
          _pos.getX() * 10,  // X in Angstroms               %8.3f
          _pos.getY() * 10,  // Y in Angstroms               %8.3f
          _pos.getZ() * 10,  // Z in Angstroms               %8.3f
          1.0,               // Occupancy                    %6.2f
          0.0,               // Temperature factor           %6.2f
          " ",               // Segment identifier           %4s
          _name.c_str(),     // Element symbol               %2s
          " ",               // Charge on the atom.          %2s
          Q00);
}

void APolarSite::WriteXyzLine(FILE *out, vec &shift, string format) {
  /*
      double int2ext = 1.0;



      if (format == "gaussian") {
          int2ext = 10.;
      }
      else {
          int2ext = 10.;
      }
  */
  vec pos = _pos + shift;
  fprintf(out, "%-2s %+4.9f %+4.9f %+4.9f \n", _name.c_str(), pos.getX() * 10,
          pos.getY() * 10, pos.getZ() * 10);
}

void APolarSite::WriteChkLine(FILE *out, vec &shift, bool split_dpl,
                              string format, double spacing) {

  vec pos = _pos + shift;

  string unit = "";

  if (format == "xyz") {
    unit = "angstrom";
  } else if (format == "gaussian") {
    unit = "angstrom";
  }

  // Take care of unit conversion
  double int2ext = 0;

  if (unit == "nanometer") {
    int2ext = 1.;
  } else if (unit == "angstrom") {
    int2ext = 10.;
  } else if (unit == "bohr") {
    assert(false);
  }

  if (format == "xyz") {
    fprintf(out, "%2s ", _name.c_str());
  }

  // Print charge line
  fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n", pos.getX() * int2ext,
          pos.getY() * int2ext, pos.getZ() * int2ext, Q00);

  // Split dipole moment onto charges (if desired)
  if (split_dpl) {

    vec tot_dpl = vec(U1x, U1y, U1z);

    if (_rank > 0) {
      tot_dpl += vec(Q1x, Q1y, Q1z);
    }

    matrix::eigensystem_t EIGEN;

    if (_rank == 2) {
      // cout << endl
      //     << "WARNING: Quadrupoles are not split onto point charges."
      //     << endl;

      int state = 0;

      double Qzz = _Qs[state + 1][4];
      double Qxx = -0.5 * _Qs[state + 1][4] + 0.5 * sqrt(3) * _Qs[state + 1][7];
      double Qyy = -0.5 * _Qs[state + 1][4] - 0.5 * sqrt(3) * _Qs[state + 1][7];

      double Qxy = 0.5 * sqrt(3) * _Qs[state + 1][8];
      double Qxz = 0.5 * sqrt(3) * _Qs[state + 1][5];
      double Qyz = 0.5 * sqrt(3) * _Qs[state + 1][6];

      matrix Q =
          matrix(vec(Qxx, Qxy, Qxz), vec(Qxy, Qyy, Qyz), vec(Qxz, Qyz, Qzz));

      Q.SolveEigensystem(EIGEN);
    }

    double a = spacing;
    double mag_d = abs(tot_dpl);
    vec dir_d_0 = tot_dpl.normalize();
    vec dir_d = dir_d_0.normalize();
    vec A = pos + 0.5 * a * dir_d;
    vec B = pos - 0.5 * a * dir_d;
    double qA = mag_d / a;
    double qB = -qA;

    if (this->eigendamp == 0 || mag_d < 1e-9) {
      A = pos + 0.1 * a * vec(1, 0, 0);  // != pos since self-energy may diverge
      B = pos - 0.1 * a * vec(1, 0, 0);
      qA = 0;
      qB = 0;
    }

    if (format == "xyz") {
      fprintf(out, " A ");
    }
    fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n", A.getX() * int2ext,
            A.getY() * int2ext, A.getZ() * int2ext, qA);

    if (format == "xyz") {
      fprintf(out, " B ");
    }
    fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n", B.getX() * int2ext,
            B.getY() * int2ext, B.getZ() * int2ext, qB);

    if (format == "xyz" && _rank == 2) {
      vec D1 = pos + 0.5 * a * EIGEN.eigenvecs[0];
      vec D2 = pos - 0.5 * a * EIGEN.eigenvecs[0];
      vec E1 = pos + 0.5 * a * EIGEN.eigenvecs[1];
      vec E2 = pos - 0.5 * a * EIGEN.eigenvecs[1];
      vec F1 = pos + 0.5 * a * EIGEN.eigenvecs[2];
      vec F2 = pos - 0.5 * a * EIGEN.eigenvecs[2];
      fprintf(out, " D %+4.9f %+4.9f %+4.9f \n", D1.getX() * int2ext,
              D1.getY() * int2ext, D1.getZ() * int2ext);
      fprintf(out, " D %+4.9f %+4.9f %+4.9f \n", D2.getX() * int2ext,
              D2.getY() * int2ext, D2.getZ() * int2ext);
      fprintf(out, " E %+4.9f %+4.9f %+4.9f \n", E1.getX() * int2ext,
              E1.getY() * int2ext, E1.getZ() * int2ext);
      fprintf(out, " E %+4.9f %+4.9f %+4.9f \n", E2.getX() * int2ext,
              E2.getY() * int2ext, E2.getZ() * int2ext);
      fprintf(out, " F %+4.9f %+4.9f %+4.9f \n", F1.getX() * int2ext,
              F1.getY() * int2ext, F1.getZ() * int2ext);
      fprintf(out, " F %+4.9f %+4.9f %+4.9f \n", F2.getX() * int2ext,
              F2.getY() * int2ext, F2.getZ() * int2ext);
    }
  }
}

void APolarSite::WriteMpsLine(std::ostream &out, string unit = "angstrom") {

  // Set conversion factor for higher-rank moments (e*nm**k to e*a0**k)
  double conv_dpl = 1. / 0.0529189379;
  double conv_qdr = conv_dpl * conv_dpl;
  // Set conversion factor for polarizabilities (nm**3 to A**3)
  double conv_pol = 1000;
  // Set conversion factor for positions (nm to ??)
  double conv_pos = 1.;
  if (unit == "angstrom") {
    conv_pos = 10.;
  } else if (unit == "nanometer") {
    conv_pos = 1.;
  } else
    assert(false);  // Units error

  out << (boost::format(" %1$2s %2$+1.7f %3$+1.7f %4$+1.7f Rank %5$d\n") %
          _name % (_pos.getX() * conv_pos) % (_pos.getY() * conv_pos) %
          (_pos.getZ() * conv_pos) % _rank);
  // Charged
  out << (boost::format("    %1$+1.7f\n") % Q00);
  if (_rank > 0) {
    // Dipole z x y
    out << (boost::format("    %1$+1.7f %2$+1.7f %3$+1.7f\n") %
            (Q1z * conv_dpl) % (Q1x * conv_dpl) % (Q1y * conv_dpl));
    if (_rank > 1) {
      // Quadrupole 20 21c 21s 22c 22s
      out << (boost::format(
                  "    %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f\n") %
              (Q20 * conv_qdr) % (Q21c * conv_qdr) % (Q21s * conv_qdr) %
              (Q22c * conv_qdr) % (Q22s * conv_qdr));
    }
  }
  // Polarizability
  out << (boost::format("     P %1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f %5$+1.7f "
                        "%6$+1.7f \n") %
          (Pxx * conv_pol) % (Pxy * conv_pol) % (Pxz * conv_pol) %
          (Pyy * conv_pol) % (Pyz * conv_pol) % (Pzz * conv_pol));
}

void APolarSite::WriteXmlLine(std::ostream &out) {
  out << "<aps>" << endl;
  out << _id << endl;
  out << _name << endl;
  out << ((_isVirtual) ? "virtual" : "notvirtual") << endl;
  out << "resolution" << (int)_resolution << endl;
  out << _pos << endl;
  out << _locX << " ";
  out << _locY << " ";
  out << _locZ << " ";
  out << endl;
  for (int state = -1; state < 2; ++state) {
    out << "<state>" << endl;
    out << state << endl;
    for (unsigned i = 0; i < _Qs[state + 1].size(); ++i) {
      out << _Qs[state + 1][i] << " ";
    }
    out << endl;
    out << "</state>" << endl;
  }
  out << _rank;
  out << endl;

  out << "<pol>" << endl;
  for (int state = -1; state < 2; ++state) {
    out << "<state>" << endl;
    out << _Ps[state + 1].get(0, 0) << " " << _Ps[state + 1].get(0, 1) << " "
        << _Ps[state + 1].get(0, 2) << endl;
    out << _Ps[state + 1].get(1, 0) << " " << _Ps[state + 1].get(1, 1) << " "
        << _Ps[state + 1].get(1, 2) << endl;
    out << _Ps[state + 1].get(2, 0) << " " << _Ps[state + 1].get(2, 1) << " "
        << _Ps[state + 1].get(2, 2) << endl;
    out << "</state>" << endl;
  }
  out << "</pol>" << endl;

  out << "<config>" << endl;
  out << Pxx << " ";
  out << Pxy << " ";
  out << Pxz << " ";
  out << endl;
  out << Pyy << " ";
  out << Pyz << " ";
  out << endl;
  out << Pzz << " ";
  out << endl;

  out << pax << " ";
  out << eigenpxx;
  out << endl;
  out << pay << " ";
  out << eigenpyy;
  out << endl;
  out << paz << " ";
  out << eigenpzz;
  out << endl;

  out << eigendamp;
  out << endl;

  out << Q00 << " ";
  out << endl;
  out << Q1x << " ";
  out << Q1y << " ";
  out << Q1z << " ";
  out << endl;
  out << Q20 << " ";
  out << Q21c << " ";
  out << Q21s << " ";
  out << Q22c << " ";
  out << Q22s << " ";
  out << endl;
  out << Qxx << " ";
  out << Qxy << " ";
  out << Qxz << " ";
  out << Qyy << " ";
  out << Qyz << " ";
  out << Qzz << " ";
  out << endl;

  out << U1x << " ";
  out << U1y << " ";
  out << U1z << " ";
  out << endl;
  out << FPx << " ";
  out << FPy << " ";
  out << FPz << " ";
  out << endl;
  out << FUx << " ";
  out << FUy << " ";
  out << FUz << " ";
  out << endl;
  out << "</config>" << endl;
  out << "</aps>" << endl;
  return;
}

vector<APolarSite *> APS_FROM_MPS(string filename, int state,
                                  QMThread *thread) {

  int poleCount = 1;
  double Q0_total = 0.0;
  string units = "";
  bool useDefaultPs = true;

  vector<APolarSite *> poles;
  APolarSite *thisPole = NULL;

  vector<double> Qs;  // <- multipole moments
  matrix P1;          // <- dipole polarizability

  std::string line;
  std::ifstream intt;
  intt.open(filename.c_str());

  if (intt.is_open()) {
    while (intt.good()) {

      std::getline(intt, line);
      vector<string> split;
      Tokenizer toker(line, " \t");
      toker.ToVector(split);

      if (!split.size() || split[0] == "!" || split[0].substr(0, 1) == "!") {
        continue;
      }

      // ! Interesting information here, e.g.
      // ! DCV2T opt
      // ! SP        RB3LYP          6-311+G(d,p)
      // Units bohr
      //
      // C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
      //  -0.3853409355
      //  -0.0002321905   0.2401559510   0.6602334308
      //  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813 -0.50937399
      //  P 1.75

      // Units used
      if (split[0] == "Units") {
        units = split[1];
        if (units != "bohr" && units != "angstrom") {
          throw std::runtime_error("Unit " + units + " in file " + filename +
                                   " not supported.");
        }
      }
      // element,  position,  rank limit
      else if (split.size() == 6) {
        Qs.clear();
        int id = poleCount++;  // <- starts from 1
        string name = split[0];

        double BOHR2NM = 0.0529189379;
        double ANGSTROM2NM = 0.1;
        double x, y, z;

        if (units == "bohr") {
          x = BOHR2NM * boost::lexical_cast<double>(split[1]);
          y = BOHR2NM * boost::lexical_cast<double>(split[2]);
          z = BOHR2NM * boost::lexical_cast<double>(split[3]);
        } else if (units == "angstrom") {
          x = ANGSTROM2NM * boost::lexical_cast<double>(split[1]);
          y = ANGSTROM2NM * boost::lexical_cast<double>(split[2]);
          z = ANGSTROM2NM * boost::lexical_cast<double>(split[3]);
        } else {
          throw std::runtime_error("Unit " + units + " in file " + filename +
                                   " not supported.");
        }

        vec pos = vec(x, y, z);
        int rank = boost::lexical_cast<int>(split[5]);
        APolarSite *newPole = new APolarSite(id, name);
        newPole->setRank(rank);
        newPole->setPos(pos);
        poles.push_back(newPole);
        thisPole = newPole;
      }
      // 'P', dipole polarizability
      else if (split[0] == "P") {
        double pxx, pxy, pxz;
        double pyy, pyz;
        double pzz;
        if (split.size() == 7) {
          pxx = 1e-3 * boost::lexical_cast<double>(split[1]);
          pxy = 1e-3 * boost::lexical_cast<double>(split[2]);
          pxz = 1e-3 * boost::lexical_cast<double>(split[3]);
          pyy = 1e-3 * boost::lexical_cast<double>(split[4]);
          pyz = 1e-3 * boost::lexical_cast<double>(split[5]);
          pzz = 1e-3 * boost::lexical_cast<double>(split[6]);
          P1 = matrix(vec(pxx, pxy, pxz), vec(pxy, pyy, pyz),
                      vec(pxz, pyz, pzz));
        } else if (split.size() == 2) {
          pxx = 1e-3 * boost::lexical_cast<double>(split[1]);
          pxy = 0.0;
          pxz = 0.0;
          pyy = pxx;
          pyz = 0.0;
          pzz = pxx;
          P1 = matrix(vec(pxx, pxy, pxz), vec(pxy, pyy, pyz),
                      vec(pxz, pyz, pzz));
        } else {
          throw std::runtime_error("Invalid line in " + filename + ": " + line);
        }
        thisPole->setPs(P1, state);
        useDefaultPs = false;
      }
      // Multipole line
      else {
        int lineRank = int(sqrt(Qs.size()) + 0.5);
        if (lineRank == 0) {
          Q0_total += boost::lexical_cast<double>(split[0]);
        }
        for (unsigned i = 0; i < split.size(); i++) {
          double qXYZ = boost::lexical_cast<double>(split[i]);
          // Convert e*(a_0)^k to e*(nm)^k where k = rank
          double BOHR2NM = 0.0529189379;
          qXYZ *= pow(BOHR2NM, lineRank);  // OVERRIDE
          Qs.push_back(qXYZ);
        }
        thisPole->setQs(Qs, state);
      }
    } /* Exit loop over lines */
  } else {
    cout << endl << "ERROR: No such file " << filename << endl;
    throw runtime_error("Please supply input file.");
  }

  if (thread == NULL)
    printf("\n... ... ... Reading %-25s -> N = %2d Q0(Sum) = %+1.7f ",
           filename.c_str(), (int)poles.size(), Q0_total);

  // Apply charge correction: Sum to closest integer
  int Q_integer = boost::math::iround(Q0_total);
  double dQ = (double(Q_integer) - Q0_total) / poles.size();

  double Q0_total_corr = 0.0;
  vector<APolarSite *>::iterator pit;
  for (pit = poles.begin(); pit < poles.end(); ++pit) {
    double Q_uncorr = (*pit)->getQs(state)[0];
    double Q_corr = Q_uncorr + dQ;
    (*pit)->setQ00(Q_corr, state);
    Q0_total_corr += (*pit)->getQs(state)[0];
  }

  if (thread == NULL)
    printf("=> dQ0 = %+1.1e, Q0(corr.) = %+1.0f", dQ, Q0_total_corr);

  if (useDefaultPs) {
    if (thread == NULL)
      cout << endl
           << "... ... ... NOTE Using default Thole polarizabilities "
           << "for charge state " << state << ". ";

    std::map<string, double> polar_table = POLAR_TABLE();

    vector<APolarSite *>::iterator pol;
    for (pol = poles.begin(); pol < poles.end(); ++pol) {
      string elem = (*pol)->getName();
      double alpha = 0.0;

      try {
        alpha = polar_table.at(elem);
      } catch (const std::exception &out_of_range) {
        cout << endl
             << "WARNING No default polarizability given for "
             << "polar site type '" << elem << "'. Defaulting to 1 A**3. "
             << flush;
        alpha = 1e-3;
      }
      //            // Original set of Thole polarizabilites
      //            if      (elem == "C") { alpha = 1.75e-3;  } // <- conversion
      //            from else if (elem == "H") { alpha = 0.696e-3; } //    A³ to
      //            nm³ = 10⁻³ else if (elem == "N") { alpha = 1.073e-3; } else
      //            if (elem == "O") { alpha = 0.837e-3; } else if (elem == "S")
      //            { alpha = 2.926e-3; }
      //            // Different set of Thole polarizabilities
      //            //if      (elem == "C") { alpha = 1.334e-3; } // <-
      //            conversion from
      //            //else if (elem == "H") { alpha = 0.496e-3; } //    A³ to
      //            nm³ = 10⁻³
      //            //else if (elem == "N") { alpha = 1.073e-3; }
      //            //else if (elem == "O") { alpha = 0.837e-3; }
      //            //else if (elem == "S") { alpha = 3.300e-3; }
      //            else { throw runtime_error("No polarizability given "
      //                                       "for polar site type " + elem +
      //                                       ". "); }

      P1 = matrix(vec(alpha, 0, 0), vec(0, alpha, 0), vec(0, 0, alpha));

      (*pol)->setPs(P1, state);
    }
  }

  vector<APolarSite *>::iterator pol;
  for (pol = poles.begin(); pol < poles.end(); ++pol) {
    (*pol)->Charge(state);
  }

  return poles;
}

map<string, double> POLAR_TABLE() {
  map<string, double> polar_table;
  polar_table["H"] = 0.496e-3;
  polar_table["C"] = 1.334e-3;
  polar_table["N"] = 1.073e-3;
  polar_table["O"] = 0.837e-3;
  polar_table["S"] = 2.926e-3;
  polar_table["F"] = 0.440e-3;
  polar_table["Si"] = 3.962e-3;  // B3LYP/6-311+g(2d,2p)
  polar_table["Zn"] = 5.962e-3;  // B3LYP/6-311+g(2d,2p)
  polar_table["Al"] =
      5.80e-3;  //[1]P. Fuentealba, “The static dipole polarizability of
                // aluminium atom: discrepancy between theory and experiment,”
                // Chemical physics letters, vol. 397, no. 4, pp. 459–461, 2004.

  return polar_table;
}

}  // namespace xtp
}  // namespace votca
