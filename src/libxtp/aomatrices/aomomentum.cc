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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>

namespace votca {
namespace xtp {

void AOMomentum::FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd> >& matrix,
                           const AOShell& shell_row,
                           const AOShell& shell_col) const {

  /* Calculating the AO matrix of the gradient operator requires
   * the raw overlap matrix (i.e. in unnormalized cartesians)
   * with lmax of the column shell increased by one:
   *
   *         phi_(ijk) = x^i y^j z^k exp(-beta r^2)
   * => d/dx phi_(ijk) = (i*x^(i-1) - 2beta x^(i+1)) y^j z^k exp(-beta r^2)
   *    d/dy phi_(ijk) = x^i (j*y^(j-1) - 2beta y^(j+1)) z^k exp(-beta r^2)
   *    d/dz phi_(ijk) = x^i y^j (k*z^(k-1) - 2beta z^(k+1)) exp(-beta r^2)
   *
   * i.e.:   d/dx phi_s  = d/dx phi_(000) = -2beta phi_(100) = -2beta phi_px
   *         d/dy phi_px = d/dy phi_(100) = -2beta phi_(110) = -2beta phi_dxy
   *         d/dz phi_pz = d/dz phi_(001) = phi_(000) - 2beta phi_(002)
   *                                      = phi_s     - 2beta phi_dxx
   *
   * and with that
   *         <s|d/dx|s>  = -2beta <s|px>
   *         <s|d/dy|px> = -2beta <s|dxy>
   *         <s|d/dz|pz> = <s|s> - 2beta <s|dxx>
   *         ...
   *
   */

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row.getLmax();
  int lmax_col = shell_col.getLmax();

  if (lmax_col > 4) {
    throw std::runtime_error(
        "Momentum transition dipoles only implemented for S,P,D,F,G functions "
        "in DFT basis!");
  }

  // set size of internal block for recursion
  int nrows = this->getBlockSize(lmax_row);
  int ncols = this->getBlockSize(lmax_col);

  // initialize local matrix block for unnormalized cartesians
  std::vector<Eigen::MatrixXd> mom;
  for (int i_comp = 0; i_comp < 3; i_comp++) {
    mom.push_back(Eigen::MatrixXd ::Zero(nrows, ncols));
  }

  std::vector<Eigen::MatrixXd> scd_mom;
  for (int i_comp = 0; i_comp < 6; i_comp++) {
    scd_mom.push_back(Eigen::MatrixXd ::Zero(nrows, ncols));
  }

  // initialize local matrix block for unnormalized cartesians of overlap
  int nrows_ol = this->getBlockSize(lmax_row + 1);
  int ncols_ol = this->getBlockSize(lmax_col + 1);

  Eigen::MatrixXd ol = Eigen::MatrixXd::Zero(nrows_ol, ncols_ol);

  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;

  double distsq = diff.squaredNorm();

  int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

  int nx[] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0,
              0, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 5, 4, 4,
              3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};

  int ny[] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1,
              0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 0, 1, 0,
              2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0};

  int nz[] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2,
              3, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 0, 1,
              0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5};

  int i_less_x[] = {0,  0,  0,  0,  1,  2,  3,  0,  0,  0,  4,  5,  6,  7,
                    8,  9,  0,  0,  0,  0,  10, 11, 12, 13, 14, 15, 16, 17,
                    18, 19, 0,  0,  0,  0,  0,  20, 21, 22, 23, 24, 25, 26,
                    27, 28, 29, 30, 31, 32, 33, 34, 0,  0,  0,  0,  0,  0};

  int i_less_y[] = {0,  0,  0,  0,  0,  1,  0,  2,  3,  0,  0,  4,  0,  5,
                    6,  0,  7,  8,  9,  0,  0,  10, 0,  11, 12, 0,  13, 14,
                    15, 0,  16, 17, 18, 19, 0,  0,  20, 0,  21, 22, 0,  23,
                    24, 25, 0,  26, 27, 28, 29, 0,  30, 31, 32, 33, 34, 0};

  int i_less_z[] = {0,  0,  0,  0,  0,  0,  1,  0,  2,  3,  0,  0,  4,  0,
                    5,  6,  0,  7,  8,  9,  0,  0,  10, 0,  11, 12, 0,  13,
                    14, 15, 0,  16, 17, 18, 19, 0,  0,  20, 0,  21, 22, 0,
                    23, 24, 25, 0,  26, 27, 28, 29, 0,  30, 31, 32, 33, 34};

  int i_more_x[] = {1,  4,  5,  6,  10, 11, 12, 13, 14, 15, 20, 21,
                    22, 23, 24, 25, 26, 27, 28, 29, 35, 36, 37, 38,
                    39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49};

  int i_more_y[] = {2,  5,  7,  8,  11, 13, 14, 16, 17, 18, 21, 23,
                    24, 26, 27, 28, 30, 31, 32, 33, 36, 38, 39, 41,
                    42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 54};

  int i_more_z[] = {3,  6,  8,  9,  12, 14, 15, 17, 18, 19, 22, 24,
                    25, 27, 28, 29, 31, 32, 33, 34, 37, 39, 40, 42,
                    43, 44, 46, 47, 48, 49, 51, 52, 53, 54, 55};

  for (const auto& gaussian_row : shell_row) {
    const double decay_row = gaussian_row.getDecay();

    for (const auto& gaussian_col : shell_col) {

      const double decay_col = gaussian_col.getDecay();

      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;
      double _exparg = fak2 * decay_row * decay_col * distsq;

      /// check if distance between postions is big, then skip step

      if (_exparg > 30.0) {
        continue;
      }

      const Eigen::Vector3d PmA =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_row;
      const Eigen::Vector3d PmB =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_col;

      // calculate s-s- overlap matrix element
      ol(0, 0) = pow(4.0 * decay_row * decay_col, 0.75) * pow(fak2, 1.5) *
                 exp(-fak2 * decay_row * decay_col * distsq);  // s-s element

      // Integrals     p - s
      ol(Cart::x, 0) = PmA(0) * ol(0, 0);
      ol(Cart::y, 0) = PmA(1) * ol(0, 0);
      ol(Cart::z, 0) = PmA(2) * ol(0, 0);
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 0) {
        double term = fak * ol(0, 0);
        ol(Cart::xx, 0) = PmA(0) * ol(Cart::x, 0) + term;
        ol(Cart::xy, 0) = PmA(0) * ol(Cart::y, 0);
        ol(Cart::xz, 0) = PmA(0) * ol(Cart::z, 0);
        ol(Cart::yy, 0) = PmA(1) * ol(Cart::y, 0) + term;
        ol(Cart::yz, 0) = PmA(1) * ol(Cart::z, 0);
        ol(Cart::zz, 0) = PmA(2) * ol(Cart::z, 0) + term;
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 1) {
        ol(Cart::xxx, 0) = PmA(0) * ol(Cart::xx, 0) + 2 * fak * ol(Cart::x, 0);
        ol(Cart::xxy, 0) = PmA(1) * ol(Cart::xx, 0);
        ol(Cart::xxz, 0) = PmA(2) * ol(Cart::xx, 0);
        ol(Cart::xyy, 0) = PmA(0) * ol(Cart::yy, 0);
        ol(Cart::xyz, 0) = PmA(0) * ol(Cart::yz, 0);
        ol(Cart::xzz, 0) = PmA(0) * ol(Cart::zz, 0);
        ol(Cart::yyy, 0) = PmA(1) * ol(Cart::yy, 0) + 2 * fak * ol(Cart::y, 0);
        ol(Cart::yyz, 0) = PmA(2) * ol(Cart::yy, 0);
        ol(Cart::yzz, 0) = PmA(1) * ol(Cart::zz, 0);
        ol(Cart::zzz, 0) = PmA(2) * ol(Cart::zz, 0) + 2 * fak * ol(Cart::z, 0);
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 2) {
        double term_xx = fak * ol(Cart::xx, 0);
        double term_yy = fak * ol(Cart::yy, 0);
        double term_zz = fak * ol(Cart::zz, 0);
        ol(Cart::xxxx, 0) = PmA(0) * ol(Cart::xxx, 0) + 3 * term_xx;
        ol(Cart::xxxy, 0) = PmA(1) * ol(Cart::xxx, 0);
        ol(Cart::xxxz, 0) = PmA(2) * ol(Cart::xxx, 0);
        ol(Cart::xxyy, 0) = PmA(0) * ol(Cart::xyy, 0) + term_yy;
        ol(Cart::xxyz, 0) = PmA(1) * ol(Cart::xxz, 0);
        ol(Cart::xxzz, 0) = PmA(0) * ol(Cart::xzz, 0) + term_zz;
        ol(Cart::xyyy, 0) = PmA(0) * ol(Cart::yyy, 0);
        ol(Cart::xyyz, 0) = PmA(0) * ol(Cart::yyz, 0);
        ol(Cart::xyzz, 0) = PmA(0) * ol(Cart::yzz, 0);
        ol(Cart::xzzz, 0) = PmA(0) * ol(Cart::zzz, 0);
        ol(Cart::yyyy, 0) = PmA(1) * ol(Cart::yyy, 0) + 3 * term_yy;
        ol(Cart::yyyz, 0) = PmA(2) * ol(Cart::yyy, 0);
        ol(Cart::yyzz, 0) = PmA(1) * ol(Cart::yzz, 0) + term_zz;
        ol(Cart::yzzz, 0) = PmA(1) * ol(Cart::zzz, 0);
        ol(Cart::zzzz, 0) = PmA(2) * ol(Cart::zzz, 0) + 3 * term_zz;
      }
      //------------------------------------------------------

      // Integrals     h - s
      if (lmax_row > 3) {
        double term_xxx = fak * ol(Cart::xxx, 0);
        double term_yyy = fak * ol(Cart::yyy, 0);
        double term_zzz = fak * ol(Cart::zzz, 0);
        ol(Cart::xxxxx, 0) = PmA(0) * ol(Cart::xxxx, 0) + 4 * term_xxx;
        ol(Cart::xxxxy, 0) = PmA(1) * ol(Cart::xxxx, 0);
        ol(Cart::xxxxz, 0) = PmA(2) * ol(Cart::xxxx, 0);
        ol(Cart::xxxyy, 0) = PmA(1) * ol(Cart::xxxy, 0) + term_xxx;
        ol(Cart::xxxyz, 0) = PmA(1) * ol(Cart::xxxz, 0);
        ol(Cart::xxxzz, 0) = PmA(2) * ol(Cart::xxxz, 0) + term_xxx;
        ol(Cart::xxyyy, 0) = PmA(0) * ol(Cart::xyyy, 0) + term_yyy;
        ol(Cart::xxyyz, 0) = PmA(2) * ol(Cart::xxyy, 0);
        ol(Cart::xxyzz, 0) = PmA(1) * ol(Cart::xxzz, 0);
        ol(Cart::xxzzz, 0) = PmA(0) * ol(Cart::xzzz, 0) + term_zzz;
        ol(Cart::xyyyy, 0) = PmA(0) * ol(Cart::yyyy, 0);
        ol(Cart::xyyyz, 0) = PmA(0) * ol(Cart::yyyz, 0);
        ol(Cart::xyyzz, 0) = PmA(0) * ol(Cart::yyzz, 0);
        ol(Cart::xyzzz, 0) = PmA(0) * ol(Cart::yzzz, 0);
        ol(Cart::xzzzz, 0) = PmA(0) * ol(Cart::zzzz, 0);
        ol(Cart::yyyyy, 0) = PmA(1) * ol(Cart::yyyy, 0) + 4 * term_yyy;
        ol(Cart::yyyyz, 0) = PmA(2) * ol(Cart::yyyy, 0);
        ol(Cart::yyyzz, 0) = PmA(2) * ol(Cart::yyyz, 0) + term_yyy;
        ol(Cart::yyzzz, 0) = PmA(1) * ol(Cart::yzzz, 0) + term_zzz;
        ol(Cart::yzzzz, 0) = PmA(1) * ol(Cart::zzzz, 0);
        ol(Cart::zzzzz, 0) = PmA(2) * ol(Cart::zzzz, 0) + 4 * term_zzz;
      }
      //------------------------------------------------------

      // Integrals     s - p
      ol(0, Cart::x) = PmB(0) * ol(0, 0);
      ol(0, Cart::y) = PmB(1) * ol(0, 0);
      ol(0, Cart::z) = PmB(2) * ol(0, 0);
      //------------------------------------------------------

      // Integrals     p - p     d - p     f - p     g - p
      for (int i = 1; i < n_orbitals[lmax_row]; i++) {
        ol(i, Cart::x) = PmB(0) * ol(i, 0) + nx[i] * fak * ol(i_less_x[i], 0);
        ol(i, Cart::y) = PmB(1) * ol(i, 0) + ny[i] * fak * ol(i_less_y[i], 0);
        ol(i, Cart::z) = PmB(2) * ol(i, 0) + nz[i] * fak * ol(i_less_z[i], 0);
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - d
        double term = fak * ol(0, 0);
        ol(0, Cart::xx) = PmB(0) * ol(0, Cart::x) + term;
        ol(0, Cart::xy) = PmB(0) * ol(0, Cart::y);
        ol(0, Cart::xz) = PmB(0) * ol(0, Cart::z);
        ol(0, Cart::yy) = PmB(1) * ol(0, Cart::y) + term;
        ol(0, Cart::yz) = PmB(1) * ol(0, Cart::z);
        ol(0, Cart::zz) = PmB(2) * ol(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term = fak * ol(i, 0);
          ol(i, Cart::xx) = PmB(0) * ol(i, Cart::x) +
                            nx[i] * fak * ol(i_less_x[i], Cart::x) + term;
          ol(i, Cart::xy) =
              PmB(0) * ol(i, Cart::y) + nx[i] * fak * ol(i_less_x[i], Cart::y);
          ol(i, Cart::xz) =
              PmB(0) * ol(i, Cart::z) + nx[i] * fak * ol(i_less_x[i], Cart::z);
          ol(i, Cart::yy) = PmB(1) * ol(i, Cart::y) +
                            ny[i] * fak * ol(i_less_y[i], Cart::y) + term;
          ol(i, Cart::yz) =
              PmB(1) * ol(i, Cart::z) + ny[i] * fak * ol(i_less_y[i], Cart::z);
          ol(i, Cart::zz) = PmB(2) * ol(i, Cart::z) +
                            nz[i] * fak * ol(i_less_z[i], Cart::z) + term;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - f
        ol(0, Cart::xxx) = PmB(0) * ol(0, Cart::xx) + 2 * fak * ol(0, Cart::x);
        ol(0, Cart::xxy) = PmB(1) * ol(0, Cart::xx);
        ol(0, Cart::xxz) = PmB(2) * ol(0, Cart::xx);
        ol(0, Cart::xyy) = PmB(0) * ol(0, Cart::yy);
        ol(0, Cart::xyz) = PmB(0) * ol(0, Cart::yz);
        ol(0, Cart::xzz) = PmB(0) * ol(0, Cart::zz);
        ol(0, Cart::yyy) = PmB(1) * ol(0, Cart::yy) + 2 * fak * ol(0, Cart::y);
        ol(0, Cart::yyz) = PmB(2) * ol(0, Cart::yy);
        ol(0, Cart::yzz) = PmB(1) * ol(0, Cart::zz);
        ol(0, Cart::zzz) = PmB(2) * ol(0, Cart::zz) + 2 * fak * ol(0, Cart::z);
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          int nx_i = nx[i];
          int ny_i = ny[i];
          int nz_i = nz[i];
          int ilx_i = i_less_x[i];
          int ily_i = i_less_y[i];
          int ilz_i = i_less_z[i];
          double term_x = 2 * fak * ol(i, Cart::x);
          double term_y = 2 * fak * ol(i, Cart::y);
          double term_z = 2 * fak * ol(i, Cart::z);
          ol(i, Cart::xxx) = PmB(0) * ol(i, Cart::xx) +
                             nx_i * fak * ol(ilx_i, Cart::xx) + term_x;
          ol(i, Cart::xxy) =
              PmB(1) * ol(i, Cart::xx) + ny_i * fak * ol(ily_i, Cart::xx);
          ol(i, Cart::xxz) =
              PmB(2) * ol(i, Cart::xx) + nz_i * fak * ol(ilz_i, Cart::xx);
          ol(i, Cart::xyy) =
              PmB(0) * ol(i, Cart::yy) + nx_i * fak * ol(ilx_i, Cart::yy);
          ol(i, Cart::xyz) =
              PmB(0) * ol(i, Cart::yz) + nx_i * fak * ol(ilx_i, Cart::yz);
          ol(i, Cart::xzz) =
              PmB(0) * ol(i, Cart::zz) + nx_i * fak * ol(ilx_i, Cart::zz);
          ol(i, Cart::yyy) = PmB(1) * ol(i, Cart::yy) +
                             ny_i * fak * ol(ily_i, Cart::yy) + term_y;
          ol(i, Cart::yyz) =
              PmB(2) * ol(i, Cart::yy) + nz_i * fak * ol(ilz_i, Cart::yy);
          ol(i, Cart::yzz) =
              PmB(1) * ol(i, Cart::zz) + ny_i * fak * ol(ily_i, Cart::zz);
          ol(i, Cart::zzz) = PmB(2) * ol(i, Cart::zz) +
                             nz_i * fak * ol(ilz_i, Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - g
        double term_xx = fak * ol(0, Cart::xx);
        double term_yy = fak * ol(0, Cart::yy);
        double term_zz = fak * ol(0, Cart::zz);
        ol(0, Cart::xxxx) = PmB(0) * ol(0, Cart::xxx) + 3 * term_xx;
        ol(0, Cart::xxxy) = PmB(1) * ol(0, Cart::xxx);
        ol(0, Cart::xxxz) = PmB(2) * ol(0, Cart::xxx);
        ol(0, Cart::xxyy) = PmB(0) * ol(0, Cart::xyy) + term_yy;
        ol(0, Cart::xxyz) = PmB(1) * ol(0, Cart::xxz);
        ol(0, Cart::xxzz) = PmB(0) * ol(0, Cart::xzz) + term_zz;
        ol(0, Cart::xyyy) = PmB(0) * ol(0, Cart::yyy);
        ol(0, Cart::xyyz) = PmB(0) * ol(0, Cart::yyz);
        ol(0, Cart::xyzz) = PmB(0) * ol(0, Cart::yzz);
        ol(0, Cart::xzzz) = PmB(0) * ol(0, Cart::zzz);
        ol(0, Cart::yyyy) = PmB(1) * ol(0, Cart::yyy) + 3 * term_yy;
        ol(0, Cart::yyyz) = PmB(2) * ol(0, Cart::yyy);
        ol(0, Cart::yyzz) = PmB(1) * ol(0, Cart::yzz) + term_zz;
        ol(0, Cart::yzzz) = PmB(1) * ol(0, Cart::zzz);
        ol(0, Cart::zzzz) = PmB(2) * ol(0, Cart::zzz) + 3 * term_zz;
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          int nx_i = nx[i];
          int ny_i = ny[i];
          int nz_i = nz[i];
          int ilx_i = i_less_x[i];
          int ily_i = i_less_y[i];
          int ilz_i = i_less_z[i];
          double term_xx = fak * ol(i, Cart::xx);
          double term_yy = fak * ol(i, Cart::yy);
          double term_zz = fak * ol(i, Cart::zz);
          ol(i, Cart::xxxx) = PmB(0) * ol(i, Cart::xxx) +
                              nx_i * fak * ol(ilx_i, Cart::xxx) + 3 * term_xx;
          ol(i, Cart::xxxy) =
              PmB(1) * ol(i, Cart::xxx) + ny_i * fak * ol(ily_i, Cart::xxx);
          ol(i, Cart::xxxz) =
              PmB(2) * ol(i, Cart::xxx) + nz_i * fak * ol(ilz_i, Cart::xxx);
          ol(i, Cart::xxyy) = PmB(0) * ol(i, Cart::xyy) +
                              nx_i * fak * ol(ilx_i, Cart::xyy) + term_yy;
          ol(i, Cart::xxyz) =
              PmB(1) * ol(i, Cart::xxz) + ny_i * fak * ol(ily_i, Cart::xxz);
          ol(i, Cart::xxzz) = PmB(0) * ol(i, Cart::xzz) +
                              nx_i * fak * ol(ilx_i, Cart::xzz) + term_zz;
          ol(i, Cart::xyyy) =
              PmB(0) * ol(i, Cart::yyy) + nx_i * fak * ol(ilx_i, Cart::yyy);
          ol(i, Cart::xyyz) =
              PmB(0) * ol(i, Cart::yyz) + nx_i * fak * ol(ilx_i, Cart::yyz);
          ol(i, Cart::xyzz) =
              PmB(0) * ol(i, Cart::yzz) + nx_i * fak * ol(ilx_i, Cart::yzz);
          ol(i, Cart::xzzz) =
              PmB(0) * ol(i, Cart::zzz) + nx_i * fak * ol(ilx_i, Cart::zzz);
          ol(i, Cart::yyyy) = PmB(1) * ol(i, Cart::yyy) +
                              ny_i * fak * ol(ily_i, Cart::yyy) + 3 * term_yy;
          ol(i, Cart::yyyz) =
              PmB(2) * ol(i, Cart::yyy) + nz_i * fak * ol(ilz_i, Cart::yyy);
          ol(i, Cart::yyzz) = PmB(1) * ol(i, Cart::yzz) +
                              ny_i * fak * ol(ily_i, Cart::yzz) + term_zz;
          ol(i, Cart::yzzz) =
              PmB(1) * ol(i, Cart::zzz) + ny_i * fak * ol(ily_i, Cart::zzz);
          ol(i, Cart::zzzz) = PmB(2) * ol(i, Cart::zzz) +
                              nz_i * fak * ol(ilz_i, Cart::zzz) + 3 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - h
        double term_xxx = fak * ol(0, Cart::xxx);
        double term_yyy = fak * ol(0, Cart::yyy);
        double term_zzz = fak * ol(0, Cart::zzz);
        ol(0, Cart::xxxxx) = PmB(0) * ol(0, Cart::xxxx) + 4 * term_xxx;
        ol(0, Cart::xxxxy) = PmB(1) * ol(0, Cart::xxxx);
        ol(0, Cart::xxxxz) = PmB(2) * ol(0, Cart::xxxx);
        ol(0, Cart::xxxyy) = PmB(1) * ol(0, Cart::xxxy) + term_xxx;
        ol(0, Cart::xxxyz) = PmB(1) * ol(0, Cart::xxxz);
        ol(0, Cart::xxxzz) = PmB(2) * ol(0, Cart::xxxz) + term_xxx;
        ol(0, Cart::xxyyy) = PmB(0) * ol(0, Cart::xyyy) + term_yyy;
        ol(0, Cart::xxyyz) = PmB(2) * ol(0, Cart::xxyy);
        ol(0, Cart::xxyzz) = PmB(1) * ol(0, Cart::xxzz);
        ol(0, Cart::xxzzz) = PmB(0) * ol(0, Cart::xzzz) + term_zzz;
        ol(0, Cart::xyyyy) = PmB(0) * ol(0, Cart::yyyy);
        ol(0, Cart::xyyyz) = PmB(0) * ol(0, Cart::yyyz);
        ol(0, Cart::xyyzz) = PmB(0) * ol(0, Cart::yyzz);
        ol(0, Cart::xyzzz) = PmB(0) * ol(0, Cart::yzzz);
        ol(0, Cart::xzzzz) = PmB(0) * ol(0, Cart::zzzz);
        ol(0, Cart::yyyyy) = PmB(1) * ol(0, Cart::yyyy) + 4 * term_yyy;
        ol(0, Cart::yyyyz) = PmB(2) * ol(0, Cart::yyyy);
        ol(0, Cart::yyyzz) = PmB(2) * ol(0, Cart::yyyz) + term_yyy;
        ol(0, Cart::yyzzz) = PmB(1) * ol(0, Cart::yzzz) + term_zzz;
        ol(0, Cart::yzzzz) = PmB(1) * ol(0, Cart::zzzz);
        ol(0, Cart::zzzzz) = PmB(2) * ol(0, Cart::zzzz) + 4 * term_zzz;
        //------------------------------------------------------

        // Integrals     p - h     d - h     f - h     g - h
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          int nx_i = nx[i];
          int ny_i = ny[i];
          int nz_i = nz[i];
          int ilx_i = i_less_x[i];
          int ily_i = i_less_y[i];
          int ilz_i = i_less_z[i];
          double term_xxx = fak * ol(i, Cart::xxx);
          double term_yyy = fak * ol(i, Cart::yyy);
          double term_zzz = fak * ol(i, Cart::zzz);
          ol(i, Cart::xxxxx) = PmB(0) * ol(i, Cart::xxxx) +
                               nx_i * fak * ol(ilx_i, Cart::xxxx) +
                               4 * term_xxx;
          ol(i, Cart::xxxxy) =
              PmB(1) * ol(i, Cart::xxxx) + ny_i * fak * ol(ily_i, Cart::xxxx);
          ol(i, Cart::xxxxz) =
              PmB(2) * ol(i, Cart::xxxx) + nz_i * fak * ol(ilz_i, Cart::xxxx);
          ol(i, Cart::xxxyy) = PmB(1) * ol(i, Cart::xxxy) +
                               ny_i * fak * ol(ily_i, Cart::xxxy) + term_xxx;
          ol(i, Cart::xxxyz) =
              PmB(1) * ol(i, Cart::xxxz) + ny_i * fak * ol(ily_i, Cart::xxxz);
          ol(i, Cart::xxxzz) = PmB(2) * ol(i, Cart::xxxz) +
                               nz_i * fak * ol(ilz_i, Cart::xxxz) + term_xxx;
          ol(i, Cart::xxyyy) = PmB(0) * ol(i, Cart::xyyy) +
                               nx_i * fak * ol(ilx_i, Cart::xyyy) + term_yyy;
          ol(i, Cart::xxyyz) =
              PmB(2) * ol(i, Cart::xxyy) + nz_i * fak * ol(ilz_i, Cart::xxyy);
          ol(i, Cart::xxyzz) =
              PmB(1) * ol(i, Cart::xxzz) + ny_i * fak * ol(ily_i, Cart::xxzz);
          ol(i, Cart::xxzzz) = PmB(0) * ol(i, Cart::xzzz) +
                               nx_i * fak * ol(ilx_i, Cart::xzzz) + term_zzz;
          ol(i, Cart::xyyyy) =
              PmB(0) * ol(i, Cart::yyyy) + nx_i * fak * ol(ilx_i, Cart::yyyy);
          ol(i, Cart::xyyyz) =
              PmB(0) * ol(i, Cart::yyyz) + nx_i * fak * ol(ilx_i, Cart::yyyz);
          ol(i, Cart::xyyzz) =
              PmB(0) * ol(i, Cart::yyzz) + nx_i * fak * ol(ilx_i, Cart::yyzz);
          ol(i, Cart::xyzzz) =
              PmB(0) * ol(i, Cart::yzzz) + nx_i * fak * ol(ilx_i, Cart::yzzz);
          ol(i, Cart::xzzzz) =
              PmB(0) * ol(i, Cart::zzzz) + nx_i * fak * ol(ilx_i, Cart::zzzz);
          ol(i, Cart::yyyyy) = PmB(1) * ol(i, Cart::yyyy) +
                               ny_i * fak * ol(ily_i, Cart::yyyy) +
                               4 * term_yyy;
          ol(i, Cart::yyyyz) =
              PmB(2) * ol(i, Cart::yyyy) + nz_i * fak * ol(ilz_i, Cart::yyyy);
          ol(i, Cart::yyyzz) = PmB(2) * ol(i, Cart::yyyz) +
                               nz_i * fak * ol(ilz_i, Cart::yyyz) + term_yyy;
          ol(i, Cart::yyzzz) = PmB(1) * ol(i, Cart::yzzz) +
                               ny_i * fak * ol(ily_i, Cart::yzzz) + term_zzz;
          ol(i, Cart::yzzzz) =
              PmB(1) * ol(i, Cart::zzzz) + ny_i * fak * ol(ily_i, Cart::zzzz);
          ol(i, Cart::zzzzz) = PmB(2) * ol(i, Cart::zzzz) +
                               nz_i * fak * ol(ilz_i, Cart::zzzz) +
                               4 * term_zzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      double alpha2 = 2.0 * decay_row;
      double beta2 = 2.0 * decay_col;
      for (int i = 0; i < ncols; i++) {

        int nx_i = nx[i];
        int ny_i = ny[i];
        int nz_i = nz[i];
        int ilx_i = i_less_x[i];
        int ily_i = i_less_y[i];
        int ilz_i = i_less_z[i];
        int imx_i = i_more_x[i];
        int imy_i = i_more_y[i];
        int imz_i = i_more_z[i];

        for (int j = 0; j < nrows; j++) {

          mom[0](j, i) = nx_i * ol(j, ilx_i) - beta2 * ol(j, imx_i);
          mom[1](j, i) = ny_i * ol(j, ily_i) - beta2 * ol(j, imy_i);
          mom[2](j, i) = nz_i * ol(j, ilz_i) - beta2 * ol(j, imz_i);

          int nx_j = nx[j];
          int ny_j = ny[j];
          int nz_j = nz[j];
          int ilx_j = i_less_x[j];
          int ily_j = i_less_y[j];
          int ilz_j = i_less_z[j];
          int imx_j = i_more_x[j];
          int imy_j = i_more_y[j];
          int imz_j = i_more_z[j];

          scd_mom[0](j, i) =
              nx_j * (beta2 * ol(ilx_j, imx_i) - nx_i * ol(ilx_j, ilx_i)) -
              alpha2 * (beta2 * ol(imx_j, imx_i) -
                        nx_i * ol(imx_j, ilx_i));  // d2/(dxdx)
          scd_mom[1](j, i) =
              nx_j * (beta2 * ol(ilx_j, imy_i) - ny_i * ol(ilx_j, ily_i)) -
              alpha2 * (beta2 * ol(imx_j, imy_i) -
                        ny_i * ol(imx_j, ily_i));  // d2/(dxdy)
          scd_mom[2](j, i) =
              nx_j * (beta2 * ol(ilx_j, imz_i) - nz_i * ol(ilx_j, ilz_i)) -
              alpha2 * (beta2 * ol(imx_j, imz_i) -
                        nz_i * ol(imx_j, ilz_i));  // d2/(dxdz)

          scd_mom[3](j, i) =
              ny_j * (beta2 * ol(ily_j, imy_i) - ny_i * ol(ily_j, ily_i)) -
              alpha2 * (beta2 * ol(imy_j, imy_i) -
                        ny_i * ol(imy_j, ily_i));  // d2/(dydy)
          scd_mom[4](j, i) =
              ny_j * (beta2 * ol(ily_j, imz_i) - nz_i * ol(ily_j, ilz_i)) -
              alpha2 * (beta2 * ol(imy_j, imz_i) -
                        nz_i * ol(imy_j, ilz_i));  // d2/(dydz)

          scd_mom[5](j, i) =
              nz_j * (beta2 * ol(ilz_j, imz_i) - nz_i * ol(ilz_j, ilz_i)) -
              alpha2 * (beta2 * ol(imz_j, imz_i) -
                        nz_i * ol(imz_j, ilz_i));  // d2/(dzdz)
        }
      }

      Eigen::MatrixXd trafo_row = getTrafo(gaussian_row);
      Eigen::MatrixXd trafo_col = getTrafo(gaussian_col);
      // cartesian -> spherical
      for (int i_comp = 0; i_comp < 3; i_comp++) {
        Eigen::MatrixXd mom_sph =
            trafo_row.transpose() * mom[i_comp] * trafo_col;
        // save to matrix
        for (unsigned i = 0; i < matrix[0].rows(); i++) {
          for (unsigned j = 0; j < matrix[0].cols(); j++) {
            matrix[i_comp](i, j) +=
                mom_sph(i + shell_row.getOffset(), j + shell_col.getOffset());
          }
        }
      }

    }  // shell_col Gaussians
  }    // shell_row Gaussians
}

}  // namespace xtp
}  // namespace votca
