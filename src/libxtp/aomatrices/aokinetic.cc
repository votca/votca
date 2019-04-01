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

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>

namespace votca {
namespace xtp {

void AOKinetic::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                          const AOShell* shell_row, const AOShell* shell_col) {

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row->getLmax();
  int lmax_col = shell_col->getLmax();

  if (lmax_col > 4 || lmax_row > 4) {
    std::cerr << "Orbitals higher than g are not yet implemented. This should "
                 "not have happened!"
              << std::flush;
    exit(1);
  }

  // set size of internal block for recursion
  int nrows = this->getBlockSize(lmax_row);
  int ncols = this->getBlockSize(lmax_col);

  // get shell positions
  const tools::vec& pos_row = shell_row->getPos();
  const tools::vec& _pos_col = shell_col->getPos();
  const tools::vec diff = pos_row - _pos_col;

  double distsq = (diff * diff);

  int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

  int nx[] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0,
              0, 0, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0};

  int ny[] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2,
              1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0};

  int nz[] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1,
              2, 3, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4};

  int i_less_x[] = {0,  0,  0,  0,  1,  2,  3, 0, 0,  0,  4,  5,
                    6,  7,  8,  9,  0,  0,  0, 0, 10, 11, 12, 13,
                    14, 15, 16, 17, 18, 19, 0, 0, 0,  0,  0};

  int i_less_y[] = {0,  0, 0,  0,  0,  1, 0,  2,  3,  0,  0, 4,
                    0,  5, 6,  0,  7,  8, 9,  0,  0,  10, 0, 11,
                    12, 0, 13, 14, 15, 0, 16, 17, 18, 19, 0};

  int i_less_z[] = {0,  0,  0, 0,  0,  0,  1, 0,  2,  3,  0,  0,
                    4,  0,  5, 6,  0,  7,  8, 9,  0,  0,  10, 0,
                    11, 12, 0, 13, 14, 15, 0, 16, 17, 18, 19};

  // cout << "row shell is " << shell_row->getSize() << " -fold contracted!" <<
  // endl;
  // cout << "col shell is " << shell_col->getSize() << " -fold contracted!" <<
  // endl;

  // iterate over Gaussians in this shell_row
  for (AOShell::GaussianIterator itr = shell_row->begin();
       itr != shell_row->end(); ++itr) {
    // iterate over Gaussians in this shell_col
    const double decay_row = itr->getDecay();

    for (AOShell::GaussianIterator itc = shell_col->begin();
         itc != shell_col->end(); ++itc) {

      // get decay constants
      const double decay_col = itc->getDecay();

      // some helpers
      const double fak = 0.5 / (decay_row + decay_col);
      const double rzeta = 2.0 * fak;

      const double xi = rzeta * decay_row * decay_col;

      // check if distance between postions is big, then skip step
      double exparg = xi * distsq;
      if (exparg > 30.0) {
        continue;
      }

      const double PmA0 =
          rzeta * (decay_row * pos_row.getX() + decay_col * _pos_col.getX()) -
          pos_row.getX();
      const double PmA1 =
          rzeta * (decay_row * pos_row.getY() + decay_col * _pos_col.getY()) -
          pos_row.getY();
      const double PmA2 =
          rzeta * (decay_row * pos_row.getZ() + decay_col * _pos_col.getZ()) -
          pos_row.getZ();

      const double PmB0 =
          rzeta * (decay_row * pos_row.getX() + decay_col * _pos_col.getX()) -
          _pos_col.getX();
      const double PmB1 =
          rzeta * (decay_row * pos_row.getY() + decay_col * _pos_col.getY()) -
          _pos_col.getY();
      const double PmB2 =
          rzeta * (decay_row * pos_row.getZ() + decay_col * _pos_col.getZ()) -
          _pos_col.getZ();

      const double xi2 = 2. * xi;              //////////////
      const double fak_a = rzeta * decay_col;  /////////////
      const double fak_b = rzeta * decay_row;  ////////////

      // matrix for kinetic energies
      Eigen::MatrixXd kin = Eigen::MatrixXd::Zero(nrows, ncols);
      // matrix for unnormalized overlap integrals
      Eigen::MatrixXd ol = Eigen::MatrixXd::Zero(nrows, ncols);

      // s-s overlap integral
      ol(Cart::s, Cart::s) = pow(rzeta, 1.5) *
                             pow(4.0 * decay_row * decay_col, 0.75) *
                             exp(-exparg);
      // s-s- kinetic energy integral
      kin(Cart::s, Cart::s) = ol(Cart::s, Cart::s) * xi * (3 - 2 * xi * distsq);

      // Integrals     p - s
      if (lmax_row > 0) {
        ol(Cart::x, 0) = PmA0 * ol(0, 0);
        ol(Cart::y, 0) = PmA1 * ol(0, 0);
        ol(Cart::z, 0) = PmA2 * ol(0, 0);
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        double term = fak * ol(0, 0);
        ol(Cart::xx, 0) = PmA0 * ol(Cart::x, 0) + term;
        ol(Cart::xy, 0) = PmA0 * ol(Cart::y, 0);
        ol(Cart::xz, 0) = PmA0 * ol(Cart::z, 0);
        ol(Cart::yy, 0) = PmA1 * ol(Cart::y, 0) + term;
        ol(Cart::yz, 0) = PmA1 * ol(Cart::z, 0);
        ol(Cart::zz, 0) = PmA2 * ol(Cart::z, 0) + term;
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        ol(Cart::xxx, 0) = PmA0 * ol(Cart::xx, 0) + 2 * fak * ol(Cart::x, 0);
        ol(Cart::xxy, 0) = PmA1 * ol(Cart::xx, 0);
        ol(Cart::xxz, 0) = PmA2 * ol(Cart::xx, 0);
        ol(Cart::xyy, 0) = PmA0 * ol(Cart::yy, 0);
        ol(Cart::xyz, 0) = PmA0 * ol(Cart::yz, 0);
        ol(Cart::xzz, 0) = PmA0 * ol(Cart::zz, 0);
        ol(Cart::yyy, 0) = PmA1 * ol(Cart::yy, 0) + 2 * fak * ol(Cart::y, 0);
        ol(Cart::yyz, 0) = PmA2 * ol(Cart::yy, 0);
        ol(Cart::yzz, 0) = PmA1 * ol(Cart::zz, 0);
        ol(Cart::zzz, 0) = PmA2 * ol(Cart::zz, 0) + 2 * fak * ol(Cart::z, 0);
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        double term_xx = fak * ol(Cart::xx, 0);
        double term_yy = fak * ol(Cart::yy, 0);
        double term_zz = fak * ol(Cart::zz, 0);
        ol(Cart::xxxx, 0) = PmA0 * ol(Cart::xxx, 0) + 3 * term_xx;
        ol(Cart::xxxy, 0) = PmA1 * ol(Cart::xxx, 0);
        ol(Cart::xxxz, 0) = PmA2 * ol(Cart::xxx, 0);
        ol(Cart::xxyy, 0) = PmA0 * ol(Cart::xyy, 0) + term_yy;
        ol(Cart::xxyz, 0) = PmA1 * ol(Cart::xxz, 0);
        ol(Cart::xxzz, 0) = PmA0 * ol(Cart::xzz, 0) + term_zz;
        ol(Cart::xyyy, 0) = PmA0 * ol(Cart::yyy, 0);
        ol(Cart::xyyz, 0) = PmA0 * ol(Cart::yyz, 0);
        ol(Cart::xyzz, 0) = PmA0 * ol(Cart::yzz, 0);
        ol(Cart::xzzz, 0) = PmA0 * ol(Cart::zzz, 0);
        ol(Cart::yyyy, 0) = PmA1 * ol(Cart::yyy, 0) + 3 * term_yy;
        ol(Cart::yyyz, 0) = PmA2 * ol(Cart::yyy, 0);
        ol(Cart::yyzz, 0) = PmA1 * ol(Cart::yzz, 0) + term_zz;
        ol(Cart::yzzz, 0) = PmA1 * ol(Cart::zzz, 0);
        ol(Cart::zzzz, 0) = PmA2 * ol(Cart::zzz, 0) + 3 * term_zz;
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        ol(0, Cart::x) = PmB0 * ol(0, 0);
        ol(0, Cart::y) = PmB1 * ol(0, 0);
        ol(0, Cart::z) = PmB2 * ol(0, 0);
        //------------------------------------------------------

        // Integrals     p - p     d - p     f - p     g - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          ol(i, Cart::x) = PmB0 * ol(i, 0) + nx[i] * fak * ol(i_less_x[i], 0);
          ol(i, Cart::y) = PmB1 * ol(i, 0) + ny[i] * fak * ol(i_less_y[i], 0);
          ol(i, Cart::z) = PmB2 * ol(i, 0) + nz[i] * fak * ol(i_less_z[i], 0);
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        double term = fak * ol(0, 0);
        ol(0, Cart::xx) = PmB0 * ol(0, Cart::x) + term;
        ol(0, Cart::xy) = PmB0 * ol(0, Cart::y);
        ol(0, Cart::xz) = PmB0 * ol(0, Cart::z);
        ol(0, Cart::yy) = PmB1 * ol(0, Cart::y) + term;
        ol(0, Cart::yz) = PmB1 * ol(0, Cart::z);
        ol(0, Cart::zz) = PmB2 * ol(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term = fak * ol(i, 0);
          ol(i, Cart::xx) = PmB0 * ol(i, Cart::x) +
                            nx[i] * fak * ol(i_less_x[i], Cart::x) + term;
          ol(i, Cart::xy) =
              PmB0 * ol(i, Cart::y) + nx[i] * fak * ol(i_less_x[i], Cart::y);
          ol(i, Cart::xz) =
              PmB0 * ol(i, Cart::z) + nx[i] * fak * ol(i_less_x[i], Cart::z);
          ol(i, Cart::yy) = PmB1 * ol(i, Cart::y) +
                            ny[i] * fak * ol(i_less_y[i], Cart::y) + term;
          ol(i, Cart::yz) =
              PmB1 * ol(i, Cart::z) + ny[i] * fak * ol(i_less_y[i], Cart::z);
          ol(i, Cart::zz) = PmB2 * ol(i, Cart::z) +
                            nz[i] * fak * ol(i_less_z[i], Cart::z) + term;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        ol(0, Cart::xxx) = PmB0 * ol(0, Cart::xx) + 2 * fak * ol(0, Cart::x);
        ol(0, Cart::xxy) = PmB1 * ol(0, Cart::xx);
        ol(0, Cart::xxz) = PmB2 * ol(0, Cart::xx);
        ol(0, Cart::xyy) = PmB0 * ol(0, Cart::yy);
        ol(0, Cart::xyz) = PmB0 * ol(0, Cart::yz);
        ol(0, Cart::xzz) = PmB0 * ol(0, Cart::zz);
        ol(0, Cart::yyy) = PmB1 * ol(0, Cart::yy) + 2 * fak * ol(0, Cart::y);
        ol(0, Cart::yyz) = PmB2 * ol(0, Cart::yy);
        ol(0, Cart::yzz) = PmB1 * ol(0, Cart::zz);
        ol(0, Cart::zzz) = PmB2 * ol(0, Cart::zz) + 2 * fak * ol(0, Cart::z);
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_x = 2 * fak * ol(i, Cart::x);
          double term_y = 2 * fak * ol(i, Cart::y);
          double term_z = 2 * fak * ol(i, Cart::z);
          ol(i, Cart::xxx) = PmB0 * ol(i, Cart::xx) +
                             nx[i] * fak * ol(i_less_x[i], Cart::xx) + term_x;
          ol(i, Cart::xxy) =
              PmB1 * ol(i, Cart::xx) + ny[i] * fak * ol(i_less_y[i], Cart::xx);
          ol(i, Cart::xxz) =
              PmB2 * ol(i, Cart::xx) + nz[i] * fak * ol(i_less_z[i], Cart::xx);
          ol(i, Cart::xyy) =
              PmB0 * ol(i, Cart::yy) + nx[i] * fak * ol(i_less_x[i], Cart::yy);
          ol(i, Cart::xyz) =
              PmB0 * ol(i, Cart::yz) + nx[i] * fak * ol(i_less_x[i], Cart::yz);
          ol(i, Cart::xzz) =
              PmB0 * ol(i, Cart::zz) + nx[i] * fak * ol(i_less_x[i], Cart::zz);
          ol(i, Cart::yyy) = PmB1 * ol(i, Cart::yy) +
                             ny[i] * fak * ol(i_less_y[i], Cart::yy) + term_y;
          ol(i, Cart::yyz) =
              PmB2 * ol(i, Cart::yy) + nz[i] * fak * ol(i_less_z[i], Cart::yy);
          ol(i, Cart::yzz) =
              PmB1 * ol(i, Cart::zz) + ny[i] * fak * ol(i_less_y[i], Cart::zz);
          ol(i, Cart::zzz) = PmB2 * ol(i, Cart::zz) +
                             nz[i] * fak * ol(i_less_z[i], Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        double term_xx = fak * ol(0, Cart::xx);
        double term_yy = fak * ol(0, Cart::yy);
        double term_zz = fak * ol(0, Cart::zz);
        ol(0, Cart::xxxx) = PmB0 * ol(0, Cart::xxx) + 3 * term_xx;
        ol(0, Cart::xxxy) = PmB1 * ol(0, Cart::xxx);
        ol(0, Cart::xxxz) = PmB2 * ol(0, Cart::xxx);
        ol(0, Cart::xxyy) = PmB0 * ol(0, Cart::xyy) + term_yy;
        ol(0, Cart::xxyz) = PmB1 * ol(0, Cart::xxz);
        ol(0, Cart::xxzz) = PmB0 * ol(0, Cart::xzz) + term_zz;
        ol(0, Cart::xyyy) = PmB0 * ol(0, Cart::yyy);
        ol(0, Cart::xyyz) = PmB0 * ol(0, Cart::yyz);
        ol(0, Cart::xyzz) = PmB0 * ol(0, Cart::yzz);
        ol(0, Cart::xzzz) = PmB0 * ol(0, Cart::zzz);
        ol(0, Cart::yyyy) = PmB1 * ol(0, Cart::yyy) + 3 * term_yy;
        ol(0, Cart::yyyz) = PmB2 * ol(0, Cart::yyy);
        ol(0, Cart::yyzz) = PmB1 * ol(0, Cart::yzz) + term_zz;
        ol(0, Cart::yzzz) = PmB1 * ol(0, Cart::zzz);
        ol(0, Cart::zzzz) = PmB2 * ol(0, Cart::zzz) + 3 * term_zz;
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xx = fak * ol(i, Cart::xx);
          double term_yy = fak * ol(i, Cart::yy);
          double term_zz = fak * ol(i, Cart::zz);
          ol(i, Cart::xxxx) = PmB0 * ol(i, Cart::xxx) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xxx) +
                              3 * term_xx;
          ol(i, Cart::xxxy) = PmB1 * ol(i, Cart::xxx) +
                              ny[i] * fak * ol(i_less_y[i], Cart::xxx);
          ol(i, Cart::xxxz) = PmB2 * ol(i, Cart::xxx) +
                              nz[i] * fak * ol(i_less_z[i], Cart::xxx);
          ol(i, Cart::xxyy) = PmB0 * ol(i, Cart::xyy) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xyy) +
                              term_yy;
          ol(i, Cart::xxyz) = PmB1 * ol(i, Cart::xxz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::xxz);
          ol(i, Cart::xxzz) = PmB0 * ol(i, Cart::xzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xzz) +
                              term_zz;
          ol(i, Cart::xyyy) = PmB0 * ol(i, Cart::yyy) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yyy);
          ol(i, Cart::xyyz) = PmB0 * ol(i, Cart::yyz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yyz);
          ol(i, Cart::xyzz) = PmB0 * ol(i, Cart::yzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yzz);
          ol(i, Cart::xzzz) = PmB0 * ol(i, Cart::zzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::zzz);
          ol(i, Cart::yyyy) = PmB1 * ol(i, Cart::yyy) +
                              ny[i] * fak * ol(i_less_y[i], Cart::yyy) +
                              3 * term_yy;
          ol(i, Cart::yyyz) = PmB2 * ol(i, Cart::yyy) +
                              nz[i] * fak * ol(i_less_z[i], Cart::yyy);
          ol(i, Cart::yyzz) = PmB1 * ol(i, Cart::yzz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::yzz) +
                              term_zz;
          ol(i, Cart::yzzz) = PmB1 * ol(i, Cart::zzz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::zzz);
          ol(i, Cart::zzzz) = PmB2 * ol(i, Cart::zzz) +
                              nz[i] * fak * ol(i_less_z[i], Cart::zzz) +
                              3 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      // Integrals     p - s
      if (lmax_row > 0) {
        kin(Cart::x, 0) = xi2 * ol(Cart::x, 0) + PmA0 * kin(0, 0);
        kin(Cart::y, 0) = xi2 * ol(Cart::y, 0) + PmA1 * kin(0, 0);
        kin(Cart::z, 0) = xi2 * ol(Cart::z, 0) + PmA2 * kin(0, 0);
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        double term = fak * kin(0, 0) - fak_a * ol(0, 0);
        kin(Cart::xx, 0) =
            xi2 * ol(Cart::xx, 0) + PmA0 * kin(Cart::x, 0) + term;
        kin(Cart::xy, 0) = xi2 * ol(Cart::xy, 0) + PmA0 * kin(Cart::y, 0);
        kin(Cart::xz, 0) = xi2 * ol(Cart::xz, 0) + PmA0 * kin(Cart::z, 0);
        kin(Cart::yy, 0) =
            xi2 * ol(Cart::yy, 0) + PmA1 * kin(Cart::y, 0) + term;
        kin(Cart::yz, 0) = xi2 * ol(Cart::yz, 0) + PmA1 * kin(Cart::z, 0);
        kin(Cart::zz, 0) =
            xi2 * ol(Cart::zz, 0) + PmA2 * kin(Cart::z, 0) + term;
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        kin(Cart::xxx, 0) =
            xi2 * ol(Cart::xxx, 0) + PmA0 * kin(Cart::xx, 0) +
            2 * (fak * kin(Cart::x, 0) - fak_a * ol(Cart::x, 0));
        kin(Cart::xxy, 0) = xi2 * ol(Cart::xxy, 0) + PmA1 * kin(Cart::xx, 0);
        kin(Cart::xxz, 0) = xi2 * ol(Cart::xxz, 0) + PmA2 * kin(Cart::xx, 0);
        kin(Cart::xyy, 0) = xi2 * ol(Cart::xyy, 0) + PmA0 * kin(Cart::yy, 0);
        kin(Cart::xyz, 0) = xi2 * ol(Cart::xyz, 0) + PmA0 * kin(Cart::yz, 0);
        kin(Cart::xzz, 0) = xi2 * ol(Cart::xzz, 0) + PmA0 * kin(Cart::zz, 0);
        kin(Cart::yyy, 0) =
            xi2 * ol(Cart::yyy, 0) + PmA1 * kin(Cart::yy, 0) +
            2 * (fak * kin(Cart::y, 0) - fak_a * ol(Cart::y, 0));
        kin(Cart::yyz, 0) = xi2 * ol(Cart::yyz, 0) + PmA2 * kin(Cart::yy, 0);
        kin(Cart::yzz, 0) = xi2 * ol(Cart::yzz, 0) + PmA1 * kin(Cart::zz, 0);
        kin(Cart::zzz, 0) =
            xi2 * ol(Cart::zzz, 0) + PmA2 * kin(Cart::zz, 0) +
            2 * (fak * kin(Cart::z, 0) - fak_a * ol(Cart::z, 0));
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        double term_xx = fak * kin(Cart::xx, 0) - fak_a * ol(Cart::xx, 0);
        double term_yy = fak * kin(Cart::yy, 0) - fak_a * ol(Cart::yy, 0);
        double term_zz = fak * kin(Cart::zz, 0) - fak_a * ol(Cart::zz, 0);
        kin(Cart::xxxx, 0) =
            xi2 * ol(Cart::xxxx, 0) + PmA0 * kin(Cart::xxx, 0) + 3 * term_xx;
        kin(Cart::xxxy, 0) = xi2 * ol(Cart::xxxy, 0) + PmA1 * kin(Cart::xxx, 0);
        kin(Cart::xxxz, 0) = xi2 * ol(Cart::xxxz, 0) + PmA2 * kin(Cart::xxx, 0);
        kin(Cart::xxyy, 0) =
            xi2 * ol(Cart::xxyy, 0) + PmA0 * kin(Cart::xyy, 0) + term_yy;
        kin(Cart::xxyz, 0) = xi2 * ol(Cart::xxyz, 0) + PmA1 * kin(Cart::xxz, 0);
        kin(Cart::xxzz, 0) =
            xi2 * ol(Cart::xxzz, 0) + PmA0 * kin(Cart::xzz, 0) + term_zz;
        kin(Cart::xyyy, 0) = xi2 * ol(Cart::xyyy, 0) + PmA0 * kin(Cart::yyy, 0);
        kin(Cart::xyyz, 0) = xi2 * ol(Cart::xyyz, 0) + PmA0 * kin(Cart::yyz, 0);
        kin(Cart::xyzz, 0) = xi2 * ol(Cart::xyzz, 0) + PmA0 * kin(Cart::yzz, 0);
        kin(Cart::xzzz, 0) = xi2 * ol(Cart::xzzz, 0) + PmA0 * kin(Cart::zzz, 0);
        kin(Cart::yyyy, 0) =
            xi2 * ol(Cart::yyyy, 0) + PmA1 * kin(Cart::yyy, 0) + 3 * term_yy;
        kin(Cart::yyyz, 0) = xi2 * ol(Cart::yyyz, 0) + PmA2 * kin(Cart::yyy, 0);
        kin(Cart::yyzz, 0) =
            xi2 * ol(Cart::yyzz, 0) + PmA1 * kin(Cart::yzz, 0) + term_zz;
        kin(Cart::yzzz, 0) = xi2 * ol(Cart::yzzz, 0) + PmA1 * kin(Cart::zzz, 0);
        kin(Cart::zzzz, 0) =
            xi2 * ol(Cart::zzzz, 0) + PmA2 * kin(Cart::zzz, 0) + 3 * term_zz;
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        kin(0, Cart::x) = xi2 * ol(0, Cart::x) + PmB0 * kin(0, 0);
        kin(0, Cart::y) = xi2 * ol(0, Cart::y) + PmB1 * kin(0, 0);
        kin(0, Cart::z) = xi2 * ol(0, Cart::z) + PmB2 * kin(0, 0);
        //------------------------------------------------------

        // Integrals     p - p     d - p     f - p     g - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          kin(i, Cart::x) = xi2 * ol(i, Cart::x) + PmB0 * kin(i, 0) +
                            nx[i] * fak * kin(i_less_x[i], 0);
          kin(i, Cart::y) = xi2 * ol(i, Cart::y) + PmB1 * kin(i, 0) +
                            ny[i] * fak * kin(i_less_y[i], 0);
          kin(i, Cart::z) = xi2 * ol(i, Cart::z) + PmB2 * kin(i, 0) +
                            nz[i] * fak * kin(i_less_z[i], 0);
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        double term = fak * kin(0, 0) - fak_b * ol(0, 0);
        kin(0, Cart::xx) =
            xi2 * ol(0, Cart::xx) + PmB0 * kin(0, Cart::x) + term;
        kin(0, Cart::xy) = xi2 * ol(0, Cart::xy) + PmB0 * kin(0, Cart::y);
        kin(0, Cart::xz) = xi2 * ol(0, Cart::xz) + PmB0 * kin(0, Cart::z);
        kin(0, Cart::yy) =
            xi2 * ol(0, Cart::yy) + PmB1 * kin(0, Cart::y) + term;
        kin(0, Cart::yz) = xi2 * ol(0, Cart::yz) + PmB1 * kin(0, Cart::z);
        kin(0, Cart::zz) =
            xi2 * ol(0, Cart::zz) + PmB2 * kin(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term = fak * kin(i, 0) - fak_b * ol(i, 0);
          kin(i, Cart::xx) = xi2 * ol(i, Cart::xx) + PmB0 * kin(i, Cart::x) +
                             nx[i] * fak * kin(i_less_x[i], Cart::x) + term;
          kin(i, Cart::xy) = xi2 * ol(i, Cart::xy) + PmB0 * kin(i, Cart::y) +
                             nx[i] * fak * kin(i_less_x[i], Cart::y);
          kin(i, Cart::xz) = xi2 * ol(i, Cart::xz) + PmB0 * kin(i, Cart::z) +
                             nx[i] * fak * kin(i_less_x[i], Cart::z);
          kin(i, Cart::yy) = xi2 * ol(i, Cart::yy) + PmB1 * kin(i, Cart::y) +
                             ny[i] * fak * kin(i_less_y[i], Cart::y) + term;
          kin(i, Cart::yz) = xi2 * ol(i, Cart::yz) + PmB1 * kin(i, Cart::z) +
                             ny[i] * fak * kin(i_less_y[i], Cart::z);
          kin(i, Cart::zz) = xi2 * ol(i, Cart::zz) + PmB2 * kin(i, Cart::z) +
                             nz[i] * fak * kin(i_less_z[i], Cart::z) + term;
          //    }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        kin(0, Cart::xxx) =
            xi2 * ol(0, Cart::xxx) + PmB0 * kin(0, Cart::xx) +
            2 * (fak * kin(0, Cart::x) - fak_b * ol(0, Cart::x));
        kin(0, Cart::xxy) = xi2 * ol(0, Cart::xxy) + PmB1 * kin(0, Cart::xx);
        kin(0, Cart::xxz) = xi2 * ol(0, Cart::xxz) + PmB2 * kin(0, Cart::xx);
        kin(0, Cart::xyy) = xi2 * ol(0, Cart::xyy) + PmB0 * kin(0, Cart::yy);
        kin(0, Cart::xyz) = xi2 * ol(0, Cart::xyz) + PmB0 * kin(0, Cart::yz);
        kin(0, Cart::xzz) = xi2 * ol(0, Cart::xzz) + PmB0 * kin(0, Cart::zz);
        kin(0, Cart::yyy) =
            xi2 * ol(0, Cart::yyy) + PmB1 * kin(0, Cart::yy) +
            2 * (fak * kin(0, Cart::y) - fak_b * ol(0, Cart::y));
        kin(0, Cart::yyz) = xi2 * ol(0, Cart::yyz) + PmB2 * kin(0, Cart::yy);
        kin(0, Cart::yzz) = xi2 * ol(0, Cart::yzz) + PmB1 * kin(0, Cart::zz);
        kin(0, Cart::zzz) =
            xi2 * ol(0, Cart::zzz) + PmB2 * kin(0, Cart::zz) +
            2 * (fak * kin(0, Cart::z) - fak_b * ol(0, Cart::z));
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_x = 2 * (fak * kin(i, Cart::x) - fak_b * ol(i, Cart::x));
          double term_y = 2 * (fak * kin(i, Cart::y) - fak_b * ol(i, Cart::y));
          double term_z = 2 * (fak * kin(i, Cart::z) - fak_b * ol(i, Cart::z));
          kin(i, Cart::xxx) = xi2 * ol(i, Cart::xxx) + PmB0 * kin(i, Cart::xx) +
                              nx[i] * fak * kin(i_less_x[i], Cart::xx) + term_x;
          kin(i, Cart::xxy) = xi2 * ol(i, Cart::xxy) + PmB1 * kin(i, Cart::xx) +
                              ny[i] * fak * kin(i_less_y[i], Cart::xx);
          kin(i, Cart::xxz) = xi2 * ol(i, Cart::xxz) + PmB2 * kin(i, Cart::xx) +
                              nz[i] * fak * kin(i_less_z[i], Cart::xx);
          kin(i, Cart::xyy) = xi2 * ol(i, Cart::xyy) + PmB0 * kin(i, Cart::yy) +
                              nx[i] * fak * kin(i_less_x[i], Cart::yy);
          kin(i, Cart::xyz) = xi2 * ol(i, Cart::xyz) + PmB0 * kin(i, Cart::yz) +
                              nx[i] * fak * kin(i_less_x[i], Cart::yz);
          kin(i, Cart::xzz) = xi2 * ol(i, Cart::xzz) + PmB0 * kin(i, Cart::zz) +
                              nx[i] * fak * kin(i_less_x[i], Cart::zz);
          kin(i, Cart::yyy) = xi2 * ol(i, Cart::yyy) + PmB1 * kin(i, Cart::yy) +
                              ny[i] * fak * kin(i_less_y[i], Cart::yy) + term_y;
          kin(i, Cart::yyz) = xi2 * ol(i, Cart::yyz) + PmB2 * kin(i, Cart::yy) +
                              nz[i] * fak * kin(i_less_z[i], Cart::yy);
          kin(i, Cart::yzz) = xi2 * ol(i, Cart::yzz) + PmB1 * kin(i, Cart::zz) +
                              ny[i] * fak * kin(i_less_y[i], Cart::zz);
          kin(i, Cart::zzz) = xi2 * ol(i, Cart::zzz) + PmB2 * kin(i, Cart::zz) +
                              nz[i] * fak * kin(i_less_z[i], Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        double term_xx = fak * kin(0, Cart::xx) - fak_b * ol(0, Cart::xx);
        double term_yy = fak * kin(0, Cart::yy) - fak_b * ol(0, Cart::yy);
        double term_zz = fak * kin(0, Cart::zz) - fak_b * ol(0, Cart::zz);
        kin(0, Cart::xxxx) =
            xi2 * ol(0, Cart::xxxx) + PmB0 * kin(0, Cart::xxx) + 3 * term_xx;
        kin(0, Cart::xxxy) = xi2 * ol(0, Cart::xxxy) + PmB1 * kin(0, Cart::xxx);
        kin(0, Cart::xxxz) = xi2 * ol(0, Cart::xxxz) + PmB2 * kin(0, Cart::xxx);
        kin(0, Cart::xxyy) =
            xi2 * ol(0, Cart::xxyy) + PmB0 * kin(0, Cart::xyy) + term_yy;
        kin(0, Cart::xxyz) = xi2 * ol(0, Cart::xxyz) + PmB1 * kin(0, Cart::xxz);
        kin(0, Cart::xxzz) =
            xi2 * ol(0, Cart::xxzz) + PmB0 * kin(0, Cart::xzz) + term_zz;
        kin(0, Cart::xyyy) = xi2 * ol(0, Cart::xyyy) + PmB0 * kin(0, Cart::yyy);
        kin(0, Cart::xyyz) = xi2 * ol(0, Cart::xyyz) + PmB0 * kin(0, Cart::yyz);
        kin(0, Cart::xyzz) = xi2 * ol(0, Cart::xyzz) + PmB0 * kin(0, Cart::yzz);
        kin(0, Cart::xzzz) = xi2 * ol(0, Cart::xzzz) + PmB0 * kin(0, Cart::zzz);
        kin(0, Cart::yyyy) =
            xi2 * ol(0, Cart::yyyy) + PmB1 * kin(0, Cart::yyy) + 3 * term_yy;
        kin(0, Cart::yyyz) = xi2 * ol(0, Cart::yyyz) + PmB2 * kin(0, Cart::yyy);
        kin(0, Cart::yyzz) =
            xi2 * ol(0, Cart::yyzz) + PmB1 * kin(0, Cart::yzz) + term_zz;
        kin(0, Cart::yzzz) = xi2 * ol(0, Cart::yzzz) + PmB1 * kin(0, Cart::zzz);
        kin(0, Cart::zzzz) =
            xi2 * ol(0, Cart::zzzz) + PmB2 * kin(0, Cart::zzz) + 3 * term_zz;
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xx = fak * kin(i, Cart::xx) - fak_b * ol(i, Cart::xx);
          double term_yy = fak * kin(i, Cart::yy) - fak_b * ol(i, Cart::yy);
          double term_zz = fak * kin(i, Cart::zz) - fak_b * ol(i, Cart::zz);
          kin(i, Cart::xxxx) =
              xi2 * ol(i, Cart::xxxx) + PmB0 * kin(i, Cart::xxx) +
              nx[i] * fak * kin(i_less_x[i], Cart::xxx) + 3 * term_xx;
          kin(i, Cart::xxxy) = xi2 * ol(i, Cart::xxxy) +
                               PmB1 * kin(i, Cart::xxx) +
                               ny[i] * fak * kin(i_less_y[i], Cart::xxx);
          kin(i, Cart::xxxz) = xi2 * ol(i, Cart::xxxz) +
                               PmB2 * kin(i, Cart::xxx) +
                               nz[i] * fak * kin(i_less_z[i], Cart::xxx);
          kin(i, Cart::xxyy) =
              xi2 * ol(i, Cart::xxyy) + PmB0 * kin(i, Cart::xyy) +
              nx[i] * fak * kin(i_less_x[i], Cart::xyy) + term_yy;
          kin(i, Cart::xxyz) = xi2 * ol(i, Cart::xxyz) +
                               PmB1 * kin(i, Cart::xxz) +
                               ny[i] * fak * kin(i_less_y[i], Cart::xxz);
          kin(i, Cart::xxzz) =
              xi2 * ol(i, Cart::xxzz) + PmB0 * kin(i, Cart::xzz) +
              nx[i] * fak * kin(i_less_x[i], Cart::xzz) + term_zz;
          kin(i, Cart::xyyy) = xi2 * ol(i, Cart::xyyy) +
                               PmB0 * kin(i, Cart::yyy) +
                               nx[i] * fak * kin(i_less_x[i], Cart::yyy);
          kin(i, Cart::xyyz) = xi2 * ol(i, Cart::xyyz) +
                               PmB0 * kin(i, Cart::yyz) +
                               nx[i] * fak * kin(i_less_x[i], Cart::yyz);
          kin(i, Cart::xyzz) = xi2 * ol(i, Cart::xyzz) +
                               PmB0 * kin(i, Cart::yzz) +
                               nx[i] * fak * kin(i_less_x[i], Cart::yzz);
          kin(i, Cart::xzzz) = xi2 * ol(i, Cart::xzzz) +
                               PmB0 * kin(i, Cart::zzz) +
                               nx[i] * fak * kin(i_less_x[i], Cart::zzz);
          kin(i, Cart::yyyy) =
              xi2 * ol(i, Cart::yyyy) + PmB1 * kin(i, Cart::yyy) +
              ny[i] * fak * kin(i_less_y[i], Cart::yyy) + 3 * term_yy;
          kin(i, Cart::yyyz) = xi2 * ol(i, Cart::yyyz) +
                               PmB2 * kin(i, Cart::yyy) +
                               nz[i] * fak * kin(i_less_z[i], Cart::yyy);
          kin(i, Cart::yyzz) =
              xi2 * ol(i, Cart::yyzz) + PmB1 * kin(i, Cart::yzz) +
              ny[i] * fak * kin(i_less_y[i], Cart::yzz) + term_zz;
          kin(i, Cart::yzzz) = xi2 * ol(i, Cart::yzzz) +
                               PmB1 * kin(i, Cart::zzz) +
                               ny[i] * fak * kin(i_less_y[i], Cart::zzz);
          kin(i, Cart::zzzz) =
              xi2 * ol(i, Cart::zzzz) + PmB2 * kin(i, Cart::zzz) +
              nz[i] * fak * kin(i_less_z[i], Cart::zzz) + 3 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      // normalization and cartesian -> spherical factors
      Eigen::MatrixXd kin_sph =
          getTrafo(*itr).transpose() * kin * getTrafo(*itc);
      // save to matrix

      for (unsigned i = 0; i < matrix.rows(); i++) {
        for (unsigned j = 0; j < matrix.cols(); j++) {
          matrix(i, j) +=
              kin_sph(i + shell_row->getOffset(), j + shell_col->getOffset());
        }
      }

    }  // col
  }    // row
  return;
}
}  // namespace xtp
}  // namespace votca
