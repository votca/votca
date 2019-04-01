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

namespace votca {
namespace xtp {

void AOPlanewave::FillBlock(Eigen::Block<Eigen::MatrixXcd>& matrix,
                            const AOShell* shell_row,
                            const AOShell* shell_col) {

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row->getLmax();
  int lmax_col = shell_col->getLmax();
  // set size of internal block for recursion
  int nrows = this->getBlockSize(lmax_row);
  int ncols = this->getBlockSize(lmax_col);
  if (lmax_col > 6 || lmax_row > 6) {
    std::cerr << "Orbitals higher than i are not yet implemented. This should "
                 "not have happened!"
              << std::flush;
    exit(1);
  }
  // get shell positions
  const tools::vec& pos_row = shell_row->getPos();  // get position R_{i}
  const tools::vec& pos_col = shell_col->getPos();  // get position R_{j}
  const tools::vec diff = pos_row - pos_col;        // get difference r_{ij}
  double distsq = (diff * diff);                    // get |R_{ij}|^2
  // get kvector modulus
  double kmodulus = (_k * _k);  // get |k|^2

  int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

  int nx[] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 4,
              3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 5, 4, 4, 3, 3, 3, 2,
              2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3,
              3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};

  int ny[] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 0,
              1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3,
              2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3,
              2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0};

  int nz[] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0,
              0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 0, 1, 0, 1, 2, 0,
              1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 0, 1, 0, 1, 2, 0,
              1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6};

  int i_less_x[] = {0,  0,  0,  0,  1,  2,  3,  0,  0,  0,  4,  5,  6,  7,
                    8,  9,  0,  0,  0,  0,  10, 11, 12, 13, 14, 15, 16, 17,
                    18, 19, 0,  0,  0,  0,  0,  20, 21, 22, 23, 24, 25, 26,
                    27, 28, 29, 30, 31, 32, 33, 34, 0,  0,  0,  0,  0,  0,
                    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                    49, 50, 51, 52, 53, 54, 55, 0,  0,  0,  0,  0,  0,  0};

  int i_less_y[] = {0,  0,  0,  0,  0,  1,  0,  2,  3,  0,  0,  4,  0,  5,
                    6,  0,  7,  8,  9,  0,  0,  10, 0,  11, 12, 0,  13, 14,
                    15, 0,  16, 17, 18, 19, 0,  0,  20, 0,  21, 22, 0,  23,
                    24, 25, 0,  26, 27, 28, 29, 0,  30, 31, 32, 33, 34, 0,
                    0,  35, 0,  36, 37, 0,  38, 39, 40, 0,  41, 42, 43, 44,
                    0,  45, 46, 47, 48, 49, 0,  50, 51, 52, 53, 54, 55, 0};

  int i_less_z[] = {0,  0,  0,  0,  0,  0,  1,  0,  2,  3,  0,  0,  4,  0,
                    5,  6,  0,  7,  8,  9,  0,  0,  10, 0,  11, 12, 0,  13,
                    14, 15, 0,  16, 17, 18, 19, 0,  0,  20, 0,  21, 22, 0,
                    23, 24, 25, 0,  26, 27, 28, 29, 0,  30, 31, 32, 33, 34,
                    0,  0,  35, 0,  36, 37, 0,  38, 39, 40, 0,  41, 42, 43,
                    44, 0,  45, 46, 47, 48, 49, 0,  50, 51, 52, 53, 54, 55};
  // iterate over Gaussians in this shell_row
  for (AOShell::GaussianIterator itr = shell_row->begin();
       itr != shell_row->end(); ++itr) {
    // iterate over Gaussians in this shell_col
    // get decay constant
    const double decay_row = itr->getDecay();

    for (AOShell::GaussianIterator itc = shell_col->begin();
         itc != shell_col->end(); ++itc) {
      // get decay constant
      const double decay_col = itc->getDecay();

      // some helpers

      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;
      double exparg = fak2 * decay_row * decay_col * distsq;

      // check if distance between postions is big, then skip step

      if (exparg > 30.0) {
        continue;
      }

      // initialize local matrix block for unnormalized cartesians
      Eigen::MatrixXcd olk = Eigen::MatrixXcd::Zero(nrows, ncols);

      typedef std::complex<double> COMPLEX;  // Define an abbreviation for
                                             // complex numbers

      const COMPLEX PmA0(
          fak2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) -
              pos_row.getX(),
          fak * _k.getX());
      const COMPLEX PmA1(
          fak2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) -
              pos_row.getY(),
          fak * _k.getY());
      const COMPLEX PmA2(
          fak2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) -
              pos_row.getZ(),
          fak * _k.getZ());

      const COMPLEX PmB0(
          fak2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) -
              pos_col.getX(),
          fak * _k.getX());
      const COMPLEX PmB1(
          fak2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) -
              pos_col.getY(),
          fak * _k.getY());
      const COMPLEX PmB2(
          fak2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) -
              pos_col.getZ(),
          fak * _k.getZ());

      const COMPLEX cfak(fak, 0.0);
      const COMPLEX cfak2(fak2, 0.0);

      // calculate s-s- overlap matrix element
      COMPLEX ssol(pow(4.0 * decay_row * decay_col, 0.75) * pow(fak2, 1.5) *
                       exp(-exparg),
                   0.0);  // s-s element

      // calculate s-W-s matrix element
      double kdotr_row = (_k * pos_row);
      double kdotr_col = (_k * pos_col);
      COMPLEX kexparg(
          fak2 * (-0.25) * (kmodulus),
          fak2 * (decay_row) * (kdotr_row) + fak2 * (decay_col) * (kdotr_col));

      olk(0, 0) = ssol * (std::exp(kexparg));  // s-W-s element

      // Integral p-W-s
      if (lmax_row > 0) {
        olk(Cart::x, 0) = PmA0 * olk(0, 0);
        olk(Cart::y, 0) = PmA1 * olk(0, 0);
        olk(Cart::z, 0) = PmA2 * olk(0, 0);
      }
      //------------------------------------------------------

      // Integrals     d - W - s
      if (lmax_row > 1) {
        COMPLEX term = (cfak) * (olk(0, 0));
        olk(Cart::xx, 0) = PmA0 * olk(Cart::x, 0) + term;
        olk(Cart::xy, 0) = PmA0 * olk(Cart::y, 0);
        olk(Cart::xz, 0) = PmA0 * olk(Cart::z, 0);
        olk(Cart::yy, 0) = PmA1 * olk(Cart::y, 0) + term;
        olk(Cart::yz, 0) = PmA1 * olk(Cart::z, 0);
        olk(Cart::zz, 0) = PmA2 * olk(Cart::z, 0) + term;
      }
      //------------------------------------------------------
      // Integrals     f - W - s
      if (lmax_row > 2) {
        olk(Cart::xxx, 0) = PmA0 * olk(Cart::xx, 0) + cfak2 * olk(Cart::x, 0);
        olk(Cart::xxy, 0) = PmA1 * olk(Cart::xx, 0);
        olk(Cart::xxz, 0) = PmA2 * olk(Cart::xx, 0);
        olk(Cart::xyy, 0) = PmA0 * olk(Cart::yy, 0);
        olk(Cart::xyz, 0) = PmA0 * olk(Cart::yz, 0);
        olk(Cart::xzz, 0) = PmA0 * olk(Cart::zz, 0);
        olk(Cart::yyy, 0) = PmA1 * olk(Cart::yy, 0) + cfak2 * olk(Cart::y, 0);
        olk(Cart::yyz, 0) = PmA2 * olk(Cart::yy, 0);
        olk(Cart::yzz, 0) = PmA1 * olk(Cart::zz, 0);
        olk(Cart::zzz, 0) = PmA2 * olk(Cart::zz, 0) + cfak2 * olk(Cart::z, 0);
      }
      //------------------------------------------------------
      // Integrals     g - W - s
      if (lmax_row > 3) {
        COMPLEX term_xx = (cfak) * (olk(Cart::xx, 0));
        COMPLEX term_yy = (cfak) * (olk(Cart::yy, 0));
        COMPLEX term_zz = (cfak) * (olk(Cart::zz, 0));
        olk(Cart::xxxx, 0) = PmA0 * olk(Cart::xxx, 0) + 3.0 * term_xx;
        olk(Cart::xxxy, 0) = PmA1 * olk(Cart::xxx, 0);
        olk(Cart::xxxz, 0) = PmA2 * olk(Cart::xxx, 0);
        olk(Cart::xxyy, 0) = PmA0 * olk(Cart::xyy, 0) + term_yy;
        olk(Cart::xxyz, 0) = PmA1 * olk(Cart::xxz, 0);
        olk(Cart::xxzz, 0) = PmA0 * olk(Cart::xzz, 0) + term_zz;
        olk(Cart::xyyy, 0) = PmA0 * olk(Cart::yyy, 0);
        olk(Cart::xyyz, 0) = PmA0 * olk(Cart::yyz, 0);
        olk(Cart::xyzz, 0) = PmA0 * olk(Cart::yzz, 0);
        olk(Cart::xzzz, 0) = PmA0 * olk(Cart::zzz, 0);
        olk(Cart::yyyy, 0) = PmA1 * olk(Cart::yyy, 0) + 3.0 * term_yy;
        olk(Cart::yyyz, 0) = PmA2 * olk(Cart::yyy, 0);
        olk(Cart::yyzz, 0) = PmA1 * olk(Cart::yzz, 0) + term_zz;
        olk(Cart::yzzz, 0) = PmA1 * olk(Cart::zzz, 0);
        olk(Cart::zzzz, 0) = PmA2 * olk(Cart::zzz, 0) + 3.0 * term_zz;
      }
      //------------------------------------------------------
      // Integrals     h - W - s
      if (lmax_row > 4) {
        COMPLEX term_xxx = (cfak) * (olk(Cart::xxx, 0));
        COMPLEX term_yyy = (cfak) * (olk(Cart::yyy, 0));
        COMPLEX term_zzz = (cfak) * (olk(Cart::zzz, 0));
        olk(Cart::xxxxx, 0) = PmA0 * olk(Cart::xxxx, 0) + 4.0 * term_xxx;
        olk(Cart::xxxxy, 0) = PmA1 * olk(Cart::xxxx, 0);
        olk(Cart::xxxxz, 0) = PmA2 * olk(Cart::xxxx, 0);
        olk(Cart::xxxyy, 0) = PmA1 * olk(Cart::xxxy, 0) + term_xxx;
        olk(Cart::xxxyz, 0) = PmA1 * olk(Cart::xxxz, 0);
        olk(Cart::xxxzz, 0) = PmA2 * olk(Cart::xxxz, 0) + term_xxx;
        olk(Cart::xxyyy, 0) = PmA0 * olk(Cart::xyyy, 0) + term_yyy;
        olk(Cart::xxyyz, 0) = PmA2 * olk(Cart::xxyy, 0);
        olk(Cart::xxyzz, 0) = PmA1 * olk(Cart::xxzz, 0);
        olk(Cart::xxzzz, 0) = PmA0 * olk(Cart::xzzz, 0) + term_zzz;
        olk(Cart::xyyyy, 0) = PmA0 * olk(Cart::yyyy, 0);
        olk(Cart::xyyyz, 0) = PmA0 * olk(Cart::yyyz, 0);
        olk(Cart::xyyzz, 0) = PmA0 * olk(Cart::yyzz, 0);
        olk(Cart::xyzzz, 0) = PmA0 * olk(Cart::yzzz, 0);
        olk(Cart::xzzzz, 0) = PmA0 * olk(Cart::zzzz, 0);
        olk(Cart::yyyyy, 0) = PmA1 * olk(Cart::yyyy, 0) + 4.0 * term_yyy;
        olk(Cart::yyyyz, 0) = PmA2 * olk(Cart::yyyy, 0);
        olk(Cart::yyyzz, 0) = PmA2 * olk(Cart::yyyz, 0) + term_yyy;
        olk(Cart::yyzzz, 0) = PmA1 * olk(Cart::yzzz, 0) + term_zzz;
        olk(Cart::yzzzz, 0) = PmA1 * olk(Cart::zzzz, 0);
        olk(Cart::zzzzz, 0) = PmA2 * olk(Cart::zzzz, 0) + 4.0 * term_zzz;
      }
      //------------------------------------------------------
      // Integrals     i -W - s
      if (lmax_row > 5) {
        COMPLEX term_xxxx = (cfak) * (olk(Cart::xxxx, 0));
        COMPLEX term_xyyy = (cfak) * (olk(Cart::xyyy, 0));
        COMPLEX term_xzzz = (cfak) * (olk(Cart::xzzz, 0));
        COMPLEX term_yyyy = (cfak) * (olk(Cart::yyyy, 0));
        COMPLEX term_yyzz = (cfak) * (olk(Cart::yyzz, 0));
        COMPLEX term_yzzz = (cfak) * (olk(Cart::yzzz, 0));
        COMPLEX term_zzzz = (cfak) * (olk(Cart::zzzz, 0));
        olk(Cart::xxxxxx, 0) = PmA0 * olk(Cart::xxxxx, 0) + 5.0 * term_xxxx;
        olk(Cart::xxxxxy, 0) = PmA1 * olk(Cart::xxxxx, 0);
        olk(Cart::xxxxxz, 0) = PmA2 * olk(Cart::xxxxx, 0);
        olk(Cart::xxxxyy, 0) = PmA1 * olk(Cart::xxxxy, 0) + term_xxxx;
        olk(Cart::xxxxyz, 0) = PmA1 * olk(Cart::xxxxz, 0);
        olk(Cart::xxxxzz, 0) = PmA2 * olk(Cart::xxxxz, 0) + term_xxxx;
        olk(Cart::xxxyyy, 0) = PmA0 * olk(Cart::xxyyy, 0) + 2.0 * term_xyyy;
        olk(Cart::xxxyyz, 0) = PmA2 * olk(Cart::xxxyy, 0);
        olk(Cart::xxxyzz, 0) = PmA1 * olk(Cart::xxxzz, 0);
        olk(Cart::xxxzzz, 0) = PmA0 * olk(Cart::xxzzz, 0) + 2.0 * term_xzzz;
        olk(Cart::xxyyyy, 0) = PmA0 * olk(Cart::xyyyy, 0) + term_yyyy;
        olk(Cart::xxyyyz, 0) = PmA2 * olk(Cart::xxyyy, 0);
        olk(Cart::xxyyzz, 0) = PmA0 * olk(Cart::xyyzz, 0) + term_yyzz;
        olk(Cart::xxyzzz, 0) = PmA1 * olk(Cart::xxzzz, 0);
        olk(Cart::xxzzzz, 0) = PmA0 * olk(Cart::xzzzz, 0) + term_zzzz;
        olk(Cart::xyyyyy, 0) = PmA0 * olk(Cart::yyyyy, 0);
        olk(Cart::xyyyyz, 0) = PmA0 * olk(Cart::yyyyz, 0);
        olk(Cart::xyyyzz, 0) = PmA0 * olk(Cart::yyyzz, 0);
        olk(Cart::xyyzzz, 0) = PmA0 * olk(Cart::yyzzz, 0);
        olk(Cart::xyzzzz, 0) = PmA0 * olk(Cart::yzzzz, 0);
        olk(Cart::xzzzzz, 0) = PmA0 * olk(Cart::zzzzz, 0);
        olk(Cart::yyyyyy, 0) = PmA1 * olk(Cart::yyyyy, 0) + 5.0 * term_yyyy;
        olk(Cart::yyyyyz, 0) = PmA2 * olk(Cart::yyyyy, 0);
        olk(Cart::yyyyzz, 0) = PmA2 * olk(Cart::yyyyz, 0) + term_yyyy;
        olk(Cart::yyyzzz, 0) = PmA1 * olk(Cart::yyzzz, 0) + 2.0 * term_yzzz;
        olk(Cart::yyzzzz, 0) = PmA1 * olk(Cart::yzzzz, 0) + term_zzzz;
        olk(Cart::yzzzzz, 0) = PmA1 * olk(Cart::zzzzz, 0);
        olk(Cart::zzzzzz, 0) = PmA2 * olk(Cart::zzzzz, 0) + 5.0 * term_zzzz;
      }
      //------------------------------------------------------
      if (lmax_col > 0) {

        // Integrals     s - W - p
        olk(0, Cart::x) = PmB0 * olk(0, 0);
        olk(0, Cart::y) = PmB1 * olk(0, 0);
        olk(0, Cart::z) = PmB2 * olk(0, 0);
        //------------------------------------------------------

        // Integrals     p - W - p     d - W - p     f - W - p     g - W - p h -
        // W - p     i - W - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          // COMPLEX cnx(nx[i] * fak, 0.0);
          olk(i, Cart::x) =
              PmB0 * olk(i, 0) + double(nx[i]) * cfak * olk(i_less_x[i], 0);
          // COMPLEX cny = (ny[i] * fak, 0.0);
          olk(i, Cart::y) =
              PmB1 * olk(i, 0) + double(ny[i]) * cfak * olk(i_less_y[i], 0);
          // COMPLEX cnz = (nz[i] * fak, 0.0);
          olk(i, Cart::z) =
              PmB2 * olk(i, 0) + double(nz[i]) * cfak * olk(i_less_z[i], 0);
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)
      if (lmax_col > 1) {

        // Integrals     s - W - d
        COMPLEX term = cfak * olk(0, 0);
        olk(0, Cart::xx) = PmB0 * olk(0, Cart::x) + term;
        olk(0, Cart::xy) = PmB0 * olk(0, Cart::y);
        olk(0, Cart::xz) = PmB0 * olk(0, Cart::z);
        olk(0, Cart::yy) = PmB1 * olk(0, Cart::y) + term;
        olk(0, Cart::yz) = PmB1 * olk(0, Cart::z);
        olk(0, Cart::zz) = PmB2 * olk(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - W - d     d - W - d     f - W - d     g - W - d h -
        // W - d     i - W - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          COMPLEX term = cfak * olk(i, 0);
          olk(i, Cart::xx) = PmB0 * olk(i, Cart::x) +
                             double(nx[i]) * cfak * olk(i_less_x[i], Cart::x) +
                             term;
          olk(i, Cart::xy) = PmB0 * olk(i, Cart::y) +
                             double(nx[i]) * cfak * olk(i_less_x[i], Cart::y);
          olk(i, Cart::xz) = PmB0 * olk(i, Cart::z) +
                             double(nx[i]) * cfak * olk(i_less_x[i], Cart::z);
          olk(i, Cart::yy) = PmB1 * olk(i, Cart::y) +
                             double(ny[i]) * cfak * olk(i_less_y[i], Cart::y) +
                             term;
          olk(i, Cart::yz) = PmB1 * olk(i, Cart::z) +
                             double(ny[i]) * cfak * olk(i_less_y[i], Cart::z);
          olk(i, Cart::zz) = PmB2 * olk(i, Cart::z) +
                             double(nz[i]) * cfak * olk(i_less_z[i], Cart::z) +
                             term;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - W - f
        olk(0, Cart::xxx) =
            PmB0 * olk(0, Cart::xx) + 2.0 * cfak * olk(0, Cart::x);
        olk(0, Cart::xxy) = PmB1 * olk(0, Cart::xx);
        olk(0, Cart::xxz) = PmB2 * olk(0, Cart::xx);
        olk(0, Cart::xyy) = PmB0 * olk(0, Cart::yy);
        olk(0, Cart::xyz) = PmB0 * olk(0, Cart::yz);
        olk(0, Cart::xzz) = PmB0 * olk(0, Cart::zz);
        olk(0, Cart::yyy) =
            PmB1 * olk(0, Cart::yy) + 2.0 * cfak * olk(0, Cart::y);
        olk(0, Cart::yyz) = PmB2 * olk(0, Cart::yy);
        olk(0, Cart::yzz) = PmB1 * olk(0, Cart::zz);
        olk(0, Cart::zzz) =
            PmB2 * olk(0, Cart::zz) + 2.0 * cfak * olk(0, Cart::z);
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f     h - f     i - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          COMPLEX term_x = 2.0 * cfak * olk(i, Cart::x);
          COMPLEX term_y = 2.0 * cfak * olk(i, Cart::y);
          COMPLEX term_z = 2.0 * cfak * olk(i, Cart::z);
          olk(i, Cart::xxx) =
              PmB0 * olk(i, Cart::xx) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xx) + term_x;
          olk(i, Cart::xxy) = PmB1 * olk(i, Cart::xx) +
                              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xx);
          olk(i, Cart::xxz) = PmB2 * olk(i, Cart::xx) +
                              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xx);
          olk(i, Cart::xyy) = PmB0 * olk(i, Cart::yy) +
                              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yy);
          olk(i, Cart::xyz) = PmB0 * olk(i, Cart::yz) +
                              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yz);
          olk(i, Cart::xzz) = PmB0 * olk(i, Cart::zz) +
                              double(nx[i]) * cfak * olk(i_less_x[i], Cart::zz);
          olk(i, Cart::yyy) =
              PmB1 * olk(i, Cart::yy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yy) + term_y;
          olk(i, Cart::yyz) = PmB2 * olk(i, Cart::yy) +
                              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yy);
          olk(i, Cart::yzz) = PmB1 * olk(i, Cart::zz) +
                              double(ny[i]) * cfak * olk(i_less_y[i], Cart::zz);
          olk(i, Cart::zzz) =
              PmB2 * olk(i, Cart::zz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - W - g
        COMPLEX term_xx = cfak * olk(0, Cart::xx);
        COMPLEX term_yy = cfak * olk(0, Cart::yy);
        COMPLEX term_zz = cfak * olk(0, Cart::zz);
        olk(0, Cart::xxxx) = PmB0 * olk(0, Cart::xxx) + 3.0 * term_xx;
        olk(0, Cart::xxxy) = PmB1 * olk(0, Cart::xxx);
        olk(0, Cart::xxxz) = PmB2 * olk(0, Cart::xxx);
        olk(0, Cart::xxyy) = PmB0 * olk(0, Cart::xyy) + term_yy;
        olk(0, Cart::xxyz) = PmB1 * olk(0, Cart::xxz);
        olk(0, Cart::xxzz) = PmB0 * olk(0, Cart::xzz) + term_zz;
        olk(0, Cart::xyyy) = PmB0 * olk(0, Cart::yyy);
        olk(0, Cart::xyyz) = PmB0 * olk(0, Cart::yyz);
        olk(0, Cart::xyzz) = PmB0 * olk(0, Cart::yzz);
        olk(0, Cart::xzzz) = PmB0 * olk(0, Cart::zzz);
        olk(0, Cart::yyyy) = PmB1 * olk(0, Cart::yyy) + 3.0 * term_yy;
        olk(0, Cart::yyyz) = PmB2 * olk(0, Cart::yyy);
        olk(0, Cart::yyzz) = PmB1 * olk(0, Cart::yzz) + term_zz;
        olk(0, Cart::yzzz) = PmB1 * olk(0, Cart::zzz);
        olk(0, Cart::zzzz) = PmB2 * olk(0, Cart::zzz) + 3.0 * term_zz;
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g     h - g     i - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          COMPLEX term_xx = cfak * olk(i, Cart::xx);
          COMPLEX term_yy = cfak * olk(i, Cart::yy);
          COMPLEX term_zz = cfak * olk(i, Cart::zz);
          olk(i, Cart::xxxx) =
              PmB0 * olk(i, Cart::xxx) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xxx) +
              3.0 * term_xx;
          olk(i, Cart::xxxy) =
              PmB1 * olk(i, Cart::xxx) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxx);
          olk(i, Cart::xxxz) =
              PmB2 * olk(i, Cart::xxx) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxx);
          olk(i, Cart::xxyy) =
              PmB0 * olk(i, Cart::xyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xyy) + term_yy;
          olk(i, Cart::xxyz) =
              PmB1 * olk(i, Cart::xxz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxz);
          olk(i, Cart::xxzz) =
              PmB0 * olk(i, Cart::xzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xzz) + term_zz;
          olk(i, Cart::xyyy) =
              PmB0 * olk(i, Cart::yyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyy);
          olk(i, Cart::xyyz) =
              PmB0 * olk(i, Cart::yyz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyz);
          olk(i, Cart::xyzz) =
              PmB0 * olk(i, Cart::yzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yzz);
          olk(i, Cart::xzzz) =
              PmB0 * olk(i, Cart::zzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::zzz);
          olk(i, Cart::yyyy) =
              PmB1 * olk(i, Cart::yyy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yyy) +
              3.0 * term_yy;
          olk(i, Cart::yyyz) =
              PmB2 * olk(i, Cart::yyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yyy);
          olk(i, Cart::yyzz) =
              PmB1 * olk(i, Cart::yzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yzz) + term_zz;
          olk(i, Cart::yzzz) =
              PmB1 * olk(i, Cart::zzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::zzz);
          olk(i, Cart::zzzz) =
              PmB2 * olk(i, Cart::zzz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::zzz) +
              3.0 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      if (lmax_col > 4) {

        // Integrals     s - h
        COMPLEX term_xxx = cfak * olk(0, Cart::xxx);
        COMPLEX term_yyy = cfak * olk(0, Cart::yyy);
        COMPLEX term_zzz = cfak * olk(0, Cart::zzz);
        olk(0, Cart::xxxxx) = PmB0 * olk(0, Cart::xxxx) + 4.0 * term_xxx;
        olk(0, Cart::xxxxy) = PmB1 * olk(0, Cart::xxxx);
        olk(0, Cart::xxxxz) = PmB2 * olk(0, Cart::xxxx);
        olk(0, Cart::xxxyy) = PmB1 * olk(0, Cart::xxxy) + term_xxx;
        olk(0, Cart::xxxyz) = PmB1 * olk(0, Cart::xxxz);
        olk(0, Cart::xxxzz) = PmB2 * olk(0, Cart::xxxz) + term_xxx;
        olk(0, Cart::xxyyy) = PmB0 * olk(0, Cart::xyyy) + term_yyy;
        olk(0, Cart::xxyyz) = PmB2 * olk(0, Cart::xxyy);
        olk(0, Cart::xxyzz) = PmB1 * olk(0, Cart::xxzz);
        olk(0, Cart::xxzzz) = PmB0 * olk(0, Cart::xzzz) + term_zzz;
        olk(0, Cart::xyyyy) = PmB0 * olk(0, Cart::yyyy);
        olk(0, Cart::xyyyz) = PmB0 * olk(0, Cart::yyyz);
        olk(0, Cart::xyyzz) = PmB0 * olk(0, Cart::yyzz);
        olk(0, Cart::xyzzz) = PmB0 * olk(0, Cart::yzzz);
        olk(0, Cart::xzzzz) = PmB0 * olk(0, Cart::zzzz);
        olk(0, Cart::yyyyy) = PmB1 * olk(0, Cart::yyyy) + 4.0 * term_yyy;
        olk(0, Cart::yyyyz) = PmB2 * olk(0, Cart::yyyy);
        olk(0, Cart::yyyzz) = PmB2 * olk(0, Cart::yyyz) + term_yyy;
        olk(0, Cart::yyzzz) = PmB1 * olk(0, Cart::yzzz) + term_zzz;
        olk(0, Cart::yzzzz) = PmB1 * olk(0, Cart::zzzz);
        olk(0, Cart::zzzzz) = PmB2 * olk(0, Cart::zzzz) + 4.0 * term_zzz;
        //------------------------------------------------------

        // Integrals     p - h     d - h     f - h     g - h     h - h     i - h
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          COMPLEX term_xxx = cfak * olk(i, Cart::xxx);
          COMPLEX term_yyy = cfak * olk(i, Cart::yyy);
          COMPLEX term_zzz = cfak * olk(i, Cart::zzz);
          olk(i, Cart::xxxxx) =
              PmB0 * olk(i, Cart::xxxx) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xxxx) +
              4.0 * term_xxx;
          olk(i, Cart::xxxxy) =
              PmB1 * olk(i, Cart::xxxx) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxx);
          olk(i, Cart::xxxxz) =
              PmB2 * olk(i, Cart::xxxx) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxxx);
          olk(i, Cart::xxxyy) =
              PmB1 * olk(i, Cart::xxxy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxy) + term_xxx;
          olk(i, Cart::xxxyz) =
              PmB1 * olk(i, Cart::xxxz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxz);
          olk(i, Cart::xxxzz) =
              PmB2 * olk(i, Cart::xxxz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxxz) + term_xxx;
          olk(i, Cart::xxyyy) =
              PmB0 * olk(i, Cart::xyyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xyyy) + term_yyy;
          olk(i, Cart::xxyyz) =
              PmB2 * olk(i, Cart::xxyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxyy);
          olk(i, Cart::xxyzz) =
              PmB1 * olk(i, Cart::xxzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxzz);
          olk(i, Cart::xxzzz) =
              PmB0 * olk(i, Cart::xzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xzzz) + term_zzz;
          olk(i, Cart::xyyyy) =
              PmB0 * olk(i, Cart::yyyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyyy);
          olk(i, Cart::xyyyz) =
              PmB0 * olk(i, Cart::yyyz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyyz);
          olk(i, Cart::xyyzz) =
              PmB0 * olk(i, Cart::yyzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyzz);
          olk(i, Cart::xyzzz) =
              PmB0 * olk(i, Cart::yzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yzzz);
          olk(i, Cart::xzzzz) =
              PmB0 * olk(i, Cart::zzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::zzzz);
          olk(i, Cart::yyyyy) =
              PmB1 * olk(i, Cart::yyyy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yyyy) +
              4.0 * term_yyy;
          olk(i, Cart::yyyyz) =
              PmB2 * olk(i, Cart::yyyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yyyy);
          olk(i, Cart::yyyzz) =
              PmB2 * olk(i, Cart::yyyz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yyyz) + term_yyy;
          olk(i, Cart::yyzzz) =
              PmB1 * olk(i, Cart::yzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yzzz) + term_zzz;
          olk(i, Cart::yzzzz) =
              PmB1 * olk(i, Cart::zzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::zzzz);
          olk(i, Cart::zzzzz) =
              PmB2 * olk(i, Cart::zzzz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::zzzz) +
              4.0 * term_zzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 4)

      if (lmax_col > 5) {

        // Integrals     s - W -i
        COMPLEX term_xxxx = cfak * olk(0, Cart::xxxx);
        COMPLEX term_xyyy = cfak * olk(0, Cart::xyyy);
        COMPLEX term_xzzz = cfak * olk(0, Cart::xzzz);
        COMPLEX term_yyyy = cfak * olk(0, Cart::yyyy);
        COMPLEX term_yyzz = cfak * olk(0, Cart::yyzz);
        COMPLEX term_yzzz = cfak * olk(0, Cart::yzzz);
        COMPLEX term_zzzz = cfak * olk(0, Cart::zzzz);
        olk(0, Cart::xxxxxx) = PmB0 * olk(0, Cart::xxxxx) + 5.0 * term_xxxx;
        olk(0, Cart::xxxxxy) = PmB1 * olk(0, Cart::xxxxx);
        olk(0, Cart::xxxxxz) = PmB2 * olk(0, Cart::xxxxx);
        olk(0, Cart::xxxxyy) = PmB1 * olk(0, Cart::xxxxy) + term_xxxx;
        olk(0, Cart::xxxxyz) = PmB1 * olk(0, Cart::xxxxz);
        olk(0, Cart::xxxxzz) = PmB2 * olk(0, Cart::xxxxz) + term_xxxx;
        olk(0, Cart::xxxyyy) = PmB0 * olk(0, Cart::xxyyy) + 2.0 * term_xyyy;
        olk(0, Cart::xxxyyz) = PmB2 * olk(0, Cart::xxxyy);
        olk(0, Cart::xxxyzz) = PmB1 * olk(0, Cart::xxxzz);
        olk(0, Cart::xxxzzz) = PmB0 * olk(0, Cart::xxzzz) + 2.0 * term_xzzz;
        olk(0, Cart::xxyyyy) = PmB0 * olk(0, Cart::xyyyy) + term_yyyy;
        olk(0, Cart::xxyyyz) = PmB2 * olk(0, Cart::xxyyy);
        olk(0, Cart::xxyyzz) = PmB0 * olk(0, Cart::xyyzz) + term_yyzz;
        olk(0, Cart::xxyzzz) = PmB1 * olk(0, Cart::xxzzz);
        olk(0, Cart::xxzzzz) = PmB0 * olk(0, Cart::xzzzz) + term_zzzz;
        olk(0, Cart::xyyyyy) = PmB0 * olk(0, Cart::yyyyy);
        olk(0, Cart::xyyyyz) = PmB0 * olk(0, Cart::yyyyz);
        olk(0, Cart::xyyyzz) = PmB0 * olk(0, Cart::yyyzz);
        olk(0, Cart::xyyzzz) = PmB0 * olk(0, Cart::yyzzz);
        olk(0, Cart::xyzzzz) = PmB0 * olk(0, Cart::yzzzz);
        olk(0, Cart::xzzzzz) = PmB0 * olk(0, Cart::zzzzz);
        olk(0, Cart::yyyyyy) = PmB1 * olk(0, Cart::yyyyy) + 5.0 * term_yyyy;
        olk(0, Cart::yyyyyz) = PmB2 * olk(0, Cart::yyyyy);
        olk(0, Cart::yyyyzz) = PmB2 * olk(0, Cart::yyyyz) + term_yyyy;
        olk(0, Cart::yyyzzz) = PmB1 * olk(0, Cart::yyzzz) + 2.0 * term_yzzz;
        olk(0, Cart::yyzzzz) = PmB1 * olk(0, Cart::yzzzz) + term_zzzz;
        olk(0, Cart::yzzzzz) = PmB1 * olk(0, Cart::zzzzz);
        olk(0, Cart::zzzzzz) = PmB2 * olk(0, Cart::zzzzz) + 5.0 * term_zzzz;
        //------------------------------------------------------

        // Integrals     p - W - i     d - W - i     f - W - i     g - W -i h -
        // W - i     i - W - i
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          COMPLEX term_xxxx = cfak * olk(i, Cart::xxxx);
          COMPLEX term_xyyy = cfak * olk(i, Cart::xyyy);
          COMPLEX term_xzzz = cfak * olk(i, Cart::xzzz);
          COMPLEX term_yyyy = cfak * olk(i, Cart::yyyy);
          COMPLEX term_yyzz = cfak * olk(i, Cart::yyzz);
          COMPLEX term_yzzz = cfak * olk(i, Cart::yzzz);
          COMPLEX term_zzzz = cfak * olk(i, Cart::zzzz);
          olk(i, Cart::xxxxxx) =
              PmB0 * olk(i, Cart::xxxxx) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xxxxx) +
              5.0 * term_xxxx;
          olk(i, Cart::xxxxxy) =
              PmB1 * olk(i, Cart::xxxxx) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxxx);
          olk(i, Cart::xxxxxz) =
              PmB2 * olk(i, Cart::xxxxx) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxxxx);
          olk(i, Cart::xxxxyy) =
              PmB1 * olk(i, Cart::xxxxy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxxy) + term_xxxx;
          olk(i, Cart::xxxxyz) =
              PmB1 * olk(i, Cart::xxxxz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxxz);
          olk(i, Cart::xxxxzz) =
              PmB2 * olk(i, Cart::xxxxz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxxxz) + term_xxxx;
          olk(i, Cart::xxxyyy) =
              PmB0 * olk(i, Cart::xxyyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xxyyy) +
              2.0 * term_xyyy;
          olk(i, Cart::xxxyyz) =
              PmB2 * olk(i, Cart::xxxyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxxyy);
          olk(i, Cart::xxxyzz) =
              PmB1 * olk(i, Cart::xxxzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxxzz);
          olk(i, Cart::xxxzzz) =
              PmB0 * olk(i, Cart::xxzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xxzzz) +
              2.0 * term_xzzz;
          olk(i, Cart::xxyyyy) =
              PmB0 * olk(i, Cart::xyyyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xyyyy) + term_yyyy;
          olk(i, Cart::xxyyyz) =
              PmB2 * olk(i, Cart::xxyyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::xxyyy);
          olk(i, Cart::xxyyzz) =
              PmB0 * olk(i, Cart::xyyzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xyyzz) + term_yyzz;
          olk(i, Cart::xxyzzz) =
              PmB1 * olk(i, Cart::xxzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::xxzzz);
          olk(i, Cart::xxzzzz) =
              PmB0 * olk(i, Cart::xzzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::xzzzz) + term_zzzz;
          olk(i, Cart::xyyyyy) =
              PmB0 * olk(i, Cart::yyyyy) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyyyy);
          olk(i, Cart::xyyyyz) =
              PmB0 * olk(i, Cart::yyyyz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyyyz);
          olk(i, Cart::xyyyzz) =
              PmB0 * olk(i, Cart::yyyzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyyzz);
          olk(i, Cart::xyyzzz) =
              PmB0 * olk(i, Cart::yyzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yyzzz);
          olk(i, Cart::xyzzzz) =
              PmB0 * olk(i, Cart::yzzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::yzzzz);
          olk(i, Cart::xzzzzz) =
              PmB0 * olk(i, Cart::zzzzz) +
              double(nx[i]) * cfak * olk(i_less_x[i], Cart::zzzzz);
          olk(i, Cart::yyyyyy) =
              PmB1 * olk(i, Cart::yyyyy) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yyyyy) +
              5.0 * term_yyyy;
          olk(i, Cart::yyyyyz) =
              PmB2 * olk(i, Cart::yyyyy) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yyyyy);
          olk(i, Cart::yyyyzz) =
              PmB2 * olk(i, Cart::yyyyz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::yyyyz) + term_yyyy;
          olk(i, Cart::yyyzzz) =
              PmB1 * olk(i, Cart::yyzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yyzzz) +
              2.0 * term_yzzz;
          olk(i, Cart::yyzzzz) =
              PmB1 * olk(i, Cart::yzzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::yzzzz) + term_zzzz;
          olk(i, Cart::yzzzzz) =
              PmB1 * olk(i, Cart::zzzzz) +
              double(ny[i]) * cfak * olk(i_less_y[i], Cart::zzzzz);
          olk(i, Cart::zzzzzz) =
              PmB2 * olk(i, Cart::zzzzz) +
              double(nz[i]) * cfak * olk(i_less_z[i], Cart::zzzzz) +
              5.0 * term_zzzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 5)

      Eigen::MatrixXd trafo_row = getTrafo(*itr);
      Eigen::MatrixXd trafo_col = getTrafo(*itc);

      // cartesian -> spherical
      Eigen::MatrixXcd olk_sph = trafo_row.transpose() * olk * trafo_col;

      // save to matrix
      for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
          matrix(i, j) +=
              olk_sph(i + shell_row->getOffset(), j + shell_col->getOffset());
        }
      }

    }  // close Gaussian shell_col

  }  // close Gaussian shell_row

}  // End AOPlanewave

void AOPlanewave::Fillextpotential(const AOBasis& aobasis,
                                   const std::vector<tools::vec>& kpoints) {

  _externalpotential =
      Eigen::MatrixXcd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

  for (const auto& kpoint : kpoints) {
    _aomatrix =
        Eigen::MatrixXcd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
    setkVector(kpoint);
    Fill(aobasis);
    _externalpotential += _aomatrix;
  }

  return;
}

}  // namespace xtp
}  // namespace votca
