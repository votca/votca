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
#include <votca/xtp/aotransform.h>

namespace votca {
namespace xtp {

void AOOverlap::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                          const AOShell& shell_row,
                          const AOShell& shell_col) const {

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row.getLmax();
  int lmax_col = shell_col.getLmax();

  // set size of internal block for recursion
  int nrows = AOTransform::getBlockSize(lmax_row);
  int ncols = AOTransform::getBlockSize(lmax_col);

  if (lmax_col > 6 || lmax_row > 6) {
    throw std::runtime_error(
        "Orbitals higher than i are not yet implemented. This should not have "
        "happened!");
  }

  /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
   * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION
   * COEFFICIENTS, AND ADD TO matrix(i,j)
   */

  // get shell positions
  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;

  double distsq = diff.squaredNorm();
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
  for (const auto& gaussian_row : shell_row) {
    // iterate over Gaussians in this shell_col
    const double decay_row = gaussian_row.getDecay();

    for (const auto& gaussian_col : shell_col) {

      const double decay_col = gaussian_col.getDecay();

      // some helpers
      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;

      // check if distance between postions is big, then skip step
      double exparg = fak2 * decay_row * decay_col * distsq;
      if (exparg > 30.0) {
        continue;
      }
      // initialize local matrix block for unnormalized cartesians
      Eigen::MatrixXd ol = Eigen::MatrixXd::Zero(nrows, ncols);

      // Definition of coefficients for recursive overlap formulas
      // A for rows (i). B for columns (j)
      const Eigen::Vector3d PmA =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_row;
      const Eigen::Vector3d PmB =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_col;

      // calculate matrix elements
      ol(0, 0) = pow(4.0 * decay_row * decay_col, 0.75) * pow(fak2, 1.5) *
                 exp(-exparg);  // s-s element

      // Integrals     p - s
      if (lmax_row > 0) {
        ol(Cart::x, 0) = PmA(0) * ol(0, 0);
        ol(Cart::y, 0) = PmA(1) * ol(0, 0);
        ol(Cart::z, 0) = PmA(2) * ol(0, 0);
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
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
      if (lmax_row > 2) {
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
      if (lmax_row > 3) {
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
      if (lmax_row > 4) {
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

      // Integrals     i - s
      if (lmax_row > 5) {
        double term_xxxx = fak * ol(Cart::xxxx, 0);
        double term_xyyy = fak * ol(Cart::xyyy, 0);
        double term_xzzz = fak * ol(Cart::xzzz, 0);
        double term_yyyy = fak * ol(Cart::yyyy, 0);
        double term_yyzz = fak * ol(Cart::yyzz, 0);
        double term_yzzz = fak * ol(Cart::yzzz, 0);
        double term_zzzz = fak * ol(Cart::zzzz, 0);
        ol(Cart::xxxxxx, 0) = PmA(0) * ol(Cart::xxxxx, 0) + 5 * term_xxxx;
        ol(Cart::xxxxxy, 0) = PmA(1) * ol(Cart::xxxxx, 0);
        ol(Cart::xxxxxz, 0) = PmA(2) * ol(Cart::xxxxx, 0);
        ol(Cart::xxxxyy, 0) = PmA(1) * ol(Cart::xxxxy, 0) + term_xxxx;
        ol(Cart::xxxxyz, 0) = PmA(1) * ol(Cart::xxxxz, 0);
        ol(Cart::xxxxzz, 0) = PmA(2) * ol(Cart::xxxxz, 0) + term_xxxx;
        ol(Cart::xxxyyy, 0) = PmA(0) * ol(Cart::xxyyy, 0) + 2 * term_xyyy;
        ol(Cart::xxxyyz, 0) = PmA(2) * ol(Cart::xxxyy, 0);
        ol(Cart::xxxyzz, 0) = PmA(1) * ol(Cart::xxxzz, 0);
        ol(Cart::xxxzzz, 0) = PmA(0) * ol(Cart::xxzzz, 0) + 2 * term_xzzz;
        ol(Cart::xxyyyy, 0) = PmA(0) * ol(Cart::xyyyy, 0) + term_yyyy;
        ol(Cart::xxyyyz, 0) = PmA(2) * ol(Cart::xxyyy, 0);
        ol(Cart::xxyyzz, 0) = PmA(0) * ol(Cart::xyyzz, 0) + term_yyzz;
        ol(Cart::xxyzzz, 0) = PmA(1) * ol(Cart::xxzzz, 0);
        ol(Cart::xxzzzz, 0) = PmA(0) * ol(Cart::xzzzz, 0) + term_zzzz;
        ol(Cart::xyyyyy, 0) = PmA(0) * ol(Cart::yyyyy, 0);
        ol(Cart::xyyyyz, 0) = PmA(0) * ol(Cart::yyyyz, 0);
        ol(Cart::xyyyzz, 0) = PmA(0) * ol(Cart::yyyzz, 0);
        ol(Cart::xyyzzz, 0) = PmA(0) * ol(Cart::yyzzz, 0);
        ol(Cart::xyzzzz, 0) = PmA(0) * ol(Cart::yzzzz, 0);
        ol(Cart::xzzzzz, 0) = PmA(0) * ol(Cart::zzzzz, 0);
        ol(Cart::yyyyyy, 0) = PmA(1) * ol(Cart::yyyyy, 0) + 5 * term_yyyy;
        ol(Cart::yyyyyz, 0) = PmA(2) * ol(Cart::yyyyy, 0);
        ol(Cart::yyyyzz, 0) = PmA(2) * ol(Cart::yyyyz, 0) + term_yyyy;
        ol(Cart::yyyzzz, 0) = PmA(1) * ol(Cart::yyzzz, 0) + 2 * term_yzzz;
        ol(Cart::yyzzzz, 0) = PmA(1) * ol(Cart::yzzzz, 0) + term_zzzz;
        ol(Cart::yzzzzz, 0) = PmA(1) * ol(Cart::zzzzz, 0);
        ol(Cart::zzzzzz, 0) = PmA(2) * ol(Cart::zzzzz, 0) + 5 * term_zzzz;
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        ol(0, Cart::x) = PmB(0) * ol(0, 0);
        ol(0, Cart::y) = PmB(1) * ol(0, 0);
        ol(0, Cart::z) = PmB(2) * ol(0, 0);
        //------------------------------------------------------

        // Integrals     p - p     d - p     f - p     g - p     h - p     i - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          ol(i, Cart::x) = PmB(0) * ol(i, 0) + nx[i] * fak * ol(i_less_x[i], 0);
          ol(i, Cart::y) = PmB(1) * ol(i, 0) + ny[i] * fak * ol(i_less_y[i], 0);
          ol(i, Cart::z) = PmB(2) * ol(i, 0) + nz[i] * fak * ol(i_less_z[i], 0);
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        double term = fak * ol(0, 0);
        ol(0, Cart::xx) = PmB(0) * ol(0, Cart::x) + term;
        ol(0, Cart::xy) = PmB(0) * ol(0, Cart::y);
        ol(0, Cart::xz) = PmB(0) * ol(0, Cart::z);
        ol(0, Cart::yy) = PmB(1) * ol(0, Cart::y) + term;
        ol(0, Cart::yz) = PmB(1) * ol(0, Cart::z);
        ol(0, Cart::zz) = PmB(2) * ol(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d     h - d     i - d
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

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

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

        // Integrals     p - f     d - f     f - f     g - f     h - f     i - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_x = 2 * fak * ol(i, Cart::x);
          double term_y = 2 * fak * ol(i, Cart::y);
          double term_z = 2 * fak * ol(i, Cart::z);
          ol(i, Cart::xxx) = PmB(0) * ol(i, Cart::xx) +
                             nx[i] * fak * ol(i_less_x[i], Cart::xx) + term_x;
          ol(i, Cart::xxy) = PmB(1) * ol(i, Cart::xx) +
                             ny[i] * fak * ol(i_less_y[i], Cart::xx);
          ol(i, Cart::xxz) = PmB(2) * ol(i, Cart::xx) +
                             nz[i] * fak * ol(i_less_z[i], Cart::xx);
          ol(i, Cart::xyy) = PmB(0) * ol(i, Cart::yy) +
                             nx[i] * fak * ol(i_less_x[i], Cart::yy);
          ol(i, Cart::xyz) = PmB(0) * ol(i, Cart::yz) +
                             nx[i] * fak * ol(i_less_x[i], Cart::yz);
          ol(i, Cart::xzz) = PmB(0) * ol(i, Cart::zz) +
                             nx[i] * fak * ol(i_less_x[i], Cart::zz);
          ol(i, Cart::yyy) = PmB(1) * ol(i, Cart::yy) +
                             ny[i] * fak * ol(i_less_y[i], Cart::yy) + term_y;
          ol(i, Cart::yyz) = PmB(2) * ol(i, Cart::yy) +
                             nz[i] * fak * ol(i_less_z[i], Cart::yy);
          ol(i, Cart::yzz) = PmB(1) * ol(i, Cart::zz) +
                             ny[i] * fak * ol(i_less_y[i], Cart::zz);
          ol(i, Cart::zzz) = PmB(2) * ol(i, Cart::zz) +
                             nz[i] * fak * ol(i_less_z[i], Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

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

        // Integrals     p - g     d - g     f - g     g - g     h - g     i - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xx = fak * ol(i, Cart::xx);
          double term_yy = fak * ol(i, Cart::yy);
          double term_zz = fak * ol(i, Cart::zz);
          ol(i, Cart::xxxx) = PmB(0) * ol(i, Cart::xxx) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xxx) +
                              3 * term_xx;
          ol(i, Cart::xxxy) = PmB(1) * ol(i, Cart::xxx) +
                              ny[i] * fak * ol(i_less_y[i], Cart::xxx);
          ol(i, Cart::xxxz) = PmB(2) * ol(i, Cart::xxx) +
                              nz[i] * fak * ol(i_less_z[i], Cart::xxx);
          ol(i, Cart::xxyy) = PmB(0) * ol(i, Cart::xyy) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xyy) +
                              term_yy;
          ol(i, Cart::xxyz) = PmB(1) * ol(i, Cart::xxz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::xxz);
          ol(i, Cart::xxzz) = PmB(0) * ol(i, Cart::xzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::xzz) +
                              term_zz;
          ol(i, Cart::xyyy) = PmB(0) * ol(i, Cart::yyy) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yyy);
          ol(i, Cart::xyyz) = PmB(0) * ol(i, Cart::yyz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yyz);
          ol(i, Cart::xyzz) = PmB(0) * ol(i, Cart::yzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::yzz);
          ol(i, Cart::xzzz) = PmB(0) * ol(i, Cart::zzz) +
                              nx[i] * fak * ol(i_less_x[i], Cart::zzz);
          ol(i, Cart::yyyy) = PmB(1) * ol(i, Cart::yyy) +
                              ny[i] * fak * ol(i_less_y[i], Cart::yyy) +
                              3 * term_yy;
          ol(i, Cart::yyyz) = PmB(2) * ol(i, Cart::yyy) +
                              nz[i] * fak * ol(i_less_z[i], Cart::yyy);
          ol(i, Cart::yyzz) = PmB(1) * ol(i, Cart::yzz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::yzz) +
                              term_zz;
          ol(i, Cart::yzzz) = PmB(1) * ol(i, Cart::zzz) +
                              ny[i] * fak * ol(i_less_y[i], Cart::zzz);
          ol(i, Cart::zzzz) = PmB(2) * ol(i, Cart::zzz) +
                              nz[i] * fak * ol(i_less_z[i], Cart::zzz) +
                              3 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      if (lmax_col > 4) {

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

        // Integrals     p - h     d - h     f - h     g - h     h - h     i - h
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xxx = fak * ol(i, Cart::xxx);
          double term_yyy = fak * ol(i, Cart::yyy);
          double term_zzz = fak * ol(i, Cart::zzz);
          ol(i, Cart::xxxxx) = PmB(0) * ol(i, Cart::xxxx) +
                               nx[i] * fak * ol(i_less_x[i], Cart::xxxx) +
                               4 * term_xxx;
          ol(i, Cart::xxxxy) = PmB(1) * ol(i, Cart::xxxx) +
                               ny[i] * fak * ol(i_less_y[i], Cart::xxxx);
          ol(i, Cart::xxxxz) = PmB(2) * ol(i, Cart::xxxx) +
                               nz[i] * fak * ol(i_less_z[i], Cart::xxxx);
          ol(i, Cart::xxxyy) = PmB(1) * ol(i, Cart::xxxy) +
                               ny[i] * fak * ol(i_less_y[i], Cart::xxxy) +
                               term_xxx;
          ol(i, Cart::xxxyz) = PmB(1) * ol(i, Cart::xxxz) +
                               ny[i] * fak * ol(i_less_y[i], Cart::xxxz);
          ol(i, Cart::xxxzz) = PmB(2) * ol(i, Cart::xxxz) +
                               nz[i] * fak * ol(i_less_z[i], Cart::xxxz) +
                               term_xxx;
          ol(i, Cart::xxyyy) = PmB(0) * ol(i, Cart::xyyy) +
                               nx[i] * fak * ol(i_less_x[i], Cart::xyyy) +
                               term_yyy;
          ol(i, Cart::xxyyz) = PmB(2) * ol(i, Cart::xxyy) +
                               nz[i] * fak * ol(i_less_z[i], Cart::xxyy);
          ol(i, Cart::xxyzz) = PmB(1) * ol(i, Cart::xxzz) +
                               ny[i] * fak * ol(i_less_y[i], Cart::xxzz);
          ol(i, Cart::xxzzz) = PmB(0) * ol(i, Cart::xzzz) +
                               nx[i] * fak * ol(i_less_x[i], Cart::xzzz) +
                               term_zzz;
          ol(i, Cart::xyyyy) = PmB(0) * ol(i, Cart::yyyy) +
                               nx[i] * fak * ol(i_less_x[i], Cart::yyyy);
          ol(i, Cart::xyyyz) = PmB(0) * ol(i, Cart::yyyz) +
                               nx[i] * fak * ol(i_less_x[i], Cart::yyyz);
          ol(i, Cart::xyyzz) = PmB(0) * ol(i, Cart::yyzz) +
                               nx[i] * fak * ol(i_less_x[i], Cart::yyzz);
          ol(i, Cart::xyzzz) = PmB(0) * ol(i, Cart::yzzz) +
                               nx[i] * fak * ol(i_less_x[i], Cart::yzzz);
          ol(i, Cart::xzzzz) = PmB(0) * ol(i, Cart::zzzz) +
                               nx[i] * fak * ol(i_less_x[i], Cart::zzzz);
          ol(i, Cart::yyyyy) = PmB(1) * ol(i, Cart::yyyy) +
                               ny[i] * fak * ol(i_less_y[i], Cart::yyyy) +
                               4 * term_yyy;
          ol(i, Cart::yyyyz) = PmB(2) * ol(i, Cart::yyyy) +
                               nz[i] * fak * ol(i_less_z[i], Cart::yyyy);
          ol(i, Cart::yyyzz) = PmB(2) * ol(i, Cart::yyyz) +
                               nz[i] * fak * ol(i_less_z[i], Cart::yyyz) +
                               term_yyy;
          ol(i, Cart::yyzzz) = PmB(1) * ol(i, Cart::yzzz) +
                               ny[i] * fak * ol(i_less_y[i], Cart::yzzz) +
                               term_zzz;
          ol(i, Cart::yzzzz) = PmB(1) * ol(i, Cart::zzzz) +
                               ny[i] * fak * ol(i_less_y[i], Cart::zzzz);
          ol(i, Cart::zzzzz) = PmB(2) * ol(i, Cart::zzzz) +
                               nz[i] * fak * ol(i_less_z[i], Cart::zzzz) +
                               4 * term_zzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 4)

      if (lmax_col > 5) {

        // Integrals     s - i
        double term_xxxx = fak * ol(0, Cart::xxxx);
        double term_xyyy = fak * ol(0, Cart::xyyy);
        double term_xzzz = fak * ol(0, Cart::xzzz);
        double term_yyyy = fak * ol(0, Cart::yyyy);
        double term_yyzz = fak * ol(0, Cart::yyzz);
        double term_yzzz = fak * ol(0, Cart::yzzz);
        double term_zzzz = fak * ol(0, Cart::zzzz);
        ol(0, Cart::xxxxxx) = PmB(0) * ol(0, Cart::xxxxx) + 5 * term_xxxx;
        ol(0, Cart::xxxxxy) = PmB(1) * ol(0, Cart::xxxxx);
        ol(0, Cart::xxxxxz) = PmB(2) * ol(0, Cart::xxxxx);
        ol(0, Cart::xxxxyy) = PmB(1) * ol(0, Cart::xxxxy) + term_xxxx;
        ol(0, Cart::xxxxyz) = PmB(1) * ol(0, Cart::xxxxz);
        ol(0, Cart::xxxxzz) = PmB(2) * ol(0, Cart::xxxxz) + term_xxxx;
        ol(0, Cart::xxxyyy) = PmB(0) * ol(0, Cart::xxyyy) + 2 * term_xyyy;
        ol(0, Cart::xxxyyz) = PmB(2) * ol(0, Cart::xxxyy);
        ol(0, Cart::xxxyzz) = PmB(1) * ol(0, Cart::xxxzz);
        ol(0, Cart::xxxzzz) = PmB(0) * ol(0, Cart::xxzzz) + 2 * term_xzzz;
        ol(0, Cart::xxyyyy) = PmB(0) * ol(0, Cart::xyyyy) + term_yyyy;
        ol(0, Cart::xxyyyz) = PmB(2) * ol(0, Cart::xxyyy);
        ol(0, Cart::xxyyzz) = PmB(0) * ol(0, Cart::xyyzz) + term_yyzz;
        ol(0, Cart::xxyzzz) = PmB(1) * ol(0, Cart::xxzzz);
        ol(0, Cart::xxzzzz) = PmB(0) * ol(0, Cart::xzzzz) + term_zzzz;
        ol(0, Cart::xyyyyy) = PmB(0) * ol(0, Cart::yyyyy);
        ol(0, Cart::xyyyyz) = PmB(0) * ol(0, Cart::yyyyz);
        ol(0, Cart::xyyyzz) = PmB(0) * ol(0, Cart::yyyzz);
        ol(0, Cart::xyyzzz) = PmB(0) * ol(0, Cart::yyzzz);
        ol(0, Cart::xyzzzz) = PmB(0) * ol(0, Cart::yzzzz);
        ol(0, Cart::xzzzzz) = PmB(0) * ol(0, Cart::zzzzz);
        ol(0, Cart::yyyyyy) = PmB(1) * ol(0, Cart::yyyyy) + 5 * term_yyyy;
        ol(0, Cart::yyyyyz) = PmB(2) * ol(0, Cart::yyyyy);
        ol(0, Cart::yyyyzz) = PmB(2) * ol(0, Cart::yyyyz) + term_yyyy;
        ol(0, Cart::yyyzzz) = PmB(1) * ol(0, Cart::yyzzz) + 2 * term_yzzz;
        ol(0, Cart::yyzzzz) = PmB(1) * ol(0, Cart::yzzzz) + term_zzzz;
        ol(0, Cart::yzzzzz) = PmB(1) * ol(0, Cart::zzzzz);
        ol(0, Cart::zzzzzz) = PmB(2) * ol(0, Cart::zzzzz) + 5 * term_zzzz;
        //------------------------------------------------------

        // Integrals     p - i     d - i     f - i     g - i     h - i     i - i
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xxxx = fak * ol(i, Cart::xxxx);
          double term_xyyy = fak * ol(i, Cart::xyyy);
          double term_xzzz = fak * ol(i, Cart::xzzz);
          double term_yyyy = fak * ol(i, Cart::yyyy);
          double term_yyzz = fak * ol(i, Cart::yyzz);
          double term_yzzz = fak * ol(i, Cart::yzzz);
          double term_zzzz = fak * ol(i, Cart::zzzz);
          ol(i, Cart::xxxxxx) = PmB(0) * ol(i, Cart::xxxxx) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xxxxx) +
                                5 * term_xxxx;
          ol(i, Cart::xxxxxy) = PmB(1) * ol(i, Cart::xxxxx) +
                                ny[i] * fak * ol(i_less_y[i], Cart::xxxxx);
          ol(i, Cart::xxxxxz) = PmB(2) * ol(i, Cart::xxxxx) +
                                nz[i] * fak * ol(i_less_z[i], Cart::xxxxx);
          ol(i, Cart::xxxxyy) = PmB(1) * ol(i, Cart::xxxxy) +
                                ny[i] * fak * ol(i_less_y[i], Cart::xxxxy) +
                                term_xxxx;
          ol(i, Cart::xxxxyz) = PmB(1) * ol(i, Cart::xxxxz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::xxxxz);
          ol(i, Cart::xxxxzz) = PmB(2) * ol(i, Cart::xxxxz) +
                                nz[i] * fak * ol(i_less_z[i], Cart::xxxxz) +
                                term_xxxx;
          ol(i, Cart::xxxyyy) = PmB(0) * ol(i, Cart::xxyyy) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xxyyy) +
                                2 * term_xyyy;
          ol(i, Cart::xxxyyz) = PmB(2) * ol(i, Cart::xxxyy) +
                                nz[i] * fak * ol(i_less_z[i], Cart::xxxyy);
          ol(i, Cart::xxxyzz) = PmB(1) * ol(i, Cart::xxxzz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::xxxzz);
          ol(i, Cart::xxxzzz) = PmB(0) * ol(i, Cart::xxzzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xxzzz) +
                                2 * term_xzzz;
          ol(i, Cart::xxyyyy) = PmB(0) * ol(i, Cart::xyyyy) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xyyyy) +
                                term_yyyy;
          ol(i, Cart::xxyyyz) = PmB(2) * ol(i, Cart::xxyyy) +
                                nz[i] * fak * ol(i_less_z[i], Cart::xxyyy);
          ol(i, Cart::xxyyzz) = PmB(0) * ol(i, Cart::xyyzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xyyzz) +
                                term_yyzz;
          ol(i, Cart::xxyzzz) = PmB(1) * ol(i, Cart::xxzzz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::xxzzz);
          ol(i, Cart::xxzzzz) = PmB(0) * ol(i, Cart::xzzzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::xzzzz) +
                                term_zzzz;
          ol(i, Cart::xyyyyy) = PmB(0) * ol(i, Cart::yyyyy) +
                                nx[i] * fak * ol(i_less_x[i], Cart::yyyyy);
          ol(i, Cart::xyyyyz) = PmB(0) * ol(i, Cart::yyyyz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::yyyyz);
          ol(i, Cart::xyyyzz) = PmB(0) * ol(i, Cart::yyyzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::yyyzz);
          ol(i, Cart::xyyzzz) = PmB(0) * ol(i, Cart::yyzzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::yyzzz);
          ol(i, Cart::xyzzzz) = PmB(0) * ol(i, Cart::yzzzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::yzzzz);
          ol(i, Cart::xzzzzz) = PmB(0) * ol(i, Cart::zzzzz) +
                                nx[i] * fak * ol(i_less_x[i], Cart::zzzzz);
          ol(i, Cart::yyyyyy) = PmB(1) * ol(i, Cart::yyyyy) +
                                ny[i] * fak * ol(i_less_y[i], Cart::yyyyy) +
                                5 * term_yyyy;
          ol(i, Cart::yyyyyz) = PmB(2) * ol(i, Cart::yyyyy) +
                                nz[i] * fak * ol(i_less_z[i], Cart::yyyyy);
          ol(i, Cart::yyyyzz) = PmB(2) * ol(i, Cart::yyyyz) +
                                nz[i] * fak * ol(i_less_z[i], Cart::yyyyz) +
                                term_yyyy;
          ol(i, Cart::yyyzzz) = PmB(1) * ol(i, Cart::yyzzz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::yyzzz) +
                                2 * term_yzzz;
          ol(i, Cart::yyzzzz) = PmB(1) * ol(i, Cart::yzzzz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::yzzzz) +
                                term_zzzz;
          ol(i, Cart::yzzzzz) = PmB(1) * ol(i, Cart::zzzzz) +
                                ny[i] * fak * ol(i_less_y[i], Cart::zzzzz);
          ol(i, Cart::zzzzzz) = PmB(2) * ol(i, Cart::zzzzz) +
                                nz[i] * fak * ol(i_less_z[i], Cart::zzzzz) +
                                5 * term_zzzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 5)

      // cout << "Done with unnormalized matrix " << endl;

      Eigen::MatrixXd ol_sph = AOTransform::getTrafo(gaussian_row).transpose() *
                               ol * AOTransform::getTrafo(gaussian_col);
      // save to matrix

      matrix += ol_sph.block(shell_row.getOffset(), shell_col.getOffset(),
                             matrix.rows(), matrix.cols());
    }  // shell_col Gaussians
  }    // shell_row Gaussians
}

Eigen::MatrixXd AOOverlap::FillShell(const AOShell& shell) const {
  Eigen::MatrixXd block =
      Eigen::MatrixXd::Zero(shell.getNumFunc(), shell.getNumFunc());
  Eigen::Block<Eigen::MatrixXd> submatrix =
      block.block(0, 0, shell.getNumFunc(), shell.getNumFunc());
  FillBlock(submatrix, shell, shell);
  return block;
}

Eigen::MatrixXd AOOverlap::Pseudo_InvSqrt(double etol) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());
  removedfunctions = 0;
  for (unsigned i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  return es.eigenvectors() * diagonal.asDiagonal() *
         es.eigenvectors().transpose();
}

Eigen::MatrixXd AOOverlap::Sqrt() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  return es.operatorSqrt();
}

}  // namespace xtp
}  // namespace votca
