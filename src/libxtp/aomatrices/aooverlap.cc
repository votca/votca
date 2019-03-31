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

void AOOverlap::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                          const AOShell* shell_row, const AOShell* shell_col) {

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

  /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
   * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION
   * COEFFICIENTS, AND ADD TO matrix(i,j)
   */

  // get shell positions
  const tools::vec& pos_row = shell_row->getPos();
  const tools::vec& pos_col = shell_col->getPos();
  const tools::vec diff = pos_row - pos_col;
  std::vector<double> _pma(3, 0.0);
  std::vector<double> _pmb(3, 0.0);

  double distsq = (diff * diff);
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
    const double decay_row = itr->getDecay();

    for (AOShell::GaussianIterator itc = shell_col->begin();
         itc != shell_col->end(); ++itc) {

      // get decay constants
      const double decay_col = itc->getDecay();

      // some helpers
      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;

      // check if distance between postions is big, then skip step
      double exparg = fak2 * decay_row * decay_col * distsq;
      if (exparg > 30.0) {
        continue;
      }
      // initialize local matrix block for unnormalized cartesians
      Eigen::MatrixXd _ol = Eigen::MatrixXd::Zero(nrows, ncols);

      // Definition of coefficients for recursive overlap formulas
      // A for rows (i). B for columns (j)
      const double PmA0 =
          fak2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) -
          pos_row.getX();
      const double PmA1 =
          fak2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) -
          pos_row.getY();
      const double PmA2 =
          fak2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) -
          pos_row.getZ();

      const double PmB0 =
          fak2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) -
          pos_col.getX();
      const double PmB1 =
          fak2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) -
          pos_col.getY();
      const double PmB2 =
          fak2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) -
          pos_col.getZ();

      // calculate matrix elements
      _ol(0, 0) = pow(4.0 * decay_row * decay_col, 0.75) * pow(fak2, 1.5) *
                  exp(-exparg);  // s-s element
      // cout << "\t setting s-s: " << _ol(0,0) << endl;

      // Integrals     p - s
      if (lmax_row > 0) {
        _ol(Cart::x, 0) = PmA0 * _ol(0, 0);
        _ol(Cart::y, 0) = PmA1 * _ol(0, 0);
        _ol(Cart::z, 0) = PmA2 * _ol(0, 0);
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        double term = fak * _ol(0, 0);
        _ol(Cart::xx, 0) = PmA0 * _ol(Cart::x, 0) + term;
        _ol(Cart::xy, 0) = PmA0 * _ol(Cart::y, 0);
        _ol(Cart::xz, 0) = PmA0 * _ol(Cart::z, 0);
        _ol(Cart::yy, 0) = PmA1 * _ol(Cart::y, 0) + term;
        _ol(Cart::yz, 0) = PmA1 * _ol(Cart::z, 0);
        _ol(Cart::zz, 0) = PmA2 * _ol(Cart::z, 0) + term;
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        _ol(Cart::xxx, 0) = PmA0 * _ol(Cart::xx, 0) + 2 * fak * _ol(Cart::x, 0);
        _ol(Cart::xxy, 0) = PmA1 * _ol(Cart::xx, 0);
        _ol(Cart::xxz, 0) = PmA2 * _ol(Cart::xx, 0);
        _ol(Cart::xyy, 0) = PmA0 * _ol(Cart::yy, 0);
        _ol(Cart::xyz, 0) = PmA0 * _ol(Cart::yz, 0);
        _ol(Cart::xzz, 0) = PmA0 * _ol(Cart::zz, 0);
        _ol(Cart::yyy, 0) = PmA1 * _ol(Cart::yy, 0) + 2 * fak * _ol(Cart::y, 0);
        _ol(Cart::yyz, 0) = PmA2 * _ol(Cart::yy, 0);
        _ol(Cart::yzz, 0) = PmA1 * _ol(Cart::zz, 0);
        _ol(Cart::zzz, 0) = PmA2 * _ol(Cart::zz, 0) + 2 * fak * _ol(Cart::z, 0);
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        double term_xx = fak * _ol(Cart::xx, 0);
        double term_yy = fak * _ol(Cart::yy, 0);
        double term_zz = fak * _ol(Cart::zz, 0);
        _ol(Cart::xxxx, 0) = PmA0 * _ol(Cart::xxx, 0) + 3 * term_xx;
        _ol(Cart::xxxy, 0) = PmA1 * _ol(Cart::xxx, 0);
        _ol(Cart::xxxz, 0) = PmA2 * _ol(Cart::xxx, 0);
        _ol(Cart::xxyy, 0) = PmA0 * _ol(Cart::xyy, 0) + term_yy;
        _ol(Cart::xxyz, 0) = PmA1 * _ol(Cart::xxz, 0);
        _ol(Cart::xxzz, 0) = PmA0 * _ol(Cart::xzz, 0) + term_zz;
        _ol(Cart::xyyy, 0) = PmA0 * _ol(Cart::yyy, 0);
        _ol(Cart::xyyz, 0) = PmA0 * _ol(Cart::yyz, 0);
        _ol(Cart::xyzz, 0) = PmA0 * _ol(Cart::yzz, 0);
        _ol(Cart::xzzz, 0) = PmA0 * _ol(Cart::zzz, 0);
        _ol(Cart::yyyy, 0) = PmA1 * _ol(Cart::yyy, 0) + 3 * term_yy;
        _ol(Cart::yyyz, 0) = PmA2 * _ol(Cart::yyy, 0);
        _ol(Cart::yyzz, 0) = PmA1 * _ol(Cart::yzz, 0) + term_zz;
        _ol(Cart::yzzz, 0) = PmA1 * _ol(Cart::zzz, 0);
        _ol(Cart::zzzz, 0) = PmA2 * _ol(Cart::zzz, 0) + 3 * term_zz;
      }
      //------------------------------------------------------

      // Integrals     h - s
      if (lmax_row > 4) {
        double term_xxx = fak * _ol(Cart::xxx, 0);
        double term_yyy = fak * _ol(Cart::yyy, 0);
        double term_zzz = fak * _ol(Cart::zzz, 0);
        _ol(Cart::xxxxx, 0) = PmA0 * _ol(Cart::xxxx, 0) + 4 * term_xxx;
        _ol(Cart::xxxxy, 0) = PmA1 * _ol(Cart::xxxx, 0);
        _ol(Cart::xxxxz, 0) = PmA2 * _ol(Cart::xxxx, 0);
        _ol(Cart::xxxyy, 0) = PmA1 * _ol(Cart::xxxy, 0) + term_xxx;
        _ol(Cart::xxxyz, 0) = PmA1 * _ol(Cart::xxxz, 0);
        _ol(Cart::xxxzz, 0) = PmA2 * _ol(Cart::xxxz, 0) + term_xxx;
        _ol(Cart::xxyyy, 0) = PmA0 * _ol(Cart::xyyy, 0) + term_yyy;
        _ol(Cart::xxyyz, 0) = PmA2 * _ol(Cart::xxyy, 0);
        _ol(Cart::xxyzz, 0) = PmA1 * _ol(Cart::xxzz, 0);
        _ol(Cart::xxzzz, 0) = PmA0 * _ol(Cart::xzzz, 0) + term_zzz;
        _ol(Cart::xyyyy, 0) = PmA0 * _ol(Cart::yyyy, 0);
        _ol(Cart::xyyyz, 0) = PmA0 * _ol(Cart::yyyz, 0);
        _ol(Cart::xyyzz, 0) = PmA0 * _ol(Cart::yyzz, 0);
        _ol(Cart::xyzzz, 0) = PmA0 * _ol(Cart::yzzz, 0);
        _ol(Cart::xzzzz, 0) = PmA0 * _ol(Cart::zzzz, 0);
        _ol(Cart::yyyyy, 0) = PmA1 * _ol(Cart::yyyy, 0) + 4 * term_yyy;
        _ol(Cart::yyyyz, 0) = PmA2 * _ol(Cart::yyyy, 0);
        _ol(Cart::yyyzz, 0) = PmA2 * _ol(Cart::yyyz, 0) + term_yyy;
        _ol(Cart::yyzzz, 0) = PmA1 * _ol(Cart::yzzz, 0) + term_zzz;
        _ol(Cart::yzzzz, 0) = PmA1 * _ol(Cart::zzzz, 0);
        _ol(Cart::zzzzz, 0) = PmA2 * _ol(Cart::zzzz, 0) + 4 * term_zzz;
      }
      //------------------------------------------------------

      // Integrals     i - s
      if (lmax_row > 5) {
        double term_xxxx = fak * _ol(Cart::xxxx, 0);
        double term_xyyy = fak * _ol(Cart::xyyy, 0);
        double term_xzzz = fak * _ol(Cart::xzzz, 0);
        double term_yyyy = fak * _ol(Cart::yyyy, 0);
        double term_yyzz = fak * _ol(Cart::yyzz, 0);
        double term_yzzz = fak * _ol(Cart::yzzz, 0);
        double term_zzzz = fak * _ol(Cart::zzzz, 0);
        _ol(Cart::xxxxxx, 0) = PmA0 * _ol(Cart::xxxxx, 0) + 5 * term_xxxx;
        _ol(Cart::xxxxxy, 0) = PmA1 * _ol(Cart::xxxxx, 0);
        _ol(Cart::xxxxxz, 0) = PmA2 * _ol(Cart::xxxxx, 0);
        _ol(Cart::xxxxyy, 0) = PmA1 * _ol(Cart::xxxxy, 0) + term_xxxx;
        _ol(Cart::xxxxyz, 0) = PmA1 * _ol(Cart::xxxxz, 0);
        _ol(Cart::xxxxzz, 0) = PmA2 * _ol(Cart::xxxxz, 0) + term_xxxx;
        _ol(Cart::xxxyyy, 0) = PmA0 * _ol(Cart::xxyyy, 0) + 2 * term_xyyy;
        _ol(Cart::xxxyyz, 0) = PmA2 * _ol(Cart::xxxyy, 0);
        _ol(Cart::xxxyzz, 0) = PmA1 * _ol(Cart::xxxzz, 0);
        _ol(Cart::xxxzzz, 0) = PmA0 * _ol(Cart::xxzzz, 0) + 2 * term_xzzz;
        _ol(Cart::xxyyyy, 0) = PmA0 * _ol(Cart::xyyyy, 0) + term_yyyy;
        _ol(Cart::xxyyyz, 0) = PmA2 * _ol(Cart::xxyyy, 0);
        _ol(Cart::xxyyzz, 0) = PmA0 * _ol(Cart::xyyzz, 0) + term_yyzz;
        _ol(Cart::xxyzzz, 0) = PmA1 * _ol(Cart::xxzzz, 0);
        _ol(Cart::xxzzzz, 0) = PmA0 * _ol(Cart::xzzzz, 0) + term_zzzz;
        _ol(Cart::xyyyyy, 0) = PmA0 * _ol(Cart::yyyyy, 0);
        _ol(Cart::xyyyyz, 0) = PmA0 * _ol(Cart::yyyyz, 0);
        _ol(Cart::xyyyzz, 0) = PmA0 * _ol(Cart::yyyzz, 0);
        _ol(Cart::xyyzzz, 0) = PmA0 * _ol(Cart::yyzzz, 0);
        _ol(Cart::xyzzzz, 0) = PmA0 * _ol(Cart::yzzzz, 0);
        _ol(Cart::xzzzzz, 0) = PmA0 * _ol(Cart::zzzzz, 0);
        _ol(Cart::yyyyyy, 0) = PmA1 * _ol(Cart::yyyyy, 0) + 5 * term_yyyy;
        _ol(Cart::yyyyyz, 0) = PmA2 * _ol(Cart::yyyyy, 0);
        _ol(Cart::yyyyzz, 0) = PmA2 * _ol(Cart::yyyyz, 0) + term_yyyy;
        _ol(Cart::yyyzzz, 0) = PmA1 * _ol(Cart::yyzzz, 0) + 2 * term_yzzz;
        _ol(Cart::yyzzzz, 0) = PmA1 * _ol(Cart::yzzzz, 0) + term_zzzz;
        _ol(Cart::yzzzzz, 0) = PmA1 * _ol(Cart::zzzzz, 0);
        _ol(Cart::zzzzzz, 0) = PmA2 * _ol(Cart::zzzzz, 0) + 5 * term_zzzz;
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        _ol(0, Cart::x) = PmB0 * _ol(0, 0);
        _ol(0, Cart::y) = PmB1 * _ol(0, 0);
        _ol(0, Cart::z) = PmB2 * _ol(0, 0);
        //------------------------------------------------------

        // Integrals     p - p     d - p     f - p     g - p     h - p     i - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          _ol(i, Cart::x) =
              PmB0 * _ol(i, 0) + nx[i] * fak * _ol(i_less_x[i], 0);
          _ol(i, Cart::y) =
              PmB1 * _ol(i, 0) + ny[i] * fak * _ol(i_less_y[i], 0);
          _ol(i, Cart::z) =
              PmB2 * _ol(i, 0) + nz[i] * fak * _ol(i_less_z[i], 0);
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        double term = fak * _ol(0, 0);
        _ol(0, Cart::xx) = PmB0 * _ol(0, Cart::x) + term;
        _ol(0, Cart::xy) = PmB0 * _ol(0, Cart::y);
        _ol(0, Cart::xz) = PmB0 * _ol(0, Cart::z);
        _ol(0, Cart::yy) = PmB1 * _ol(0, Cart::y) + term;
        _ol(0, Cart::yz) = PmB1 * _ol(0, Cart::z);
        _ol(0, Cart::zz) = PmB2 * _ol(0, Cart::z) + term;
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d     h - d     i - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term = fak * _ol(i, 0);
          _ol(i, Cart::xx) = PmB0 * _ol(i, Cart::x) +
                             nx[i] * fak * _ol(i_less_x[i], Cart::x) + term;
          _ol(i, Cart::xy) =
              PmB0 * _ol(i, Cart::y) + nx[i] * fak * _ol(i_less_x[i], Cart::y);
          _ol(i, Cart::xz) =
              PmB0 * _ol(i, Cart::z) + nx[i] * fak * _ol(i_less_x[i], Cart::z);
          _ol(i, Cart::yy) = PmB1 * _ol(i, Cart::y) +
                             ny[i] * fak * _ol(i_less_y[i], Cart::y) + term;
          _ol(i, Cart::yz) =
              PmB1 * _ol(i, Cart::z) + ny[i] * fak * _ol(i_less_y[i], Cart::z);
          _ol(i, Cart::zz) = PmB2 * _ol(i, Cart::z) +
                             nz[i] * fak * _ol(i_less_z[i], Cart::z) + term;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        _ol(0, Cart::xxx) = PmB0 * _ol(0, Cart::xx) + 2 * fak * _ol(0, Cart::x);
        _ol(0, Cart::xxy) = PmB1 * _ol(0, Cart::xx);
        _ol(0, Cart::xxz) = PmB2 * _ol(0, Cart::xx);
        _ol(0, Cart::xyy) = PmB0 * _ol(0, Cart::yy);
        _ol(0, Cart::xyz) = PmB0 * _ol(0, Cart::yz);
        _ol(0, Cart::xzz) = PmB0 * _ol(0, Cart::zz);
        _ol(0, Cart::yyy) = PmB1 * _ol(0, Cart::yy) + 2 * fak * _ol(0, Cart::y);
        _ol(0, Cart::yyz) = PmB2 * _ol(0, Cart::yy);
        _ol(0, Cart::yzz) = PmB1 * _ol(0, Cart::zz);
        _ol(0, Cart::zzz) = PmB2 * _ol(0, Cart::zz) + 2 * fak * _ol(0, Cart::z);
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f     h - f     i - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_x = 2 * fak * _ol(i, Cart::x);
          double term_y = 2 * fak * _ol(i, Cart::y);
          double term_z = 2 * fak * _ol(i, Cart::z);
          _ol(i, Cart::xxx) = PmB0 * _ol(i, Cart::xx) +
                              nx[i] * fak * _ol(i_less_x[i], Cart::xx) + term_x;
          _ol(i, Cart::xxy) = PmB1 * _ol(i, Cart::xx) +
                              ny[i] * fak * _ol(i_less_y[i], Cart::xx);
          _ol(i, Cart::xxz) = PmB2 * _ol(i, Cart::xx) +
                              nz[i] * fak * _ol(i_less_z[i], Cart::xx);
          _ol(i, Cart::xyy) = PmB0 * _ol(i, Cart::yy) +
                              nx[i] * fak * _ol(i_less_x[i], Cart::yy);
          _ol(i, Cart::xyz) = PmB0 * _ol(i, Cart::yz) +
                              nx[i] * fak * _ol(i_less_x[i], Cart::yz);
          _ol(i, Cart::xzz) = PmB0 * _ol(i, Cart::zz) +
                              nx[i] * fak * _ol(i_less_x[i], Cart::zz);
          _ol(i, Cart::yyy) = PmB1 * _ol(i, Cart::yy) +
                              ny[i] * fak * _ol(i_less_y[i], Cart::yy) + term_y;
          _ol(i, Cart::yyz) = PmB2 * _ol(i, Cart::yy) +
                              nz[i] * fak * _ol(i_less_z[i], Cart::yy);
          _ol(i, Cart::yzz) = PmB1 * _ol(i, Cart::zz) +
                              ny[i] * fak * _ol(i_less_y[i], Cart::zz);
          _ol(i, Cart::zzz) = PmB2 * _ol(i, Cart::zz) +
                              nz[i] * fak * _ol(i_less_z[i], Cart::zz) + term_z;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        double term_xx = fak * _ol(0, Cart::xx);
        double term_yy = fak * _ol(0, Cart::yy);
        double term_zz = fak * _ol(0, Cart::zz);
        _ol(0, Cart::xxxx) = PmB0 * _ol(0, Cart::xxx) + 3 * term_xx;
        _ol(0, Cart::xxxy) = PmB1 * _ol(0, Cart::xxx);
        _ol(0, Cart::xxxz) = PmB2 * _ol(0, Cart::xxx);
        _ol(0, Cart::xxyy) = PmB0 * _ol(0, Cart::xyy) + term_yy;
        _ol(0, Cart::xxyz) = PmB1 * _ol(0, Cart::xxz);
        _ol(0, Cart::xxzz) = PmB0 * _ol(0, Cart::xzz) + term_zz;
        _ol(0, Cart::xyyy) = PmB0 * _ol(0, Cart::yyy);
        _ol(0, Cart::xyyz) = PmB0 * _ol(0, Cart::yyz);
        _ol(0, Cart::xyzz) = PmB0 * _ol(0, Cart::yzz);
        _ol(0, Cart::xzzz) = PmB0 * _ol(0, Cart::zzz);
        _ol(0, Cart::yyyy) = PmB1 * _ol(0, Cart::yyy) + 3 * term_yy;
        _ol(0, Cart::yyyz) = PmB2 * _ol(0, Cart::yyy);
        _ol(0, Cart::yyzz) = PmB1 * _ol(0, Cart::yzz) + term_zz;
        _ol(0, Cart::yzzz) = PmB1 * _ol(0, Cart::zzz);
        _ol(0, Cart::zzzz) = PmB2 * _ol(0, Cart::zzz) + 3 * term_zz;
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g     h - g     i - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xx = fak * _ol(i, Cart::xx);
          double term_yy = fak * _ol(i, Cart::yy);
          double term_zz = fak * _ol(i, Cart::zz);
          _ol(i, Cart::xxxx) = PmB0 * _ol(i, Cart::xxx) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::xxx) +
                               3 * term_xx;
          _ol(i, Cart::xxxy) = PmB1 * _ol(i, Cart::xxx) +
                               ny[i] * fak * _ol(i_less_y[i], Cart::xxx);
          _ol(i, Cart::xxxz) = PmB2 * _ol(i, Cart::xxx) +
                               nz[i] * fak * _ol(i_less_z[i], Cart::xxx);
          _ol(i, Cart::xxyy) = PmB0 * _ol(i, Cart::xyy) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::xyy) +
                               term_yy;
          _ol(i, Cart::xxyz) = PmB1 * _ol(i, Cart::xxz) +
                               ny[i] * fak * _ol(i_less_y[i], Cart::xxz);
          _ol(i, Cart::xxzz) = PmB0 * _ol(i, Cart::xzz) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::xzz) +
                               term_zz;
          _ol(i, Cart::xyyy) = PmB0 * _ol(i, Cart::yyy) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::yyy);
          _ol(i, Cart::xyyz) = PmB0 * _ol(i, Cart::yyz) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::yyz);
          _ol(i, Cart::xyzz) = PmB0 * _ol(i, Cart::yzz) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::yzz);
          _ol(i, Cart::xzzz) = PmB0 * _ol(i, Cart::zzz) +
                               nx[i] * fak * _ol(i_less_x[i], Cart::zzz);
          _ol(i, Cart::yyyy) = PmB1 * _ol(i, Cart::yyy) +
                               ny[i] * fak * _ol(i_less_y[i], Cart::yyy) +
                               3 * term_yy;
          _ol(i, Cart::yyyz) = PmB2 * _ol(i, Cart::yyy) +
                               nz[i] * fak * _ol(i_less_z[i], Cart::yyy);
          _ol(i, Cart::yyzz) = PmB1 * _ol(i, Cart::yzz) +
                               ny[i] * fak * _ol(i_less_y[i], Cart::yzz) +
                               term_zz;
          _ol(i, Cart::yzzz) = PmB1 * _ol(i, Cart::zzz) +
                               ny[i] * fak * _ol(i_less_y[i], Cart::zzz);
          _ol(i, Cart::zzzz) = PmB2 * _ol(i, Cart::zzz) +
                               nz[i] * fak * _ol(i_less_z[i], Cart::zzz) +
                               3 * term_zz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      if (lmax_col > 4) {

        // Integrals     s - h
        double term_xxx = fak * _ol(0, Cart::xxx);
        double term_yyy = fak * _ol(0, Cart::yyy);
        double term_zzz = fak * _ol(0, Cart::zzz);
        _ol(0, Cart::xxxxx) = PmB0 * _ol(0, Cart::xxxx) + 4 * term_xxx;
        _ol(0, Cart::xxxxy) = PmB1 * _ol(0, Cart::xxxx);
        _ol(0, Cart::xxxxz) = PmB2 * _ol(0, Cart::xxxx);
        _ol(0, Cart::xxxyy) = PmB1 * _ol(0, Cart::xxxy) + term_xxx;
        _ol(0, Cart::xxxyz) = PmB1 * _ol(0, Cart::xxxz);
        _ol(0, Cart::xxxzz) = PmB2 * _ol(0, Cart::xxxz) + term_xxx;
        _ol(0, Cart::xxyyy) = PmB0 * _ol(0, Cart::xyyy) + term_yyy;
        _ol(0, Cart::xxyyz) = PmB2 * _ol(0, Cart::xxyy);
        _ol(0, Cart::xxyzz) = PmB1 * _ol(0, Cart::xxzz);
        _ol(0, Cart::xxzzz) = PmB0 * _ol(0, Cart::xzzz) + term_zzz;
        _ol(0, Cart::xyyyy) = PmB0 * _ol(0, Cart::yyyy);
        _ol(0, Cart::xyyyz) = PmB0 * _ol(0, Cart::yyyz);
        _ol(0, Cart::xyyzz) = PmB0 * _ol(0, Cart::yyzz);
        _ol(0, Cart::xyzzz) = PmB0 * _ol(0, Cart::yzzz);
        _ol(0, Cart::xzzzz) = PmB0 * _ol(0, Cart::zzzz);
        _ol(0, Cart::yyyyy) = PmB1 * _ol(0, Cart::yyyy) + 4 * term_yyy;
        _ol(0, Cart::yyyyz) = PmB2 * _ol(0, Cart::yyyy);
        _ol(0, Cart::yyyzz) = PmB2 * _ol(0, Cart::yyyz) + term_yyy;
        _ol(0, Cart::yyzzz) = PmB1 * _ol(0, Cart::yzzz) + term_zzz;
        _ol(0, Cart::yzzzz) = PmB1 * _ol(0, Cart::zzzz);
        _ol(0, Cart::zzzzz) = PmB2 * _ol(0, Cart::zzzz) + 4 * term_zzz;
        //------------------------------------------------------

        // Integrals     p - h     d - h     f - h     g - h     h - h     i - h
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xxx = fak * _ol(i, Cart::xxx);
          double term_yyy = fak * _ol(i, Cart::yyy);
          double term_zzz = fak * _ol(i, Cart::zzz);
          _ol(i, Cart::xxxxx) = PmB0 * _ol(i, Cart::xxxx) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::xxxx) +
                                4 * term_xxx;
          _ol(i, Cart::xxxxy) = PmB1 * _ol(i, Cart::xxxx) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::xxxx);
          _ol(i, Cart::xxxxz) = PmB2 * _ol(i, Cart::xxxx) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::xxxx);
          _ol(i, Cart::xxxyy) = PmB1 * _ol(i, Cart::xxxy) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::xxxy) +
                                term_xxx;
          _ol(i, Cart::xxxyz) = PmB1 * _ol(i, Cart::xxxz) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::xxxz);
          _ol(i, Cart::xxxzz) = PmB2 * _ol(i, Cart::xxxz) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::xxxz) +
                                term_xxx;
          _ol(i, Cart::xxyyy) = PmB0 * _ol(i, Cart::xyyy) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::xyyy) +
                                term_yyy;
          _ol(i, Cart::xxyyz) = PmB2 * _ol(i, Cart::xxyy) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::xxyy);
          _ol(i, Cart::xxyzz) = PmB1 * _ol(i, Cart::xxzz) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::xxzz);
          _ol(i, Cart::xxzzz) = PmB0 * _ol(i, Cart::xzzz) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::xzzz) +
                                term_zzz;
          _ol(i, Cart::xyyyy) = PmB0 * _ol(i, Cart::yyyy) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::yyyy);
          _ol(i, Cart::xyyyz) = PmB0 * _ol(i, Cart::yyyz) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::yyyz);
          _ol(i, Cart::xyyzz) = PmB0 * _ol(i, Cart::yyzz) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::yyzz);
          _ol(i, Cart::xyzzz) = PmB0 * _ol(i, Cart::yzzz) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::yzzz);
          _ol(i, Cart::xzzzz) = PmB0 * _ol(i, Cart::zzzz) +
                                nx[i] * fak * _ol(i_less_x[i], Cart::zzzz);
          _ol(i, Cart::yyyyy) = PmB1 * _ol(i, Cart::yyyy) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::yyyy) +
                                4 * term_yyy;
          _ol(i, Cart::yyyyz) = PmB2 * _ol(i, Cart::yyyy) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::yyyy);
          _ol(i, Cart::yyyzz) = PmB2 * _ol(i, Cart::yyyz) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::yyyz) +
                                term_yyy;
          _ol(i, Cart::yyzzz) = PmB1 * _ol(i, Cart::yzzz) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::yzzz) +
                                term_zzz;
          _ol(i, Cart::yzzzz) = PmB1 * _ol(i, Cart::zzzz) +
                                ny[i] * fak * _ol(i_less_y[i], Cart::zzzz);
          _ol(i, Cart::zzzzz) = PmB2 * _ol(i, Cart::zzzz) +
                                nz[i] * fak * _ol(i_less_z[i], Cart::zzzz) +
                                4 * term_zzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 4)

      if (lmax_col > 5) {

        // Integrals     s - i
        double term_xxxx = fak * _ol(0, Cart::xxxx);
        double term_xyyy = fak * _ol(0, Cart::xyyy);
        double term_xzzz = fak * _ol(0, Cart::xzzz);
        double term_yyyy = fak * _ol(0, Cart::yyyy);
        double term_yyzz = fak * _ol(0, Cart::yyzz);
        double term_yzzz = fak * _ol(0, Cart::yzzz);
        double term_zzzz = fak * _ol(0, Cart::zzzz);
        _ol(0, Cart::xxxxxx) = PmB0 * _ol(0, Cart::xxxxx) + 5 * term_xxxx;
        _ol(0, Cart::xxxxxy) = PmB1 * _ol(0, Cart::xxxxx);
        _ol(0, Cart::xxxxxz) = PmB2 * _ol(0, Cart::xxxxx);
        _ol(0, Cart::xxxxyy) = PmB1 * _ol(0, Cart::xxxxy) + term_xxxx;
        _ol(0, Cart::xxxxyz) = PmB1 * _ol(0, Cart::xxxxz);
        _ol(0, Cart::xxxxzz) = PmB2 * _ol(0, Cart::xxxxz) + term_xxxx;
        _ol(0, Cart::xxxyyy) = PmB0 * _ol(0, Cart::xxyyy) + 2 * term_xyyy;
        _ol(0, Cart::xxxyyz) = PmB2 * _ol(0, Cart::xxxyy);
        _ol(0, Cart::xxxyzz) = PmB1 * _ol(0, Cart::xxxzz);
        _ol(0, Cart::xxxzzz) = PmB0 * _ol(0, Cart::xxzzz) + 2 * term_xzzz;
        _ol(0, Cart::xxyyyy) = PmB0 * _ol(0, Cart::xyyyy) + term_yyyy;
        _ol(0, Cart::xxyyyz) = PmB2 * _ol(0, Cart::xxyyy);
        _ol(0, Cart::xxyyzz) = PmB0 * _ol(0, Cart::xyyzz) + term_yyzz;
        _ol(0, Cart::xxyzzz) = PmB1 * _ol(0, Cart::xxzzz);
        _ol(0, Cart::xxzzzz) = PmB0 * _ol(0, Cart::xzzzz) + term_zzzz;
        _ol(0, Cart::xyyyyy) = PmB0 * _ol(0, Cart::yyyyy);
        _ol(0, Cart::xyyyyz) = PmB0 * _ol(0, Cart::yyyyz);
        _ol(0, Cart::xyyyzz) = PmB0 * _ol(0, Cart::yyyzz);
        _ol(0, Cart::xyyzzz) = PmB0 * _ol(0, Cart::yyzzz);
        _ol(0, Cart::xyzzzz) = PmB0 * _ol(0, Cart::yzzzz);
        _ol(0, Cart::xzzzzz) = PmB0 * _ol(0, Cart::zzzzz);
        _ol(0, Cart::yyyyyy) = PmB1 * _ol(0, Cart::yyyyy) + 5 * term_yyyy;
        _ol(0, Cart::yyyyyz) = PmB2 * _ol(0, Cart::yyyyy);
        _ol(0, Cart::yyyyzz) = PmB2 * _ol(0, Cart::yyyyz) + term_yyyy;
        _ol(0, Cart::yyyzzz) = PmB1 * _ol(0, Cart::yyzzz) + 2 * term_yzzz;
        _ol(0, Cart::yyzzzz) = PmB1 * _ol(0, Cart::yzzzz) + term_zzzz;
        _ol(0, Cart::yzzzzz) = PmB1 * _ol(0, Cart::zzzzz);
        _ol(0, Cart::zzzzzz) = PmB2 * _ol(0, Cart::zzzzz) + 5 * term_zzzz;
        //------------------------------------------------------

        // Integrals     p - i     d - i     f - i     g - i     h - i     i - i
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          double term_xxxx = fak * _ol(i, Cart::xxxx);
          double term_xyyy = fak * _ol(i, Cart::xyyy);
          double term_xzzz = fak * _ol(i, Cart::xzzz);
          double term_yyyy = fak * _ol(i, Cart::yyyy);
          double term_yyzz = fak * _ol(i, Cart::yyzz);
          double term_yzzz = fak * _ol(i, Cart::yzzz);
          double term_zzzz = fak * _ol(i, Cart::zzzz);
          _ol(i, Cart::xxxxxx) = PmB0 * _ol(i, Cart::xxxxx) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xxxxx) +
                                 5 * term_xxxx;
          _ol(i, Cart::xxxxxy) = PmB1 * _ol(i, Cart::xxxxx) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::xxxxx);
          _ol(i, Cart::xxxxxz) = PmB2 * _ol(i, Cart::xxxxx) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::xxxxx);
          _ol(i, Cart::xxxxyy) = PmB1 * _ol(i, Cart::xxxxy) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::xxxxy) +
                                 term_xxxx;
          _ol(i, Cart::xxxxyz) = PmB1 * _ol(i, Cart::xxxxz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::xxxxz);
          _ol(i, Cart::xxxxzz) = PmB2 * _ol(i, Cart::xxxxz) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::xxxxz) +
                                 term_xxxx;
          _ol(i, Cart::xxxyyy) = PmB0 * _ol(i, Cart::xxyyy) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xxyyy) +
                                 2 * term_xyyy;
          _ol(i, Cart::xxxyyz) = PmB2 * _ol(i, Cart::xxxyy) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::xxxyy);
          _ol(i, Cart::xxxyzz) = PmB1 * _ol(i, Cart::xxxzz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::xxxzz);
          _ol(i, Cart::xxxzzz) = PmB0 * _ol(i, Cart::xxzzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xxzzz) +
                                 2 * term_xzzz;
          _ol(i, Cart::xxyyyy) = PmB0 * _ol(i, Cart::xyyyy) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xyyyy) +
                                 term_yyyy;
          _ol(i, Cart::xxyyyz) = PmB2 * _ol(i, Cart::xxyyy) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::xxyyy);
          _ol(i, Cart::xxyyzz) = PmB0 * _ol(i, Cart::xyyzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xyyzz) +
                                 term_yyzz;
          _ol(i, Cart::xxyzzz) = PmB1 * _ol(i, Cart::xxzzz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::xxzzz);
          _ol(i, Cart::xxzzzz) = PmB0 * _ol(i, Cart::xzzzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::xzzzz) +
                                 term_zzzz;
          _ol(i, Cart::xyyyyy) = PmB0 * _ol(i, Cart::yyyyy) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::yyyyy);
          _ol(i, Cart::xyyyyz) = PmB0 * _ol(i, Cart::yyyyz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::yyyyz);
          _ol(i, Cart::xyyyzz) = PmB0 * _ol(i, Cart::yyyzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::yyyzz);
          _ol(i, Cart::xyyzzz) = PmB0 * _ol(i, Cart::yyzzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::yyzzz);
          _ol(i, Cart::xyzzzz) = PmB0 * _ol(i, Cart::yzzzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::yzzzz);
          _ol(i, Cart::xzzzzz) = PmB0 * _ol(i, Cart::zzzzz) +
                                 nx[i] * fak * _ol(i_less_x[i], Cart::zzzzz);
          _ol(i, Cart::yyyyyy) = PmB1 * _ol(i, Cart::yyyyy) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::yyyyy) +
                                 5 * term_yyyy;
          _ol(i, Cart::yyyyyz) = PmB2 * _ol(i, Cart::yyyyy) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::yyyyy);
          _ol(i, Cart::yyyyzz) = PmB2 * _ol(i, Cart::yyyyz) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::yyyyz) +
                                 term_yyyy;
          _ol(i, Cart::yyyzzz) = PmB1 * _ol(i, Cart::yyzzz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::yyzzz) +
                                 2 * term_yzzz;
          _ol(i, Cart::yyzzzz) = PmB1 * _ol(i, Cart::yzzzz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::yzzzz) +
                                 term_zzzz;
          _ol(i, Cart::yzzzzz) = PmB1 * _ol(i, Cart::zzzzz) +
                                 ny[i] * fak * _ol(i_less_y[i], Cart::zzzzz);
          _ol(i, Cart::zzzzzz) = PmB2 * _ol(i, Cart::zzzzz) +
                                 nz[i] * fak * _ol(i_less_z[i], Cart::zzzzz) +
                                 5 * term_zzzz;
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 5)

      // cout << "Done with unnormalized matrix " << endl;

      Eigen::MatrixXd _ol_sph =
          getTrafo(*itr).transpose() * _ol * getTrafo(*itc);
      // save to matrix

      for (unsigned i = 0; i < matrix.rows(); i++) {
        for (unsigned j = 0; j < matrix.cols(); j++) {
          matrix(i, j) +=
              _ol_sph(i + shell_row->getOffset(), j + shell_col->getOffset());
        }
      }
    }  // shell_col Gaussians
  }    // shell_row Gaussians
}

Eigen::MatrixXd AOOverlap::FillShell(const AOShell* shell) {
  Eigen::MatrixXd block =
      Eigen::MatrixXd::Zero(shell->getNumFunc(), shell->getNumFunc());
  Eigen::Block<Eigen::MatrixXd> submatrix =
      block.block(0, 0, shell->getNumFunc(), shell->getNumFunc());
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
