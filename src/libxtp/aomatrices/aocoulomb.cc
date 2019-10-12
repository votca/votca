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
 * olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aotransform.h>
namespace votca {
namespace xtp {

void AOCoulomb::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                          const AOShell& shell_row,
                          const AOShell& shell_col) const {

  // shell info, only lmax tells how far to go
  const int lmax_row = shell_row.getLmax();
  const int lmax_col = shell_col.getLmax();

  // set size of internal block for recursion
  int nrows = AOTransform::getBlockSize(lmax_row);
  int ncols = AOTransform::getBlockSize(lmax_col);
  const int mmax = lmax_row + lmax_col;
  const int nextra = mmax + 1;

  // get shell positions
  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;
  double distsq = diff.squaredNorm();

  const double pi = boost::math::constants::pi<double>();

  int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};
  // for alphabetical order

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
    const double rdecay_row = 0.5 / decay_row;
    const double powfactor_row = gaussian_row.getPowfactor();
    for (const auto& gaussian_col : shell_col) {

      // get decay constants
      const double decay_col = gaussian_col.getDecay();
      const double rdecay_col = 0.5 / decay_col;
      const double powfactor_col = gaussian_col.getPowfactor();

      Eigen::Tensor<double, 3> cou(nrows, ncols, nextra);
      cou.setZero();

      const double decay = decay_row + decay_col;
      const double r_decay = 0.5 / decay;
      const double r_decay_2 = 2. * r_decay;
      const double fac_a_ac = decay_row / decay;
      const double fac_c_ac = decay_col / decay;

      const Eigen::Vector3d wmp =
          r_decay_2 * (decay_row * pos_row + decay_col * pos_col) - pos_row;
      const Eigen::Vector3d wmq =
          r_decay_2 * (decay_row * pos_row + decay_col * pos_col) - pos_col;

      const double T = fac_a_ac * decay_col * distsq;

      double fak = 2.0 * pow(pi, 2.5) /
                   (decay_row * decay_col * sqrt(decay_row + decay_col));
      fak = fak * powfactor_col * powfactor_row;

      const Eigen::VectorXd FmT = AOTransform::XIntegrate(nextra, T);

      // get initial data from FmT -> s-s element
      for (int i = 0; i != nextra; ++i) {
        cou(0, 0, i) = fak * FmT[i];
      }

      // Integrals     p - s
      if (lmax_row > 0) {
        for (int m = 0; m < mmax; m++) {
          cou(Cart::x, 0, m) = wmp(0) * cou(0, 0, m + 1);
          cou(Cart::y, 0, m) = wmp(1) * cou(0, 0, m + 1);
          cou(Cart::z, 0, m) = wmp(2) * cou(0, 0, m + 1);
        }
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        for (int m = 0; m < mmax - 1; m++) {
          double term =
              rdecay_row * (cou(0, 0, m) - fac_c_ac * cou(0, 0, m + 1));
          cou(Cart::xx, 0, m) = wmp(0) * cou(Cart::x, 0, m + 1) + term;
          cou(Cart::xy, 0, m) = wmp(0) * cou(Cart::y, 0, m + 1);
          cou(Cart::xz, 0, m) = wmp(0) * cou(Cart::z, 0, m + 1);
          cou(Cart::yy, 0, m) = wmp(1) * cou(Cart::y, 0, m + 1) + term;
          cou(Cart::yz, 0, m) = wmp(1) * cou(Cart::z, 0, m + 1);
          cou(Cart::zz, 0, m) = wmp(2) * cou(Cart::z, 0, m + 1) + term;
        }
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        for (int m = 0; m < mmax - 2; m++) {
          cou(Cart::xxx, 0, m) =
              wmp(0) * cou(Cart::xx, 0, m + 1) +
              2 * rdecay_row *
                  (cou(Cart::x, 0, m) - fac_c_ac * cou(Cart::x, 0, m + 1));
          cou(Cart::xxy, 0, m) = wmp(1) * cou(Cart::xx, 0, m + 1);
          cou(Cart::xxz, 0, m) = wmp(2) * cou(Cart::xx, 0, m + 1);
          cou(Cart::xyy, 0, m) = wmp(0) * cou(Cart::yy, 0, m + 1);
          cou(Cart::xyz, 0, m) = wmp(0) * cou(Cart::yz, 0, m + 1);
          cou(Cart::xzz, 0, m) = wmp(0) * cou(Cart::zz, 0, m + 1);
          cou(Cart::yyy, 0, m) =
              wmp(1) * cou(Cart::yy, 0, m + 1) +
              2 * rdecay_row *
                  (cou(Cart::y, 0, m) - fac_c_ac * cou(Cart::y, 0, m + 1));
          cou(Cart::yyz, 0, m) = wmp(2) * cou(Cart::yy, 0, m + 1);
          cou(Cart::yzz, 0, m) = wmp(1) * cou(Cart::zz, 0, m + 1);
          cou(Cart::zzz, 0, m) =
              wmp(2) * cou(Cart::zz, 0, m + 1) +
              2 * rdecay_row *
                  (cou(Cart::z, 0, m) - fac_c_ac * cou(Cart::z, 0, m + 1));
        }
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        for (int m = 0; m < mmax - 3; m++) {
          double term_xx = rdecay_row * (cou(Cart::xx, 0, m) -
                                         fac_c_ac * cou(Cart::xx, 0, m + 1));
          double term_yy = rdecay_row * (cou(Cart::yy, 0, m) -
                                         fac_c_ac * cou(Cart::yy, 0, m + 1));
          double term_zz = rdecay_row * (cou(Cart::zz, 0, m) -
                                         fac_c_ac * cou(Cart::zz, 0, m + 1));
          cou(Cart::xxxx, 0, m) =
              wmp(0) * cou(Cart::xxx, 0, m + 1) + 3 * term_xx;
          cou(Cart::xxxy, 0, m) = wmp(1) * cou(Cart::xxx, 0, m + 1);
          cou(Cart::xxxz, 0, m) = wmp(2) * cou(Cart::xxx, 0, m + 1);
          cou(Cart::xxyy, 0, m) = wmp(0) * cou(Cart::xyy, 0, m + 1) + term_yy;
          cou(Cart::xxyz, 0, m) = wmp(1) * cou(Cart::xxz, 0, m + 1);
          cou(Cart::xxzz, 0, m) = wmp(0) * cou(Cart::xzz, 0, m + 1) + term_zz;
          cou(Cart::xyyy, 0, m) = wmp(0) * cou(Cart::yyy, 0, m + 1);
          cou(Cart::xyyz, 0, m) = wmp(0) * cou(Cart::yyz, 0, m + 1);
          cou(Cart::xyzz, 0, m) = wmp(0) * cou(Cart::yzz, 0, m + 1);
          cou(Cart::xzzz, 0, m) = wmp(0) * cou(Cart::zzz, 0, m + 1);
          cou(Cart::yyyy, 0, m) =
              wmp(1) * cou(Cart::yyy, 0, m + 1) + 3 * term_yy;
          cou(Cart::yyyz, 0, m) = wmp(2) * cou(Cart::yyy, 0, m + 1);
          cou(Cart::yyzz, 0, m) = wmp(1) * cou(Cart::yzz, 0, m + 1) + term_zz;
          cou(Cart::yzzz, 0, m) = wmp(1) * cou(Cart::zzz, 0, m + 1);
          cou(Cart::zzzz, 0, m) =
              wmp(2) * cou(Cart::zzz, 0, m + 1) + 3 * term_zz;
        }
      }
      //------------------------------------------------------

      // Integrals     h - s
      if (lmax_row > 4) {
        for (int m = 0; m < mmax - 4; m++) {
          double term_xxx = rdecay_row * (cou(Cart::xxx, 0, m) -
                                          fac_c_ac * cou(Cart::xxx, 0, m + 1));
          double term_yyy = rdecay_row * (cou(Cart::yyy, 0, m) -
                                          fac_c_ac * cou(Cart::yyy, 0, m + 1));
          double term_zzz = rdecay_row * (cou(Cart::zzz, 0, m) -
                                          fac_c_ac * cou(Cart::zzz, 0, m + 1));
          cou(Cart::xxxxx, 0, m) =
              wmp(0) * cou(Cart::xxxx, 0, m + 1) + 4 * term_xxx;
          cou(Cart::xxxxy, 0, m) = wmp(1) * cou(Cart::xxxx, 0, m + 1);
          cou(Cart::xxxxz, 0, m) = wmp(2) * cou(Cart::xxxx, 0, m + 1);
          cou(Cart::xxxyy, 0, m) =
              wmp(1) * cou(Cart::xxxy, 0, m + 1) + term_xxx;
          cou(Cart::xxxyz, 0, m) = wmp(1) * cou(Cart::xxxz, 0, m + 1);
          cou(Cart::xxxzz, 0, m) =
              wmp(2) * cou(Cart::xxxz, 0, m + 1) + term_xxx;
          cou(Cart::xxyyy, 0, m) =
              wmp(0) * cou(Cart::xyyy, 0, m + 1) + term_yyy;
          cou(Cart::xxyyz, 0, m) = wmp(2) * cou(Cart::xxyy, 0, m + 1);
          cou(Cart::xxyzz, 0, m) = wmp(1) * cou(Cart::xxzz, 0, m + 1);
          cou(Cart::xxzzz, 0, m) =
              wmp(0) * cou(Cart::xzzz, 0, m + 1) + term_zzz;
          cou(Cart::xyyyy, 0, m) = wmp(0) * cou(Cart::yyyy, 0, m + 1);
          cou(Cart::xyyyz, 0, m) = wmp(0) * cou(Cart::yyyz, 0, m + 1);
          cou(Cart::xyyzz, 0, m) = wmp(0) * cou(Cart::yyzz, 0, m + 1);
          cou(Cart::xyzzz, 0, m) = wmp(0) * cou(Cart::yzzz, 0, m + 1);
          cou(Cart::xzzzz, 0, m) = wmp(0) * cou(Cart::zzzz, 0, m + 1);
          cou(Cart::yyyyy, 0, m) =
              wmp(1) * cou(Cart::yyyy, 0, m + 1) + 4 * term_yyy;
          cou(Cart::yyyyz, 0, m) = wmp(2) * cou(Cart::yyyy, 0, m + 1);
          cou(Cart::yyyzz, 0, m) =
              wmp(2) * cou(Cart::yyyz, 0, m + 1) + term_yyy;
          cou(Cart::yyzzz, 0, m) =
              wmp(1) * cou(Cart::yzzz, 0, m + 1) + term_zzz;
          cou(Cart::yzzzz, 0, m) = wmp(1) * cou(Cart::zzzz, 0, m + 1);
          cou(Cart::zzzzz, 0, m) =
              wmp(2) * cou(Cart::zzzz, 0, m + 1) + 4 * term_zzz;
        }
      }
      //------------------------------------------------------

      // Integrals     i - s
      if (lmax_row > 5) {
        for (int m = 0; m < mmax - 5; m++) {
          double term_xxxx =
              rdecay_row *
              (cou(Cart::xxxx, 0, m) - fac_c_ac * cou(Cart::xxxx, 0, m + 1));
          double term_xyyy =
              rdecay_row *
              (cou(Cart::xyyy, 0, m) - fac_c_ac * cou(Cart::xyyy, 0, m + 1));
          double term_xzzz =
              rdecay_row *
              (cou(Cart::xzzz, 0, m) - fac_c_ac * cou(Cart::xzzz, 0, m + 1));
          double term_yyyy =
              rdecay_row *
              (cou(Cart::yyyy, 0, m) - fac_c_ac * cou(Cart::yyyy, 0, m + 1));
          double term_yyzz =
              rdecay_row *
              (cou(Cart::yyzz, 0, m) - fac_c_ac * cou(Cart::yyzz, 0, m + 1));
          double term_yzzz =
              rdecay_row *
              (cou(Cart::yzzz, 0, m) - fac_c_ac * cou(Cart::yzzz, 0, m + 1));
          double term_zzzz =
              rdecay_row *
              (cou(Cart::zzzz, 0, m) - fac_c_ac * cou(Cart::zzzz, 0, m + 1));
          cou(Cart::xxxxxx, 0, m) =
              wmp(0) * cou(Cart::xxxxx, 0, m + 1) + 5 * term_xxxx;
          cou(Cart::xxxxxy, 0, m) = wmp(1) * cou(Cart::xxxxx, 0, m + 1);
          cou(Cart::xxxxxz, 0, m) = wmp(2) * cou(Cart::xxxxx, 0, m + 1);
          cou(Cart::xxxxyy, 0, m) =
              wmp(1) * cou(Cart::xxxxy, 0, m + 1) + term_xxxx;
          cou(Cart::xxxxyz, 0, m) = wmp(1) * cou(Cart::xxxxz, 0, m + 1);
          cou(Cart::xxxxzz, 0, m) =
              wmp(2) * cou(Cart::xxxxz, 0, m + 1) + term_xxxx;
          cou(Cart::xxxyyy, 0, m) =
              wmp(0) * cou(Cart::xxyyy, 0, m + 1) + 2 * term_xyyy;
          cou(Cart::xxxyyz, 0, m) = wmp(2) * cou(Cart::xxxyy, 0, m + 1);
          cou(Cart::xxxyzz, 0, m) = wmp(1) * cou(Cart::xxxzz, 0, m + 1);
          cou(Cart::xxxzzz, 0, m) =
              wmp(0) * cou(Cart::xxzzz, 0, m + 1) + 2 * term_xzzz;
          cou(Cart::xxyyyy, 0, m) =
              wmp(0) * cou(Cart::xyyyy, 0, m + 1) + term_yyyy;
          cou(Cart::xxyyyz, 0, m) = wmp(2) * cou(Cart::xxyyy, 0, m + 1);
          cou(Cart::xxyyzz, 0, m) =
              wmp(0) * cou(Cart::xyyzz, 0, m + 1) + term_yyzz;
          cou(Cart::xxyzzz, 0, m) = wmp(1) * cou(Cart::xxzzz, 0, m + 1);
          cou(Cart::xxzzzz, 0, m) =
              wmp(0) * cou(Cart::xzzzz, 0, m + 1) + term_zzzz;
          cou(Cart::xyyyyy, 0, m) = wmp(0) * cou(Cart::yyyyy, 0, m + 1);
          cou(Cart::xyyyyz, 0, m) = wmp(0) * cou(Cart::yyyyz, 0, m + 1);
          cou(Cart::xyyyzz, 0, m) = wmp(0) * cou(Cart::yyyzz, 0, m + 1);
          cou(Cart::xyyzzz, 0, m) = wmp(0) * cou(Cart::yyzzz, 0, m + 1);
          cou(Cart::xyzzzz, 0, m) = wmp(0) * cou(Cart::yzzzz, 0, m + 1);
          cou(Cart::xzzzzz, 0, m) = wmp(0) * cou(Cart::zzzzz, 0, m + 1);
          cou(Cart::yyyyyy, 0, m) =
              wmp(1) * cou(Cart::yyyyy, 0, m + 1) + 5 * term_yyyy;
          cou(Cart::yyyyyz, 0, m) = wmp(2) * cou(Cart::yyyyy, 0, m + 1);
          cou(Cart::yyyyzz, 0, m) =
              wmp(2) * cou(Cart::yyyyz, 0, m + 1) + term_yyyy;
          cou(Cart::yyyzzz, 0, m) =
              wmp(1) * cou(Cart::yyzzz, 0, m + 1) + 2 * term_yzzz;
          cou(Cart::yyzzzz, 0, m) =
              wmp(1) * cou(Cart::yzzzz, 0, m + 1) + term_zzzz;
          cou(Cart::yzzzzz, 0, m) = wmp(1) * cou(Cart::zzzzz, 0, m + 1);
          cou(Cart::zzzzzz, 0, m) =
              wmp(2) * cou(Cart::zzzzz, 0, m + 1) + 5 * term_zzzz;
        }
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        for (int m = 0; m < lmax_col; m++) {
          cou(0, Cart::x, m) = wmq(0) * cou(0, 0, m + 1);
          cou(0, Cart::y, m) = wmq(1) * cou(0, 0, m + 1);
          cou(0, Cart::z, m) = wmq(2) * cou(0, 0, m + 1);
        }
        //------------------------------------------------------

        // Integrals     p - p
        if (lmax_row > 0) {
          for (int m = 0; m < lmax_col; m++) {
            double term = r_decay * cou(0, 0, m + 1);
            for (int i = 1; i < 4; i++) {
              cou(i, Cart::x, m) = wmq(0) * cou(i, 0, m + 1) + nx[i] * term;
              cou(i, Cart::y, m) = wmq(1) * cou(i, 0, m + 1) + ny[i] * term;
              cou(i, Cart::z, m) = wmq(2) * cou(i, 0, m + 1) + nz[i] * term;
            }
          }
        }
        //------------------------------------------------------

        // Integrals     d - p     f - p     g - p     h - p     i - p
        for (int i_row = 2; i_row < lmax_row + 1; i_row++) {
          for (int m = 0; m < lmax_col; m++) {
            for (int i = 4; i < n_orbitals[lmax_row]; i++) {
              cou(i, Cart::x, m) = wmq(0) * cou(i, 0, m + 1) +
                                   nx[i] * r_decay * cou(i_less_x[i], 0, m + 1);
              cou(i, Cart::y, m) = wmq(1) * cou(i, 0, m + 1) +
                                   ny[i] * r_decay * cou(i_less_y[i], 0, m + 1);
              cou(i, Cart::z, m) = wmq(2) * cou(i, 0, m + 1) +
                                   nz[i] * r_decay * cou(i_less_z[i], 0, m + 1);
            }
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        for (int m = 0; m < lmax_col - 1; m++) {
          double term =
              rdecay_col * (cou(0, 0, m) - fac_a_ac * cou(0, 0, m + 1));
          cou(0, Cart::xx, m) = wmq(0) * cou(0, Cart::x, m + 1) + term;
          cou(0, Cart::xy, m) = wmq(0) * cou(0, Cart::y, m + 1);
          cou(0, Cart::xz, m) = wmq(0) * cou(0, Cart::z, m + 1);
          cou(0, Cart::yy, m) = wmq(1) * cou(0, Cart::y, m + 1) + term;
          cou(0, Cart::yz, m) = wmq(1) * cou(0, Cart::z, m + 1);
          cou(0, Cart::zz, m) = wmq(2) * cou(0, Cart::z, m + 1) + term;
        }
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d     h - d     i - d
        for (int m = 0; m < lmax_col - 1; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term =
                rdecay_col * (cou(i, 0, m) - fac_a_ac * cou(i, 0, m + 1));
            cou(i, Cart::xx, m) =
                wmq(0) * cou(i, Cart::x, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::x, m + 1) + term;
            cou(i, Cart::xy, m) =
                wmq(0) * cou(i, Cart::y, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::y, m + 1);
            cou(i, Cart::xz, m) =
                wmq(0) * cou(i, Cart::z, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::z, m + 1);
            cou(i, Cart::yy, m) =
                wmq(1) * cou(i, Cart::y, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::y, m + 1) + term;
            cou(i, Cart::yz, m) =
                wmq(1) * cou(i, Cart::z, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::z, m + 1);
            cou(i, Cart::zz, m) =
                wmq(2) * cou(i, Cart::z, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::z, m + 1) + term;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        for (int m = 0; m < lmax_col - 2; m++) {
          cou(0, Cart::xxx, m) =
              wmq(0) * cou(0, Cart::xx, m + 1) +
              2 * rdecay_col *
                  (cou(0, Cart::x, m) - fac_a_ac * cou(0, Cart::x, m + 1));
          cou(0, Cart::xxy, m) = wmq(1) * cou(0, Cart::xx, m + 1);
          cou(0, Cart::xxz, m) = wmq(2) * cou(0, Cart::xx, m + 1);
          cou(0, Cart::xyy, m) = wmq(0) * cou(0, Cart::yy, m + 1);
          cou(0, Cart::xyz, m) = wmq(0) * cou(0, Cart::yz, m + 1);
          cou(0, Cart::xzz, m) = wmq(0) * cou(0, Cart::zz, m + 1);
          cou(0, Cart::yyy, m) =
              wmq(1) * cou(0, Cart::yy, m + 1) +
              2 * rdecay_col *
                  (cou(0, Cart::y, m) - fac_a_ac * cou(0, Cart::y, m + 1));
          cou(0, Cart::yyz, m) = wmq(2) * cou(0, Cart::yy, m + 1);
          cou(0, Cart::yzz, m) = wmq(1) * cou(0, Cart::zz, m + 1);
          cou(0, Cart::zzz, m) =
              wmq(2) * cou(0, Cart::zz, m + 1) +
              2 * rdecay_col *
                  (cou(0, Cart::z, m) - fac_a_ac * cou(0, Cart::z, m + 1));
        }
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f     h - f     i - f
        for (int m = 0; m < lmax_col - 2; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_x =
                2 * rdecay_col *
                (cou(i, Cart::x, m) - fac_a_ac * cou(i, Cart::x, m + 1));
            double term_y =
                2 * rdecay_col *
                (cou(i, Cart::y, m) - fac_a_ac * cou(i, Cart::y, m + 1));
            double term_z =
                2 * rdecay_col *
                (cou(i, Cart::z, m) - fac_a_ac * cou(i, Cart::z, m + 1));
            cou(i, Cart::xxx, m) =
                wmq(0) * cou(i, Cart::xx, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xx, m + 1) + term_x;
            cou(i, Cart::xxy, m) =
                wmq(1) * cou(i, Cart::xx, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xx, m + 1);
            cou(i, Cart::xxz, m) =
                wmq(2) * cou(i, Cart::xx, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xx, m + 1);
            cou(i, Cart::xyy, m) =
                wmq(0) * cou(i, Cart::yy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yy, m + 1);
            cou(i, Cart::xyz, m) =
                wmq(0) * cou(i, Cart::yz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yz, m + 1);
            cou(i, Cart::xzz, m) =
                wmq(0) * cou(i, Cart::zz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::zz, m + 1);
            cou(i, Cart::yyy, m) =
                wmq(1) * cou(i, Cart::yy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yy, m + 1) + term_y;
            cou(i, Cart::yyz, m) =
                wmq(2) * cou(i, Cart::yy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yy, m + 1);
            cou(i, Cart::yzz, m) =
                wmq(1) * cou(i, Cart::zz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::zz, m + 1);
            cou(i, Cart::zzz, m) =
                wmq(2) * cou(i, Cart::zz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::zz, m + 1) + term_z;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        for (int m = 0; m < lmax_col - 3; m++) {
          double term_xx = rdecay_col * (cou(0, Cart::xx, m) -
                                         fac_a_ac * cou(0, Cart::xx, m + 1));
          double term_yy = rdecay_col * (cou(0, Cart::yy, m) -
                                         fac_a_ac * cou(0, Cart::yy, m + 1));
          double term_zz = rdecay_col * (cou(0, Cart::zz, m) -
                                         fac_a_ac * cou(0, Cart::zz, m + 1));
          cou(0, Cart::xxxx, m) =
              wmq(0) * cou(0, Cart::xxx, m + 1) + 3 * term_xx;
          cou(0, Cart::xxxy, m) = wmq(1) * cou(0, Cart::xxx, m + 1);
          cou(0, Cart::xxxz, m) = wmq(2) * cou(0, Cart::xxx, m + 1);
          cou(0, Cart::xxyy, m) = wmq(0) * cou(0, Cart::xyy, m + 1) + term_yy;
          cou(0, Cart::xxyz, m) = wmq(1) * cou(0, Cart::xxz, m + 1);
          cou(0, Cart::xxzz, m) = wmq(0) * cou(0, Cart::xzz, m + 1) + term_zz;
          cou(0, Cart::xyyy, m) = wmq(0) * cou(0, Cart::yyy, m + 1);
          cou(0, Cart::xyyz, m) = wmq(0) * cou(0, Cart::yyz, m + 1);
          cou(0, Cart::xyzz, m) = wmq(0) * cou(0, Cart::yzz, m + 1);
          cou(0, Cart::xzzz, m) = wmq(0) * cou(0, Cart::zzz, m + 1);
          cou(0, Cart::yyyy, m) =
              wmq(1) * cou(0, Cart::yyy, m + 1) + 3 * term_yy;
          cou(0, Cart::yyyz, m) = wmq(2) * cou(0, Cart::yyy, m + 1);
          cou(0, Cart::yyzz, m) = wmq(1) * cou(0, Cart::yzz, m + 1) + term_zz;
          cou(0, Cart::yzzz, m) = wmq(1) * cou(0, Cart::zzz, m + 1);
          cou(0, Cart::zzzz, m) =
              wmq(2) * cou(0, Cart::zzz, m + 1) + 3 * term_zz;
        }
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g     h - g     i - g
        for (int m = 0; m < lmax_col - 3; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_xx = rdecay_col * (cou(i, Cart::xx, m) -
                                           fac_a_ac * cou(i, Cart::xx, m + 1));
            double term_yy = rdecay_col * (cou(i, Cart::yy, m) -
                                           fac_a_ac * cou(i, Cart::yy, m + 1));
            double term_zz = rdecay_col * (cou(i, Cart::zz, m) -
                                           fac_a_ac * cou(i, Cart::zz, m + 1));
            cou(i, Cart::xxxx, m) =
                wmq(0) * cou(i, Cart::xxx, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xxx, m + 1) +
                3 * term_xx;
            cou(i, Cart::xxxy, m) =
                wmq(1) * cou(i, Cart::xxx, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxx, m + 1);
            cou(i, Cart::xxxz, m) =
                wmq(2) * cou(i, Cart::xxx, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxx, m + 1);
            cou(i, Cart::xxyy, m) =
                wmq(0) * cou(i, Cart::xyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xyy, m + 1) + term_yy;
            cou(i, Cart::xxyz, m) =
                wmq(1) * cou(i, Cart::xxz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxz, m + 1);
            cou(i, Cart::xxzz, m) =
                wmq(0) * cou(i, Cart::xzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xzz, m + 1) + term_zz;
            cou(i, Cart::xyyy, m) =
                wmq(0) * cou(i, Cart::yyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyy, m + 1);
            cou(i, Cart::xyyz, m) =
                wmq(0) * cou(i, Cart::yyz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyz, m + 1);
            cou(i, Cart::xyzz, m) =
                wmq(0) * cou(i, Cart::yzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yzz, m + 1);
            cou(i, Cart::xzzz, m) =
                wmq(0) * cou(i, Cart::zzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::zzz, m + 1);
            cou(i, Cart::yyyy, m) =
                wmq(1) * cou(i, Cart::yyy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yyy, m + 1) +
                3 * term_yy;
            cou(i, Cart::yyyz, m) =
                wmq(2) * cou(i, Cart::yyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yyy, m + 1);
            cou(i, Cart::yyzz, m) =
                wmq(1) * cou(i, Cart::yzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yzz, m + 1) + term_zz;
            cou(i, Cart::yzzz, m) =
                wmq(1) * cou(i, Cart::zzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::zzz, m + 1);
            cou(i, Cart::zzzz, m) =
                wmq(2) * cou(i, Cart::zzz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::zzz, m + 1) +
                3 * term_zz;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      if (lmax_col > 4) {

        // Integrals     s - h
        for (int m = 0; m < lmax_col - 4; m++) {
          double term_xxx = rdecay_col * (cou(0, Cart::xxx, m) -
                                          fac_a_ac * cou(0, Cart::xxx, m + 1));
          double term_yyy = rdecay_col * (cou(0, Cart::yyy, m) -
                                          fac_a_ac * cou(0, Cart::yyy, m + 1));
          double term_zzz = rdecay_col * (cou(0, Cart::zzz, m) -
                                          fac_a_ac * cou(0, Cart::zzz, m + 1));
          cou(0, Cart::xxxxx, m) =
              wmq(0) * cou(0, Cart::xxxx, m + 1) + 4 * term_xxx;
          cou(0, Cart::xxxxy, m) = wmq(1) * cou(0, Cart::xxxx, m + 1);
          cou(0, Cart::xxxxz, m) = wmq(2) * cou(0, Cart::xxxx, m + 1);
          cou(0, Cart::xxxyy, m) =
              wmq(1) * cou(0, Cart::xxxy, m + 1) + term_xxx;
          cou(0, Cart::xxxyz, m) = wmq(1) * cou(0, Cart::xxxz, m + 1);
          cou(0, Cart::xxxzz, m) =
              wmq(2) * cou(0, Cart::xxxz, m + 1) + term_xxx;
          cou(0, Cart::xxyyy, m) =
              wmq(0) * cou(0, Cart::xyyy, m + 1) + term_yyy;
          cou(0, Cart::xxyyz, m) = wmq(2) * cou(0, Cart::xxyy, m + 1);
          cou(0, Cart::xxyzz, m) = wmq(1) * cou(0, Cart::xxzz, m + 1);
          cou(0, Cart::xxzzz, m) =
              wmq(0) * cou(0, Cart::xzzz, m + 1) + term_zzz;
          cou(0, Cart::xyyyy, m) = wmq(0) * cou(0, Cart::yyyy, m + 1);
          cou(0, Cart::xyyyz, m) = wmq(0) * cou(0, Cart::yyyz, m + 1);
          cou(0, Cart::xyyzz, m) = wmq(0) * cou(0, Cart::yyzz, m + 1);
          cou(0, Cart::xyzzz, m) = wmq(0) * cou(0, Cart::yzzz, m + 1);
          cou(0, Cart::xzzzz, m) = wmq(0) * cou(0, Cart::zzzz, m + 1);
          cou(0, Cart::yyyyy, m) =
              wmq(1) * cou(0, Cart::yyyy, m + 1) + 4 * term_yyy;
          cou(0, Cart::yyyyz, m) = wmq(2) * cou(0, Cart::yyyy, m + 1);
          cou(0, Cart::yyyzz, m) =
              wmq(2) * cou(0, Cart::yyyz, m + 1) + term_yyy;
          cou(0, Cart::yyzzz, m) =
              wmq(1) * cou(0, Cart::yzzz, m + 1) + term_zzz;
          cou(0, Cart::yzzzz, m) = wmq(1) * cou(0, Cart::zzzz, m + 1);
          cou(0, Cart::zzzzz, m) =
              wmq(2) * cou(0, Cart::zzzz, m + 1) + 4 * term_zzz;
        }
        //------------------------------------------------------

        // Integrals     p - h     d - h     f - h     g - h     h - h     i - h
        for (int m = 0; m < lmax_col - 4; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_xxx =
                rdecay_col *
                (cou(i, Cart::xxx, m) - fac_a_ac * cou(i, Cart::xxx, m + 1));
            double term_yyy =
                rdecay_col *
                (cou(i, Cart::yyy, m) - fac_a_ac * cou(i, Cart::yyy, m + 1));
            double term_zzz =
                rdecay_col *
                (cou(i, Cart::zzz, m) - fac_a_ac * cou(i, Cart::zzz, m + 1));
            cou(i, Cart::xxxxx, m) =
                wmq(0) * cou(i, Cart::xxxx, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xxxx, m + 1) +
                4 * term_xxx;
            cou(i, Cart::xxxxy, m) =
                wmq(1) * cou(i, Cart::xxxx, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxx, m + 1);
            cou(i, Cart::xxxxz, m) =
                wmq(2) * cou(i, Cart::xxxx, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxxx, m + 1);
            cou(i, Cart::xxxyy, m) =
                wmq(1) * cou(i, Cart::xxxy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxy, m + 1) +
                term_xxx;
            cou(i, Cart::xxxyz, m) =
                wmq(1) * cou(i, Cart::xxxz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxz, m + 1);
            cou(i, Cart::xxxzz, m) =
                wmq(2) * cou(i, Cart::xxxz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxxz, m + 1) +
                term_xxx;
            cou(i, Cart::xxyyy, m) =
                wmq(0) * cou(i, Cart::xyyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xyyy, m + 1) +
                term_yyy;
            cou(i, Cart::xxyyz, m) =
                wmq(2) * cou(i, Cart::xxyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxyy, m + 1);
            cou(i, Cart::xxyzz, m) =
                wmq(1) * cou(i, Cart::xxzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxzz, m + 1);
            cou(i, Cart::xxzzz, m) =
                wmq(0) * cou(i, Cart::xzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xzzz, m + 1) +
                term_zzz;
            cou(i, Cart::xyyyy, m) =
                wmq(0) * cou(i, Cart::yyyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyyy, m + 1);
            cou(i, Cart::xyyyz, m) =
                wmq(0) * cou(i, Cart::yyyz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyyz, m + 1);
            cou(i, Cart::xyyzz, m) =
                wmq(0) * cou(i, Cart::yyzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyzz, m + 1);
            cou(i, Cart::xyzzz, m) =
                wmq(0) * cou(i, Cart::yzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yzzz, m + 1);
            cou(i, Cart::xzzzz, m) =
                wmq(0) * cou(i, Cart::zzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::zzzz, m + 1);
            cou(i, Cart::yyyyy, m) =
                wmq(1) * cou(i, Cart::yyyy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yyyy, m + 1) +
                4 * term_yyy;
            cou(i, Cart::yyyyz, m) =
                wmq(2) * cou(i, Cart::yyyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yyyy, m + 1);
            cou(i, Cart::yyyzz, m) =
                wmq(2) * cou(i, Cart::yyyz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yyyz, m + 1) +
                term_yyy;
            cou(i, Cart::yyzzz, m) =
                wmq(1) * cou(i, Cart::yzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yzzz, m + 1) +
                term_zzz;
            cou(i, Cart::yzzzz, m) =
                wmq(1) * cou(i, Cart::zzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::zzzz, m + 1);
            cou(i, Cart::zzzzz, m) =
                wmq(2) * cou(i, Cart::zzzz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::zzzz, m + 1) +
                4 * term_zzz;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 4)

      if (lmax_col > 5) {

        // Integrals     s - i
        for (int m = 0; m < lmax_col - 5; m++) {
          double term_xxxx =
              rdecay_col *
              (cou(0, Cart::xxxx, m) - fac_a_ac * cou(0, Cart::xxxx, m + 1));
          double term_xyyy =
              rdecay_col *
              (cou(0, Cart::xyyy, m) - fac_a_ac * cou(0, Cart::xyyy, m + 1));
          double term_xzzz =
              rdecay_col *
              (cou(0, Cart::xzzz, m) - fac_a_ac * cou(0, Cart::xzzz, m + 1));
          double term_yyyy =
              rdecay_col *
              (cou(0, Cart::yyyy, m) - fac_a_ac * cou(0, Cart::yyyy, m + 1));
          double term_yyzz =
              rdecay_col *
              (cou(0, Cart::yyzz, m) - fac_a_ac * cou(0, Cart::yyzz, m + 1));
          double term_yzzz =
              rdecay_col *
              (cou(0, Cart::yzzz, m) - fac_a_ac * cou(0, Cart::yzzz, m + 1));
          double term_zzzz =
              rdecay_col *
              (cou(0, Cart::zzzz, m) - fac_a_ac * cou(0, Cart::zzzz, m + 1));
          cou(0, Cart::xxxxxx, m) =
              wmq(0) * cou(0, Cart::xxxxx, m + 1) + 5 * term_xxxx;
          cou(0, Cart::xxxxxy, m) = wmq(1) * cou(0, Cart::xxxxx, m + 1);
          cou(0, Cart::xxxxxz, m) = wmq(2) * cou(0, Cart::xxxxx, m + 1);
          cou(0, Cart::xxxxyy, m) =
              wmq(1) * cou(0, Cart::xxxxy, m + 1) + term_xxxx;
          cou(0, Cart::xxxxyz, m) = wmq(1) * cou(0, Cart::xxxxz, m + 1);
          cou(0, Cart::xxxxzz, m) =
              wmq(2) * cou(0, Cart::xxxxz, m + 1) + term_xxxx;
          cou(0, Cart::xxxyyy, m) =
              wmq(0) * cou(0, Cart::xxyyy, m + 1) + 2 * term_xyyy;
          cou(0, Cart::xxxyyz, m) = wmq(2) * cou(0, Cart::xxxyy, m + 1);
          cou(0, Cart::xxxyzz, m) = wmq(1) * cou(0, Cart::xxxzz, m + 1);
          cou(0, Cart::xxxzzz, m) =
              wmq(0) * cou(0, Cart::xxzzz, m + 1) + 2 * term_xzzz;
          cou(0, Cart::xxyyyy, m) =
              wmq(0) * cou(0, Cart::xyyyy, m + 1) + term_yyyy;
          cou(0, Cart::xxyyyz, m) = wmq(2) * cou(0, Cart::xxyyy, m + 1);
          cou(0, Cart::xxyyzz, m) =
              wmq(0) * cou(0, Cart::xyyzz, m + 1) + term_yyzz;
          cou(0, Cart::xxyzzz, m) = wmq(1) * cou(0, Cart::xxzzz, m + 1);
          cou(0, Cart::xxzzzz, m) =
              wmq(0) * cou(0, Cart::xzzzz, m + 1) + term_zzzz;
          cou(0, Cart::xyyyyy, m) = wmq(0) * cou(0, Cart::yyyyy, m + 1);
          cou(0, Cart::xyyyyz, m) = wmq(0) * cou(0, Cart::yyyyz, m + 1);
          cou(0, Cart::xyyyzz, m) = wmq(0) * cou(0, Cart::yyyzz, m + 1);
          cou(0, Cart::xyyzzz, m) = wmq(0) * cou(0, Cart::yyzzz, m + 1);
          cou(0, Cart::xyzzzz, m) = wmq(0) * cou(0, Cart::yzzzz, m + 1);
          cou(0, Cart::xzzzzz, m) = wmq(0) * cou(0, Cart::zzzzz, m + 1);
          cou(0, Cart::yyyyyy, m) =
              wmq(1) * cou(0, Cart::yyyyy, m + 1) + 5 * term_yyyy;
          cou(0, Cart::yyyyyz, m) = wmq(2) * cou(0, Cart::yyyyy, m + 1);
          cou(0, Cart::yyyyzz, m) =
              wmq(2) * cou(0, Cart::yyyyz, m + 1) + term_yyyy;
          cou(0, Cart::yyyzzz, m) =
              wmq(1) * cou(0, Cart::yyzzz, m + 1) + 2 * term_yzzz;
          cou(0, Cart::yyzzzz, m) =
              wmq(1) * cou(0, Cart::yzzzz, m + 1) + term_zzzz;
          cou(0, Cart::yzzzzz, m) = wmq(1) * cou(0, Cart::zzzzz, m + 1);
          cou(0, Cart::zzzzzz, m) =
              wmq(2) * cou(0, Cart::zzzzz, m + 1) + 5 * term_zzzz;
        }
        //------------------------------------------------------

        // Integrals     p - i     d - i     f - i     g - i     h - i     i - i
        for (int m = 0; m < lmax_col - 5; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_xxxx =
                rdecay_col *
                (cou(i, Cart::xxxx, m) - fac_a_ac * cou(i, Cart::xxxx, m + 1));
            double term_xyyy =
                rdecay_col *
                (cou(i, Cart::xyyy, m) - fac_a_ac * cou(i, Cart::xyyy, m + 1));
            double term_xzzz =
                rdecay_col *
                (cou(i, Cart::xzzz, m) - fac_a_ac * cou(i, Cart::xzzz, m + 1));
            double term_yyyy =
                rdecay_col *
                (cou(i, Cart::yyyy, m) - fac_a_ac * cou(i, Cart::yyyy, m + 1));
            double term_yyzz =
                rdecay_col *
                (cou(i, Cart::yyzz, m) - fac_a_ac * cou(i, Cart::yyzz, m + 1));
            double term_yzzz =
                rdecay_col *
                (cou(i, Cart::yzzz, m) - fac_a_ac * cou(i, Cart::yzzz, m + 1));
            double term_zzzz =
                rdecay_col *
                (cou(i, Cart::zzzz, m) - fac_a_ac * cou(i, Cart::zzzz, m + 1));
            cou(i, Cart::xxxxxx, m) =
                wmq(0) * cou(i, Cart::xxxxx, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xxxxx, m + 1) +
                5 * term_xxxx;
            cou(i, Cart::xxxxxy, m) =
                wmq(1) * cou(i, Cart::xxxxx, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxxx, m + 1);
            cou(i, Cart::xxxxxz, m) =
                wmq(2) * cou(i, Cart::xxxxx, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxxxx, m + 1);
            cou(i, Cart::xxxxyy, m) =
                wmq(1) * cou(i, Cart::xxxxy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxxy, m + 1) +
                term_xxxx;
            cou(i, Cart::xxxxyz, m) =
                wmq(1) * cou(i, Cart::xxxxz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxxz, m + 1);
            cou(i, Cart::xxxxzz, m) =
                wmq(2) * cou(i, Cart::xxxxz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxxxz, m + 1) +
                term_xxxx;
            cou(i, Cart::xxxyyy, m) =
                wmq(0) * cou(i, Cart::xxyyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xxyyy, m + 1) +
                2 * term_xyyy;
            cou(i, Cart::xxxyyz, m) =
                wmq(2) * cou(i, Cart::xxxyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxxyy, m + 1);
            cou(i, Cart::xxxyzz, m) =
                wmq(1) * cou(i, Cart::xxxzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxxzz, m + 1);
            cou(i, Cart::xxxzzz, m) =
                wmq(0) * cou(i, Cart::xxzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xxzzz, m + 1) +
                2 * term_xzzz;
            cou(i, Cart::xxyyyy, m) =
                wmq(0) * cou(i, Cart::xyyyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xyyyy, m + 1) +
                term_yyyy;
            cou(i, Cart::xxyyyz, m) =
                wmq(2) * cou(i, Cart::xxyyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::xxyyy, m + 1);
            cou(i, Cart::xxyyzz, m) =
                wmq(0) * cou(i, Cart::xyyzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xyyzz, m + 1) +
                term_yyzz;
            cou(i, Cart::xxyzzz, m) =
                wmq(1) * cou(i, Cart::xxzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::xxzzz, m + 1);
            cou(i, Cart::xxzzzz, m) =
                wmq(0) * cou(i, Cart::xzzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::xzzzz, m + 1) +
                term_zzzz;
            cou(i, Cart::xyyyyy, m) =
                wmq(0) * cou(i, Cart::yyyyy, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyyyy, m + 1);
            cou(i, Cart::xyyyyz, m) =
                wmq(0) * cou(i, Cart::yyyyz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyyyz, m + 1);
            cou(i, Cart::xyyyzz, m) =
                wmq(0) * cou(i, Cart::yyyzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyyzz, m + 1);
            cou(i, Cart::xyyzzz, m) =
                wmq(0) * cou(i, Cart::yyzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yyzzz, m + 1);
            cou(i, Cart::xyzzzz, m) =
                wmq(0) * cou(i, Cart::yzzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::yzzzz, m + 1);
            cou(i, Cart::xzzzzz, m) =
                wmq(0) * cou(i, Cart::zzzzz, m + 1) +
                nx[i] * r_decay * cou(i_less_x[i], Cart::zzzzz, m + 1);
            cou(i, Cart::yyyyyy, m) =
                wmq(1) * cou(i, Cart::yyyyy, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yyyyy, m + 1) +
                5 * term_yyyy;
            cou(i, Cart::yyyyyz, m) =
                wmq(2) * cou(i, Cart::yyyyy, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yyyyy, m + 1);
            cou(i, Cart::yyyyzz, m) =
                wmq(2) * cou(i, Cart::yyyyz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::yyyyz, m + 1) +
                term_yyyy;
            cou(i, Cart::yyyzzz, m) =
                wmq(1) * cou(i, Cart::yyzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yyzzz, m + 1) +
                2 * term_yzzz;
            cou(i, Cart::yyzzzz, m) =
                wmq(1) * cou(i, Cart::yzzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::yzzzz, m + 1) +
                term_zzzz;
            cou(i, Cart::yzzzzz, m) =
                wmq(1) * cou(i, Cart::zzzzz, m + 1) +
                ny[i] * r_decay * cou(i_less_y[i], Cart::zzzzz, m + 1);
            cou(i, Cart::zzzzzz, m) =
                wmq(2) * cou(i, Cart::zzzzz, m + 1) +
                nz[i] * r_decay * cou(i_less_z[i], Cart::zzzzz, m + 1) +
                5 * term_zzzz;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 5)

      // put cou(i,j,0) into eigen map
      Eigen::Map<Eigen::MatrixXd> coumat =
          Eigen::Map<Eigen::MatrixXd>(cou.data(), nrows, ncols);

      Eigen::MatrixXd cou_sph =
          AOTransform::getTrafo(gaussian_row).transpose() *
          coumat.bottomRightCorner(shell_row.getCartesianNumFunc(),
                                   shell_col.getCartesianNumFunc()) *
          AOTransform::getTrafo(gaussian_col);
      // save to matrix
      matrix += cou_sph;

    }  // shell_col Gaussians
  }    // shell_row Gaussians
  return;
}

// This converts V into ((S-1/2 V S-1/2)-1/2 S-1/2)T, which is needed to
// construct 4c integrals,
Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                                double etol) {

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(auxoverlap.Matrix());
  removedfunctions = 0;
  Eigen::VectorXd diagonal_overlap =
      Eigen::VectorXd::Zero(eo.eigenvalues().size());
  for (unsigned i = 0; i < diagonal_overlap.size(); ++i) {
    if (eo.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal_overlap(i) = 1.0 / std::sqrt(eo.eigenvalues()(i));
    }
  }
  Eigen::MatrixXd Ssqrt = eo.eigenvectors() * diagonal_overlap.asDiagonal() *
                          eo.eigenvectors().transpose();

  Eigen::MatrixXd ortho = Ssqrt * _aomatrix * Ssqrt;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());

  for (unsigned i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  Eigen::MatrixXd Vm1 =
      es.eigenvectors() * diagonal.asDiagonal() * es.eigenvectors().transpose();
  Eigen::MatrixXd result = (Vm1 * Ssqrt).transpose();
  return result;
}

Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt(double etol) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
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

}  // namespace xtp
}  // namespace votca
