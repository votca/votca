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

#include <map>
#include <string>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>

#include "votca/xtp/qmatom.h"

namespace votca {
namespace xtp {

void AOESP::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                      const AOShell& shell_row, const AOShell& shell_col) {

  const double pi = boost::math::constants::pi<double>();

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row.getLmax();
  int lmax_col = shell_col.getLmax();
  int lsum = lmax_row + lmax_col;
  // set size of internal block for recursion
  int nrows = this->getBlockSize(lmax_row);
  int ncols = this->getBlockSize(lmax_col);

  // initialize local matrix block for unnormalized cartesians
  Eigen::MatrixXd nuc = Eigen::MatrixXd::Zero(nrows, ncols);

  // cout << nuc.size1() << ":" << nuc.size2() << endl;

  /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
   * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION
   * COEFFICIENTS, AND ADD TO matrix(i,j)
   */

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

  // get shell positions
  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;
  // initialize some helper

  double distsq = diff.squaredNorm();

  // iterate over Gaussians in this shell_row
  for (const auto& gaussian_row : shell_row) {
    // iterate over Gaussians in this shell_col
    const double decay_row = gaussian_row.getDecay();

    for (const auto& gaussian_col : shell_col) {
      const double decay_col = gaussian_col.getDecay();

      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;

      double exparg = fak2 * decay_row * decay_col * distsq;

      // check if distance between postions is big, then skip step

      if (exparg > 30.0) {
        continue;
      }

      // some helpers

      const double zeta = decay_row + decay_col;

      const Eigen::Vector3d PmA =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_row;
      const Eigen::Vector3d PmB =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_col;
      const Eigen::Vector3d PmC =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - _r;

      const double U =
          zeta * (PmC(0) * PmC(0) + PmC(1) * PmC(1) + PmC(2) * PmC(2));

      const std::vector<double> _FmU = XIntegrate(lsum + 1, U);
      // cout << endl;

      // (s-s element normiert )
      double prefactor = 2 * sqrt(1.0 / pi) *
                         pow(4.0 * decay_row * decay_col, 0.75) * fak2 *
                         exp(-exparg);
      nuc(Cart::s, Cart::s) = prefactor * _FmU[0];

      typedef boost::multi_array<double, 3> ma_type;
      ma_type nuc3(boost::extents[nrows][ncols][lsum + 1]);
      typedef ma_type::index index;

      for (index i = 0; i < nrows; ++i) {
        for (index j = 0; j < ncols; ++j) {
          for (index k = 0; k < lsum + 1; ++k) {
            nuc3[i][j][k] = 0.;
          }
        }
      }

      for (int i = 0; i < lsum + 1; i++) {  //////////////////////
        nuc3[0][0][i] = prefactor * _FmU[i];
      }

      // Integrals     p - s
      if (lmax_row > 0) {
        for (int m = 0; m < lsum; m++) {
          nuc3[Cart::x][0][m] =
              PmA(0) * nuc3[0][0][m] - PmC(0) * nuc3[0][0][m + 1];
          nuc3[Cart::y][0][m] =
              PmA(1) * nuc3[0][0][m] - PmC(1) * nuc3[0][0][m + 1];
          nuc3[Cart::z][0][m] =
              PmA(2) * nuc3[0][0][m] - PmC(2) * nuc3[0][0][m + 1];
        }
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        for (int m = 0; m < lsum - 1; m++) {
          double term = fak * (nuc3[0][0][m] - nuc3[0][0][m + 1]);
          nuc3[Cart::xx][0][m] = PmA(0) * nuc3[Cart::x][0][m] -
                                 PmC(0) * nuc3[Cart::x][0][m + 1] + term;
          nuc3[Cart::xy][0][m] =
              PmA(0) * nuc3[Cart::y][0][m] - PmC(0) * nuc3[Cart::y][0][m + 1];
          nuc3[Cart::xz][0][m] =
              PmA(0) * nuc3[Cart::z][0][m] - PmC(0) * nuc3[Cart::z][0][m + 1];
          nuc3[Cart::yy][0][m] = PmA(1) * nuc3[Cart::y][0][m] -
                                 PmC(1) * nuc3[Cart::y][0][m + 1] + term;
          nuc3[Cart::yz][0][m] =
              PmA(1) * nuc3[Cart::z][0][m] - PmC(1) * nuc3[Cart::z][0][m + 1];
          nuc3[Cart::zz][0][m] = PmA(2) * nuc3[Cart::z][0][m] -
                                 PmC(2) * nuc3[Cart::z][0][m + 1] + term;
        }
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        for (int m = 0; m < lsum - 2; m++) {
          nuc3[Cart::xxx][0][m] =
              PmA(0) * nuc3[Cart::xx][0][m] -
              PmC(0) * nuc3[Cart::xx][0][m + 1] +
              2 * fak * (nuc3[Cart::x][0][m] - nuc3[Cart::x][0][m + 1]);
          nuc3[Cart::xxy][0][m] =
              PmA(1) * nuc3[Cart::xx][0][m] - PmC(1) * nuc3[Cart::xx][0][m + 1];
          nuc3[Cart::xxz][0][m] =
              PmA(2) * nuc3[Cart::xx][0][m] - PmC(2) * nuc3[Cart::xx][0][m + 1];
          nuc3[Cart::xyy][0][m] =
              PmA(0) * nuc3[Cart::yy][0][m] - PmC(0) * nuc3[Cart::yy][0][m + 1];
          nuc3[Cart::xyz][0][m] =
              PmA(0) * nuc3[Cart::yz][0][m] - PmC(0) * nuc3[Cart::yz][0][m + 1];
          nuc3[Cart::xzz][0][m] =
              PmA(0) * nuc3[Cart::zz][0][m] - PmC(0) * nuc3[Cart::zz][0][m + 1];
          nuc3[Cart::yyy][0][m] =
              PmA(1) * nuc3[Cart::yy][0][m] -
              PmC(1) * nuc3[Cart::yy][0][m + 1] +
              2 * fak * (nuc3[Cart::y][0][m] - nuc3[Cart::y][0][m + 1]);
          nuc3[Cart::yyz][0][m] =
              PmA(2) * nuc3[Cart::yy][0][m] - PmC(2) * nuc3[Cart::yy][0][m + 1];
          nuc3[Cart::yzz][0][m] =
              PmA(1) * nuc3[Cart::zz][0][m] - PmC(1) * nuc3[Cart::zz][0][m + 1];
          nuc3[Cart::zzz][0][m] =
              PmA(2) * nuc3[Cart::zz][0][m] -
              PmC(2) * nuc3[Cart::zz][0][m + 1] +
              2 * fak * (nuc3[Cart::z][0][m] - nuc3[Cart::z][0][m + 1]);
        }
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        for (int m = 0; m < lsum - 3; m++) {
          double term_xx =
              fak * (nuc3[Cart::xx][0][m] - nuc3[Cart::xx][0][m + 1]);
          double term_yy =
              fak * (nuc3[Cart::yy][0][m] - nuc3[Cart::yy][0][m + 1]);
          double term_zz =
              fak * (nuc3[Cart::zz][0][m] - nuc3[Cart::zz][0][m + 1]);
          nuc3[Cart::xxxx][0][m] = PmA(0) * nuc3[Cart::xxx][0][m] -
                                   PmC(0) * nuc3[Cart::xxx][0][m + 1] +
                                   3 * term_xx;
          nuc3[Cart::xxxy][0][m] = PmA(1) * nuc3[Cart::xxx][0][m] -
                                   PmC(1) * nuc3[Cart::xxx][0][m + 1];
          nuc3[Cart::xxxz][0][m] = PmA(2) * nuc3[Cart::xxx][0][m] -
                                   PmC(2) * nuc3[Cart::xxx][0][m + 1];
          nuc3[Cart::xxyy][0][m] = PmA(0) * nuc3[Cart::xyy][0][m] -
                                   PmC(0) * nuc3[Cart::xyy][0][m + 1] + term_yy;
          nuc3[Cart::xxyz][0][m] = PmA(1) * nuc3[Cart::xxz][0][m] -
                                   PmC(1) * nuc3[Cart::xxz][0][m + 1];
          nuc3[Cart::xxzz][0][m] = PmA(0) * nuc3[Cart::xzz][0][m] -
                                   PmC(0) * nuc3[Cart::xzz][0][m + 1] + term_zz;
          nuc3[Cart::xyyy][0][m] = PmA(0) * nuc3[Cart::yyy][0][m] -
                                   PmC(0) * nuc3[Cart::yyy][0][m + 1];
          nuc3[Cart::xyyz][0][m] = PmA(0) * nuc3[Cart::yyz][0][m] -
                                   PmC(0) * nuc3[Cart::yyz][0][m + 1];
          nuc3[Cart::xyzz][0][m] = PmA(0) * nuc3[Cart::yzz][0][m] -
                                   PmC(0) * nuc3[Cart::yzz][0][m + 1];
          nuc3[Cart::xzzz][0][m] = PmA(0) * nuc3[Cart::zzz][0][m] -
                                   PmC(0) * nuc3[Cart::zzz][0][m + 1];
          nuc3[Cart::yyyy][0][m] = PmA(1) * nuc3[Cart::yyy][0][m] -
                                   PmC(1) * nuc3[Cart::yyy][0][m + 1] +
                                   3 * term_yy;
          nuc3[Cart::yyyz][0][m] = PmA(2) * nuc3[Cart::yyy][0][m] -
                                   PmC(2) * nuc3[Cart::yyy][0][m + 1];
          nuc3[Cart::yyzz][0][m] = PmA(1) * nuc3[Cart::yzz][0][m] -
                                   PmC(1) * nuc3[Cart::yzz][0][m + 1] + term_zz;
          nuc3[Cart::yzzz][0][m] = PmA(1) * nuc3[Cart::zzz][0][m] -
                                   PmC(1) * nuc3[Cart::zzz][0][m + 1];
          nuc3[Cart::zzzz][0][m] = PmA(2) * nuc3[Cart::zzz][0][m] -
                                   PmC(2) * nuc3[Cart::zzz][0][m + 1] +
                                   3 * term_zz;
        }
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        for (int m = 0; m < lmax_col; m++) {
          nuc3[0][Cart::x][m] =
              PmB(0) * nuc3[0][0][m] - PmC(0) * nuc3[0][0][m + 1];
          nuc3[0][Cart::y][m] =
              PmB(1) * nuc3[0][0][m] - PmC(1) * nuc3[0][0][m + 1];
          nuc3[0][Cart::z][m] =
              PmB(2) * nuc3[0][0][m] - PmC(2) * nuc3[0][0][m + 1];
        }
        //------------------------------------------------------

        // Integrals     p - p
        if (lmax_row > 0) {
          for (int m = 0; m < lmax_col; m++) {
            double term = fak * (nuc3[0][0][m] - nuc3[0][0][m + 1]);
            for (int i = 1; i < 4; i++) {
              nuc3[i][Cart::x][m] = PmB(0) * nuc3[i][0][m] -
                                    PmC(0) * nuc3[i][0][m + 1] + nx[i] * term;
              nuc3[i][Cart::y][m] = PmB(1) * nuc3[i][0][m] -
                                    PmC(1) * nuc3[i][0][m + 1] + ny[i] * term;
              nuc3[i][Cart::z][m] = PmB(2) * nuc3[i][0][m] -
                                    PmC(2) * nuc3[i][0][m + 1] + nz[i] * term;
            }
          }
        }
        //------------------------------------------------------

        // Integrals     d - p     f - p     g - p
        for (int m = 0; m < lmax_col; m++) {
          for (int i = 4; i < n_orbitals[lmax_row]; i++) {
            nuc3[i][Cart::x][m] =
                PmB(0) * nuc3[i][0][m] - PmC(0) * nuc3[i][0][m + 1] +
                nx[i] * fak *
                    (nuc3[i_less_x[i]][0][m] - nuc3[i_less_x[i]][0][m + 1]);
            nuc3[i][Cart::y][m] =
                PmB(1) * nuc3[i][0][m] - PmC(1) * nuc3[i][0][m + 1] +
                ny[i] * fak *
                    (nuc3[i_less_y[i]][0][m] - nuc3[i_less_y[i]][0][m + 1]);
            nuc3[i][Cart::z][m] =
                PmB(2) * nuc3[i][0][m] - PmC(2) * nuc3[i][0][m + 1] +
                nz[i] * fak *
                    (nuc3[i_less_z[i]][0][m] - nuc3[i_less_z[i]][0][m + 1]);
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        for (int m = 0; m < lmax_col - 1; m++) {
          double term = fak * (nuc3[0][0][m] - nuc3[0][0][m + 1]);
          nuc3[0][Cart::xx][m] = PmB(0) * nuc3[0][Cart::x][m] -
                                 PmC(0) * nuc3[0][Cart::x][m + 1] + term;
          nuc3[0][Cart::xy][m] =
              PmB(0) * nuc3[0][Cart::y][m] - PmC(0) * nuc3[0][Cart::y][m + 1];
          nuc3[0][Cart::xz][m] =
              PmB(0) * nuc3[0][Cart::z][m] - PmC(0) * nuc3[0][Cart::z][m + 1];
          nuc3[0][Cart::yy][m] = PmB(1) * nuc3[0][Cart::y][m] -
                                 PmC(1) * nuc3[0][Cart::y][m + 1] + term;
          nuc3[0][Cart::yz][m] =
              PmB(1) * nuc3[0][Cart::z][m] - PmC(1) * nuc3[0][Cart::z][m + 1];
          nuc3[0][Cart::zz][m] = PmB(2) * nuc3[0][Cart::z][m] -
                                 PmC(2) * nuc3[0][Cart::z][m + 1] + term;
        }
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d
        for (int m = 0; m < lmax_col - 1; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term = fak * (nuc3[i][0][m] - nuc3[i][0][m + 1]);
            nuc3[i][Cart::xx][m] = PmB(0) * nuc3[i][Cart::x][m] -
                                   PmC(0) * nuc3[i][Cart::x][m + 1] +
                                   nx[i] * fak *
                                       (nuc3[i_less_x[i]][Cart::x][m] -
                                        nuc3[i_less_x[i]][Cart::x][m + 1]) +
                                   term;
            nuc3[i][Cart::xy][m] = PmB(0) * nuc3[i][Cart::y][m] -
                                   PmC(0) * nuc3[i][Cart::y][m + 1] +
                                   nx[i] * fak *
                                       (nuc3[i_less_x[i]][Cart::y][m] -
                                        nuc3[i_less_x[i]][Cart::y][m + 1]);
            nuc3[i][Cart::xz][m] = PmB(0) * nuc3[i][Cart::z][m] -
                                   PmC(0) * nuc3[i][Cart::z][m + 1] +
                                   nx[i] * fak *
                                       (nuc3[i_less_x[i]][Cart::z][m] -
                                        nuc3[i_less_x[i]][Cart::z][m + 1]);
            nuc3[i][Cart::yy][m] = PmB(1) * nuc3[i][Cart::y][m] -
                                   PmC(1) * nuc3[i][Cart::y][m + 1] +
                                   ny[i] * fak *
                                       (nuc3[i_less_y[i]][Cart::y][m] -
                                        nuc3[i_less_y[i]][Cart::y][m + 1]) +
                                   term;
            nuc3[i][Cart::yz][m] = PmB(1) * nuc3[i][Cart::z][m] -
                                   PmC(1) * nuc3[i][Cart::z][m + 1] +
                                   ny[i] * fak *
                                       (nuc3[i_less_y[i]][Cart::z][m] -
                                        nuc3[i_less_y[i]][Cart::z][m + 1]);
            nuc3[i][Cart::zz][m] = PmB(2) * nuc3[i][Cart::z][m] -
                                   PmC(2) * nuc3[i][Cart::z][m + 1] +
                                   nz[i] * fak *
                                       (nuc3[i_less_z[i]][Cart::z][m] -
                                        nuc3[i_less_z[i]][Cart::z][m + 1]) +
                                   term;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        for (int m = 0; m < lmax_col - 2; m++) {
          nuc3[0][Cart::xxx][m] =
              PmB(0) * nuc3[0][Cart::xx][m] -
              PmC(0) * nuc3[0][Cart::xx][m + 1] +
              2 * fak * (nuc3[0][Cart::x][m] - nuc3[0][Cart::x][m + 1]);
          nuc3[0][Cart::xxy][m] =
              PmB(1) * nuc3[0][Cart::xx][m] - PmC(1) * nuc3[0][Cart::xx][m + 1];
          nuc3[0][Cart::xxz][m] =
              PmB(2) * nuc3[0][Cart::xx][m] - PmC(2) * nuc3[0][Cart::xx][m + 1];
          nuc3[0][Cart::xyy][m] =
              PmB(0) * nuc3[0][Cart::yy][m] - PmC(0) * nuc3[0][Cart::yy][m + 1];
          nuc3[0][Cart::xyz][m] =
              PmB(0) * nuc3[0][Cart::yz][m] - PmC(0) * nuc3[0][Cart::yz][m + 1];
          nuc3[0][Cart::xzz][m] =
              PmB(0) * nuc3[0][Cart::zz][m] - PmC(0) * nuc3[0][Cart::zz][m + 1];
          nuc3[0][Cart::yyy][m] =
              PmB(1) * nuc3[0][Cart::yy][m] -
              PmC(1) * nuc3[0][Cart::yy][m + 1] +
              2 * fak * (nuc3[0][Cart::y][m] - nuc3[0][Cart::y][m + 1]);
          nuc3[0][Cart::yyz][m] =
              PmB(2) * nuc3[0][Cart::yy][m] - PmC(2) * nuc3[0][Cart::yy][m + 1];
          nuc3[0][Cart::yzz][m] =
              PmB(1) * nuc3[0][Cart::zz][m] - PmC(1) * nuc3[0][Cart::zz][m + 1];
          nuc3[0][Cart::zzz][m] =
              PmB(2) * nuc3[0][Cart::zz][m] -
              PmC(2) * nuc3[0][Cart::zz][m + 1] +
              2 * fak * (nuc3[0][Cart::z][m] - nuc3[0][Cart::z][m + 1]);
        }
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f
        for (int m = 0; m < lmax_col - 2; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_x =
                2 * fak * (nuc3[i][Cart::x][m] - nuc3[i][Cart::x][m + 1]);
            double term_y =
                2 * fak * (nuc3[i][Cart::y][m] - nuc3[i][Cart::y][m + 1]);
            double term_z =
                2 * fak * (nuc3[i][Cart::z][m] - nuc3[i][Cart::z][m + 1]);
            nuc3[i][Cart::xxx][m] = PmB(0) * nuc3[i][Cart::xx][m] -
                                    PmC(0) * nuc3[i][Cart::xx][m + 1] +
                                    nx[i] * fak *
                                        (nuc3[i_less_x[i]][Cart::xx][m] -
                                         nuc3[i_less_x[i]][Cart::xx][m + 1]) +
                                    term_x;
            nuc3[i][Cart::xxy][m] = PmB(1) * nuc3[i][Cart::xx][m] -
                                    PmC(1) * nuc3[i][Cart::xx][m + 1] +
                                    ny[i] * fak *
                                        (nuc3[i_less_y[i]][Cart::xx][m] -
                                         nuc3[i_less_y[i]][Cart::xx][m + 1]);
            nuc3[i][Cart::xxz][m] = PmB(2) * nuc3[i][Cart::xx][m] -
                                    PmC(2) * nuc3[i][Cart::xx][m + 1] +
                                    nz[i] * fak *
                                        (nuc3[i_less_z[i]][Cart::xx][m] -
                                         nuc3[i_less_z[i]][Cart::xx][m + 1]);
            nuc3[i][Cart::xyy][m] = PmB(0) * nuc3[i][Cart::yy][m] -
                                    PmC(0) * nuc3[i][Cart::yy][m + 1] +
                                    nx[i] * fak *
                                        (nuc3[i_less_x[i]][Cart::yy][m] -
                                         nuc3[i_less_x[i]][Cart::yy][m + 1]);
            nuc3[i][Cart::xyz][m] = PmB(0) * nuc3[i][Cart::yz][m] -
                                    PmC(0) * nuc3[i][Cart::yz][m + 1] +
                                    nx[i] * fak *
                                        (nuc3[i_less_x[i]][Cart::yz][m] -
                                         nuc3[i_less_x[i]][Cart::yz][m + 1]);
            nuc3[i][Cart::xzz][m] = PmB(0) * nuc3[i][Cart::zz][m] -
                                    PmC(0) * nuc3[i][Cart::zz][m + 1] +
                                    nx[i] * fak *
                                        (nuc3[i_less_x[i]][Cart::zz][m] -
                                         nuc3[i_less_x[i]][Cart::zz][m + 1]);
            nuc3[i][Cart::yyy][m] = PmB(1) * nuc3[i][Cart::yy][m] -
                                    PmC(1) * nuc3[i][Cart::yy][m + 1] +
                                    ny[i] * fak *
                                        (nuc3[i_less_y[i]][Cart::yy][m] -
                                         nuc3[i_less_y[i]][Cart::yy][m + 1]) +
                                    term_y;
            nuc3[i][Cart::yyz][m] = PmB(2) * nuc3[i][Cart::yy][m] -
                                    PmC(2) * nuc3[i][Cart::yy][m + 1] +
                                    nz[i] * fak *
                                        (nuc3[i_less_z[i]][Cart::yy][m] -
                                         nuc3[i_less_z[i]][Cart::yy][m + 1]);
            nuc3[i][Cart::yzz][m] = PmB(1) * nuc3[i][Cart::zz][m] -
                                    PmC(1) * nuc3[i][Cart::zz][m + 1] +
                                    ny[i] * fak *
                                        (nuc3[i_less_y[i]][Cart::zz][m] -
                                         nuc3[i_less_y[i]][Cart::zz][m + 1]);
            nuc3[i][Cart::zzz][m] = PmB(2) * nuc3[i][Cart::zz][m] -
                                    PmC(2) * nuc3[i][Cart::zz][m + 1] +
                                    nz[i] * fak *
                                        (nuc3[i_less_z[i]][Cart::zz][m] -
                                         nuc3[i_less_z[i]][Cart::zz][m + 1]) +
                                    term_z;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        for (int m = 0; m < lmax_col - 3; m++) {
          double term_xx =
              fak * (nuc3[0][Cart::xx][m] - nuc3[0][Cart::xx][m + 1]);
          double term_yy =
              fak * (nuc3[0][Cart::yy][m] - nuc3[0][Cart::yy][m + 1]);
          double term_zz =
              fak * (nuc3[0][Cart::zz][m] - nuc3[0][Cart::zz][m + 1]);
          nuc3[0][Cart::xxxx][m] = PmB(0) * nuc3[0][Cart::xxx][m] -
                                   PmC(0) * nuc3[0][Cart::xxx][m + 1] +
                                   3 * term_xx;
          nuc3[0][Cart::xxxy][m] = PmB(1) * nuc3[0][Cart::xxx][m] -
                                   PmC(1) * nuc3[0][Cart::xxx][m + 1];
          nuc3[0][Cart::xxxz][m] = PmB(2) * nuc3[0][Cart::xxx][m] -
                                   PmC(2) * nuc3[0][Cart::xxx][m + 1];
          nuc3[0][Cart::xxyy][m] = PmB(0) * nuc3[0][Cart::xyy][m] -
                                   PmC(0) * nuc3[0][Cart::xyy][m + 1] + term_yy;
          nuc3[0][Cart::xxyz][m] = PmB(1) * nuc3[0][Cart::xxz][m] -
                                   PmC(1) * nuc3[0][Cart::xxz][m + 1];
          nuc3[0][Cart::xxzz][m] = PmB(0) * nuc3[0][Cart::xzz][m] -
                                   PmC(0) * nuc3[0][Cart::xzz][m + 1] + term_zz;
          nuc3[0][Cart::xyyy][m] = PmB(0) * nuc3[0][Cart::yyy][m] -
                                   PmC(0) * nuc3[0][Cart::yyy][m + 1];
          nuc3[0][Cart::xyyz][m] = PmB(0) * nuc3[0][Cart::yyz][m] -
                                   PmC(0) * nuc3[0][Cart::yyz][m + 1];
          nuc3[0][Cart::xyzz][m] = PmB(0) * nuc3[0][Cart::yzz][m] -
                                   PmC(0) * nuc3[0][Cart::yzz][m + 1];
          nuc3[0][Cart::xzzz][m] = PmB(0) * nuc3[0][Cart::zzz][m] -
                                   PmC(0) * nuc3[0][Cart::zzz][m + 1];
          nuc3[0][Cart::yyyy][m] = PmB(1) * nuc3[0][Cart::yyy][m] -
                                   PmC(1) * nuc3[0][Cart::yyy][m + 1] +
                                   3 * term_yy;
          nuc3[0][Cart::yyyz][m] = PmB(2) * nuc3[0][Cart::yyy][m] -
                                   PmC(2) * nuc3[0][Cart::yyy][m + 1];
          nuc3[0][Cart::yyzz][m] = PmB(1) * nuc3[0][Cart::yzz][m] -
                                   PmC(1) * nuc3[0][Cart::yzz][m + 1] + term_zz;
          nuc3[0][Cart::yzzz][m] = PmB(1) * nuc3[0][Cart::zzz][m] -
                                   PmC(1) * nuc3[0][Cart::zzz][m + 1];
          nuc3[0][Cart::zzzz][m] = PmB(2) * nuc3[0][Cart::zzz][m] -
                                   PmC(2) * nuc3[0][Cart::zzz][m + 1] +
                                   3 * term_zz;
        }
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g
        for (int m = 0; m < lmax_col - 3; m++) {
          for (int i = 1; i < n_orbitals[lmax_row]; i++) {
            double term_xx =
                fak * (nuc3[i][Cart::xx][m] - nuc3[i][Cart::xx][m + 1]);
            double term_yy =
                fak * (nuc3[i][Cart::yy][m] - nuc3[i][Cart::yy][m + 1]);
            double term_zz =
                fak * (nuc3[i][Cart::zz][m] - nuc3[i][Cart::zz][m + 1]);
            nuc3[i][Cart::xxxx][m] = PmB(0) * nuc3[i][Cart::xxx][m] -
                                     PmC(0) * nuc3[i][Cart::xxx][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::xxx][m] -
                                          nuc3[i_less_x[i]][Cart::xxx][m + 1]) +
                                     3 * term_xx;
            nuc3[i][Cart::xxxy][m] = PmB(1) * nuc3[i][Cart::xxx][m] -
                                     PmC(1) * nuc3[i][Cart::xxx][m + 1] +
                                     ny[i] * fak *
                                         (nuc3[i_less_y[i]][Cart::xxx][m] -
                                          nuc3[i_less_y[i]][Cart::xxx][m + 1]);
            nuc3[i][Cart::xxxz][m] = PmB(2) * nuc3[i][Cart::xxx][m] -
                                     PmC(2) * nuc3[i][Cart::xxx][m + 1] +
                                     nz[i] * fak *
                                         (nuc3[i_less_z[i]][Cart::xxx][m] -
                                          nuc3[i_less_z[i]][Cart::xxx][m + 1]);
            nuc3[i][Cart::xxyy][m] = PmB(0) * nuc3[i][Cart::xyy][m] -
                                     PmC(0) * nuc3[i][Cart::xyy][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::xyy][m] -
                                          nuc3[i_less_x[i]][Cart::xyy][m + 1]) +
                                     term_yy;
            nuc3[i][Cart::xxyz][m] = PmB(1) * nuc3[i][Cart::xxz][m] -
                                     PmC(1) * nuc3[i][Cart::xxz][m + 1] +
                                     ny[i] * fak *
                                         (nuc3[i_less_y[i]][Cart::xxz][m] -
                                          nuc3[i_less_y[i]][Cart::xxz][m + 1]);
            nuc3[i][Cart::xxzz][m] = PmB(0) * nuc3[i][Cart::xzz][m] -
                                     PmC(0) * nuc3[i][Cart::xzz][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::xzz][m] -
                                          nuc3[i_less_x[i]][Cart::xzz][m + 1]) +
                                     term_zz;
            nuc3[i][Cart::xyyy][m] = PmB(0) * nuc3[i][Cart::yyy][m] -
                                     PmC(0) * nuc3[i][Cart::yyy][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::yyy][m] -
                                          nuc3[i_less_x[i]][Cart::yyy][m + 1]);
            nuc3[i][Cart::xyyz][m] = PmB(0) * nuc3[i][Cart::yyz][m] -
                                     PmC(0) * nuc3[i][Cart::yyz][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::yyz][m] -
                                          nuc3[i_less_x[i]][Cart::yyz][m + 1]);
            nuc3[i][Cart::xyzz][m] = PmB(0) * nuc3[i][Cart::yzz][m] -
                                     PmC(0) * nuc3[i][Cart::yzz][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::yzz][m] -
                                          nuc3[i_less_x[i]][Cart::yzz][m + 1]);
            nuc3[i][Cart::xzzz][m] = PmB(0) * nuc3[i][Cart::zzz][m] -
                                     PmC(0) * nuc3[i][Cart::zzz][m + 1] +
                                     nx[i] * fak *
                                         (nuc3[i_less_x[i]][Cart::zzz][m] -
                                          nuc3[i_less_x[i]][Cart::zzz][m + 1]);
            nuc3[i][Cart::yyyy][m] = PmB(1) * nuc3[i][Cart::yyy][m] -
                                     PmC(1) * nuc3[i][Cart::yyy][m + 1] +
                                     ny[i] * fak *
                                         (nuc3[i_less_y[i]][Cart::yyy][m] -
                                          nuc3[i_less_y[i]][Cart::yyy][m + 1]) +
                                     3 * term_yy;
            nuc3[i][Cart::yyyz][m] = PmB(2) * nuc3[i][Cart::yyy][m] -
                                     PmC(2) * nuc3[i][Cart::yyy][m + 1] +
                                     nz[i] * fak *
                                         (nuc3[i_less_z[i]][Cart::yyy][m] -
                                          nuc3[i_less_z[i]][Cart::yyy][m + 1]);
            nuc3[i][Cart::yyzz][m] = PmB(1) * nuc3[i][Cart::yzz][m] -
                                     PmC(1) * nuc3[i][Cart::yzz][m + 1] +
                                     ny[i] * fak *
                                         (nuc3[i_less_y[i]][Cart::yzz][m] -
                                          nuc3[i_less_y[i]][Cart::yzz][m + 1]) +
                                     term_zz;
            nuc3[i][Cart::yzzz][m] = PmB(1) * nuc3[i][Cart::zzz][m] -
                                     PmC(1) * nuc3[i][Cart::zzz][m + 1] +
                                     ny[i] * fak *
                                         (nuc3[i_less_y[i]][Cart::zzz][m] -
                                          nuc3[i_less_y[i]][Cart::zzz][m + 1]);
            nuc3[i][Cart::zzzz][m] = PmB(2) * nuc3[i][Cart::zzz][m] -
                                     PmC(2) * nuc3[i][Cart::zzz][m + 1] +
                                     nz[i] * fak *
                                         (nuc3[i_less_z[i]][Cart::zzz][m] -
                                          nuc3[i_less_z[i]][Cart::zzz][m + 1]) +
                                     3 * term_zz;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
          nuc(i, j) = nuc3[i][j][0];
        }
      }

      Eigen::MatrixXd nuc_sph =
          getTrafo(gaussian_row).transpose() * nuc * getTrafo(gaussian_col);
      // save to matrix

      for (unsigned i = 0; i < matrix.rows(); i++) {
        for (unsigned j = 0; j < matrix.cols(); j++) {
          matrix(i, j) +=
              nuc_sph(i + shell_row.getOffset(), j + shell_col.getOffset());
        }
      }

    }  // shell_col Gaussians
  }    // shell_row Gaussians
  return;
}

// Calculates the electrostatic potential matrix element between two basis
// functions, for an array of atomcores.

void AOESP::Fillnucpotential(const AOBasis& aobasis, const QMMolecule& atoms) {
  _nuclearpotential =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

  for (const auto& atom : atoms) {
    double Znuc = atom.getNuccharge();
    _aomatrix =
        Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
    setPosition(atom.getPos());
    Fill(aobasis);
    _nuclearpotential -= (Znuc)*_aomatrix;
  }
  return;
}

void AOESP::Fillextpotential(
    const AOBasis& aobasis,
    const std::vector<std::unique_ptr<StaticSite> >& externalsites) {

  _externalpotential =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

  for (const std::unique_ptr<StaticSite>& site : externalsites) {
    _aomatrix =
        Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
    setPosition(site->getPos());
    Fill(aobasis);
    _externalpotential -= site->getCharge() * _aomatrix;
  }
  return;
}

}  // namespace xtp
}  // namespace votca
