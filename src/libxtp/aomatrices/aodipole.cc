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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/aotransform.h>
namespace votca {
namespace xtp {

void AODipole::FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd> >& matrix,
                         const AOShell& shell_row,
                         const AOShell& shell_col) const {

  /* Calculating the AO matrix of the gradient operator requires
   * the raw overlap matrix (i.e. in unnormalized cartesians)

   */

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row.getLmax();
  int lmax_col = shell_col.getLmax();

  if (std::max(lmax_col, lmax_row) > 4) {
    throw std::runtime_error(
        "Dipole transition dipoles only implemented for S,P,D,F,G functions in "
        "DFT basis!");
  }

  // set size of internal block for recursion
  int nrows = AOTransform::getBlockSize(lmax_row);
  int ncols = AOTransform::getBlockSize(lmax_col);

  // initialize local matrix block for unnormalized cartesians
  std::array<Eigen::MatrixXd, 3> dip;
  for (int i_comp = 0; i_comp < 3; i_comp++) {
    dip[i_comp] = Eigen::MatrixXd::Zero(nrows, ncols);
  }

  // get shell positions
  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;
  double distsq = diff.squaredNorm();

  std::array<int, 9> n_orbitals = AOTransform::n_orbitals();
  std::array<int, 165> nx = AOTransform::nx();
  std::array<int, 165> ny = AOTransform::ny();
  std::array<int, 165> nz = AOTransform::nz();
  std::array<int, 165> i_less_x = AOTransform::i_less_x();
  std::array<int, 165> i_less_y = AOTransform::i_less_y();
  std::array<int, 165> i_less_z = AOTransform::i_less_z();

  // iterate over Gaussians in this shell_row
  for (const auto& gaussian_row : shell_row) {
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

      const Eigen::Vector3d PmA =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_row;
      const Eigen::Vector3d PmB =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - pos_col;
      const Eigen::Vector3d pmc =
          fak2 * (decay_row * pos_row + decay_col * pos_col) - _r;

      AOOverlap overlap;
      Eigen::MatrixXd ol =
          overlap.Primitive_Overlap(gaussian_row, gaussian_col);

      // s-s dipole moment integrals
      for (int i_comp = 0; i_comp < 3; i_comp++) {
        dip[i_comp](0, 0) = pmc[i_comp] * ol(0, 0);
      }

      // Integrals     p - s
      if (lmax_row > 0) {
        for (int k = 0; k < 3; k++) {
          dip[k](Cart::x, 0) =
              PmA(0) * dip[k](0, 0) + (k == 0) * fak * ol(0, 0);
          dip[k](Cart::y, 0) =
              PmA(1) * dip[k](0, 0) + (k == 1) * fak * ol(0, 0);
          dip[k](Cart::z, 0) =
              PmA(2) * dip[k](0, 0) + (k == 2) * fak * ol(0, 0);
        }
      }
      //------------------------------------------------------

      // Integrals     d - s
      if (lmax_row > 1) {
        for (int k = 0; k < 3; k++) {
          double term = fak * dip[k](0, 0);
          dip[k](Cart::xx, 0) = PmA(0) * dip[k](Cart::x, 0) +
                                (k == 0) * fak * ol(Cart::x, 0) + term;
          dip[k](Cart::xy, 0) =
              PmA(0) * dip[k](Cart::y, 0) + (k == 0) * fak * ol(Cart::y, 0);
          dip[k](Cart::xz, 0) =
              PmA(0) * dip[k](Cart::z, 0) + (k == 0) * fak * ol(Cart::z, 0);
          dip[k](Cart::yy, 0) = PmA(1) * dip[k](Cart::y, 0) +
                                (k == 1) * fak * ol(Cart::y, 0) + term;
          dip[k](Cart::yz, 0) =
              PmA(1) * dip[k](Cart::z, 0) + (k == 1) * fak * ol(Cart::z, 0);
          dip[k](Cart::zz, 0) = PmA(2) * dip[k](Cart::z, 0) +
                                (k == 2) * fak * ol(Cart::z, 0) + term;
        }
      }
      //------------------------------------------------------

      // Integrals     f - s
      if (lmax_row > 2) {
        for (int k = 0; k < 3; k++) {
          dip[k](Cart::xxx, 0) = PmA(0) * dip[k](Cart::xx, 0) +
                                 (k == 0) * fak * ol(Cart::xx, 0) +
                                 2 * fak * dip[k](Cart::x, 0);
          dip[k](Cart::xxy, 0) =
              PmA(1) * dip[k](Cart::xx, 0) + (k == 1) * fak * ol(Cart::xx, 0);
          dip[k](Cart::xxz, 0) =
              PmA(2) * dip[k](Cart::xx, 0) + (k == 2) * fak * ol(Cart::xx, 0);
          dip[k](Cart::xyy, 0) =
              PmA(0) * dip[k](Cart::yy, 0) + (k == 0) * fak * ol(Cart::yy, 0);
          dip[k](Cart::xyz, 0) =
              PmA(0) * dip[k](Cart::yz, 0) + (k == 0) * fak * ol(Cart::yz, 0);
          dip[k](Cart::xzz, 0) =
              PmA(0) * dip[k](Cart::zz, 0) + (k == 0) * fak * ol(Cart::zz, 0);
          dip[k](Cart::yyy, 0) = PmA(1) * dip[k](Cart::yy, 0) +
                                 (k == 1) * fak * ol(Cart::yy, 0) +
                                 2 * fak * dip[k](Cart::y, 0);
          dip[k](Cart::yyz, 0) =
              PmA(2) * dip[k](Cart::yy, 0) + (k == 2) * fak * ol(Cart::yy, 0);
          dip[k](Cart::yzz, 0) =
              PmA(1) * dip[k](Cart::zz, 0) + (k == 1) * fak * ol(Cart::zz, 0);
          dip[k](Cart::zzz, 0) = PmA(2) * dip[k](Cart::zz, 0) +
                                 (k == 2) * fak * ol(Cart::zz, 0) +
                                 2 * fak * dip[k](Cart::z, 0);
        }
      }
      //------------------------------------------------------

      // Integrals     g - s
      if (lmax_row > 3) {
        for (int k = 0; k < 3; k++) {
          double term_xx = fak * dip[k](Cart::xx, 0);
          double term_yy = fak * dip[k](Cart::yy, 0);
          double term_zz = fak * dip[k](Cart::zz, 0);
          dip[k](Cart::xxxx, 0) = PmA(0) * dip[k](Cart::xxx, 0) +
                                  (k == 0) * fak * ol(Cart::xxx, 0) +
                                  3 * term_xx;
          dip[k](Cart::xxxy, 0) =
              PmA(1) * dip[k](Cart::xxx, 0) + (k == 1) * fak * ol(Cart::xxx, 0);
          dip[k](Cart::xxxz, 0) =
              PmA(2) * dip[k](Cart::xxx, 0) + (k == 2) * fak * ol(Cart::xxx, 0);
          dip[k](Cart::xxyy, 0) = PmA(0) * dip[k](Cart::xyy, 0) +
                                  (k == 0) * fak * ol(Cart::xyy, 0) + term_yy;
          dip[k](Cart::xxyz, 0) =
              PmA(1) * dip[k](Cart::xxz, 0) + (k == 1) * fak * ol(Cart::xxz, 0);
          dip[k](Cart::xxzz, 0) = PmA(0) * dip[k](Cart::xzz, 0) +
                                  (k == 0) * fak * ol(Cart::xzz, 0) + term_zz;
          dip[k](Cart::xyyy, 0) =
              PmA(0) * dip[k](Cart::yyy, 0) + (k == 0) * fak * ol(Cart::yyy, 0);
          dip[k](Cart::xyyz, 0) =
              PmA(0) * dip[k](Cart::yyz, 0) + (k == 0) * fak * ol(Cart::yyz, 0);
          dip[k](Cart::xyzz, 0) =
              PmA(0) * dip[k](Cart::yzz, 0) + (k == 0) * fak * ol(Cart::yzz, 0);
          dip[k](Cart::xzzz, 0) =
              PmA(0) * dip[k](Cart::zzz, 0) + (k == 0) * fak * ol(Cart::zzz, 0);
          dip[k](Cart::yyyy, 0) = PmA(1) * dip[k](Cart::yyy, 0) +
                                  (k == 1) * fak * ol(Cart::yyy, 0) +
                                  3 * term_yy;
          dip[k](Cart::yyyz, 0) =
              PmA(2) * dip[k](Cart::yyy, 0) + (k == 2) * fak * ol(Cart::yyy, 0);
          dip[k](Cart::yyzz, 0) = PmA(1) * dip[k](Cart::yzz, 0) +
                                  (k == 1) * fak * ol(Cart::yzz, 0) + term_zz;
          dip[k](Cart::yzzz, 0) =
              PmA(1) * dip[k](Cart::zzz, 0) + (k == 1) * fak * ol(Cart::zzz, 0);
          dip[k](Cart::zzzz, 0) = PmA(2) * dip[k](Cart::zzz, 0) +
                                  (k == 2) * fak * ol(Cart::zzz, 0) +
                                  3 * term_zz;
        }
      }
      //------------------------------------------------------

      if (lmax_col > 0) {

        // Integrals     s - p
        for (int k = 0; k < 3; k++) {
          dip[k](0, Cart::x) =
              PmB(0) * dip[k](0, 0) + (k == 0) * fak * ol(0, 0);
          dip[k](0, Cart::y) =
              PmB(1) * dip[k](0, 0) + (k == 1) * fak * ol(0, 0);
          dip[k](0, Cart::z) =
              PmB(2) * dip[k](0, 0) + (k == 2) * fak * ol(0, 0);
        }
        //------------------------------------------------------

        // Integrals     p - p     d - p     f - p     g - p
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          for (int k = 0; k < 3; k++) {
            dip[k](i, Cart::x) = PmB(0) * dip[k](i, 0) +
                                 (k == 0) * fak * ol(i, 0) +
                                 nx[i] * fak * dip[k](i_less_x[i], 0);
            dip[k](i, Cart::y) = PmB(1) * dip[k](i, 0) +
                                 (k == 1) * fak * ol(i, 0) +
                                 ny[i] * fak * dip[k](i_less_y[i], 0);
            dip[k](i, Cart::z) = PmB(2) * dip[k](i, 0) +
                                 (k == 2) * fak * ol(i, 0) +
                                 nz[i] * fak * dip[k](i_less_z[i], 0);
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 0)

      if (lmax_col > 1) {

        // Integrals     s - d
        for (int k = 0; k < 3; k++) {
          double term = fak * dip[k](0, 0);
          dip[k](0, Cart::xx) = PmB(0) * dip[k](0, Cart::x) +
                                (k == 0) * fak * ol(0, Cart::x) + term;
          dip[k](0, Cart::xy) =
              PmB(0) * dip[k](0, Cart::y) + (k == 0) * fak * ol(0, Cart::y);
          dip[k](0, Cart::xz) =
              PmB(0) * dip[k](0, Cart::z) + (k == 0) * fak * ol(0, Cart::z);
          dip[k](0, Cart::yy) = PmB(1) * dip[k](0, Cart::y) +
                                (k == 1) * fak * ol(0, Cart::y) + term;
          dip[k](0, Cart::yz) =
              PmB(1) * dip[k](0, Cart::z) + (k == 1) * fak * ol(0, Cart::z);
          dip[k](0, Cart::zz) = PmB(2) * dip[k](0, Cart::z) +
                                (k == 2) * fak * ol(0, Cart::z) + term;
        }
        //------------------------------------------------------

        // Integrals     p - d     d - d     f - d     g - d
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          for (int k = 0; k < 3; k++) {
            double term = fak * dip[k](i, 0);
            dip[k](i, Cart::xx) =
                PmB(0) * dip[k](i, Cart::x) + (k == 0) * fak * ol(i, Cart::x) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::x) + term;
            dip[k](i, Cart::xy) = PmB(0) * dip[k](i, Cart::y) +
                                  (k == 0) * fak * ol(i, Cart::y) +
                                  nx[i] * fak * dip[k](i_less_x[i], Cart::y);
            dip[k](i, Cart::xz) = PmB(0) * dip[k](i, Cart::z) +
                                  (k == 0) * fak * ol(i, Cart::z) +
                                  nx[i] * fak * dip[k](i_less_x[i], Cart::z);
            dip[k](i, Cart::yy) =
                PmB(1) * dip[k](i, Cart::y) + (k == 1) * fak * ol(i, Cart::y) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::y) + term;
            dip[k](i, Cart::yz) = PmB(1) * dip[k](i, Cart::z) +
                                  (k == 1) * fak * ol(i, Cart::z) +
                                  ny[i] * fak * dip[k](i_less_y[i], Cart::z);
            dip[k](i, Cart::zz) =
                PmB(2) * dip[k](i, Cart::z) + (k == 2) * fak * ol(i, Cart::z) +
                nz[i] * fak * dip[k](i_less_z[i], Cart::z) + term;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 1)

      if (lmax_col > 2) {

        // Integrals     s - f
        for (int k = 0; k < 3; k++) {
          dip[k](0, Cart::xxx) = PmB(0) * dip[k](0, Cart::xx) +
                                 (k == 0) * fak * ol(0, Cart::xx) +
                                 2 * fak * dip[k](0, Cart::x);
          dip[k](0, Cart::xxy) =
              PmB(1) * dip[k](0, Cart::xx) + (k == 1) * fak * ol(0, Cart::xx);
          dip[k](0, Cart::xxz) =
              PmB(2) * dip[k](0, Cart::xx) + (k == 2) * fak * ol(0, Cart::xx);
          dip[k](0, Cart::xyy) =
              PmB(0) * dip[k](0, Cart::yy) + (k == 0) * fak * ol(0, Cart::yy);
          dip[k](0, Cart::xyz) =
              PmB(0) * dip[k](0, Cart::yz) + (k == 0) * fak * ol(0, Cart::yz);
          dip[k](0, Cart::xzz) =
              PmB(0) * dip[k](0, Cart::zz) + (k == 0) * fak * ol(0, Cart::zz);
          dip[k](0, Cart::yyy) = PmB(1) * dip[k](0, Cart::yy) +
                                 (k == 1) * fak * ol(0, Cart::yy) +
                                 2 * fak * dip[k](0, Cart::y);
          dip[k](0, Cart::yyz) =
              PmB(2) * dip[k](0, Cart::yy) + (k == 2) * fak * ol(0, Cart::yy);
          dip[k](0, Cart::yzz) =
              PmB(1) * dip[k](0, Cart::zz) + (k == 1) * fak * ol(0, Cart::zz);
          dip[k](0, Cart::zzz) = PmB(2) * dip[k](0, Cart::zz) +
                                 (k == 2) * fak * ol(0, Cart::zz) +
                                 2 * fak * dip[k](0, Cart::z);
        }
        //------------------------------------------------------

        // Integrals     p - f     d - f     f - f     g - f
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          for (int k = 0; k < 3; k++) {
            double term_x = 2 * fak * dip[k](i, Cart::x);
            double term_y = 2 * fak * dip[k](i, Cart::y);
            double term_z = 2 * fak * dip[k](i, Cart::z);
            dip[k](i, Cart::xxx) = PmB(0) * dip[k](i, Cart::xx) +
                                   (k == 0) * fak * ol(i, Cart::xx) +
                                   nx[i] * fak * dip[k](i_less_x[i], Cart::xx) +
                                   term_x;
            dip[k](i, Cart::xxy) = PmB(1) * dip[k](i, Cart::xx) +
                                   (k == 1) * fak * ol(i, Cart::xx) +
                                   ny[i] * fak * dip[k](i_less_y[i], Cart::xx);
            dip[k](i, Cart::xxz) = PmB(2) * dip[k](i, Cart::xx) +
                                   (k == 2) * fak * ol(i, Cart::xx) +
                                   nz[i] * fak * dip[k](i_less_z[i], Cart::xx);
            dip[k](i, Cart::xyy) = PmB(0) * dip[k](i, Cart::yy) +
                                   (k == 0) * fak * ol(i, Cart::yy) +
                                   nx[i] * fak * dip[k](i_less_x[i], Cart::yy);
            dip[k](i, Cart::xyz) = PmB(0) * dip[k](i, Cart::yz) +
                                   (k == 0) * fak * ol(i, Cart::yz) +
                                   nx[i] * fak * dip[k](i_less_x[i], Cart::yz);
            dip[k](i, Cart::xzz) = PmB(0) * dip[k](i, Cart::zz) +
                                   (k == 0) * fak * ol(i, Cart::zz) +
                                   nx[i] * fak * dip[k](i_less_x[i], Cart::zz);
            dip[k](i, Cart::yyy) = PmB(1) * dip[k](i, Cart::yy) +
                                   (k == 1) * fak * ol(i, Cart::yy) +
                                   ny[i] * fak * dip[k](i_less_y[i], Cart::yy) +
                                   term_y;
            dip[k](i, Cart::yyz) = PmB(2) * dip[k](i, Cart::yy) +
                                   (k == 2) * fak * ol(i, Cart::yy) +
                                   nz[i] * fak * dip[k](i_less_z[i], Cart::yy);
            dip[k](i, Cart::yzz) = PmB(1) * dip[k](i, Cart::zz) +
                                   (k == 1) * fak * ol(i, Cart::zz) +
                                   ny[i] * fak * dip[k](i_less_y[i], Cart::zz);
            dip[k](i, Cart::zzz) = PmB(2) * dip[k](i, Cart::zz) +
                                   (k == 2) * fak * ol(i, Cart::zz) +
                                   nz[i] * fak * dip[k](i_less_z[i], Cart::zz) +
                                   term_z;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 2)

      if (lmax_col > 3) {

        // Integrals     s - g
        for (int k = 0; k < 3; k++) {
          double term_xx = fak * dip[k](0, Cart::xx);
          double term_yy = fak * dip[k](0, Cart::yy);
          double term_zz = fak * dip[k](0, Cart::zz);
          dip[k](0, Cart::xxxx) = PmB(0) * dip[k](0, Cart::xxx) +
                                  (k == 0) * fak * ol(0, Cart::xxx) +
                                  3 * term_xx;
          dip[k](0, Cart::xxxy) =
              PmB(1) * dip[k](0, Cart::xxx) + (k == 1) * fak * ol(0, Cart::xxx);
          dip[k](0, Cart::xxxz) =
              PmB(2) * dip[k](0, Cart::xxx) + (k == 2) * fak * ol(0, Cart::xxx);
          dip[k](0, Cart::xxyy) = PmB(0) * dip[k](0, Cart::xyy) +
                                  (k == 0) * fak * ol(0, Cart::xyy) + term_yy;
          dip[k](0, Cart::xxyz) =
              PmB(1) * dip[k](0, Cart::xxz) + (k == 1) * fak * ol(0, Cart::xxz);
          dip[k](0, Cart::xxzz) = PmB(0) * dip[k](0, Cart::xzz) +
                                  (k == 0) * fak * ol(0, Cart::xzz) + term_zz;
          dip[k](0, Cart::xyyy) =
              PmB(0) * dip[k](0, Cart::yyy) + (k == 0) * fak * ol(0, Cart::yyy);
          dip[k](0, Cart::xyyz) =
              PmB(0) * dip[k](0, Cart::yyz) + (k == 0) * fak * ol(0, Cart::yyz);
          dip[k](0, Cart::xyzz) =
              PmB(0) * dip[k](0, Cart::yzz) + (k == 0) * fak * ol(0, Cart::yzz);
          dip[k](0, Cart::xzzz) =
              PmB(0) * dip[k](0, Cart::zzz) + (k == 0) * fak * ol(0, Cart::zzz);
          dip[k](0, Cart::yyyy) = PmB(1) * dip[k](0, Cart::yyy) +
                                  (k == 1) * fak * ol(0, Cart::yyy) +
                                  3 * term_yy;
          dip[k](0, Cart::yyyz) =
              PmB(2) * dip[k](0, Cart::yyy) + (k == 2) * fak * ol(0, Cart::yyy);
          dip[k](0, Cart::yyzz) = PmB(1) * dip[k](0, Cart::yzz) +
                                  (k == 1) * fak * ol(0, Cart::yzz) + term_zz;
          dip[k](0, Cart::yzzz) =
              PmB(1) * dip[k](0, Cart::zzz) + (k == 1) * fak * ol(0, Cart::zzz);
          dip[k](0, Cart::zzzz) = PmB(2) * dip[k](0, Cart::zzz) +
                                  (k == 2) * fak * ol(0, Cart::zzz) +
                                  3 * term_zz;
        }
        //------------------------------------------------------

        // Integrals     p - g     d - g     f - g     g - g
        for (int i = 1; i < n_orbitals[lmax_row]; i++) {
          for (int k = 0; k < 3; k++) {
            double term_xx = fak * dip[k](i, Cart::xx);
            double term_yy = fak * dip[k](i, Cart::yy);
            double term_zz = fak * dip[k](i, Cart::zz);
            dip[k](i, Cart::xxxx) =
                PmB(0) * dip[k](i, Cart::xxx) +
                (k == 0) * fak * ol(i, Cart::xxx) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::xxx) + 3 * term_xx;
            dip[k](i, Cart::xxxy) =
                PmB(1) * dip[k](i, Cart::xxx) +
                (k == 1) * fak * ol(i, Cart::xxx) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::xxx);
            dip[k](i, Cart::xxxz) =
                PmB(2) * dip[k](i, Cart::xxx) +
                (k == 2) * fak * ol(i, Cart::xxx) +
                nz[i] * fak * dip[k](i_less_z[i], Cart::xxx);
            dip[k](i, Cart::xxyy) =
                PmB(0) * dip[k](i, Cart::xyy) +
                (k == 0) * fak * ol(i, Cart::xyy) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::xyy) + term_yy;
            dip[k](i, Cart::xxyz) =
                PmB(1) * dip[k](i, Cart::xxz) +
                (k == 1) * fak * ol(i, Cart::xxz) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::xxz);
            dip[k](i, Cart::xxzz) =
                PmB(0) * dip[k](i, Cart::xzz) +
                (k == 0) * fak * ol(i, Cart::xzz) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::xzz) + term_zz;
            dip[k](i, Cart::xyyy) =
                PmB(0) * dip[k](i, Cart::yyy) +
                (k == 0) * fak * ol(i, Cart::yyy) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::yyy);
            dip[k](i, Cart::xyyz) =
                PmB(0) * dip[k](i, Cart::yyz) +
                (k == 0) * fak * ol(i, Cart::yyz) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::yyz);
            dip[k](i, Cart::xyzz) =
                PmB(0) * dip[k](i, Cart::yzz) +
                (k == 0) * fak * ol(i, Cart::yzz) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::yzz);
            dip[k](i, Cart::xzzz) =
                PmB(0) * dip[k](i, Cart::zzz) +
                (k == 0) * fak * ol(i, Cart::zzz) +
                nx[i] * fak * dip[k](i_less_x[i], Cart::zzz);
            dip[k](i, Cart::yyyy) =
                PmB(1) * dip[k](i, Cart::yyy) +
                (k == 1) * fak * ol(i, Cart::yyy) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::yyy) + 3 * term_yy;
            dip[k](i, Cart::yyyz) =
                PmB(2) * dip[k](i, Cart::yyy) +
                (k == 2) * fak * ol(i, Cart::yyy) +
                nz[i] * fak * dip[k](i_less_z[i], Cart::yyy);
            dip[k](i, Cart::yyzz) =
                PmB(1) * dip[k](i, Cart::yzz) +
                (k == 1) * fak * ol(i, Cart::yzz) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::yzz) + term_zz;
            dip[k](i, Cart::yzzz) =
                PmB(1) * dip[k](i, Cart::zzz) +
                (k == 1) * fak * ol(i, Cart::zzz) +
                ny[i] * fak * dip[k](i_less_y[i], Cart::zzz);
            dip[k](i, Cart::zzzz) =
                PmB(2) * dip[k](i, Cart::zzz) +
                (k == 2) * fak * ol(i, Cart::zzz) +
                nz[i] * fak * dip[k](i_less_z[i], Cart::zzz) + 3 * term_zz;
          }
        }
        //------------------------------------------------------

      }  // end if (lmax_col > 3)

      Eigen::MatrixXd trafo_row = AOTransform::getTrafo(gaussian_row);
      Eigen::MatrixXd trafo_col = AOTransform::getTrafo(gaussian_col);

      // cartesian -> spherical

      for (int i = 0; i < 3; i++) {
        Eigen::MatrixXd dip_sph =
            trafo_row.transpose() *
            dip[i].bottomRightCorner(shell_row.getCartesianNumFunc(),
                                     shell_col.getCartesianNumFunc()) *
            trafo_col;
        matrix[i] += dip_sph;
      }

    }  // shell_col Gaussians
  }    // shell_row Gaussians
}

}  // namespace xtp
}  // namespace votca
