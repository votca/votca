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
 * distributed under the License is distributed on an "AS iS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aotransform.h>
#include <votca/xtp/fourcenter.h>

namespace votca {
namespace xtp {

/*
 * Calculate 4-center electron repulsion integrals
 *    R_{abcd} = int{ phi_a(r) phi_b(r) phi_c(r') phi_d(r') /(r-r') d3rd3r' }
 * for a given set of a b c d as in
 *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
 * section iI.B for cartesian Gaussians, then transforming
 * to spherical (angular momentum) Gaussians ("complete" shells
 * from S to Lmax, and finally cutting out those angular momentum
 * components actually present in shell-shell-shell combination.
 * Currently supported for
 *      S,P,D,F,G  functions in DFT basis
 *
 */

bool FCMatrix::FillFourCenterRepBlock(Eigen::Tensor<double, 4>& block,
                                      const AOShell& shell_1,
                                      const AOShell& shell_2,
                                      const AOShell& shell_3,
                                      const AOShell& shell_4) const {

  const double pi = boost::math::constants::pi<double>();

  bool does_contribute = true;

  // shell info, only lmax tells how far to go

  Index lmax_1 = shell_1.getLmax();
  Index lmax_2 = shell_2.getLmax();
  Index lmax_3 = shell_3.getLmax();
  Index lmax_4 = shell_4.getLmax();

  Index mmax = lmax_1 + lmax_2 + lmax_3 + lmax_4;

  // set size of internal block for recursion

  const AOShell* shell_alpha;
  const AOShell* shell_beta;
  const AOShell* shell_gamma;
  const AOShell* shell_delta;
  bool alphabetaswitch = false;
  bool gammadeltaswitch = false;
  bool ab_cd_switch = false;

  // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we
  // calculate (ps,s), due to symmetry they are the same.

  if (lmax_1 < lmax_2) {
    alphabetaswitch = true;
  }
  if (lmax_3 < lmax_4) {
    gammadeltaswitch = true;
  }
  if ((lmax_1 + lmax_2) < (lmax_3 + lmax_4)) {
    ab_cd_switch = true;
  }

  if (ab_cd_switch == true) {

    if (alphabetaswitch == true) {
      shell_alpha = &shell_4;
      shell_beta = &shell_3;
    } else {
      shell_alpha = &shell_3;
      shell_beta = &shell_4;
    }

    if (gammadeltaswitch == true) {
      shell_gamma = &shell_2;
      shell_delta = &shell_1;
    } else {
      shell_gamma = &shell_1;
      shell_delta = &shell_2;
    }

  } else {

    if (alphabetaswitch == true) {
      shell_alpha = &shell_2;
      shell_beta = &shell_1;
    } else {
      shell_alpha = &shell_1;
      shell_beta = &shell_2;
    }

    if (gammadeltaswitch == true) {
      shell_gamma = &shell_4;
      shell_delta = &shell_3;
    } else {
      shell_gamma = &shell_3;
      shell_delta = &shell_4;
    }
  }

  const Eigen::Vector3d& pos_alpha = shell_alpha->getPos();
  const Eigen::Vector3d& pos_beta = shell_beta->getPos();
  const Eigen::Vector3d& pos_gamma = shell_gamma->getPos();
  const Eigen::Vector3d& pos_delta = shell_delta->getPos();

  Index lmax_alpha = shell_alpha->getLmax();
  Index lmax_beta = shell_beta->getLmax();
  Index lmax_gamma = shell_gamma->getLmax();
  Index lmax_delta = shell_delta->getLmax();

  std::array<int, 9> n_orbitals = AOTransform::n_orbitals();
  std::array<int, 165> nx = AOTransform::nx();
  std::array<int, 165> ny = AOTransform::ny();
  std::array<int, 165> nz = AOTransform::nz();
  std::array<int, 165> i_less_x = AOTransform::i_less_x();
  std::array<int, 165> i_less_y = AOTransform::i_less_y();
  std::array<int, 165> i_less_z = AOTransform::i_less_z();
  std::array<int, 120> i_more_x = AOTransform::i_more_x();
  std::array<int, 120> i_more_y = AOTransform::i_more_y();
  std::array<int, 120> i_more_z = AOTransform::i_more_z();

  Index nbeta = AOTransform::getBlockSize(lmax_beta);
  Index ndelta = AOTransform::getBlockSize(lmax_delta);
  Index ncombined_ab = AOTransform::getBlockSize(lmax_alpha + lmax_beta);
  Index ncombined_cd = AOTransform::getBlockSize(lmax_gamma + lmax_delta);

  double dist_AB = (pos_alpha - pos_beta).squaredNorm();
  double dist_CD = (pos_gamma - pos_delta).squaredNorm();

  Eigen::Vector3d amb = pos_alpha - pos_beta;

  Eigen::Vector3d cmd = pos_gamma - pos_delta;

  for (const auto& gaussian_alpha : *shell_alpha) {
    const double decay_alpha = gaussian_alpha.getDecay();

    for (const auto& gaussian_beta : *shell_beta) {
      const double decay_beta = gaussian_beta.getDecay();

      for (const auto& gaussian_gamma : *shell_gamma) {
        const double decay_gamma = gaussian_gamma.getDecay();

        for (const auto& gaussian_delta : *shell_delta) {
          const double decay_delta = gaussian_delta.getDecay();

          double zeta = decay_alpha + decay_beta;
          double eta = decay_gamma + decay_delta;
          double decay = zeta + eta;
          double rho = (zeta * eta) / decay;
          double rzeta = 0.5 / zeta;
          double reta = 0.5 / eta;
          double rdecay = 0.5 / decay;
          double gfak = eta / decay;
          double cfak = zeta / decay;
          Eigen::Vector3d P =
              (decay_alpha * pos_alpha + decay_beta * pos_beta) / zeta;
          Eigen::Vector3d Q =
              (decay_gamma * pos_gamma + decay_delta * pos_delta) / eta;
          Eigen::Vector3d W = (zeta * P + eta * Q) / decay;
          double U = rho * (P - Q).squaredNorm();

          Eigen::Vector3d pma = P - pos_alpha;
          Eigen::Vector3d qmc = Q - pos_gamma;
          Eigen::Vector3d wmp = W - P;
          Eigen::Vector3d wmq = W - Q;

          Eigen::Tensor<double, 3> R_temp(ncombined_ab, ncombined_cd, mmax + 1);
          R_temp.setZero();

          const Eigen::VectorXd FmT = AOTransform::XIntegrate(mmax + 1, U);

          double exp_AB = exp(-2. * decay_alpha * decay_beta * rzeta * dist_AB);
          double exp_CD = exp(-2. * decay_gamma * decay_delta * reta * dist_CD);
          double ssss =
              (16. *
               pow(decay_alpha * decay_beta * decay_gamma * decay_delta, .75) *
               exp_AB * exp_CD) /
              (zeta * eta * sqrt(pi * decay));

          // ss integrals
          for (Index i = 0; i < mmax + 1; i++) {
            R_temp(0, 0, i) = ssss * FmT[i];
          }

          Index lmax_alpha_beta = lmax_alpha + lmax_beta;
          Index lmax_gamma_delta = lmax_gamma + lmax_delta;

          // Integrals     p-s - s-s
          if (lmax_alpha_beta > 0) {
            for (Index m = 0; m < mmax; m++) {
              R_temp(Cart::x, 0, m) =
                  pma(0) * R_temp(0, 0, m) + wmp(0) * R_temp(0, 0, m + 1);
              R_temp(Cart::y, 0, m) =
                  pma(1) * R_temp(0, 0, m) + wmp(1) * R_temp(0, 0, m + 1);
              R_temp(Cart::z, 0, m) =
                  pma(2) * R_temp(0, 0, m) + wmp(2) * R_temp(0, 0, m + 1);
            }
          }
          //------------------------------------------------------

          // Integrals     d-s - s-s
          if (lmax_alpha_beta > 1) {
            for (Index m = 0; m < mmax - 1; m++) {
              double term =
                  rzeta * (R_temp(0, 0, m) - gfak * R_temp(0, 0, m + 1));
              R_temp(Cart::xx, 0, m) = pma(0) * R_temp(Cart::x, 0, m) +
                                       wmp(0) * R_temp(Cart::x, 0, m + 1) +
                                       term;
              R_temp(Cart::xy, 0, m) = pma(0) * R_temp(Cart::y, 0, m) +
                                       wmp(0) * R_temp(Cart::y, 0, m + 1);
              R_temp(Cart::xz, 0, m) = pma(0) * R_temp(Cart::z, 0, m) +
                                       wmp(0) * R_temp(Cart::z, 0, m + 1);
              R_temp(Cart::yy, 0, m) = pma(1) * R_temp(Cart::y, 0, m) +
                                       wmp(1) * R_temp(Cart::y, 0, m + 1) +
                                       term;
              R_temp(Cart::yz, 0, m) = pma(1) * R_temp(Cart::z, 0, m) +
                                       wmp(1) * R_temp(Cart::z, 0, m + 1);
              R_temp(Cart::zz, 0, m) = pma(2) * R_temp(Cart::z, 0, m) +
                                       wmp(2) * R_temp(Cart::z, 0, m + 1) +
                                       term;
            }
          }
          //------------------------------------------------------

          // Integrals     f-s - s-s
          if (lmax_alpha_beta > 2) {
            for (Index m = 0; m < mmax - 2; m++) {
              R_temp(Cart::xxx, 0, m) = pma(0) * R_temp(Cart::xx, 0, m) +
                                        wmp(0) * R_temp(Cart::xx, 0, m + 1) +
                                        2 * rzeta *
                                            (R_temp(Cart::x, 0, m) -
                                             gfak * R_temp(Cart::x, 0, m + 1));
              R_temp(Cart::xxy, 0, m) = pma(1) * R_temp(Cart::xx, 0, m) +
                                        wmp(1) * R_temp(Cart::xx, 0, m + 1);
              R_temp(Cart::xxz, 0, m) = pma(2) * R_temp(Cart::xx, 0, m) +
                                        wmp(2) * R_temp(Cart::xx, 0, m + 1);
              R_temp(Cart::xyy, 0, m) = pma(0) * R_temp(Cart::yy, 0, m) +
                                        wmp(0) * R_temp(Cart::yy, 0, m + 1);
              R_temp(Cart::xyz, 0, m) = pma(0) * R_temp(Cart::yz, 0, m) +
                                        wmp(0) * R_temp(Cart::yz, 0, m + 1);
              R_temp(Cart::xzz, 0, m) = pma(0) * R_temp(Cart::zz, 0, m) +
                                        wmp(0) * R_temp(Cart::zz, 0, m + 1);
              R_temp(Cart::yyy, 0, m) = pma(1) * R_temp(Cart::yy, 0, m) +
                                        wmp(1) * R_temp(Cart::yy, 0, m + 1) +
                                        2 * rzeta *
                                            (R_temp(Cart::y, 0, m) -
                                             gfak * R_temp(Cart::y, 0, m + 1));
              R_temp(Cart::yyz, 0, m) = pma(2) * R_temp(Cart::yy, 0, m) +
                                        wmp(2) * R_temp(Cart::yy, 0, m + 1);
              R_temp(Cart::yzz, 0, m) = pma(1) * R_temp(Cart::zz, 0, m) +
                                        wmp(1) * R_temp(Cart::zz, 0, m + 1);
              R_temp(Cart::zzz, 0, m) = pma(2) * R_temp(Cart::zz, 0, m) +
                                        wmp(2) * R_temp(Cart::zz, 0, m + 1) +
                                        2 * rzeta *
                                            (R_temp(Cart::z, 0, m) -
                                             gfak * R_temp(Cart::z, 0, m + 1));
            }
          }
          //------------------------------------------------------

          // Integrals     g-s - s-s     h-s - s-s     i-s - s-s     j-s - s-s
          // k-s - s-s     . . .
          for (Index l = 4; l < lmax_alpha_beta + 1; l++) {
            int norb = n_orbitals[l];
            int norb_1 = n_orbitals[l - 1];
            int norb_2 = n_orbitals[l - 2];
            int norb_3 = n_orbitals[l - 3];
            Index ncart_1 = (l * (l + 1)) / 2;
            Index ncart_2 = ncart_1 - l;
            Index ncart_3 = ncart_2 + 1 - l;
            for (Index m = 0; m < mmax + 1 - l; m++) {
              R_temp(norb_1, 0, m) =
                  pma(0) * R_temp(norb_2, 0, m) +
                  wmp(0) * R_temp(norb_2, 0, m + 1) +
                  double(l - 1) * rzeta *
                      (R_temp(norb_3, 0, m) - gfak * R_temp(norb_3, 0, m + 1));
              R_temp(norb_1 + 1, 0, m) = pma(1) * R_temp(norb_2, 0, m) +
                                         wmp(1) * R_temp(norb_2, 0, m + 1);
              R_temp(norb_1 + 2, 0, m) = pma(2) * R_temp(norb_2, 0, m) +
                                         wmp(2) * R_temp(norb_2, 0, m + 1);
              Index ntimes = 3;
              Index itimes = 3;
              for (Index k = 3; k < ncart_2; k++) {
                R_temp(norb_1 + k, 0, m) =
                    pma(0) * R_temp(norb_2 + k, 0, m) +
                    wmp(0) * R_temp(norb_2 + k, 0, m + 1) +
                    double(l - ntimes) * rzeta *
                        (R_temp(norb_3 + k, 0, m) -
                         gfak * R_temp(norb_3 + k, 0, m + 1));
                itimes--;
                if (itimes == 0) {
                  ntimes++;
                  itimes = ntimes;
                }
              }
              for (Index k = 0; k < l - 1; k++) {
                Index k2 = norb_2 + ncart_2 + k;
                R_temp(norb_1 + ncart_2 + k, 0, m) =
                    pma(0) * R_temp(k2, 0, m) + wmp(0) * R_temp(k2, 0, m + 1);
                R_temp(norb_1 + ncart_1 + k, 0, m) =
                    pma(1) * R_temp(k2, 0, m) + wmp(1) * R_temp(k2, 0, m + 1) +
                    double(l - 1 - k) * rzeta *
                        (R_temp(norb_3 + ncart_3 + k, 0, m) -
                         gfak * R_temp(norb_3 + ncart_3 + k, 0, m + 1));
              }
              R_temp(norb_1 + ncart_2 + l - 1, 0, m) =
                  pma(0) * R_temp(norb_2 + ncart_2 + l - 1, 0, m) +
                  wmp(0) * R_temp(norb_2 + ncart_2 + l - 1, 0, m + 1);
              R_temp(norb - 2, 0, m) = pma(1) * R_temp(norb_1 - 1, 0, m) +
                                       wmp(1) * R_temp(norb_1 - 1, 0, m + 1);
              R_temp(norb - 1, 0, m) =
                  pma(2) * R_temp(norb_1 - 1, 0, m) +
                  wmp(2) * R_temp(norb_1 - 1, 0, m + 1) +
                  double(l - 1) * rzeta *
                      (R_temp(norb_2 - 1, 0, m) -
                       gfak * R_temp(norb_2 - 1, 0, m + 1));
            }
          }
          //------------------------------------------------------

          if (lmax_gamma_delta > 0) {

            // Integrals     s-s - p-s
            for (Index m = 0; m < lmax_gamma_delta; m++) {
              R_temp(0, Cart::x, m) =
                  qmc(0) * R_temp(0, 0, m) + wmq(0) * R_temp(0, 0, m + 1);
              R_temp(0, Cart::y, m) =
                  qmc(1) * R_temp(0, 0, m) + wmq(1) * R_temp(0, 0, m + 1);
              R_temp(0, Cart::z, m) =
                  qmc(2) * R_temp(0, 0, m) + wmq(2) * R_temp(0, 0, m + 1);
            }
            //------------------------------------------------------

            // Integrals     p-s - p-s
            if (lmax_alpha_beta > 0) {
              for (Index m = 0; m < lmax_gamma_delta; m++) {
                double term = rdecay * R_temp(0, 0, m + 1);
                for (Index i = 1; i < 4; i++) {
                  R_temp(i, Cart::x, m) = qmc(0) * R_temp(i, 0, m) +
                                          wmq(0) * R_temp(i, 0, m + 1) +
                                          nx[i] * term;
                  R_temp(i, Cart::y, m) = qmc(1) * R_temp(i, 0, m) +
                                          wmq(1) * R_temp(i, 0, m + 1) +
                                          ny[i] * term;
                  R_temp(i, Cart::z, m) = qmc(2) * R_temp(i, 0, m) +
                                          wmq(2) * R_temp(i, 0, m + 1) +
                                          nz[i] * term;
                }
              }
            }
            //------------------------------------------------------

            // Integrals     d-s - p-s     f-s - p-s     g-s - p-s     h-s - p-s
            // i-s - p-s     j-s - p-s     k-s - p-s     . . .
            for (Index m = 0; m < lmax_gamma_delta; m++) {
              for (Index i = 4; i < ncombined_ab; i++) {
                R_temp(i, Cart::x, m) =
                    qmc(0) * R_temp(i, 0, m) + wmq(0) * R_temp(i, 0, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], 0, m + 1);
                R_temp(i, Cart::y, m) =
                    qmc(1) * R_temp(i, 0, m) + wmq(1) * R_temp(i, 0, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], 0, m + 1);
                R_temp(i, Cart::z, m) =
                    qmc(2) * R_temp(i, 0, m) + wmq(2) * R_temp(i, 0, m + 1) +
                    nz[i] * rdecay * R_temp(i_less_z[i], 0, m + 1);
              }
            }
            //------------------------------------------------------

          }  // end if (lmax_gamma_delta > 0)

          if (lmax_gamma_delta > 1) {

            // Integrals     s-s - d-s
            for (Index m = 0; m < lmax_gamma_delta - 1; m++) {
              double term =
                  reta * (R_temp(0, 0, m) - cfak * R_temp(0, 0, m + 1));
              R_temp(0, Cart::xx, m) = qmc(0) * R_temp(0, Cart::x, m) +
                                       wmq(0) * R_temp(0, Cart::x, m + 1) +
                                       term;
              R_temp(0, Cart::xy, m) = qmc(0) * R_temp(0, Cart::y, m) +
                                       wmq(0) * R_temp(0, Cart::y, m + 1);
              R_temp(0, Cart::xz, m) = qmc(0) * R_temp(0, Cart::z, m) +
                                       wmq(0) * R_temp(0, Cart::z, m + 1);
              R_temp(0, Cart::yy, m) = qmc(1) * R_temp(0, Cart::y, m) +
                                       wmq(1) * R_temp(0, Cart::y, m + 1) +
                                       term;
              R_temp(0, Cart::yz, m) = qmc(1) * R_temp(0, Cart::z, m) +
                                       wmq(1) * R_temp(0, Cart::z, m + 1);
              R_temp(0, Cart::zz, m) = qmc(2) * R_temp(0, Cart::z, m) +
                                       wmq(2) * R_temp(0, Cart::z, m + 1) +
                                       term;
            }
            //------------------------------------------------------

            // Integrals     p-s - d-s     d-s - d-s     f-s - d-s     g-s - d-s
            // h-s - d-s     i-s - d-s     j-s - d-s     k-s - d-s     . . .
            for (Index m = 0; m < lmax_gamma_delta - 1; m++) {
              for (Index i = 1; i < ncombined_ab; i++) {
                double term =
                    reta * (R_temp(i, 0, m) - cfak * R_temp(i, 0, m + 1));
                R_temp(i, Cart::xx, m) =
                    qmc(0) * R_temp(i, Cart::x, m) +
                    wmq(0) * R_temp(i, Cart::x, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::x, m + 1) + term;
                R_temp(i, Cart::xy, m) =
                    qmc(0) * R_temp(i, Cart::y, m) +
                    wmq(0) * R_temp(i, Cart::y, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::y, m + 1);
                R_temp(i, Cart::xz, m) =
                    qmc(0) * R_temp(i, Cart::z, m) +
                    wmq(0) * R_temp(i, Cart::z, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::z, m + 1);
                R_temp(i, Cart::yy, m) =
                    qmc(1) * R_temp(i, Cart::y, m) +
                    wmq(1) * R_temp(i, Cart::y, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], Cart::y, m + 1) + term;
                R_temp(i, Cart::yz, m) =
                    qmc(1) * R_temp(i, Cart::z, m) +
                    wmq(1) * R_temp(i, Cart::z, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], Cart::z, m + 1);
                R_temp(i, Cart::zz, m) =
                    qmc(2) * R_temp(i, Cart::z, m) +
                    wmq(2) * R_temp(i, Cart::z, m + 1) +
                    nz[i] * rdecay * R_temp(i_less_z[i], Cart::z, m + 1) + term;
              }
            }
            //------------------------------------------------------
          }  // end if (lmax_gamma_delta > 1)

          if (lmax_gamma_delta > 2) {

            // Integrals     s-s - f-s
            for (Index m = 0; m < lmax_gamma_delta - 2; m++) {
              R_temp(0, Cart::xxx, m) = qmc(0) * R_temp(0, Cart::xx, m) +
                                        wmq(0) * R_temp(0, Cart::xx, m + 1) +
                                        2 * reta *
                                            (R_temp(0, Cart::x, m) -
                                             cfak * R_temp(0, Cart::x, m + 1));
              R_temp(0, Cart::xxy, m) = qmc(1) * R_temp(0, Cart::xx, m) +
                                        wmq(1) * R_temp(0, Cart::xx, m + 1);
              R_temp(0, Cart::xxz, m) = qmc(2) * R_temp(0, Cart::xx, m) +
                                        wmq(2) * R_temp(0, Cart::xx, m + 1);
              R_temp(0, Cart::xyy, m) = qmc(0) * R_temp(0, Cart::yy, m) +
                                        wmq(0) * R_temp(0, Cart::yy, m + 1);
              R_temp(0, Cart::xyz, m) = qmc(0) * R_temp(0, Cart::yz, m) +
                                        wmq(0) * R_temp(0, Cart::yz, m + 1);
              R_temp(0, Cart::xzz, m) = qmc(0) * R_temp(0, Cart::zz, m) +
                                        wmq(0) * R_temp(0, Cart::zz, m + 1);
              R_temp(0, Cart::yyy, m) = qmc(1) * R_temp(0, Cart::yy, m) +
                                        wmq(1) * R_temp(0, Cart::yy, m + 1) +
                                        2 * reta *
                                            (R_temp(0, Cart::y, m) -
                                             cfak * R_temp(0, Cart::y, m + 1));
              R_temp(0, Cart::yyz, m) = qmc(2) * R_temp(0, Cart::yy, m) +
                                        wmq(2) * R_temp(0, Cart::yy, m + 1);
              R_temp(0, Cart::yzz, m) = qmc(1) * R_temp(0, Cart::zz, m) +
                                        wmq(1) * R_temp(0, Cart::zz, m + 1);
              R_temp(0, Cart::zzz, m) = qmc(2) * R_temp(0, Cart::zz, m) +
                                        wmq(2) * R_temp(0, Cart::zz, m + 1) +
                                        2 * reta *
                                            (R_temp(0, Cart::z, m) -
                                             cfak * R_temp(0, Cart::z, m + 1));
            }
            //------------------------------------------------------

            // Integrals     p-s - f-s     d-s - f-s     f-s - f-s     g-s - f-s
            // h-s - f-s     i-s - f-s     j-s - f-s     k-s - f-s     . . .
            for (Index m = 0; m < lmax_gamma_delta - 2; m++) {
              for (Index i = 1; i < ncombined_ab; i++) {
                double term_x =
                    2 * reta *
                    (R_temp(i, Cart::x, m) - cfak * R_temp(i, Cart::x, m + 1));
                double term_y =
                    2 * reta *
                    (R_temp(i, Cart::y, m) - cfak * R_temp(i, Cart::y, m + 1));
                double term_z =
                    2 * reta *
                    (R_temp(i, Cart::z, m) - cfak * R_temp(i, Cart::z, m + 1));
                R_temp(i, Cart::xxx, m) =
                    qmc(0) * R_temp(i, Cart::xx, m) +
                    wmq(0) * R_temp(i, Cart::xx, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::xx, m + 1) +
                    term_x;
                R_temp(i, Cart::xxy, m) =
                    qmc(1) * R_temp(i, Cart::xx, m) +
                    wmq(1) * R_temp(i, Cart::xx, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], Cart::xx, m + 1);
                R_temp(i, Cart::xxz, m) =
                    qmc(2) * R_temp(i, Cart::xx, m) +
                    wmq(2) * R_temp(i, Cart::xx, m + 1) +
                    nz[i] * rdecay * R_temp(i_less_z[i], Cart::xx, m + 1);
                R_temp(i, Cart::xyy, m) =
                    qmc(0) * R_temp(i, Cart::yy, m) +
                    wmq(0) * R_temp(i, Cart::yy, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::yy, m + 1);
                R_temp(i, Cart::xyz, m) =
                    qmc(0) * R_temp(i, Cart::yz, m) +
                    wmq(0) * R_temp(i, Cart::yz, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::yz, m + 1);
                R_temp(i, Cart::xzz, m) =
                    qmc(0) * R_temp(i, Cart::zz, m) +
                    wmq(0) * R_temp(i, Cart::zz, m + 1) +
                    nx[i] * rdecay * R_temp(i_less_x[i], Cart::zz, m + 1);
                R_temp(i, Cart::yyy, m) =
                    qmc(1) * R_temp(i, Cart::yy, m) +
                    wmq(1) * R_temp(i, Cart::yy, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], Cart::yy, m + 1) +
                    term_y;
                R_temp(i, Cart::yyz, m) =
                    qmc(2) * R_temp(i, Cart::yy, m) +
                    wmq(2) * R_temp(i, Cart::yy, m + 1) +
                    nz[i] * rdecay * R_temp(i_less_z[i], Cart::yy, m + 1);
                R_temp(i, Cart::yzz, m) =
                    qmc(1) * R_temp(i, Cart::zz, m) +
                    wmq(1) * R_temp(i, Cart::zz, m + 1) +
                    ny[i] * rdecay * R_temp(i_less_y[i], Cart::zz, m + 1);
                R_temp(i, Cart::zzz, m) =
                    qmc(2) * R_temp(i, Cart::zz, m) +
                    wmq(2) * R_temp(i, Cart::zz, m + 1) +
                    nz[i] * rdecay * R_temp(i_less_z[i], Cart::zz, m + 1) +
                    term_z;
              }
            }
            //------------------------------------------------------
          }  // end if (lmax_gamma_delta > 2)

          // Integrals     s-s - g-s     p-s - g-s     d-s - g-s     f-s - g-s
          // g-s - g-s     h-s - g-s     i-s - g-s     j-s - g-s     k-s - g-s .
          // . .
          //              s-s - h-s     p-s - h-s     d-s - h-s     f-s - h-s
          //              g-s - h-s     h-s - h-s     i-s - h-s     j-s - h-s
          //              k-s - h-s     . . . s-s - i-s     p-s - i-s     d-s -
          //              i-s     f-s - i-s     g-s - i-s     h-s - i-s     i-s
          //              - i-s     j-s - i-s     k-s - i-s     . . .
          //                    j             j             j             j j j
          //                    j             j             j       . . . . . .
          //                    .             .             .             . . .
          //                    . . . .             .             . . . . . . .
          //                    . . .
          for (Index l = 4; l < lmax_gamma_delta + 1; l++) {
            int norb = n_orbitals[l];
            int norb_1 = n_orbitals[l - 1];
            int norb_2 = n_orbitals[l - 2];
            int norb_3 = n_orbitals[l - 3];
            Index ncart_1 = (l * (l + 1)) / 2;
            Index ncart_2 = ncart_1 - l;
            Index ncart_3 = ncart_2 + 1 - l;

            for (Index m = 0; m < lmax_gamma_delta + 1 - l; m++) {
              R_temp(0, norb_1, m) =
                  qmc(0) * R_temp(0, norb_2, m) +
                  wmq(0) * R_temp(0, norb_2, m + 1) +
                  double(l - 1) * reta *
                      (R_temp(0, norb_3, m) - cfak * R_temp(0, norb_3, m + 1));
              R_temp(0, norb_1 + 1, m) = qmc(1) * R_temp(0, norb_2, m) +
                                         wmq(1) * R_temp(0, norb_2, m + 1);
              R_temp(0, norb_1 + 2, m) = qmc(2) * R_temp(0, norb_2, m) +
                                         wmq(2) * R_temp(0, norb_2, m + 1);
              Index ntimes = 3;
              Index itimes = 3;
              for (Index k = 3; k < ncart_2; k++) {
                R_temp(0, norb_1 + k, m) =
                    qmc(0) * R_temp(0, norb_2 + k, m) +
                    wmq(0) * R_temp(0, norb_2 + k, m + 1) +
                    double(l - ntimes) * reta *
                        (R_temp(0, norb_3 + k, m) -
                         cfak * R_temp(0, norb_3 + k, m + 1));
                itimes--;
                if (itimes == 0) {
                  ntimes++;
                  itimes = ntimes;
                }
              }
              for (Index k = 0; k < l - 1; k++) {
                R_temp(0, norb_1 + ncart_2 + k, m) =
                    qmc(0) * R_temp(0, norb_2 + ncart_2 + k, m) +
                    wmq(0) * R_temp(0, norb_2 + ncart_2 + k, m + 1);
                R_temp(0, norb_1 + ncart_1 + k, m) =
                    qmc(1) * R_temp(0, norb_2 + ncart_2 + k, m) +
                    wmq(1) * R_temp(0, norb_2 + ncart_2 + k, m + 1) +
                    double(l - 1 - k) * reta *
                        (R_temp(0, norb_3 + ncart_3 + k, m) -
                         cfak * R_temp(0, norb_3 + ncart_3 + k, m + 1));
              }
              R_temp(0, norb_1 + ncart_2 + l - 1, m) =
                  qmc(0) * R_temp(0, norb_2 + ncart_2 + l - 1, m) +
                  wmq(0) * R_temp(0, norb_2 + ncart_2 + l - 1, m + 1);
              R_temp(0, norb - 2, m) = qmc(1) * R_temp(0, norb_1 - 1, m) +
                                       wmq(1) * R_temp(0, norb_1 - 1, m + 1);
              R_temp(0, norb - 1, m) =
                  qmc(2) * R_temp(0, norb_1 - 1, m) +
                  wmq(2) * R_temp(0, norb_1 - 1, m + 1) +
                  double(l - 1) * reta *
                      (R_temp(0, norb_2 - 1, m) -
                       cfak * R_temp(0, norb_2 - 1, m + 1));
            }

            for (Index m = 0; m < lmax_gamma_delta + 1 - l; m++) {
              for (Index i = 1; i < ncombined_ab; i++) {
                int nx_i = nx[i];
                int ny_i = ny[i];
                int nz_i = nz[i];
                int ilx_i = i_less_x[i];
                int ily_i = i_less_y[i];
                int ilz_i = i_less_z[i];

                R_temp(i, norb_1, m) =
                    qmc(0) * R_temp(i, norb_2, m) +
                    wmq(0) * R_temp(i, norb_2, m + 1) +
                    nx_i * rdecay * R_temp(ilx_i, norb_2, m + 1) +
                    double(l - 1) * reta *
                        (R_temp(i, norb_3, m) -
                         cfak * R_temp(i, norb_3, m + 1));
                R_temp(i, norb_1 + 1, m) =
                    qmc(1) * R_temp(i, norb_2, m) +
                    wmq(1) * R_temp(i, norb_2, m + 1) +
                    ny_i * rdecay * R_temp(ily_i, norb_2, m + 1);
                R_temp(i, norb_1 + 2, m) =
                    qmc(2) * R_temp(i, norb_2, m) +
                    wmq(2) * R_temp(i, norb_2, m + 1) +
                    nz_i * rdecay * R_temp(ilz_i, norb_2, m + 1);
                Index ntimes = 3;
                Index itimes = 3;
                for (Index k = 3; k < ncart_2; k++) {
                  R_temp(i, norb_1 + k, m) =
                      qmc(0) * R_temp(i, norb_2 + k, m) +
                      wmq(0) * R_temp(i, norb_2 + k, m + 1) +
                      nx_i * rdecay * R_temp(ilx_i, norb_2 + k, m + 1) +
                      double(l - ntimes) * reta *
                          (R_temp(i, norb_3 + k, m) -
                           cfak * R_temp(i, norb_3 + k, m + 1));
                  itimes--;
                  if (itimes == 0) {
                    ntimes++;
                    itimes = ntimes;
                  }
                }
                for (Index k = 0; k < l - 1; k++) {
                  Index k2 = norb_2 + ncart_2 + k;
                  R_temp(i, norb_1 + ncart_2 + k, m) =
                      qmc(0) * R_temp(i, k2, m) +
                      wmq(0) * R_temp(i, k2, m + 1) +
                      nx_i * rdecay *
                          R_temp(ilx_i, norb_2 + ncart_2 + k, m + 1);
                  R_temp(i, norb_1 + ncart_1 + k, m) =
                      qmc(1) * R_temp(i, k2, m) +
                      wmq(1) * R_temp(i, k2, m + 1) +
                      ny_i * rdecay *
                          R_temp(ily_i, norb_2 + ncart_2 + k, m + 1) +
                      double(l - 1 - k) * reta *
                          (R_temp(i, norb_3 + ncart_3 + k, m) -
                           cfak * R_temp(i, norb_3 + ncart_3 + k, m + 1));
                }
                R_temp(i, norb_1 + ncart_2 + l - 1, m) =
                    qmc(0) * R_temp(i, norb_2 + ncart_2 + l - 1, m) +
                    wmq(0) * R_temp(i, norb_2 + ncart_2 + l - 1, m + 1) +
                    nx_i * rdecay *
                        R_temp(ilx_i, norb_2 + ncart_2 + l - 1, m + 1);
                R_temp(i, norb - 2, m) =
                    qmc(1) * R_temp(i, norb_1 - 1, m) +
                    wmq(1) * R_temp(i, norb_1 - 1, m + 1) +
                    ny_i * rdecay * R_temp(ily_i, norb_1 - 1, m + 1);
                R_temp(i, norb - 1, m) =
                    qmc(2) * R_temp(i, norb_1 - 1, m) +
                    wmq(2) * R_temp(i, norb_1 - 1, m + 1) +
                    nz_i * rdecay * R_temp(ilz_i, norb_1 - 1, m + 1) +
                    double(l - 1) * reta *
                        (R_temp(i, norb_2 - 1, m) -
                         cfak * R_temp(i, norb_2 - 1, m + 1));
              }
            }
          }
          //------------------------------------------------------

          // copy into new array for 3D use.
          Eigen::Tensor<double, 3> R(ncombined_ab, nbeta, ncombined_cd);
          R.setZero();
          for (Index k = 0; k < ncombined_cd; ++k) {
            for (Index i = 0; i < ncombined_ab; ++i) {
              R(i, 0, k) = R_temp(i, k, 0);
            }
          }

          if (lmax_beta > 0) {
            // Integrals     s-p - *-s     p-p - *-s     d-p - *-s     f-p - *-s
            // g-p - *-s     h-p - *-s     i-p - *-s     j-p - *-s     . . .
            for (Index i = 0; i < n_orbitals[lmax_alpha_beta - 1]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index j = 0; j < ncombined_cd; j++) {
                R(i, Cart::x, j) = R(imx_i, 0, j) + amb(0) * R(i, 0, j);
                R(i, Cart::y, j) = R(imy_i, 0, j) + amb(1) * R(i, 0, j);
                R(i, Cart::z, j) = R(imz_i, 0, j) + amb(2) * R(i, 0, j);
              }
            }
            //------------------------------------------------------
          }

          if (lmax_beta > 1) {
            // Integrals     s-d - *-s     p-d - *-s     d-d - *-s     f-d - *-s
            // g-d - *-s     h-d - *-s     i-d - *-s     . . .
            for (Index i = 0; i < n_orbitals[lmax_alpha_beta - 2]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index j = 0; j < ncombined_cd; j++) {
                R(i, Cart::xx, j) =
                    R(imx_i, Cart::x, j) + amb(0) * R(i, Cart::x, j);
                R(i, Cart::xy, j) =
                    R(imx_i, Cart::y, j) + amb(0) * R(i, Cart::y, j);
                R(i, Cart::xz, j) =
                    R(imx_i, Cart::z, j) + amb(0) * R(i, Cart::z, j);
                R(i, Cart::yy, j) =
                    R(imy_i, Cart::y, j) + amb(1) * R(i, Cart::y, j);
                R(i, Cart::yz, j) =
                    R(imy_i, Cart::z, j) + amb(1) * R(i, Cart::z, j);
                R(i, Cart::zz, j) =
                    R(imz_i, Cart::z, j) + amb(2) * R(i, Cart::z, j);
              }
            }
            //------------------------------------------------------
          }

          if (lmax_beta > 2) {
            // Integrals     s-f - *-s     p-f - *-s     d-f - *-s     f-f - *-s
            // g-f - *-s     h-f - *-s     . . .
            for (Index i = 0; i < n_orbitals[lmax_alpha_beta - 3]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index j = 0; j < ncombined_cd; j++) {
                R(i, Cart::xxx, j) =
                    R(imx_i, Cart::xx, j) + amb(0) * R(i, Cart::xx, j);
                R(i, Cart::xxy, j) =
                    R(imx_i, Cart::xy, j) + amb(0) * R(i, Cart::xy, j);
                R(i, Cart::xxz, j) =
                    R(imx_i, Cart::xz, j) + amb(0) * R(i, Cart::xz, j);
                R(i, Cart::xyy, j) =
                    R(imx_i, Cart::yy, j) + amb(0) * R(i, Cart::yy, j);
                R(i, Cart::xyz, j) =
                    R(imx_i, Cart::yz, j) + amb(0) * R(i, Cart::yz, j);
                R(i, Cart::xzz, j) =
                    R(imx_i, Cart::zz, j) + amb(0) * R(i, Cart::zz, j);
                R(i, Cart::yyy, j) =
                    R(imy_i, Cart::yy, j) + amb(1) * R(i, Cart::yy, j);
                R(i, Cart::yyz, j) =
                    R(imy_i, Cart::yz, j) + amb(1) * R(i, Cart::yz, j);
                R(i, Cart::yzz, j) =
                    R(imy_i, Cart::zz, j) + amb(1) * R(i, Cart::zz, j);
                R(i, Cart::zzz, j) =
                    R(imz_i, Cart::zz, j) + amb(2) * R(i, Cart::zz, j);
              }
            }
            //------------------------------------------------------
          }

          // Integrals     s-g - *-s     p-g - *-s     d-g - *-s     f-g - *-s
          // g-g - *-s     . . .
          //              s-h - *-s     p-h - *-s     d-h - *-s     f-h - *-s .
          //              . . s-i - *-s     p-i - *-s     d-i - *-s     . . .
          for (Index l = 4; l < lmax_beta + 1; l++) {
            int norb = n_orbitals[l];
            int norb_1 = n_orbitals[l - 1];
            int norb_2 = n_orbitals[l - 2];
            Index ncart_1 = (l * (l + 1)) / 2;
            Index ncart_2 = ncart_1 - l;
            for (Index i = 0; i < n_orbitals[lmax_alpha_beta - l]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index j = 0; j < ncombined_cd; j++) {
                for (Index k = 0; k < ncart_2; k++) {
                  R(i, norb_1 + k, j) =
                      R(imx_i, norb_2 + k, j) + amb(0) * R(i, norb_2 + k, j);
                }
                for (Index k = 0; k < l; k++) {
                  Index k2 = norb_2 + ncart_2 + k;
                  R(i, norb_1 + ncart_2 + k, j) =
                      R(imx_i, k2, j) + amb(0) * R(i, k2, j);
                  R(i, norb_1 + ncart_1 + k, j) =
                      R(imy_i, k2, j) + amb(1) * R(i, k2, j);
                }
                R(i, norb - 1, j) =
                    R(imz_i, norb_1 - 1, j) + amb(2) * R(i, norb_1 - 1, j);
              }
            }
          }
          //------------------------------------------------------

          // prepare transformation matrices
          Index ntrafo_alpha = shell_alpha->getNumFunc();
          Index ntrafo_beta = shell_beta->getNumFunc();

          // get transformation matrices
          const Eigen::MatrixXd trafo_alpha =
              AOTransform::getTrafo(gaussian_alpha);
          const Eigen::MatrixXd trafo_beta =
              AOTransform::getTrafo(gaussian_beta);

          Eigen::Tensor<double, 4> R4_ab_sph(ntrafo_alpha, ntrafo_beta,
                                             ncombined_cd, ndelta);

          Index numfunc_cart_alpha = shell_alpha->getCartesianNumFunc();
          Index numfunc_cart_beta = shell_beta->getCartesianNumFunc();
          Index cartoffset_alpha = shell_alpha->getCartesianOffset();
          Index cartoffset_beta = shell_beta->getCartesianOffset();
          R4_ab_sph.setZero();

          for (Index j = 0; j < ncombined_cd; j++) {
            for (Index i_beta = 0; i_beta < ntrafo_beta; i_beta++) {
              for (Index i_alpha = 0; i_alpha < ntrafo_alpha; i_alpha++) {

                for (Index i_beta_c = 0; i_beta_c < numfunc_cart_beta;
                     i_beta_c++) {
                  for (Index i_alpha_c = 0; i_alpha_c < numfunc_cart_alpha;
                       i_alpha_c++) {

                    R4_ab_sph(i_alpha, i_beta, j, 0) +=
                        R(i_alpha_c + cartoffset_alpha,
                          i_beta_c + cartoffset_beta, j) *
                        trafo_alpha(i_alpha_c, i_alpha) *
                        trafo_beta(i_beta_c, i_beta);
                  }
                }
              }
            }
          }

          if (lmax_delta > 0) {
            // Integrals     *-* - s-p     *-* - p-p     *-* - d-p     *-* - f-p
            // *-* - g-p     *-* - h-p     *-* - i-p     *-* - j-p     . . .
            for (Index i = 0; i < n_orbitals[lmax_gamma_delta - 1]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index k = 0; k < ntrafo_beta; k++) {
                for (Index j = 0; j < ntrafo_alpha; j++) {
                  R4_ab_sph(j, k, i, Cart::x) = R4_ab_sph(j, k, imx_i, 0) +
                                                cmd(0) * R4_ab_sph(j, k, i, 0);
                  R4_ab_sph(j, k, i, Cart::y) = R4_ab_sph(j, k, imy_i, 0) +
                                                cmd(1) * R4_ab_sph(j, k, i, 0);
                  R4_ab_sph(j, k, i, Cart::z) = R4_ab_sph(j, k, imz_i, 0) +
                                                cmd(2) * R4_ab_sph(j, k, i, 0);
                }
              }
            }
          }

          if (lmax_delta > 1) {
            // Integrals     *-* - s-d     *-* - p-d     *-* - d-d     *-* - f-d
            // *-* - g-d     *-* - h-d     *-* - i-d     . . .
            for (Index i = 0; i < n_orbitals[lmax_gamma_delta - 2]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index k = 0; k < ntrafo_beta; k++) {
                for (Index j = 0; j < ntrafo_alpha; j++) {
                  R4_ab_sph(j, k, i, Cart::xx) =
                      R4_ab_sph(j, k, imx_i, Cart::x) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::x);
                  R4_ab_sph(j, k, i, Cart::xy) =
                      R4_ab_sph(j, k, imx_i, Cart::y) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::y);
                  R4_ab_sph(j, k, i, Cart::xz) =
                      R4_ab_sph(j, k, imx_i, Cart::z) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::z);
                  R4_ab_sph(j, k, i, Cart::yy) =
                      R4_ab_sph(j, k, imy_i, Cart::y) +
                      cmd(1) * R4_ab_sph(j, k, i, Cart::y);
                  R4_ab_sph(j, k, i, Cart::yz) =
                      R4_ab_sph(j, k, imy_i, Cart::z) +
                      cmd(1) * R4_ab_sph(j, k, i, Cart::z);
                  R4_ab_sph(j, k, i, Cart::zz) =
                      R4_ab_sph(j, k, imz_i, Cart::z) +
                      cmd(2) * R4_ab_sph(j, k, i, Cart::z);
                }
              }
            }
          }

          if (lmax_delta > 2) {
            // Integrals     *-* - s-f     *-* - p-f     *-* - d-f     *-* - f-f
            // *-* - g-f     *-* - h-f     . . .
            for (Index i = 0; i < n_orbitals[lmax_gamma_delta - 3]; i++) {
              int imx_i = i_more_x[i];
              int imy_i = i_more_y[i];
              int imz_i = i_more_z[i];
              for (Index k = 0; k < ntrafo_beta; k++) {
                for (Index j = 0; j < ntrafo_alpha; j++) {
                  R4_ab_sph(j, k, i, Cart::xxx) =
                      R4_ab_sph(j, k, imx_i, Cart::xx) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::xx);
                  R4_ab_sph(j, k, i, Cart::xxy) =
                      R4_ab_sph(j, k, imx_i, Cart::xy) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::xy);
                  R4_ab_sph(j, k, i, Cart::xxz) =
                      R4_ab_sph(j, k, imx_i, Cart::xz) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::xz);
                  R4_ab_sph(j, k, i, Cart::xyy) =
                      R4_ab_sph(j, k, imx_i, Cart::yy) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::yy);
                  R4_ab_sph(j, k, i, Cart::xyz) =
                      R4_ab_sph(j, k, imx_i, Cart::yz) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::yz);
                  R4_ab_sph(j, k, i, Cart::xzz) =
                      R4_ab_sph(j, k, imx_i, Cart::zz) +
                      cmd(0) * R4_ab_sph(j, k, i, Cart::zz);
                  R4_ab_sph(j, k, i, Cart::yyy) =
                      R4_ab_sph(j, k, imy_i, Cart::yy) +
                      cmd(1) * R4_ab_sph(j, k, i, Cart::yy);
                  R4_ab_sph(j, k, i, Cart::yyz) =
                      R4_ab_sph(j, k, imy_i, Cart::yz) +
                      cmd(1) * R4_ab_sph(j, k, i, Cart::yz);
                  R4_ab_sph(j, k, i, Cart::yzz) =
                      R4_ab_sph(j, k, imy_i, Cart::zz) +
                      cmd(1) * R4_ab_sph(j, k, i, Cart::zz);
                  R4_ab_sph(j, k, i, Cart::zzz) =
                      R4_ab_sph(j, k, imz_i, Cart::zz) +
                      cmd(2) * R4_ab_sph(j, k, i, Cart::zz);
                }
              }
            }
          }

          // Integrals     *-* - s-g     *-* - p-g     *-* - d-g     *-* - f-g
          // *-* - g-g     . . .
          //              *-* - s-h     *-* - p-h     *-* - d-h     *-* - f-h .
          //              . .
          //              *-* - s-i     *-* - p-i     *-* - d-i     . . .
          for (Index l = 4; l < lmax_delta + 1; l++) {
            int norb = n_orbitals[l];
            int norb_1 = n_orbitals[l - 1];
            int norb_2 = n_orbitals[l - 2];
            Index ncart_1 = (l * (l + 1)) / 2;
            Index ncart_2 = ncart_1 - l;
            for (Index i = 0; i < n_orbitals[lmax_gamma_delta - l]; i++) {
              Index imx_i = i_more_x[i];
              Index imy_i = i_more_y[i];
              Index imz_i = i_more_z[i];
              for (Index k = 0; k < ntrafo_beta; k++) {
                for (Index j = 0; j < ntrafo_alpha; j++) {
                  for (Index m = 0; m < ncart_2; m++) {
                    R4_ab_sph(j, k, i, norb_1 + m) =
                        R4_ab_sph(j, k, imx_i, norb_2 + m) +
                        cmd(0) * R4_ab_sph(j, k, i, norb_2 + m);
                  }
                  for (Index m = 0; m < l; m++) {
                    Index n = norb_2 + ncart_2 + m;
                    R4_ab_sph(j, k, i, norb_1 + ncart_2 + m) =
                        R4_ab_sph(j, k, imx_i, n) +
                        cmd(0) * R4_ab_sph(j, k, i, n);
                    R4_ab_sph(j, k, i, norb_1 + ncart_1 + m) =
                        R4_ab_sph(j, k, imy_i, n) +
                        cmd(1) * R4_ab_sph(j, k, i, n);
                  }
                  R4_ab_sph(j, k, i, norb - 1) =
                      R4_ab_sph(j, k, imz_i, norb_1 - 1) +
                      cmd(2) * R4_ab_sph(j, k, i, norb_1 - 1);
                }
              }
            }
          }

          // prepare transformation matrices
          Index ntrafo_gamma = shell_gamma->getNumFunc();
          Index ntrafo_delta = shell_delta->getNumFunc();

          const Eigen::MatrixXd trafo_gamma =
              AOTransform::getTrafo(gaussian_gamma);
          const Eigen::MatrixXd trafo_delta =
              AOTransform::getTrafo(gaussian_delta);

          Eigen::Tensor<double, 4> R4_sph(ntrafo_alpha, ntrafo_beta,
                                          ntrafo_gamma, ntrafo_delta);
          R4_sph.setZero();

          Index numfunc_cart_gamma = shell_gamma->getCartesianNumFunc();
          Index numfunc_cart_delta = shell_delta->getCartesianNumFunc();
          Index cartoffset_gamma = shell_gamma->getCartesianOffset();
          Index cartoffset_delta = shell_delta->getCartesianOffset();

          for (Index i_delta = 0; i_delta < ntrafo_delta; i_delta++) {
            for (Index i_gamma = 0; i_gamma < ntrafo_gamma; i_gamma++) {

              for (Index k = 0; k < ntrafo_beta; k++) {
                for (Index j = 0; j < ntrafo_alpha; j++) {

                  for (Index i_delta_t = 0; i_delta_t < numfunc_cart_delta;
                       i_delta_t++) {
                    for (Index i_gamma_t = 0; i_gamma_t < numfunc_cart_gamma;
                         i_gamma_t++) {

                      R4_sph(j, k, i_gamma, i_delta) +=
                          R4_ab_sph(j, k, i_gamma_t + cartoffset_gamma,
                                    i_delta_t + cartoffset_delta) *
                          trafo_gamma(i_gamma_t, i_gamma) *
                          trafo_delta(i_delta_t, i_delta);
                    }
                  }
                }
              }
            }
          }

          Index NumFunc_alpha = shell_alpha->getNumFunc();
          Index NumFunc_beta = shell_beta->getNumFunc();
          Index NumFunc_gamma = shell_gamma->getNumFunc();
          Index NumFunc_delta = shell_delta->getNumFunc();

          for (Index i_delta = 0; i_delta < NumFunc_delta; i_delta++) {
            for (Index i_gamma = 0; i_gamma < NumFunc_gamma; i_gamma++) {

              Index c = i_gamma;
              Index d = i_delta;
              if (gammadeltaswitch) {
                c = i_delta;
                d = i_gamma;
              }

              for (Index i_beta = 0; i_beta < NumFunc_beta; i_beta++) {
                for (Index i_alpha = 0; i_alpha < NumFunc_alpha; i_alpha++) {
                  Index a = i_alpha;
                  Index b = i_beta;
                  if (alphabetaswitch) {
                    a = i_beta;
                    b = i_alpha;
                  }
                  if (ab_cd_switch) {
                    block(c, d, a, b) +=
                        R4_sph(i_alpha, i_beta, i_gamma, i_delta);
                  } else {
                    block(a, b, c, d) +=
                        R4_sph(i_alpha, i_beta, i_gamma, i_delta);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return does_contribute;
}  // TCrawMatrix::FillFourCenterRepBlock

}  // namespace xtp
}  // namespace votca
