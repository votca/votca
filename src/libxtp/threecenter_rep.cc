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

#include <votca/xtp/aotransform.h>
#include <votca/xtp/threecenter.h>

using namespace std;

namespace votca {
namespace xtp {

/*
 * Calculate 3-center electron repulsion integrals
 *    R_{abc} = int{ phi_a(r)^DFT phi_b(r)^DFT phi_c(r')^AUX/(r-r') d3rd3r' }
 * for a given set of a b c as in
 *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
 * section II.B for cartesian Gaussians, then transforming
 * to spherical (angular momentum) Gaussians ("complete" shells
 * from S to Lmax, and finally cutting out those angular momentum
 * components actually present in shell-shell-shell combination.
 * Currently supported for
 *      S,P,D   functions in DFT basis and
 *      S,P,D,F functions in AUX  basis
 *
 */

bool TCMatrix::FillThreeCenterRepBlock(Eigen::Tensor<double, 3>& threec_block,
                                       const AOShell& shell_3,
                                       const AOShell& shell_1,
                                       const AOShell& shell_2) const {

  const double pi = boost::math::constants::pi<double>();
  const double gwaccuracy = 1.e-11;

  bool does_contribute = false;

  // shell info, only lmax tells how far to go

  int lmax_1 = shell_1.getLmax();
  int lmax_2 = shell_2.getLmax();
  int lmax_3 = shell_3.getLmax();

  int mmax = lmax_1 + lmax_2 + lmax_3;

  // set size of internal block for recursion

  const AOShell* shell_alpha = &shell_1;
  const AOShell* shell_beta = &shell_2;
  const AOShell* shell_gamma = &shell_3;
  bool alphabetaswitch = false;

  // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we
  // calculate (ps,s), due to symmetry they are the same.

  if (lmax_1 < lmax_2) {
    shell_alpha = &shell_2;
    shell_beta = &shell_1;
    alphabetaswitch = true;
  }

  const Eigen::Vector3d& pos_alpha = shell_alpha->getPos();
  const Eigen::Vector3d& pos_beta = shell_beta->getPos();
  const Eigen::Vector3d& pos_gamma = shell_gamma->getPos();

  int lmax_alpha = shell_alpha->getLmax();
  int lmax_beta = shell_beta->getLmax();
  int lmax_gamma = shell_gamma->getLmax();

  int ngamma = AOTransform::getBlockSize(lmax_gamma);
  int nbeta = AOTransform::getBlockSize(lmax_beta);
  int ncombined = AOTransform::getBlockSize(lmax_alpha + lmax_beta);

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

  const Eigen::Vector3d amb = pos_alpha - pos_beta;
  double dist3 = amb.squaredNorm();

  for (const auto& gaussian_alpha : *shell_alpha) {
    const double decay_alpha = gaussian_alpha.getDecay();

    for (const auto& gaussian_beta : *shell_beta) {
      const double decay_beta = gaussian_beta.getDecay();
      double rzeta = 0.5 / (decay_alpha + decay_beta);
      const Eigen::Vector3d P =
          2.0 * (decay_alpha * pos_alpha + decay_beta * pos_beta) * rzeta;
      const Eigen::Vector3d pma = P - pos_alpha;
      double xi = 2.0 * decay_alpha * decay_beta * rzeta;
      double fact_alpha_beta = 16.0 * xi *
                               pow(pi / (decay_alpha * decay_beta), 0.25) *
                               exp(-xi * dist3);

      for (const auto& gaussian_gamma : *shell_gamma) {
        const double decay_gamma = gaussian_gamma.getDecay();

        double decay = decay_alpha + decay_beta + decay_gamma;
        double rgamma = 0.5 / decay_gamma;
        double rdecay = 0.5 / decay;

        double sss = fact_alpha_beta * pow(rdecay * rdecay * rgamma, 0.25);

        if (sss < gwaccuracy) {
          continue;
        }

        does_contribute = true;

        double gfak = decay_gamma / decay;
        double cfak = (decay_alpha + decay_beta) / decay;
        const Eigen::Vector3d W =
            (decay_alpha * pos_alpha + decay_beta * pos_beta +
             decay_gamma * pos_gamma) /
            decay;
        double U = (decay_alpha + decay_beta) * decay_gamma / decay *
                   (P - pos_gamma).squaredNorm();

        const Eigen::Vector3d wmp = W - P;
        const Eigen::Vector3d wmc = W - pos_gamma;

        Eigen::Tensor<double, 3> R_temp(ncombined, ngamma,
                                        std::max(2, mmax + 1));
        R_temp.setZero();

        const Eigen::VectorXd FmT = AOTransform::XIntegrate(mmax + 1, U);

        // ss integrals

        for (int i = 0; i < mmax + 1; i++) {
          R_temp(0, 0, i) = sss * FmT[i];
        }

        int lmax_alpha_beta = lmax_alpha + lmax_beta;

        // Integral  p - s - s
        if (lmax_alpha_beta > 0) {
          for (int m = 0; m < mmax; m++) {
            R_temp(Cart::x, 0, m) =
                pma(0) * R_temp(0, 0, m) + wmp(0) * R_temp(0, 0, m + 1);
            R_temp(Cart::y, 0, m) =
                pma(1) * R_temp(0, 0, m) + wmp(1) * R_temp(0, 0, m + 1);
            R_temp(Cart::z, 0, m) =
                pma(2) * R_temp(0, 0, m) + wmp(2) * R_temp(0, 0, m + 1);
          }
        }
        //------------------------------------------------------

        // Integral  d - s - s
        if (lmax_alpha_beta > 1) {
          for (int m = 0; m < mmax - 1; m++) {
            double term =
                rzeta * (R_temp(0, 0, m) - gfak * R_temp(0, 0, m + 1));
            R_temp(Cart::xx, 0, m) = pma(0) * R_temp(Cart::x, 0, m) +
                                     wmp(0) * R_temp(Cart::x, 0, m + 1) + term;
            R_temp(Cart::xy, 0, m) = pma(0) * R_temp(Cart::y, 0, m) +
                                     wmp(0) * R_temp(Cart::y, 0, m + 1);
            R_temp(Cart::xz, 0, m) = pma(0) * R_temp(Cart::z, 0, m) +
                                     wmp(0) * R_temp(Cart::z, 0, m + 1);
            R_temp(Cart::yy, 0, m) = pma(1) * R_temp(Cart::y, 0, m) +
                                     wmp(1) * R_temp(Cart::y, 0, m + 1) + term;
            R_temp(Cart::yz, 0, m) = pma(1) * R_temp(Cart::z, 0, m) +
                                     wmp(1) * R_temp(Cart::z, 0, m + 1);
            R_temp(Cart::zz, 0, m) = pma(2) * R_temp(Cart::z, 0, m) +
                                     wmp(2) * R_temp(Cart::z, 0, m + 1) + term;
          }
        }
        //------------------------------------------------------

        // Integral  f - s - s
        if (lmax_alpha_beta > 2) {
          for (int m = 0; m < mmax - 2; m++) {
            R_temp(Cart::xxx, 0, m) =
                pma(0) * R_temp(Cart::xx, 0, m) +
                wmp(0) * R_temp(Cart::xx, 0, m + 1) +
                2 * rzeta *
                    (R_temp(Cart::x, 0, m) - gfak * R_temp(Cart::x, 0, m + 1));
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
            R_temp(Cart::yyy, 0, m) =
                pma(1) * R_temp(Cart::yy, 0, m) +
                wmp(1) * R_temp(Cart::yy, 0, m + 1) +
                2 * rzeta *
                    (R_temp(Cart::y, 0, m) - gfak * R_temp(Cart::y, 0, m + 1));
            R_temp(Cart::yyz, 0, m) = pma(2) * R_temp(Cart::yy, 0, m) +
                                      wmp(2) * R_temp(Cart::yy, 0, m + 1);
            R_temp(Cart::yzz, 0, m) = pma(1) * R_temp(Cart::zz, 0, m) +
                                      wmp(1) * R_temp(Cart::zz, 0, m + 1);
            R_temp(Cart::zzz, 0, m) =
                pma(2) * R_temp(Cart::zz, 0, m) +
                wmp(2) * R_temp(Cart::zz, 0, m + 1) +
                2 * rzeta *
                    (R_temp(Cart::z, 0, m) - gfak * R_temp(Cart::z, 0, m + 1));
          }
        }
        //------------------------------------------------------

        // Integral  g - s - s
        if (lmax_alpha_beta > 3) {
          for (int m = 0; m < mmax - 3; m++) {
            double term_xx = rzeta * (R_temp(Cart::xx, 0, m) -
                                      gfak * R_temp(Cart::xx, 0, m + 1));
            double term_yy = rzeta * (R_temp(Cart::yy, 0, m) -
                                      gfak * R_temp(Cart::yy, 0, m + 1));
            double term_zz = rzeta * (R_temp(Cart::zz, 0, m) -
                                      gfak * R_temp(Cart::zz, 0, m + 1));
            R_temp(Cart::xxxx, 0, m) = pma(0) * R_temp(Cart::xxx, 0, m) +
                                       wmp(0) * R_temp(Cart::xxx, 0, m + 1) +
                                       3 * term_xx;
            R_temp(Cart::xxxy, 0, m) = pma(1) * R_temp(Cart::xxx, 0, m) +
                                       wmp(1) * R_temp(Cart::xxx, 0, m + 1);
            R_temp(Cart::xxxz, 0, m) = pma(2) * R_temp(Cart::xxx, 0, m) +
                                       wmp(2) * R_temp(Cart::xxx, 0, m + 1);
            R_temp(Cart::xxyy, 0, m) = pma(0) * R_temp(Cart::xyy, 0, m) +
                                       wmp(0) * R_temp(Cart::xyy, 0, m + 1) +
                                       term_yy;
            R_temp(Cart::xxyz, 0, m) = pma(1) * R_temp(Cart::xxz, 0, m) +
                                       wmp(1) * R_temp(Cart::xxz, 0, m + 1);
            R_temp(Cart::xxzz, 0, m) = pma(0) * R_temp(Cart::xzz, 0, m) +
                                       wmp(0) * R_temp(Cart::xzz, 0, m + 1) +
                                       term_zz;
            R_temp(Cart::xyyy, 0, m) = pma(0) * R_temp(Cart::yyy, 0, m) +
                                       wmp(0) * R_temp(Cart::yyy, 0, m + 1);
            R_temp(Cart::xyyz, 0, m) = pma(0) * R_temp(Cart::yyz, 0, m) +
                                       wmp(0) * R_temp(Cart::yyz, 0, m + 1);
            R_temp(Cart::xyzz, 0, m) = pma(0) * R_temp(Cart::yzz, 0, m) +
                                       wmp(0) * R_temp(Cart::yzz, 0, m + 1);
            R_temp(Cart::xzzz, 0, m) = pma(0) * R_temp(Cart::zzz, 0, m) +
                                       wmp(0) * R_temp(Cart::zzz, 0, m + 1);
            R_temp(Cart::yyyy, 0, m) = pma(1) * R_temp(Cart::yyy, 0, m) +
                                       wmp(1) * R_temp(Cart::yyy, 0, m + 1) +
                                       3 * term_yy;
            R_temp(Cart::yyyz, 0, m) = pma(2) * R_temp(Cart::yyy, 0, m) +
                                       wmp(2) * R_temp(Cart::yyy, 0, m + 1);
            R_temp(Cart::yyzz, 0, m) = pma(1) * R_temp(Cart::yzz, 0, m) +
                                       wmp(1) * R_temp(Cart::yzz, 0, m + 1) +
                                       term_zz;
            R_temp(Cart::yzzz, 0, m) = pma(1) * R_temp(Cart::zzz, 0, m) +
                                       wmp(1) * R_temp(Cart::zzz, 0, m + 1);
            R_temp(Cart::zzzz, 0, m) = pma(2) * R_temp(Cart::zzz, 0, m) +
                                       wmp(2) * R_temp(Cart::zzz, 0, m + 1) +
                                       3 * term_zz;
          }
        }
        //------------------------------------------------------

        // Integral  h - s - s
        if (lmax_alpha_beta > 4) {
          for (int m = 0; m < mmax - 4; m++) {
            double term_xxx = rzeta * (R_temp(Cart::xxx, 0, m) -
                                       gfak * R_temp(Cart::xxx, 0, m + 1));
            double term_yyy = rzeta * (R_temp(Cart::yyy, 0, m) -
                                       gfak * R_temp(Cart::yyy, 0, m + 1));
            double term_zzz = rzeta * (R_temp(Cart::zzz, 0, m) -
                                       gfak * R_temp(Cart::zzz, 0, m + 1));
            R_temp(Cart::xxxxx, 0, m) = pma(0) * R_temp(Cart::xxxx, 0, m) +
                                        wmp(0) * R_temp(Cart::xxxx, 0, m + 1) +
                                        4 * term_xxx;
            R_temp(Cart::xxxxy, 0, m) = pma(1) * R_temp(Cart::xxxx, 0, m) +
                                        wmp(1) * R_temp(Cart::xxxx, 0, m + 1);
            R_temp(Cart::xxxxz, 0, m) = pma(2) * R_temp(Cart::xxxx, 0, m) +
                                        wmp(2) * R_temp(Cart::xxxx, 0, m + 1);
            R_temp(Cart::xxxyy, 0, m) = pma(1) * R_temp(Cart::xxxy, 0, m) +
                                        wmp(1) * R_temp(Cart::xxxy, 0, m + 1) +
                                        term_xxx;
            R_temp(Cart::xxxyz, 0, m) = pma(1) * R_temp(Cart::xxxz, 0, m) +
                                        wmp(1) * R_temp(Cart::xxxz, 0, m + 1);
            R_temp(Cart::xxxzz, 0, m) = pma(2) * R_temp(Cart::xxxz, 0, m) +
                                        wmp(2) * R_temp(Cart::xxxz, 0, m + 1) +
                                        term_xxx;
            R_temp(Cart::xxyyy, 0, m) = pma(0) * R_temp(Cart::xyyy, 0, m) +
                                        wmp(0) * R_temp(Cart::xyyy, 0, m + 1) +
                                        term_yyy;
            R_temp(Cart::xxyyz, 0, m) = pma(2) * R_temp(Cart::xxyy, 0, m) +
                                        wmp(2) * R_temp(Cart::xxyy, 0, m + 1);
            R_temp(Cart::xxyzz, 0, m) = pma(1) * R_temp(Cart::xxzz, 0, m) +
                                        wmp(1) * R_temp(Cart::xxzz, 0, m + 1);
            R_temp(Cart::xxzzz, 0, m) = pma(0) * R_temp(Cart::xzzz, 0, m) +
                                        wmp(0) * R_temp(Cart::xzzz, 0, m + 1) +
                                        term_zzz;
            R_temp(Cart::xyyyy, 0, m) = pma(0) * R_temp(Cart::yyyy, 0, m) +
                                        wmp(0) * R_temp(Cart::yyyy, 0, m + 1);
            R_temp(Cart::xyyyz, 0, m) = pma(0) * R_temp(Cart::yyyz, 0, m) +
                                        wmp(0) * R_temp(Cart::yyyz, 0, m + 1);
            R_temp(Cart::xyyzz, 0, m) = pma(0) * R_temp(Cart::yyzz, 0, m) +
                                        wmp(0) * R_temp(Cart::yyzz, 0, m + 1);
            R_temp(Cart::xyzzz, 0, m) = pma(0) * R_temp(Cart::yzzz, 0, m) +
                                        wmp(0) * R_temp(Cart::yzzz, 0, m + 1);
            R_temp(Cart::xzzzz, 0, m) = pma(0) * R_temp(Cart::zzzz, 0, m) +
                                        wmp(0) * R_temp(Cart::zzzz, 0, m + 1);
            R_temp(Cart::yyyyy, 0, m) = pma(1) * R_temp(Cart::yyyy, 0, m) +
                                        wmp(1) * R_temp(Cart::yyyy, 0, m + 1) +
                                        4 * term_yyy;
            R_temp(Cart::yyyyz, 0, m) = pma(2) * R_temp(Cart::yyyy, 0, m) +
                                        wmp(2) * R_temp(Cart::yyyy, 0, m + 1);
            R_temp(Cart::yyyzz, 0, m) = pma(2) * R_temp(Cart::yyyz, 0, m) +
                                        wmp(2) * R_temp(Cart::yyyz, 0, m + 1) +
                                        term_yyy;
            R_temp(Cart::yyzzz, 0, m) = pma(1) * R_temp(Cart::yzzz, 0, m) +
                                        wmp(1) * R_temp(Cart::yzzz, 0, m + 1) +
                                        term_zzz;
            R_temp(Cart::yzzzz, 0, m) = pma(1) * R_temp(Cart::zzzz, 0, m) +
                                        wmp(1) * R_temp(Cart::zzzz, 0, m + 1);
            R_temp(Cart::zzzzz, 0, m) = pma(2) * R_temp(Cart::zzzz, 0, m) +
                                        wmp(2) * R_temp(Cart::zzzz, 0, m + 1) +
                                        4 * term_zzz;
          }
        }
        //------------------------------------------------------

        // Integral  i - s - s
        if (lmax_alpha_beta > 5) {
          for (int m = 0; m < mmax - 5; m++) {
            double term_xxxx = rzeta * (R_temp(Cart::xxxx, 0, m) -
                                        gfak * R_temp(Cart::xxxx, 0, m + 1));
            double term_xyyy = rzeta * (R_temp(Cart::xyyy, 0, m) -
                                        gfak * R_temp(Cart::xyyy, 0, m + 1));
            double term_xzzz = rzeta * (R_temp(Cart::xzzz, 0, m) -
                                        gfak * R_temp(Cart::xzzz, 0, m + 1));
            double term_yyyy = rzeta * (R_temp(Cart::yyyy, 0, m) -
                                        gfak * R_temp(Cart::yyyy, 0, m + 1));
            double term_yyzz = rzeta * (R_temp(Cart::yyzz, 0, m) -
                                        gfak * R_temp(Cart::yyzz, 0, m + 1));
            double term_yzzz = rzeta * (R_temp(Cart::yzzz, 0, m) -
                                        gfak * R_temp(Cart::yzzz, 0, m + 1));
            double term_zzzz = rzeta * (R_temp(Cart::zzzz, 0, m) -
                                        gfak * R_temp(Cart::zzzz, 0, m + 1));
            R_temp(Cart::xxxxxx, 0, m) =
                pma(0) * R_temp(Cart::xxxxx, 0, m) +
                wmp(0) * R_temp(Cart::xxxxx, 0, m + 1) + 5 * term_xxxx;
            R_temp(Cart::xxxxxy, 0, m) = pma(1) * R_temp(Cart::xxxxx, 0, m) +
                                         wmp(1) * R_temp(Cart::xxxxx, 0, m + 1);
            R_temp(Cart::xxxxxz, 0, m) = pma(2) * R_temp(Cart::xxxxx, 0, m) +
                                         wmp(2) * R_temp(Cart::xxxxx, 0, m + 1);
            R_temp(Cart::xxxxyy, 0, m) =
                pma(1) * R_temp(Cart::xxxxy, 0, m) +
                wmp(1) * R_temp(Cart::xxxxy, 0, m + 1) + term_xxxx;
            R_temp(Cart::xxxxyz, 0, m) = pma(1) * R_temp(Cart::xxxxz, 0, m) +
                                         wmp(1) * R_temp(Cart::xxxxz, 0, m + 1);
            R_temp(Cart::xxxxzz, 0, m) =
                pma(2) * R_temp(Cart::xxxxz, 0, m) +
                wmp(2) * R_temp(Cart::xxxxz, 0, m + 1) + term_xxxx;
            R_temp(Cart::xxxyyy, 0, m) =
                pma(0) * R_temp(Cart::xxyyy, 0, m) +
                wmp(0) * R_temp(Cart::xxyyy, 0, m + 1) + 2 * term_xyyy;
            R_temp(Cart::xxxyyz, 0, m) = pma(2) * R_temp(Cart::xxxyy, 0, m) +
                                         wmp(2) * R_temp(Cart::xxxyy, 0, m + 1);
            R_temp(Cart::xxxyzz, 0, m) = pma(1) * R_temp(Cart::xxxzz, 0, m) +
                                         wmp(1) * R_temp(Cart::xxxzz, 0, m + 1);
            R_temp(Cart::xxxzzz, 0, m) =
                pma(0) * R_temp(Cart::xxzzz, 0, m) +
                wmp(0) * R_temp(Cart::xxzzz, 0, m + 1) + 2 * term_xzzz;
            R_temp(Cart::xxyyyy, 0, m) =
                pma(0) * R_temp(Cart::xyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xyyyy, 0, m + 1) + term_yyyy;
            R_temp(Cart::xxyyyz, 0, m) = pma(2) * R_temp(Cart::xxyyy, 0, m) +
                                         wmp(2) * R_temp(Cart::xxyyy, 0, m + 1);
            R_temp(Cart::xxyyzz, 0, m) =
                pma(0) * R_temp(Cart::xyyzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyzz, 0, m + 1) + term_yyzz;
            R_temp(Cart::xxyzzz, 0, m) = pma(1) * R_temp(Cart::xxzzz, 0, m) +
                                         wmp(1) * R_temp(Cart::xxzzz, 0, m + 1);
            R_temp(Cart::xxzzzz, 0, m) =
                pma(0) * R_temp(Cart::xzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xzzzz, 0, m + 1) + term_zzzz;
            R_temp(Cart::xyyyyy, 0, m) = pma(0) * R_temp(Cart::yyyyy, 0, m) +
                                         wmp(0) * R_temp(Cart::yyyyy, 0, m + 1);
            R_temp(Cart::xyyyyz, 0, m) = pma(0) * R_temp(Cart::yyyyz, 0, m) +
                                         wmp(0) * R_temp(Cart::yyyyz, 0, m + 1);
            R_temp(Cart::xyyyzz, 0, m) = pma(0) * R_temp(Cart::yyyzz, 0, m) +
                                         wmp(0) * R_temp(Cart::yyyzz, 0, m + 1);
            R_temp(Cart::xyyzzz, 0, m) = pma(0) * R_temp(Cart::yyzzz, 0, m) +
                                         wmp(0) * R_temp(Cart::yyzzz, 0, m + 1);
            R_temp(Cart::xyzzzz, 0, m) = pma(0) * R_temp(Cart::yzzzz, 0, m) +
                                         wmp(0) * R_temp(Cart::yzzzz, 0, m + 1);
            R_temp(Cart::xzzzzz, 0, m) = pma(0) * R_temp(Cart::zzzzz, 0, m) +
                                         wmp(0) * R_temp(Cart::zzzzz, 0, m + 1);
            R_temp(Cart::yyyyyy, 0, m) =
                pma(1) * R_temp(Cart::yyyyy, 0, m) +
                wmp(1) * R_temp(Cart::yyyyy, 0, m + 1) + 5 * term_yyyy;
            R_temp(Cart::yyyyyz, 0, m) = pma(2) * R_temp(Cart::yyyyy, 0, m) +
                                         wmp(2) * R_temp(Cart::yyyyy, 0, m + 1);
            R_temp(Cart::yyyyzz, 0, m) =
                pma(2) * R_temp(Cart::yyyyz, 0, m) +
                wmp(2) * R_temp(Cart::yyyyz, 0, m + 1) + term_yyyy;
            R_temp(Cart::yyyzzz, 0, m) =
                pma(1) * R_temp(Cart::yyzzz, 0, m) +
                wmp(1) * R_temp(Cart::yyzzz, 0, m + 1) + 2 * term_yzzz;
            R_temp(Cart::yyzzzz, 0, m) =
                pma(1) * R_temp(Cart::yzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yzzzz, 0, m + 1) + term_zzzz;
            R_temp(Cart::yzzzzz, 0, m) = pma(1) * R_temp(Cart::zzzzz, 0, m) +
                                         wmp(1) * R_temp(Cart::zzzzz, 0, m + 1);
            R_temp(Cart::zzzzzz, 0, m) =
                pma(2) * R_temp(Cart::zzzzz, 0, m) +
                wmp(2) * R_temp(Cart::zzzzz, 0, m + 1) + 5 * term_zzzz;
          }
        }
        //------------------------------------------------------

        // Integral  j - s - s
        if (lmax_alpha_beta > 6) {
          for (int m = 0; m < mmax - 6; m++) {
            double term_xxxxx = rzeta * (R_temp(Cart::xxxxx, 0, m) -
                                         gfak * R_temp(Cart::xxxxx, 0, m + 1));
            double term_xxxxy = rzeta * (R_temp(Cart::xxxxy, 0, m) -
                                         gfak * R_temp(Cart::xxxxy, 0, m + 1));
            double term_xxxxz = rzeta * (R_temp(Cart::xxxxz, 0, m) -
                                         gfak * R_temp(Cart::xxxxz, 0, m + 1));
            double term_xxxzz = rzeta * (R_temp(Cart::xxxzz, 0, m) -
                                         gfak * R_temp(Cart::xxxzz, 0, m + 1));
            double term_xyyyy = rzeta * (R_temp(Cart::xyyyy, 0, m) -
                                         gfak * R_temp(Cart::xyyyy, 0, m + 1));
            double term_xzzzz = rzeta * (R_temp(Cart::xzzzz, 0, m) -
                                         gfak * R_temp(Cart::xzzzz, 0, m + 1));
            double term_yyyyy = rzeta * (R_temp(Cart::yyyyy, 0, m) -
                                         gfak * R_temp(Cart::yyyyy, 0, m + 1));
            double term_yyyyz = rzeta * (R_temp(Cart::yyyyz, 0, m) -
                                         gfak * R_temp(Cart::yyyyz, 0, m + 1));
            double term_yyyzz = rzeta * (R_temp(Cart::yyyzz, 0, m) -
                                         gfak * R_temp(Cart::yyyzz, 0, m + 1));
            double term_yyzzz = rzeta * (R_temp(Cart::yyzzz, 0, m) -
                                         gfak * R_temp(Cart::yyzzz, 0, m + 1));
            double term_yzzzz = rzeta * (R_temp(Cart::yzzzz, 0, m) -
                                         gfak * R_temp(Cart::yzzzz, 0, m + 1));
            double term_zzzzz = rzeta * (R_temp(Cart::zzzzz, 0, m) -
                                         gfak * R_temp(Cart::zzzzz, 0, m + 1));
            R_temp(Cart::xxxxxxx, 0, m) =
                pma(0) * R_temp(Cart::xxxxxx, 0, m) +
                wmp(0) * R_temp(Cart::xxxxxx, 0, m + 1) + 6 * term_xxxxx;
            R_temp(Cart::xxxxxxy, 0, m) =
                pma(1) * R_temp(Cart::xxxxxx, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxx, 0, m + 1);
            R_temp(Cart::xxxxxxz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxx, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxx, 0, m + 1);
            R_temp(Cart::xxxxxyy, 0, m) =
                pma(1) * R_temp(Cart::xxxxxy, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxy, 0, m + 1) + term_xxxxx;
            R_temp(Cart::xxxxxyz, 0, m) =
                pma(1) * R_temp(Cart::xxxxxz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxz, 0, m + 1);
            R_temp(Cart::xxxxxzz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxz, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxz, 0, m + 1) + term_xxxxx;
            R_temp(Cart::xxxxyyy, 0, m) =
                pma(1) * R_temp(Cart::xxxxyy, 0, m) +
                wmp(1) * R_temp(Cart::xxxxyy, 0, m + 1) + 2 * term_xxxxy;
            R_temp(Cart::xxxxyyz, 0, m) =
                pma(2) * R_temp(Cart::xxxxyy, 0, m) +
                wmp(2) * R_temp(Cart::xxxxyy, 0, m + 1);
            R_temp(Cart::xxxxyzz, 0, m) =
                pma(1) * R_temp(Cart::xxxxzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxzz, 0, m + 1);
            R_temp(Cart::xxxxzzz, 0, m) =
                pma(2) * R_temp(Cart::xxxxzz, 0, m) +
                wmp(2) * R_temp(Cart::xxxxzz, 0, m + 1) + 2 * term_xxxxz;
            R_temp(Cart::xxxyyyy, 0, m) =
                pma(0) * R_temp(Cart::xxyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xxyyyy, 0, m + 1) + 2 * term_xyyyy;
            R_temp(Cart::xxxyyyz, 0, m) =
                pma(2) * R_temp(Cart::xxxyyy, 0, m) +
                wmp(2) * R_temp(Cart::xxxyyy, 0, m + 1);
            R_temp(Cart::xxxyyzz, 0, m) =
                pma(1) * R_temp(Cart::xxxyzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxyzz, 0, m + 1) + term_xxxzz;
            R_temp(Cart::xxxyzzz, 0, m) =
                pma(1) * R_temp(Cart::xxxzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxzzz, 0, m + 1);
            R_temp(Cart::xxxzzzz, 0, m) =
                pma(0) * R_temp(Cart::xxzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xxzzzz, 0, m + 1) + 2 * term_xzzzz;
            R_temp(Cart::xxyyyyy, 0, m) =
                pma(0) * R_temp(Cart::xyyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xyyyyy, 0, m + 1) + term_yyyyy;
            R_temp(Cart::xxyyyyz, 0, m) =
                pma(2) * R_temp(Cart::xxyyyy, 0, m) +
                wmp(2) * R_temp(Cart::xxyyyy, 0, m + 1);
            R_temp(Cart::xxyyyzz, 0, m) =
                pma(0) * R_temp(Cart::xyyyzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyyzz, 0, m + 1) + term_yyyzz;
            R_temp(Cart::xxyyzzz, 0, m) =
                pma(0) * R_temp(Cart::xyyzzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyzzz, 0, m + 1) + term_yyzzz;
            R_temp(Cart::xxyzzzz, 0, m) =
                pma(1) * R_temp(Cart::xxzzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxzzzz, 0, m + 1);
            R_temp(Cart::xxzzzzz, 0, m) =
                pma(0) * R_temp(Cart::xzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xzzzzz, 0, m + 1) + term_zzzzz;
            R_temp(Cart::xyyyyyy, 0, m) =
                pma(0) * R_temp(Cart::yyyyyy, 0, m) +
                wmp(0) * R_temp(Cart::yyyyyy, 0, m + 1);
            R_temp(Cart::xyyyyyz, 0, m) =
                pma(0) * R_temp(Cart::yyyyyz, 0, m) +
                wmp(0) * R_temp(Cart::yyyyyz, 0, m + 1);
            R_temp(Cart::xyyyyzz, 0, m) =
                pma(0) * R_temp(Cart::yyyyzz, 0, m) +
                wmp(0) * R_temp(Cart::yyyyzz, 0, m + 1);
            R_temp(Cart::xyyyzzz, 0, m) =
                pma(0) * R_temp(Cart::yyyzzz, 0, m) +
                wmp(0) * R_temp(Cart::yyyzzz, 0, m + 1);
            R_temp(Cart::xyyzzzz, 0, m) =
                pma(0) * R_temp(Cart::yyzzzz, 0, m) +
                wmp(0) * R_temp(Cart::yyzzzz, 0, m + 1);
            R_temp(Cart::xyzzzzz, 0, m) =
                pma(0) * R_temp(Cart::yzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::yzzzzz, 0, m + 1);
            R_temp(Cart::xzzzzzz, 0, m) =
                pma(0) * R_temp(Cart::zzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::zzzzzz, 0, m + 1);
            R_temp(Cart::yyyyyyy, 0, m) =
                pma(1) * R_temp(Cart::yyyyyy, 0, m) +
                wmp(1) * R_temp(Cart::yyyyyy, 0, m + 1) + 6 * term_yyyyy;
            R_temp(Cart::yyyyyyz, 0, m) =
                pma(2) * R_temp(Cart::yyyyyy, 0, m) +
                wmp(2) * R_temp(Cart::yyyyyy, 0, m + 1);
            R_temp(Cart::yyyyyzz, 0, m) =
                pma(2) * R_temp(Cart::yyyyyz, 0, m) +
                wmp(2) * R_temp(Cart::yyyyyz, 0, m + 1) + term_yyyyy;
            R_temp(Cart::yyyyzzz, 0, m) =
                pma(2) * R_temp(Cart::yyyyzz, 0, m) +
                wmp(2) * R_temp(Cart::yyyyzz, 0, m + 1) + 2 * term_yyyyz;
            R_temp(Cart::yyyzzzz, 0, m) =
                pma(1) * R_temp(Cart::yyzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yyzzzz, 0, m + 1) + 2 * term_yzzzz;
            R_temp(Cart::yyzzzzz, 0, m) =
                pma(1) * R_temp(Cart::yzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yzzzzz, 0, m + 1) + term_zzzzz;
            R_temp(Cart::yzzzzzz, 0, m) =
                pma(1) * R_temp(Cart::zzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::zzzzzz, 0, m + 1);
            R_temp(Cart::zzzzzzz, 0, m) =
                pma(2) * R_temp(Cart::zzzzzz, 0, m) +
                wmp(2) * R_temp(Cart::zzzzzz, 0, m + 1) + 6 * term_zzzzz;
          }
        }
        //------------------------------------------------------

        // Integral  k - s - s
        if (lmax_alpha_beta > 7) {
          for (int m = 0; m < mmax - 7; m++) {
            double term_xxxxxx =
                rzeta * (R_temp(Cart::xxxxxx, 0, m) -
                         gfak * R_temp(Cart::xxxxxx, 0, m + 1));
            double term_xxxxxy =
                rzeta * (R_temp(Cart::xxxxxy, 0, m) -
                         gfak * R_temp(Cart::xxxxxy, 0, m + 1));
            double term_xxxxxz =
                rzeta * (R_temp(Cart::xxxxxz, 0, m) -
                         gfak * R_temp(Cart::xxxxxz, 0, m + 1));
            double term_xxxxzz =
                rzeta * (R_temp(Cart::xxxxzz, 0, m) -
                         gfak * R_temp(Cart::xxxxzz, 0, m + 1));
            double term_xxxyyy =
                rzeta * (R_temp(Cart::xxxyyy, 0, m) -
                         gfak * R_temp(Cart::xxxyyy, 0, m + 1));
            double term_xxxzzz =
                rzeta * (R_temp(Cart::xxxzzz, 0, m) -
                         gfak * R_temp(Cart::xxxzzz, 0, m + 1));
            double term_xxyyyy =
                rzeta * (R_temp(Cart::xxyyyy, 0, m) -
                         gfak * R_temp(Cart::xxyyyy, 0, m + 1));
            double term_xxzzzz =
                rzeta * (R_temp(Cart::xxzzzz, 0, m) -
                         gfak * R_temp(Cart::xxzzzz, 0, m + 1));
            double term_xyyyyy =
                rzeta * (R_temp(Cart::xyyyyy, 0, m) -
                         gfak * R_temp(Cart::xyyyyy, 0, m + 1));
            double term_xzzzzz =
                rzeta * (R_temp(Cart::xzzzzz, 0, m) -
                         gfak * R_temp(Cart::xzzzzz, 0, m + 1));
            double term_yyyyyy =
                rzeta * (R_temp(Cart::yyyyyy, 0, m) -
                         gfak * R_temp(Cart::yyyyyy, 0, m + 1));
            double term_yyyyyz =
                rzeta * (R_temp(Cart::yyyyyz, 0, m) -
                         gfak * R_temp(Cart::yyyyyz, 0, m + 1));
            double term_yyyyzz =
                rzeta * (R_temp(Cart::yyyyzz, 0, m) -
                         gfak * R_temp(Cart::yyyyzz, 0, m + 1));
            double term_yyyzzz =
                rzeta * (R_temp(Cart::yyyzzz, 0, m) -
                         gfak * R_temp(Cart::yyyzzz, 0, m + 1));
            double term_yyzzzz =
                rzeta * (R_temp(Cart::yyzzzz, 0, m) -
                         gfak * R_temp(Cart::yyzzzz, 0, m + 1));
            double term_yzzzzz =
                rzeta * (R_temp(Cart::yzzzzz, 0, m) -
                         gfak * R_temp(Cart::yzzzzz, 0, m + 1));
            double term_zzzzzz =
                rzeta * (R_temp(Cart::zzzzzz, 0, m) -
                         gfak * R_temp(Cart::zzzzzz, 0, m + 1));
            R_temp(Cart::xxxxxxxx, 0, m) =
                pma(0) * R_temp(Cart::xxxxxxx, 0, m) +
                wmp(0) * R_temp(Cart::xxxxxxx, 0, m + 1) + 7 * term_xxxxxx;
            R_temp(Cart::xxxxxxxy, 0, m) =
                pma(1) * R_temp(Cart::xxxxxxx, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxxx, 0, m + 1);
            R_temp(Cart::xxxxxxxz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxxx, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxxx, 0, m + 1);
            R_temp(Cart::xxxxxxyy, 0, m) =
                pma(1) * R_temp(Cart::xxxxxxy, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxxy, 0, m + 1) + term_xxxxxx;
            R_temp(Cart::xxxxxxyz, 0, m) =
                pma(1) * R_temp(Cart::xxxxxxz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxxz, 0, m + 1);
            R_temp(Cart::xxxxxxzz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxxz, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxxz, 0, m + 1) + term_xxxxxx;
            R_temp(Cart::xxxxxyyy, 0, m) =
                pma(1) * R_temp(Cart::xxxxxyy, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxyy, 0, m + 1) + 2 * term_xxxxxy;
            R_temp(Cart::xxxxxyyz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxyy, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxyy, 0, m + 1);
            R_temp(Cart::xxxxxyzz, 0, m) =
                pma(1) * R_temp(Cart::xxxxxzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxxzz, 0, m + 1);
            R_temp(Cart::xxxxxzzz, 0, m) =
                pma(2) * R_temp(Cart::xxxxxzz, 0, m) +
                wmp(2) * R_temp(Cart::xxxxxzz, 0, m + 1) + 2 * term_xxxxxz;
            R_temp(Cart::xxxxyyyy, 0, m) =
                pma(0) * R_temp(Cart::xxxyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xxxyyyy, 0, m + 1) + 3 * term_xxyyyy;
            R_temp(Cart::xxxxyyyz, 0, m) =
                pma(2) * R_temp(Cart::xxxxyyy, 0, m) +
                wmp(2) * R_temp(Cart::xxxxyyy, 0, m + 1);
            R_temp(Cart::xxxxyyzz, 0, m) =
                pma(1) * R_temp(Cart::xxxxyzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxyzz, 0, m + 1) + term_xxxxzz;
            R_temp(Cart::xxxxyzzz, 0, m) =
                pma(1) * R_temp(Cart::xxxxzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxxzzz, 0, m + 1);
            R_temp(Cart::xxxxzzzz, 0, m) =
                pma(0) * R_temp(Cart::xxxzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xxxzzzz, 0, m + 1) + 3 * term_xxzzzz;
            R_temp(Cart::xxxyyyyy, 0, m) =
                pma(0) * R_temp(Cart::xxyyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xxyyyyy, 0, m + 1) + 2 * term_xyyyyy;
            R_temp(Cart::xxxyyyyz, 0, m) =
                pma(2) * R_temp(Cart::xxxyyyy, 0, m) +
                wmp(2) * R_temp(Cart::xxxyyyy, 0, m + 1);
            R_temp(Cart::xxxyyyzz, 0, m) =
                pma(2) * R_temp(Cart::xxxyyyz, 0, m) +
                wmp(2) * R_temp(Cart::xxxyyyz, 0, m + 1) + term_xxxyyy;
            R_temp(Cart::xxxyyzzz, 0, m) =
                pma(1) * R_temp(Cart::xxxyzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxyzzz, 0, m + 1) + term_xxxzzz;
            R_temp(Cart::xxxyzzzz, 0, m) =
                pma(1) * R_temp(Cart::xxxzzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxxzzzz, 0, m + 1);
            R_temp(Cart::xxxzzzzz, 0, m) =
                pma(0) * R_temp(Cart::xxzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xxzzzzz, 0, m + 1) + 2 * term_xzzzzz;
            R_temp(Cart::xxyyyyyy, 0, m) =
                pma(0) * R_temp(Cart::xyyyyyy, 0, m) +
                wmp(0) * R_temp(Cart::xyyyyyy, 0, m + 1) + term_yyyyyy;
            R_temp(Cart::xxyyyyyz, 0, m) =
                pma(2) * R_temp(Cart::xxyyyyy, 0, m) +
                wmp(2) * R_temp(Cart::xxyyyyy, 0, m + 1);
            R_temp(Cart::xxyyyyzz, 0, m) =
                pma(0) * R_temp(Cart::xyyyyzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyyyzz, 0, m + 1) + term_yyyyzz;
            R_temp(Cart::xxyyyzzz, 0, m) =
                pma(0) * R_temp(Cart::xyyyzzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyyzzz, 0, m + 1) + term_yyyzzz;
            R_temp(Cart::xxyyzzzz, 0, m) =
                pma(0) * R_temp(Cart::xyyzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xyyzzzz, 0, m + 1) + term_yyzzzz;
            R_temp(Cart::xxyzzzzz, 0, m) =
                pma(1) * R_temp(Cart::xxzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::xxzzzzz, 0, m + 1);
            R_temp(Cart::xxzzzzzz, 0, m) =
                pma(0) * R_temp(Cart::xzzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::xzzzzzz, 0, m + 1) + term_zzzzzz;
            R_temp(Cart::xyyyyyyy, 0, m) =
                pma(0) * R_temp(Cart::yyyyyyy, 0, m) +
                wmp(0) * R_temp(Cart::yyyyyyy, 0, m + 1);
            R_temp(Cart::xyyyyyyz, 0, m) =
                pma(0) * R_temp(Cart::yyyyyyz, 0, m) +
                wmp(0) * R_temp(Cart::yyyyyyz, 0, m + 1);
            R_temp(Cart::xyyyyyzz, 0, m) =
                pma(0) * R_temp(Cart::yyyyyzz, 0, m) +
                wmp(0) * R_temp(Cart::yyyyyzz, 0, m + 1);
            R_temp(Cart::xyyyyzzz, 0, m) =
                pma(0) * R_temp(Cart::yyyyzzz, 0, m) +
                wmp(0) * R_temp(Cart::yyyyzzz, 0, m + 1);
            R_temp(Cart::xyyyzzzz, 0, m) =
                pma(0) * R_temp(Cart::yyyzzzz, 0, m) +
                wmp(0) * R_temp(Cart::yyyzzzz, 0, m + 1);
            R_temp(Cart::xyyzzzzz, 0, m) =
                pma(0) * R_temp(Cart::yyzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::yyzzzzz, 0, m + 1);
            R_temp(Cart::xyzzzzzz, 0, m) =
                pma(0) * R_temp(Cart::yzzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::yzzzzzz, 0, m + 1);
            R_temp(Cart::xzzzzzzz, 0, m) =
                pma(0) * R_temp(Cart::zzzzzzz, 0, m) +
                wmp(0) * R_temp(Cart::zzzzzzz, 0, m + 1);
            R_temp(Cart::yyyyyyyy, 0, m) =
                pma(1) * R_temp(Cart::yyyyyyy, 0, m) +
                wmp(1) * R_temp(Cart::yyyyyyy, 0, m + 1) + 7 * term_yyyyyy;
            R_temp(Cart::yyyyyyyz, 0, m) =
                pma(2) * R_temp(Cart::yyyyyyy, 0, m) +
                wmp(2) * R_temp(Cart::yyyyyyy, 0, m + 1);
            R_temp(Cart::yyyyyyzz, 0, m) =
                pma(2) * R_temp(Cart::yyyyyyz, 0, m) +
                wmp(2) * R_temp(Cart::yyyyyyz, 0, m + 1) + term_yyyyyy;
            R_temp(Cart::yyyyyzzz, 0, m) =
                pma(2) * R_temp(Cart::yyyyyzz, 0, m) +
                wmp(2) * R_temp(Cart::yyyyyzz, 0, m + 1) + 2 * term_yyyyyz;
            R_temp(Cart::yyyyzzzz, 0, m) =
                pma(1) * R_temp(Cart::yyyzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yyyzzzz, 0, m + 1) + 3 * term_yyzzzz;
            R_temp(Cart::yyyzzzzz, 0, m) =
                pma(1) * R_temp(Cart::yyzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yyzzzzz, 0, m + 1) + 2 * term_yzzzzz;
            R_temp(Cart::yyzzzzzz, 0, m) =
                pma(1) * R_temp(Cart::yzzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::yzzzzzz, 0, m + 1) + term_zzzzzz;
            R_temp(Cart::yzzzzzzz, 0, m) =
                pma(1) * R_temp(Cart::zzzzzzz, 0, m) +
                wmp(1) * R_temp(Cart::zzzzzzz, 0, m + 1);
            R_temp(Cart::zzzzzzzz, 0, m) =
                pma(2) * R_temp(Cart::zzzzzzz, 0, m) +
                wmp(2) * R_temp(Cart::zzzzzzz, 0, m + 1) + 7 * term_zzzzzz;
          }
        }
        //------------------------------------------------------

        if (lmax_gamma > 0) {

          // Integral  s - s - p
          for (int m = 0; m < lmax_gamma; m++) {
            R_temp(0, Cart::x, m) = wmc(0) * R_temp(0, 0, m + 1);
            R_temp(0, Cart::y, m) = wmc(1) * R_temp(0, 0, m + 1);
            R_temp(0, Cart::z, m) = wmc(2) * R_temp(0, 0, m + 1);
          }
          //------------------------------------------------------

          // Integral  p - s - p
          if (lmax_alpha_beta > 0) {
            for (int m = 0; m < lmax_gamma; m++) {
              double term = rdecay * R_temp(0, 0, m + 1);
              for (int i = 1; i < 4; i++) {
                R_temp(i, Cart::x, m) =
                    wmc(0) * R_temp(i, 0, m + 1) + nx[i] * term;
                R_temp(i, Cart::y, m) =
                    wmc(1) * R_temp(i, 0, m + 1) + ny[i] * term;
                R_temp(i, Cart::z, m) =
                    wmc(2) * R_temp(i, 0, m + 1) + nz[i] * term;
              }
            }
          }
          //------------------------------------------------------

          // Integrals     d - s - p     f - s - p     g - s - p     h - s - p
          // i - s - p     j - s - p     k - s - p
          for (int m = 0; m < lmax_gamma; m++) {
            for (int i = 4; i < n_orbitals[lmax_alpha_beta]; i++) {
              R_temp(i, Cart::x, m) =
                  wmc(0) * R_temp(i, 0, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], 0, m + 1);
              R_temp(i, Cart::y, m) =
                  wmc(1) * R_temp(i, 0, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], 0, m + 1);
              R_temp(i, Cart::z, m) =
                  wmc(2) * R_temp(i, 0, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], 0, m + 1);
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 0)

        if (lmax_gamma > 1) {

          // Integral  s - s - d
          for (int m = 0; m < lmax_gamma - 1; m++) {
            double term =
                rgamma * (R_temp(0, 0, m) - cfak * R_temp(0, 0, m + 1));
            R_temp(0, Cart::xx, m) = wmc(0) * R_temp(0, Cart::x, m + 1) + term;
            R_temp(0, Cart::xy, m) = wmc(0) * R_temp(0, Cart::y, m + 1);
            R_temp(0, Cart::xz, m) = wmc(0) * R_temp(0, Cart::z, m + 1);
            R_temp(0, Cart::yy, m) = wmc(1) * R_temp(0, Cart::y, m + 1) + term;
            R_temp(0, Cart::yz, m) = wmc(1) * R_temp(0, Cart::z, m + 1);
            R_temp(0, Cart::zz, m) = wmc(2) * R_temp(0, Cart::z, m + 1) + term;
          }
          //------------------------------------------------------

          // Integrals     p - s - d     d - s - d     f - s - d     g - s - d
          // h - s - d     i - s - d     j - s - d     k - s - d
          for (int m = 0; m < lmax_gamma - 1; m++) {
            for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
              double term =
                  rgamma * (R_temp(i, 0, m) - cfak * R_temp(i, 0, m + 1));
              R_temp(i, Cart::xx, m) =
                  wmc(0) * R_temp(i, Cart::x, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::x, m + 1) + term;
              R_temp(i, Cart::xy, m) =
                  wmc(0) * R_temp(i, Cart::y, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::y, m + 1);
              R_temp(i, Cart::xz, m) =
                  wmc(0) * R_temp(i, Cart::z, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::z, m + 1);
              R_temp(i, Cart::yy, m) =
                  wmc(1) * R_temp(i, Cart::y, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::y, m + 1) + term;
              R_temp(i, Cart::yz, m) =
                  wmc(1) * R_temp(i, Cart::z, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::z, m + 1);
              R_temp(i, Cart::zz, m) =
                  wmc(2) * R_temp(i, Cart::z, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::z, m + 1) + term;
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 1)

        if (lmax_gamma > 2) {

          // Integral  s - s - f
          for (int m = 0; m < lmax_gamma - 2; m++) {
            R_temp(0, Cart::xxx, m) =
                wmc(0) * R_temp(0, Cart::xx, m + 1) +
                2 * rgamma *
                    (R_temp(0, Cart::x, m) - cfak * R_temp(0, Cart::x, m + 1));
            R_temp(0, Cart::xxy, m) = wmc(1) * R_temp(0, Cart::xx, m + 1);
            R_temp(0, Cart::xxz, m) = wmc(2) * R_temp(0, Cart::xx, m + 1);
            R_temp(0, Cart::xyy, m) = wmc(0) * R_temp(0, Cart::yy, m + 1);
            R_temp(0, Cart::xyz, m) = wmc(0) * R_temp(0, Cart::yz, m + 1);
            R_temp(0, Cart::xzz, m) = wmc(0) * R_temp(0, Cart::zz, m + 1);
            R_temp(0, Cart::yyy, m) =
                wmc(1) * R_temp(0, Cart::yy, m + 1) +
                2 * rgamma *
                    (R_temp(0, Cart::y, m) - cfak * R_temp(0, Cart::y, m + 1));
            R_temp(0, Cart::yyz, m) = wmc(2) * R_temp(0, Cart::yy, m + 1);
            R_temp(0, Cart::yzz, m) = wmc(1) * R_temp(0, Cart::zz, m + 1);
            R_temp(0, Cart::zzz, m) =
                wmc(2) * R_temp(0, Cart::zz, m + 1) +
                2 * rgamma *
                    (R_temp(0, Cart::z, m) - cfak * R_temp(0, Cart::z, m + 1));
          }
          //------------------------------------------------------

          // Integrals     p - s - f     d - s - f     f - s - f     g - s - f
          // h - s - f     i - s - f     j - s - f     k - s - f
          for (int m = 0; m < lmax_gamma - 2; m++) {
            for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
              double term_x =
                  2 * rgamma *
                  (R_temp(i, Cart::x, m) - cfak * R_temp(i, Cart::x, m + 1));
              double term_y =
                  2 * rgamma *
                  (R_temp(i, Cart::y, m) - cfak * R_temp(i, Cart::y, m + 1));
              double term_z =
                  2 * rgamma *
                  (R_temp(i, Cart::z, m) - cfak * R_temp(i, Cart::z, m + 1));
              R_temp(i, Cart::xxx, m) =
                  wmc(0) * R_temp(i, Cart::xx, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xx, m + 1) +
                  term_x;
              R_temp(i, Cart::xxy, m) =
                  wmc(1) * R_temp(i, Cart::xx, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xx, m + 1);
              R_temp(i, Cart::xxz, m) =
                  wmc(2) * R_temp(i, Cart::xx, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xx, m + 1);
              R_temp(i, Cart::xyy, m) =
                  wmc(0) * R_temp(i, Cart::yy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yy, m + 1);
              R_temp(i, Cart::xyz, m) =
                  wmc(0) * R_temp(i, Cart::yz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yz, m + 1);
              R_temp(i, Cart::xzz, m) =
                  wmc(0) * R_temp(i, Cart::zz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::zz, m + 1);
              R_temp(i, Cart::yyy, m) =
                  wmc(1) * R_temp(i, Cart::yy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yy, m + 1) +
                  term_y;
              R_temp(i, Cart::yyz, m) =
                  wmc(2) * R_temp(i, Cart::yy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yy, m + 1);
              R_temp(i, Cart::yzz, m) =
                  wmc(1) * R_temp(i, Cart::zz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::zz, m + 1);
              R_temp(i, Cart::zzz, m) =
                  wmc(2) * R_temp(i, Cart::zz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::zz, m + 1) +
                  term_z;
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 2)

        if (lmax_gamma > 3) {

          // Integral  s - s - g
          for (int m = 0; m < lmax_gamma - 3; m++) {
            double term_xx = rgamma * (R_temp(0, Cart::xx, m) -
                                       cfak * R_temp(0, Cart::xx, m + 1));
            double term_yy = rgamma * (R_temp(0, Cart::yy, m) -
                                       cfak * R_temp(0, Cart::yy, m + 1));
            double term_zz = rgamma * (R_temp(0, Cart::zz, m) -
                                       cfak * R_temp(0, Cart::zz, m + 1));
            R_temp(0, Cart::xxxx, m) =
                wmc(0) * R_temp(0, Cart::xxx, m + 1) + 3 * term_xx;
            R_temp(0, Cart::xxxy, m) = wmc(1) * R_temp(0, Cart::xxx, m + 1);
            R_temp(0, Cart::xxxz, m) = wmc(2) * R_temp(0, Cart::xxx, m + 1);
            R_temp(0, Cart::xxyy, m) =
                wmc(0) * R_temp(0, Cart::xyy, m + 1) + term_yy;
            R_temp(0, Cart::xxyz, m) = wmc(1) * R_temp(0, Cart::xxz, m + 1);
            R_temp(0, Cart::xxzz, m) =
                wmc(0) * R_temp(0, Cart::xzz, m + 1) + term_zz;
            R_temp(0, Cart::xyyy, m) = wmc(0) * R_temp(0, Cart::yyy, m + 1);
            R_temp(0, Cart::xyyz, m) = wmc(0) * R_temp(0, Cart::yyz, m + 1);
            R_temp(0, Cart::xyzz, m) = wmc(0) * R_temp(0, Cart::yzz, m + 1);
            R_temp(0, Cart::xzzz, m) = wmc(0) * R_temp(0, Cart::zzz, m + 1);
            R_temp(0, Cart::yyyy, m) =
                wmc(1) * R_temp(0, Cart::yyy, m + 1) + 3 * term_yy;
            R_temp(0, Cart::yyyz, m) = wmc(2) * R_temp(0, Cart::yyy, m + 1);
            R_temp(0, Cart::yyzz, m) =
                wmc(1) * R_temp(0, Cart::yzz, m + 1) + term_zz;
            R_temp(0, Cart::yzzz, m) = wmc(1) * R_temp(0, Cart::zzz, m + 1);
            R_temp(0, Cart::zzzz, m) =
                wmc(2) * R_temp(0, Cart::zzz, m + 1) + 3 * term_zz;
          }
          //------------------------------------------------------

          // Integrals     p - s - g     d - s - g     f - s - g     g - s - g
          // h - s - g     i - s - g     j - s - g     k - s - g
          for (int m = 0; m < lmax_gamma - 3; m++) {
            for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
              double term_xx = rgamma * (R_temp(i, Cart::xx, m) -
                                         cfak * R_temp(i, Cart::xx, m + 1));
              double term_yy = rgamma * (R_temp(i, Cart::yy, m) -
                                         cfak * R_temp(i, Cart::yy, m + 1));
              double term_zz = rgamma * (R_temp(i, Cart::zz, m) -
                                         cfak * R_temp(i, Cart::zz, m + 1));
              R_temp(i, Cart::xxxx, m) =
                  wmc(0) * R_temp(i, Cart::xxx, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xxx, m + 1) +
                  3 * term_xx;
              R_temp(i, Cart::xxxy, m) =
                  wmc(1) * R_temp(i, Cart::xxx, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxx, m + 1);
              R_temp(i, Cart::xxxz, m) =
                  wmc(2) * R_temp(i, Cart::xxx, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxx, m + 1);
              R_temp(i, Cart::xxyy, m) =
                  wmc(0) * R_temp(i, Cart::xyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xyy, m + 1) +
                  term_yy;
              R_temp(i, Cart::xxyz, m) =
                  wmc(1) * R_temp(i, Cart::xxz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxz, m + 1);
              R_temp(i, Cart::xxzz, m) =
                  wmc(0) * R_temp(i, Cart::xzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xzz, m + 1) +
                  term_zz;
              R_temp(i, Cart::xyyy, m) =
                  wmc(0) * R_temp(i, Cart::yyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyy, m + 1);
              R_temp(i, Cart::xyyz, m) =
                  wmc(0) * R_temp(i, Cart::yyz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyz, m + 1);
              R_temp(i, Cart::xyzz, m) =
                  wmc(0) * R_temp(i, Cart::yzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yzz, m + 1);
              R_temp(i, Cart::xzzz, m) =
                  wmc(0) * R_temp(i, Cart::zzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::zzz, m + 1);
              R_temp(i, Cart::yyyy, m) =
                  wmc(1) * R_temp(i, Cart::yyy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yyy, m + 1) +
                  3 * term_yy;
              R_temp(i, Cart::yyyz, m) =
                  wmc(2) * R_temp(i, Cart::yyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yyy, m + 1);
              R_temp(i, Cart::yyzz, m) =
                  wmc(1) * R_temp(i, Cart::yzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yzz, m + 1) +
                  term_zz;
              R_temp(i, Cart::yzzz, m) =
                  wmc(1) * R_temp(i, Cart::zzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::zzz, m + 1);
              R_temp(i, Cart::zzzz, m) =
                  wmc(2) * R_temp(i, Cart::zzz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::zzz, m + 1) +
                  3 * term_zz;
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 3)

        if (lmax_gamma > 4) {

          // Integral  s - s - h
          for (int m = 0; m < lmax_gamma - 4; m++) {
            double term_xxx = rgamma * (R_temp(0, Cart::xxx, m) -
                                        cfak * R_temp(0, Cart::xxx, m + 1));
            double term_yyy = rgamma * (R_temp(0, Cart::yyy, m) -
                                        cfak * R_temp(0, Cart::yyy, m + 1));
            double term_zzz = rgamma * (R_temp(0, Cart::zzz, m) -
                                        cfak * R_temp(0, Cart::zzz, m + 1));
            R_temp(0, Cart::xxxxx, m) =
                wmc(0) * R_temp(0, Cart::xxxx, m + 1) + 4 * term_xxx;
            R_temp(0, Cart::xxxxy, m) = wmc(1) * R_temp(0, Cart::xxxx, m + 1);
            R_temp(0, Cart::xxxxz, m) = wmc(2) * R_temp(0, Cart::xxxx, m + 1);
            R_temp(0, Cart::xxxyy, m) =
                wmc(1) * R_temp(0, Cart::xxxy, m + 1) + term_xxx;
            R_temp(0, Cart::xxxyz, m) = wmc(1) * R_temp(0, Cart::xxxz, m + 1);
            R_temp(0, Cart::xxxzz, m) =
                wmc(2) * R_temp(0, Cart::xxxz, m + 1) + term_xxx;
            R_temp(0, Cart::xxyyy, m) =
                wmc(0) * R_temp(0, Cart::xyyy, m + 1) + term_yyy;
            R_temp(0, Cart::xxyyz, m) = wmc(2) * R_temp(0, Cart::xxyy, m + 1);
            R_temp(0, Cart::xxyzz, m) = wmc(1) * R_temp(0, Cart::xxzz, m + 1);
            R_temp(0, Cart::xxzzz, m) =
                wmc(0) * R_temp(0, Cart::xzzz, m + 1) + term_zzz;
            R_temp(0, Cart::xyyyy, m) = wmc(0) * R_temp(0, Cart::yyyy, m + 1);
            R_temp(0, Cart::xyyyz, m) = wmc(0) * R_temp(0, Cart::yyyz, m + 1);
            R_temp(0, Cart::xyyzz, m) = wmc(0) * R_temp(0, Cart::yyzz, m + 1);
            R_temp(0, Cart::xyzzz, m) = wmc(0) * R_temp(0, Cart::yzzz, m + 1);
            R_temp(0, Cart::xzzzz, m) = wmc(0) * R_temp(0, Cart::zzzz, m + 1);
            R_temp(0, Cart::yyyyy, m) =
                wmc(1) * R_temp(0, Cart::yyyy, m + 1) + 4 * term_yyy;
            R_temp(0, Cart::yyyyz, m) = wmc(2) * R_temp(0, Cart::yyyy, m + 1);
            R_temp(0, Cart::yyyzz, m) =
                wmc(2) * R_temp(0, Cart::yyyz, m + 1) + term_yyy;
            R_temp(0, Cart::yyzzz, m) =
                wmc(1) * R_temp(0, Cart::yzzz, m + 1) + term_zzz;
            R_temp(0, Cart::yzzzz, m) = wmc(1) * R_temp(0, Cart::zzzz, m + 1);
            R_temp(0, Cart::zzzzz, m) =
                wmc(2) * R_temp(0, Cart::zzzz, m + 1) + 4 * term_zzz;
          }
          //------------------------------------------------------

          // Integrals     p - s - h     d - s - h     f - s - h     g - s - h
          // h - s - h     i - s - h     j - s - h     k - s - h
          for (int m = 0; m < lmax_gamma - 4; m++) {
            for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
              double term_xxx = rgamma * (R_temp(i, Cart::xxx, m) -
                                          cfak * R_temp(i, Cart::xxx, m + 1));
              double term_yyy = rgamma * (R_temp(i, Cart::yyy, m) -
                                          cfak * R_temp(i, Cart::yyy, m + 1));
              double term_zzz = rgamma * (R_temp(i, Cart::zzz, m) -
                                          cfak * R_temp(i, Cart::zzz, m + 1));
              R_temp(i, Cart::xxxxx, m) =
                  wmc(0) * R_temp(i, Cart::xxxx, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xxxx, m + 1) +
                  4 * term_xxx;
              R_temp(i, Cart::xxxxy, m) =
                  wmc(1) * R_temp(i, Cart::xxxx, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxx, m + 1);
              R_temp(i, Cart::xxxxz, m) =
                  wmc(2) * R_temp(i, Cart::xxxx, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxxx, m + 1);
              R_temp(i, Cart::xxxyy, m) =
                  wmc(1) * R_temp(i, Cart::xxxy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxy, m + 1) +
                  term_xxx;
              R_temp(i, Cart::xxxyz, m) =
                  wmc(1) * R_temp(i, Cart::xxxz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxz, m + 1);
              R_temp(i, Cart::xxxzz, m) =
                  wmc(2) * R_temp(i, Cart::xxxz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxxz, m + 1) +
                  term_xxx;
              R_temp(i, Cart::xxyyy, m) =
                  wmc(0) * R_temp(i, Cart::xyyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xyyy, m + 1) +
                  term_yyy;
              R_temp(i, Cart::xxyyz, m) =
                  wmc(2) * R_temp(i, Cart::xxyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxyy, m + 1);
              R_temp(i, Cart::xxyzz, m) =
                  wmc(1) * R_temp(i, Cart::xxzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxzz, m + 1);
              R_temp(i, Cart::xxzzz, m) =
                  wmc(0) * R_temp(i, Cart::xzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xzzz, m + 1) +
                  term_zzz;
              R_temp(i, Cart::xyyyy, m) =
                  wmc(0) * R_temp(i, Cart::yyyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyyy, m + 1);
              R_temp(i, Cart::xyyyz, m) =
                  wmc(0) * R_temp(i, Cart::yyyz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyyz, m + 1);
              R_temp(i, Cart::xyyzz, m) =
                  wmc(0) * R_temp(i, Cart::yyzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyzz, m + 1);
              R_temp(i, Cart::xyzzz, m) =
                  wmc(0) * R_temp(i, Cart::yzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yzzz, m + 1);
              R_temp(i, Cart::xzzzz, m) =
                  wmc(0) * R_temp(i, Cart::zzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::zzzz, m + 1);
              R_temp(i, Cart::yyyyy, m) =
                  wmc(1) * R_temp(i, Cart::yyyy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yyyy, m + 1) +
                  4 * term_yyy;
              R_temp(i, Cart::yyyyz, m) =
                  wmc(2) * R_temp(i, Cart::yyyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yyyy, m + 1);
              R_temp(i, Cart::yyyzz, m) =
                  wmc(2) * R_temp(i, Cart::yyyz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yyyz, m + 1) +
                  term_yyy;
              R_temp(i, Cart::yyzzz, m) =
                  wmc(1) * R_temp(i, Cart::yzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yzzz, m + 1) +
                  term_zzz;
              R_temp(i, Cart::yzzzz, m) =
                  wmc(1) * R_temp(i, Cart::zzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::zzzz, m + 1);
              R_temp(i, Cart::zzzzz, m) =
                  wmc(2) * R_temp(i, Cart::zzzz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::zzzz, m + 1) +
                  4 * term_zzz;
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 4)

        if (lmax_gamma > 5) {

          // Integral  s - s - i
          for (int m = 0; m < lmax_gamma - 5; m++) {
            double term_xxxx = rgamma * (R_temp(0, Cart::xxxx, m) -
                                         cfak * R_temp(0, Cart::xxxx, m + 1));
            double term_xyyy = rgamma * (R_temp(0, Cart::xyyy, m) -
                                         cfak * R_temp(0, Cart::xyyy, m + 1));
            double term_xzzz = rgamma * (R_temp(0, Cart::xzzz, m) -
                                         cfak * R_temp(0, Cart::xzzz, m + 1));
            double term_yyyy = rgamma * (R_temp(0, Cart::yyyy, m) -
                                         cfak * R_temp(0, Cart::yyyy, m + 1));
            double term_yyzz = rgamma * (R_temp(0, Cart::yyzz, m) -
                                         cfak * R_temp(0, Cart::yyzz, m + 1));
            double term_yzzz = rgamma * (R_temp(0, Cart::yzzz, m) -
                                         cfak * R_temp(0, Cart::yzzz, m + 1));
            double term_zzzz = rgamma * (R_temp(0, Cart::zzzz, m) -
                                         cfak * R_temp(0, Cart::zzzz, m + 1));
            R_temp(0, Cart::xxxxxx, m) =
                wmc(0) * R_temp(0, Cart::xxxxx, m + 1) + 5 * term_xxxx;
            R_temp(0, Cart::xxxxxy, m) = wmc(1) * R_temp(0, Cart::xxxxx, m + 1);
            R_temp(0, Cart::xxxxxz, m) = wmc(2) * R_temp(0, Cart::xxxxx, m + 1);
            R_temp(0, Cart::xxxxyy, m) =
                wmc(1) * R_temp(0, Cart::xxxxy, m + 1) + term_xxxx;
            R_temp(0, Cart::xxxxyz, m) = wmc(1) * R_temp(0, Cart::xxxxz, m + 1);
            R_temp(0, Cart::xxxxzz, m) =
                wmc(2) * R_temp(0, Cart::xxxxz, m + 1) + term_xxxx;
            R_temp(0, Cart::xxxyyy, m) =
                wmc(0) * R_temp(0, Cart::xxyyy, m + 1) + 2 * term_xyyy;
            R_temp(0, Cart::xxxyyz, m) = wmc(2) * R_temp(0, Cart::xxxyy, m + 1);
            R_temp(0, Cart::xxxyzz, m) = wmc(1) * R_temp(0, Cart::xxxzz, m + 1);
            R_temp(0, Cart::xxxzzz, m) =
                wmc(0) * R_temp(0, Cart::xxzzz, m + 1) + 2 * term_xzzz;
            R_temp(0, Cart::xxyyyy, m) =
                wmc(0) * R_temp(0, Cart::xyyyy, m + 1) + term_yyyy;
            R_temp(0, Cart::xxyyyz, m) = wmc(2) * R_temp(0, Cart::xxyyy, m + 1);
            R_temp(0, Cart::xxyyzz, m) =
                wmc(0) * R_temp(0, Cart::xyyzz, m + 1) + term_yyzz;
            R_temp(0, Cart::xxyzzz, m) = wmc(1) * R_temp(0, Cart::xxzzz, m + 1);
            R_temp(0, Cart::xxzzzz, m) =
                wmc(0) * R_temp(0, Cart::xzzzz, m + 1) + term_zzzz;
            R_temp(0, Cart::xyyyyy, m) = wmc(0) * R_temp(0, Cart::yyyyy, m + 1);
            R_temp(0, Cart::xyyyyz, m) = wmc(0) * R_temp(0, Cart::yyyyz, m + 1);
            R_temp(0, Cart::xyyyzz, m) = wmc(0) * R_temp(0, Cart::yyyzz, m + 1);
            R_temp(0, Cart::xyyzzz, m) = wmc(0) * R_temp(0, Cart::yyzzz, m + 1);
            R_temp(0, Cart::xyzzzz, m) = wmc(0) * R_temp(0, Cart::yzzzz, m + 1);
            R_temp(0, Cart::xzzzzz, m) = wmc(0) * R_temp(0, Cart::zzzzz, m + 1);
            R_temp(0, Cart::yyyyyy, m) =
                wmc(1) * R_temp(0, Cart::yyyyy, m + 1) + 5 * term_yyyy;
            R_temp(0, Cart::yyyyyz, m) = wmc(2) * R_temp(0, Cart::yyyyy, m + 1);
            R_temp(0, Cart::yyyyzz, m) =
                wmc(2) * R_temp(0, Cart::yyyyz, m + 1) + term_yyyy;
            R_temp(0, Cart::yyyzzz, m) =
                wmc(1) * R_temp(0, Cart::yyzzz, m + 1) + 2 * term_yzzz;
            R_temp(0, Cart::yyzzzz, m) =
                wmc(1) * R_temp(0, Cart::yzzzz, m + 1) + term_zzzz;
            R_temp(0, Cart::yzzzzz, m) = wmc(1) * R_temp(0, Cart::zzzzz, m + 1);
            R_temp(0, Cart::zzzzzz, m) =
                wmc(2) * R_temp(0, Cart::zzzzz, m + 1) + 5 * term_zzzz;
          }
          //------------------------------------------------------

          // Integrals     p - s - i     d - s - i     f - s - i     g - s - i
          // h - s - i     i - s - i     j - s - i     k - s - i
          for (int m = 0; m < lmax_gamma - 5; m++) {
            for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
              double term_xxxx = rgamma * (R_temp(i, Cart::xxxx, m) -
                                           cfak * R_temp(i, Cart::xxxx, m + 1));
              double term_xyyy = rgamma * (R_temp(i, Cart::xyyy, m) -
                                           cfak * R_temp(i, Cart::xyyy, m + 1));
              double term_xzzz = rgamma * (R_temp(i, Cart::xzzz, m) -
                                           cfak * R_temp(i, Cart::xzzz, m + 1));
              double term_yyyy = rgamma * (R_temp(i, Cart::yyyy, m) -
                                           cfak * R_temp(i, Cart::yyyy, m + 1));
              double term_yyzz = rgamma * (R_temp(i, Cart::yyzz, m) -
                                           cfak * R_temp(i, Cart::yyzz, m + 1));
              double term_yzzz = rgamma * (R_temp(i, Cart::yzzz, m) -
                                           cfak * R_temp(i, Cart::yzzz, m + 1));
              double term_zzzz = rgamma * (R_temp(i, Cart::zzzz, m) -
                                           cfak * R_temp(i, Cart::zzzz, m + 1));
              R_temp(i, Cart::xxxxxx, m) =
                  wmc(0) * R_temp(i, Cart::xxxxx, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xxxxx, m + 1) +
                  5 * term_xxxx;
              R_temp(i, Cart::xxxxxy, m) =
                  wmc(1) * R_temp(i, Cart::xxxxx, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxxx, m + 1);
              R_temp(i, Cart::xxxxxz, m) =
                  wmc(2) * R_temp(i, Cart::xxxxx, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxxxx, m + 1);
              R_temp(i, Cart::xxxxyy, m) =
                  wmc(1) * R_temp(i, Cart::xxxxy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxxy, m + 1) +
                  term_xxxx;
              R_temp(i, Cart::xxxxyz, m) =
                  wmc(1) * R_temp(i, Cart::xxxxz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxxz, m + 1);
              R_temp(i, Cart::xxxxzz, m) =
                  wmc(2) * R_temp(i, Cart::xxxxz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxxxz, m + 1) +
                  term_xxxx;
              R_temp(i, Cart::xxxyyy, m) =
                  wmc(0) * R_temp(i, Cart::xxyyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xxyyy, m + 1) +
                  2 * term_xyyy;
              R_temp(i, Cart::xxxyyz, m) =
                  wmc(2) * R_temp(i, Cart::xxxyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxxyy, m + 1);
              R_temp(i, Cart::xxxyzz, m) =
                  wmc(1) * R_temp(i, Cart::xxxzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxxzz, m + 1);
              R_temp(i, Cart::xxxzzz, m) =
                  wmc(0) * R_temp(i, Cart::xxzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xxzzz, m + 1) +
                  2 * term_xzzz;
              R_temp(i, Cart::xxyyyy, m) =
                  wmc(0) * R_temp(i, Cart::xyyyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xyyyy, m + 1) +
                  term_yyyy;
              R_temp(i, Cart::xxyyyz, m) =
                  wmc(2) * R_temp(i, Cart::xxyyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::xxyyy, m + 1);
              R_temp(i, Cart::xxyyzz, m) =
                  wmc(0) * R_temp(i, Cart::xyyzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xyyzz, m + 1) +
                  term_yyzz;
              R_temp(i, Cart::xxyzzz, m) =
                  wmc(1) * R_temp(i, Cart::xxzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::xxzzz, m + 1);
              R_temp(i, Cart::xxzzzz, m) =
                  wmc(0) * R_temp(i, Cart::xzzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::xzzzz, m + 1) +
                  term_zzzz;
              R_temp(i, Cart::xyyyyy, m) =
                  wmc(0) * R_temp(i, Cart::yyyyy, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyyyy, m + 1);
              R_temp(i, Cart::xyyyyz, m) =
                  wmc(0) * R_temp(i, Cart::yyyyz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyyyz, m + 1);
              R_temp(i, Cart::xyyyzz, m) =
                  wmc(0) * R_temp(i, Cart::yyyzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyyzz, m + 1);
              R_temp(i, Cart::xyyzzz, m) =
                  wmc(0) * R_temp(i, Cart::yyzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yyzzz, m + 1);
              R_temp(i, Cart::xyzzzz, m) =
                  wmc(0) * R_temp(i, Cart::yzzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::yzzzz, m + 1);
              R_temp(i, Cart::xzzzzz, m) =
                  wmc(0) * R_temp(i, Cart::zzzzz, m + 1) +
                  nx[i] * rdecay * R_temp(i_less_x[i], Cart::zzzzz, m + 1);
              R_temp(i, Cart::yyyyyy, m) =
                  wmc(1) * R_temp(i, Cart::yyyyy, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yyyyy, m + 1) +
                  5 * term_yyyy;
              R_temp(i, Cart::yyyyyz, m) =
                  wmc(2) * R_temp(i, Cart::yyyyy, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yyyyy, m + 1);
              R_temp(i, Cart::yyyyzz, m) =
                  wmc(2) * R_temp(i, Cart::yyyyz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::yyyyz, m + 1) +
                  term_yyyy;
              R_temp(i, Cart::yyyzzz, m) =
                  wmc(1) * R_temp(i, Cart::yyzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yyzzz, m + 1) +
                  2 * term_yzzz;
              R_temp(i, Cart::yyzzzz, m) =
                  wmc(1) * R_temp(i, Cart::yzzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::yzzzz, m + 1) +
                  term_zzzz;
              R_temp(i, Cart::yzzzzz, m) =
                  wmc(1) * R_temp(i, Cart::zzzzz, m + 1) +
                  ny[i] * rdecay * R_temp(i_less_y[i], Cart::zzzzz, m + 1);
              R_temp(i, Cart::zzzzzz, m) =
                  wmc(2) * R_temp(i, Cart::zzzzz, m + 1) +
                  nz[i] * rdecay * R_temp(i_less_z[i], Cart::zzzzz, m + 1) +
                  5 * term_zzzz;
            }
          }
          //------------------------------------------------------

        }  // end if (lmax_gamma > 5)

        const Eigen::VectorXd& contractions_gamma =
            gaussian_gamma.getContraction();

        // s-functions
        double factor = contractions_gamma[0];
        for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
          R_temp(i, 0, 1) = factor * R_temp(i, 0, 0);  /// Y 0,0
        }

        if (lmax_gamma > 0) {
          // p-functions
          factor = 2. * sqrt(decay_gamma) * contractions_gamma[1];
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 1, 1) = factor * R_temp(i, 3, 0);  /// Y 1,0
            R_temp(i, 2, 1) = factor * R_temp(i, 2, 0);  /// Y 1,-1
            R_temp(i, 3, 1) = factor * R_temp(i, 1, 0);  /// Y 1,1
          }
        }

        if (lmax_gamma > 1) {
          // d-functions
          factor = 2.0 * decay_gamma * contractions_gamma[2];
          double factor_1 = factor / sqrt(3.0);
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 4, 1) =
                factor_1 *
                (2.0 * R_temp(i, Cart::zz, 0) - R_temp(i, Cart::xx, 0) -
                 R_temp(i, Cart::yy, 0));  /// d3z2-r2  Y 2,0
            R_temp(i, 5, 1) =
                2. * factor * R_temp(i, Cart::yz, 0);  /// dyz  Y 2,-1
            R_temp(i, 6, 1) =
                2. * factor * R_temp(i, Cart::xz, 0);  /// dxz  Y 2,1
            R_temp(i, 7, 1) =
                2. * factor * R_temp(i, Cart::xy, 0);  /// dxy  Y 2,-2
            R_temp(i, 8, 1) =
                factor * (R_temp(i, Cart::xx, 0) -
                          R_temp(i, Cart::yy, 0));  /// dx2-y2  Y 2,2
          }
        }

        if (lmax_gamma > 2) {
          // f-functions
          factor = 2.0 * pow(decay_gamma, 1.5) * contractions_gamma[3];
          double factor_1 = factor * 2. / sqrt(15.);
          double factor_2 = factor * sqrt(2.) / sqrt(5.);
          double factor_3 = factor * sqrt(2.) / sqrt(3.);
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 9, 1) =
                factor_1 *
                (2. * R_temp(i, Cart::zzz, 0) - 3. * R_temp(i, Cart::xxz, 0) -
                 3. * R_temp(i, Cart::yyz, 0));  /// Y 3,0
            R_temp(i, 10, 1) = factor_2 * (4. * R_temp(i, Cart::yzz, 0) -
                                           R_temp(i, Cart::xxy, 0) -
                                           R_temp(i, Cart::yyy, 0));  /// Y 3,-1
            R_temp(i, 11, 1) = factor_2 * (4. * R_temp(i, Cart::xzz, 0) -
                                           R_temp(i, Cart::xxx, 0) -
                                           R_temp(i, Cart::xyy, 0));  /// Y 3,1
            R_temp(i, 12, 1) =
                4. * factor * R_temp(i, Cart::xyz, 0);  /// Y 3,-2
            R_temp(i, 13, 1) =
                2. * factor *
                (R_temp(i, Cart::xxz, 0) - R_temp(i, Cart::yyz, 0));  /// Y 3,2
            R_temp(i, 14, 1) = factor_3 * (3. * R_temp(i, Cart::xxy, 0) -
                                           R_temp(i, Cart::yyy, 0));  /// Y 3,-3
            R_temp(i, 15, 1) =
                factor_3 * (R_temp(i, Cart::xxx, 0) -
                            3. * R_temp(i, Cart::xyy, 0));  /// Y 3,3
          }
        }

        if (lmax_gamma > 3) {
          // g-functions
          factor =
              2. / sqrt(3.) * decay_gamma * decay_gamma * contractions_gamma[4];
          double factor_1 = factor / sqrt(35.);
          double factor_2 = factor * 4. / sqrt(14.);
          double factor_3 = factor * 2. / sqrt(7.);
          double factor_4 = factor * 2. * sqrt(2.);
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 16, 1) =
                factor_1 *
                (3. * (R_temp(i, Cart::xxxx, 0) + R_temp(i, Cart::yyyy, 0)) +
                 6. * R_temp(i, Cart::xxyy, 0) -
                 24. * (R_temp(i, Cart::xxzz, 0) + R_temp(i, Cart::yyzz, 0)) +
                 8. * R_temp(i, Cart::zzzz, 0));  /// Y 4,0
            R_temp(i, 17, 1) =
                factor_2 *
                (-3. * (R_temp(i, Cart::xxyz, 0) + R_temp(i, Cart::yyyz, 0)) +
                 4. * R_temp(i, Cart::yzzz, 0));  /// Y 4,-1
            R_temp(i, 18, 1) =
                factor_2 *
                (-3. * (R_temp(i, Cart::xxxz, 0) + R_temp(i, Cart::xyyz, 0)) +
                 4. * R_temp(i, Cart::xzzz, 0));  /// Y 4,1
            R_temp(i, 19, 1) =
                2. * factor_3 *
                (-R_temp(i, Cart::xxxy, 0) - R_temp(i, Cart::xyyy, 0) +
                 6. * R_temp(i, Cart::xyzz, 0));  /// Y 4,-2
            R_temp(i, 20, 1) =
                factor_3 *
                (-R_temp(i, Cart::xxxx, 0) +
                 6. * (R_temp(i, Cart::xxzz, 0) - R_temp(i, Cart::yyzz, 0)) +
                 R_temp(i, Cart::yyyy, 0));  /// Y 4,2
            R_temp(i, 21, 1) =
                factor_4 * (3. * R_temp(i, Cart::xxyz, 0) -
                            R_temp(i, Cart::yyyz, 0));  /// Y 4,-3
            R_temp(i, 22, 1) =
                factor_4 * (R_temp(i, Cart::xxxz, 0) -
                            3. * R_temp(i, Cart::xyyz, 0));  /// Y 4,3
            R_temp(i, 23, 1) = 4. * factor *
                               (R_temp(i, Cart::xxxy, 0) -
                                R_temp(i, Cart::xyyy, 0));  /// Y 4,-4
            R_temp(i, 24, 1) = factor * (R_temp(i, Cart::xxxx, 0) -
                                         6. * R_temp(i, Cart::xxyy, 0) +
                                         R_temp(i, Cart::yyyy, 0));  /// Y 4,4
          }
        }

        if (lmax_gamma > 4) {
          // h-functions
          factor = (2. / 3.) * pow(decay_gamma, 2.5) * contractions_gamma[5];
          double factor_1 = factor * 2. / sqrt(105.);
          double factor_2 = factor * 2. / sqrt(7.);
          double factor_3 = factor * sqrt(6.) / 3.;
          double factor_4 = factor * 2. * sqrt(3.);
          double factor_5 = factor * .2 * sqrt(30.);
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 25, 1) =
                factor_1 *
                (15. * (R_temp(i, Cart::xxxxz, 0) + R_temp(i, Cart::yyyyz, 0)) +
                 30. * R_temp(i, Cart::xxyyz, 0) -
                 40. * (R_temp(i, Cart::xxzzz, 0) + R_temp(i, Cart::yyzzz, 0)) +
                 8. * R_temp(i, Cart::zzzzz, 0));  /// Y 5,0

            R_temp(i, 26, 1) =
                factor_2 *
                (R_temp(i, Cart::xxxxy, 0) + 2. * R_temp(i, Cart::xxyyy, 0) -
                 12. * (R_temp(i, Cart::xxyzz, 0) + R_temp(i, Cart::yyyzz, 0)) +
                 R_temp(i, Cart::yyyyy, 0) +
                 8. * R_temp(i, Cart::yzzzz, 0));  /// Y 5,-1

            R_temp(i, 27, 1) =
                factor_2 *
                (R_temp(i, Cart::xxxxx, 0) + 2. * R_temp(i, Cart::xxxyy, 0) -
                 12. * (R_temp(i, Cart::xxxzz, 0) + R_temp(i, Cart::xyyzz, 0)) +
                 R_temp(i, Cart::xyyyy, 0) +
                 8. * R_temp(i, Cart::xzzzz, 0));  /// Y 5,1

            R_temp(i, 28, 1) =
                8. * factor *
                (-R_temp(i, Cart::xxxyz, 0) - R_temp(i, Cart::xyyyz, 0) +
                 2. * R_temp(i, Cart::xyzzz, 0));  /// Y 5,-2

            R_temp(i, 29, 1) =
                4. * factor *
                (-R_temp(i, Cart::xxxxz, 0) +
                 2. * (R_temp(i, Cart::xxzzz, 0) - R_temp(i, Cart::yyzzz, 0)) +
                 R_temp(i, Cart::yyyyz, 0));  /// Y 5,2

            R_temp(i, 30, 1) =
                factor_3 *
                (-3. * R_temp(i, Cart::xxxxy, 0) -
                 2. * R_temp(i, Cart::xxyyy, 0) +
                 24. * R_temp(i, Cart::xxyzz, 0) + R_temp(i, Cart::yyyyy, 0) -
                 8. * R_temp(i, Cart::yyyzz, 0));  /// Y 5,-3

            R_temp(i, 31, 1) =
                factor_3 *
                (-R_temp(i, Cart::xxxxx, 0) + 2. * R_temp(i, Cart::xxxyy, 0) +
                 8. * R_temp(i, Cart::xxxzz, 0) +
                 3. * R_temp(i, Cart::xyyyy, 0) -
                 24. * R_temp(i, Cart::xyyzz, 0));  /// Y 5,3

            R_temp(i, 32, 1) = 4. * factor_4 *
                               (R_temp(i, Cart::xxxyz, 0) -
                                R_temp(i, Cart::xyyyz, 0));  /// Y 5,-4

            R_temp(i, 33, 1) =
                factor_4 *
                (R_temp(i, Cart::xxxxz, 0) - 6. * R_temp(i, Cart::xxyyz, 0) +
                 R_temp(i, Cart::yyyyz, 0));  /// Y 5,4

            R_temp(i, 34, 1) =
                factor_5 * (5. * R_temp(i, Cart::xxxxy, 0) -
                            10. * R_temp(i, Cart::xxyyy, 0) +
                            R_temp(i, Cart::yyyyy, 0));  /// Y 5,-5

            R_temp(i, 35, 1) =
                factor_5 *
                (R_temp(i, Cart::xxxxx, 0) - 10. * R_temp(i, Cart::xxxyy, 0) +
                 5. * R_temp(i, Cart::xyyyy, 0));  /// Y 5,5
          }
        }

        if (lmax_gamma > 5) {
          // i-functions
          factor = (2. / 3.) * decay_gamma * decay_gamma * decay_gamma *
                   contractions_gamma[6];
          double factor_1 = factor * 2. / sqrt(1155.);
          double factor_2 = factor * 4. / sqrt(55.);
          double factor_3 = factor * sqrt(22.) / 11.;
          double factor_4 = factor * 2. * sqrt(165.) / 55.;
          double factor_5 = factor * .4 * sqrt(30.);
          double factor_6 = factor * .2 * sqrt(10.);
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; i++) {
            R_temp(i, 36, 1) =
                factor_1 *
                (-5. *
                     (R_temp(i, Cart::xxxxxx, 0) + R_temp(i, Cart::yyyyyy, 0)) -
                 15. *
                     (R_temp(i, Cart::xxxxyy, 0) + R_temp(i, Cart::xxyyyy, 0)) +
                 90. *
                     (R_temp(i, Cart::xxxxzz, 0) + R_temp(i, Cart::yyyyzz, 0)) +
                 180. * R_temp(i, Cart::xxyyzz, 0) -
                 120. *
                     (R_temp(i, Cart::xxzzzz, 0) + R_temp(i, Cart::yyzzzz, 0)) +
                 16. * R_temp(i, Cart::zzzzzz, 0));  /// Y 6,0

            R_temp(i, 37, 1) =
                factor_2 * (5. * (R_temp(i, Cart::xxxxyz, 0) +
                                  R_temp(i, Cart::yyyyyz, 0)) +
                            10. * R_temp(i, Cart::xxyyyz, 0) -
                            20. * (R_temp(i, Cart::xxyzzz, 0) +
                                   R_temp(i, Cart::yyyzzz, 0)) +
                            8. * R_temp(i, Cart::yzzzzz, 0));  /// Y 6,-1

            R_temp(i, 38, 1) =
                factor_2 * (5. * (R_temp(i, Cart::xxxxxz, 0) +
                                  R_temp(i, Cart::xyyyyz, 0)) +
                            10. * R_temp(i, Cart::xxxyyz, 0) -
                            20. * (R_temp(i, Cart::xxxzzz, 0) +
                                   R_temp(i, Cart::xyyzzz, 0)) +
                            8. * R_temp(i, Cart::xzzzzz, 0));  /// Y 6,1

            R_temp(i, 39, 1) =
                2. * factor_3 *
                (R_temp(i, Cart::xxxxxy, 0) + 2. * R_temp(i, Cart::xxxyyy, 0) -
                 16. *
                     (R_temp(i, Cart::xxxyzz, 0) + R_temp(i, Cart::xyyyzz, 0) -
                      R_temp(i, Cart::xyzzzz, 0)) +
                 R_temp(i, Cart::xyyyyy, 0));  /// Y 6,-2

            R_temp(i, 40, 1) =
                factor_3 *
                (R_temp(i, Cart::xxxxxy, 0) + R_temp(i, Cart::xxxxyy, 0) -
                 16. *
                     (R_temp(i, Cart::xxxxzz, 0) - R_temp(i, Cart::xxzzzz, 0) -
                      R_temp(i, Cart::yyyyzz, 0) + R_temp(i, Cart::yyzzzz, 0)) -
                 R_temp(i, Cart::xxyyyy, 0) -
                 R_temp(i, Cart::yyyyyy, 0));  /// Y 6,2

            R_temp(i, 41, 1) = 2. * factor_3 *
                               (-9. * R_temp(i, Cart::xxxxyz, 0) -
                                6. * R_temp(i, Cart::xxyyyz, 0) +
                                24. * R_temp(i, Cart::xxyzzz, 0) +
                                3. * R_temp(i, Cart::yyyyyz, 0) -
                                8. * R_temp(i, Cart::yyyzzz, 0));  /// Y 6,-3

            R_temp(i, 42, 1) = 2. * factor_3 *
                               (-3. * R_temp(i, Cart::xxxxxz, 0) +
                                6. * R_temp(i, Cart::xxxyyz, 0) +
                                8. * R_temp(i, Cart::xxxzzz, 0) +
                                9. * R_temp(i, Cart::xyyyyz, 0) -
                                24. * R_temp(i, Cart::xyyzzz, 0));  /// Y 6,3

            R_temp(i, 43, 1) = 4. * factor_4 *
                               (-R_temp(i, Cart::xxxxxy, 0) +
                                10. * (R_temp(i, Cart::xxxyzz, 0) -
                                       R_temp(i, Cart::xyyyzz, 0)) +
                                R_temp(i, Cart::xyyyyy, 0));  /// Y 6,-4

            R_temp(i, 44, 1) =
                factor_4 * (-R_temp(i, Cart::xxxxxx, 0) +
                            5. * (R_temp(i, Cart::xxxxyy, 0) +
                                  R_temp(i, Cart::xxyyyy, 0)) +
                            10. * (R_temp(i, Cart::xxxxzz, 0) +
                                   R_temp(i, Cart::yyyyzz, 0)) -
                            60. * R_temp(i, Cart::xxyyzz, 0) -
                            R_temp(i, Cart::yyyyyy, 0));  /// Y 6,4

            R_temp(i, 45, 1) =
                factor_5 * (5. * R_temp(i, Cart::xxxxyz, 0) -
                            10. * R_temp(i, Cart::xxyyyz, 0) +
                            R_temp(i, Cart::yyyyyz, 0));  /// Y 6,-5

            R_temp(i, 46, 1) =
                factor_5 *
                (R_temp(i, Cart::xxxxxz, 0) - 10. * R_temp(i, Cart::xxxyyz, 0) +
                 5. * R_temp(i, Cart::xyyyyz, 0));  /// Y 6,5

            R_temp(i, 47, 1) = 2. * factor_6 *
                               (3. * R_temp(i, Cart::xxxxxy, 0) -
                                10. * R_temp(i, Cart::xxxyyy, 0) +
                                3. * R_temp(i, Cart::xyyyyy, 0));  /// Y 6,-6

            R_temp(i, 48, 1) =
                factor_6 * (R_temp(i, Cart::xxxxxx, 0) -
                            15. * (R_temp(i, Cart::xxxxyy, 0) -
                                   R_temp(i, Cart::xxyyyy, 0)) -
                            R_temp(i, Cart::yyyyyy, 0));  /// Y 6,6
          }
        }

        // copy into new array for 3D use.

        int gamma_num_func = shell_gamma->getNumFunc();
        Eigen::Tensor<double, 3> R(ncombined, nbeta, gamma_num_func);
        R.setZero();
        for (int k = 0; k < gamma_num_func; ++k) {
          for (int i = 0; i < n_orbitals[lmax_alpha_beta]; ++i) {
            R(i, 0, k) = R_temp(i, k + shell_gamma->getOffset(), 1);
          }
        }

        if (lmax_beta > 0) {
          // Integrals    s - p - *    p - p - *    d - p - *    f - p - *    g
          // - p - *    h - p - *    i - p - *    j - p - *
          for (int i = 0; i < gamma_num_func; i++) {
            for (int j = 0; j < n_orbitals[lmax_alpha_beta - 1]; j++) {
              R(j, Cart::x, i) = R(i_more_x[j], 0, i) + amb(0) * R(j, 0, i);
              R(j, Cart::y, i) = R(i_more_y[j], 0, i) + amb(1) * R(j, 0, i);
              R(j, Cart::z, i) = R(i_more_z[j], 0, i) + amb(2) * R(j, 0, i);
            }
          }
          //------------------------------------------------------
        }

        if (lmax_beta > 1) {
          // Integrals    s - d - *    p - d - *    d - d - *    f - d - *    g
          // - d - *    h - d - *    i - d - *
          for (int i = 0; i < gamma_num_func; i++) {
            for (int j = 0; j < n_orbitals[lmax_alpha_beta - 2]; j++) {
              R(j, Cart::xx, i) =
                  R(i_more_x[j], Cart::x, i) + amb(0) * R(j, Cart::x, i);
              R(j, Cart::xy, i) =
                  R(i_more_x[j], Cart::y, i) + amb(0) * R(j, Cart::y, i);
              R(j, Cart::xz, i) =
                  R(i_more_x[j], Cart::z, i) + amb(0) * R(j, Cart::z, i);
              R(j, Cart::yy, i) =
                  R(i_more_y[j], Cart::y, i) + amb(1) * R(j, Cart::y, i);
              R(j, Cart::yz, i) =
                  R(i_more_y[j], Cart::z, i) + amb(1) * R(j, Cart::z, i);
              R(j, Cart::zz, i) =
                  R(i_more_z[j], Cart::z, i) + amb(2) * R(j, Cart::z, i);
            }
          }
          //------------------------------------------------------
        }

        if (lmax_beta > 2) {
          // Integrals    s - f - *    p - f - *    d - f - *    f - f - *    g
          // - f - *    h - f - *
          for (int i = 0; i < gamma_num_func; i++) {
            for (int j = 0; j < n_orbitals[lmax_alpha_beta - 3]; j++) {
              R(j, Cart::xxx, i) =
                  R(i_more_x[j], Cart::xx, i) + amb(0) * R(j, Cart::xx, i);
              R(j, Cart::xxy, i) =
                  R(i_more_x[j], Cart::xy, i) + amb(0) * R(j, Cart::xy, i);
              R(j, Cart::xxz, i) =
                  R(i_more_x[j], Cart::xz, i) + amb(0) * R(j, Cart::xz, i);
              R(j, Cart::xyy, i) =
                  R(i_more_x[j], Cart::yy, i) + amb(0) * R(j, Cart::yy, i);
              R(j, Cart::xyz, i) =
                  R(i_more_x[j], Cart::yz, i) + amb(0) * R(j, Cart::yz, i);
              R(j, Cart::xzz, i) =
                  R(i_more_x[j], Cart::zz, i) + amb(0) * R(j, Cart::zz, i);
              R(j, Cart::yyy, i) =
                  R(i_more_y[j], Cart::yy, i) + amb(1) * R(j, Cart::yy, i);
              R(j, Cart::yyz, i) =
                  R(i_more_y[j], Cart::yz, i) + amb(1) * R(j, Cart::yz, i);
              R(j, Cart::yzz, i) =
                  R(i_more_y[j], Cart::zz, i) + amb(1) * R(j, Cart::zz, i);
              R(j, Cart::zzz, i) =
                  R(i_more_z[j], Cart::zz, i) + amb(2) * R(j, Cart::zz, i);
            }
          }
          //------------------------------------------------------
        }

        if (lmax_beta > 3) {
          // Integrals    s - g - *    p - g - *    d - g - *    f - g - *    g
          // - g - *
          for (int i = 0; i < gamma_num_func; i++) {
            for (int j = 0; j < n_orbitals[lmax_alpha_beta - 4]; j++) {
              R(j, Cart::xxxx, i) =
                  R(i_more_x[j], Cart::xxx, i) + amb(0) * R(j, Cart::xxx, i);
              R(j, Cart::xxxy, i) =
                  R(i_more_x[j], Cart::xxy, i) + amb(0) * R(j, Cart::xxy, i);
              R(j, Cart::xxxz, i) =
                  R(i_more_x[j], Cart::xxz, i) + amb(0) * R(j, Cart::xxz, i);
              R(j, Cart::xxyy, i) =
                  R(i_more_x[j], Cart::xyy, i) + amb(0) * R(j, Cart::xyy, i);
              R(j, Cart::xxyz, i) =
                  R(i_more_x[j], Cart::xyz, i) + amb(0) * R(j, Cart::xyz, i);
              R(j, Cart::xxzz, i) =
                  R(i_more_x[j], Cart::xzz, i) + amb(0) * R(j, Cart::xzz, i);
              R(j, Cart::xyyy, i) =
                  R(i_more_x[j], Cart::yyy, i) + amb(0) * R(j, Cart::yyy, i);
              R(j, Cart::xyyz, i) =
                  R(i_more_x[j], Cart::yyz, i) + amb(0) * R(j, Cart::yyz, i);
              R(j, Cart::xyzz, i) =
                  R(i_more_x[j], Cart::yzz, i) + amb(0) * R(j, Cart::yzz, i);
              R(j, Cart::xzzz, i) =
                  R(i_more_x[j], Cart::zzz, i) + amb(0) * R(j, Cart::zzz, i);
              R(j, Cart::yyyy, i) =
                  R(i_more_y[j], Cart::yyy, i) + amb(1) * R(j, Cart::yyy, i);
              R(j, Cart::yyyz, i) =
                  R(i_more_y[j], Cart::yyz, i) + amb(1) * R(j, Cart::yyz, i);
              R(j, Cart::yyzz, i) =
                  R(i_more_y[j], Cart::yzz, i) + amb(1) * R(j, Cart::yzz, i);
              R(j, Cart::yzzz, i) =
                  R(i_more_y[j], Cart::zzz, i) + amb(1) * R(j, Cart::zzz, i);
              R(j, Cart::zzzz, i) =
                  R(i_more_z[j], Cart::zzz, i) + amb(2) * R(j, Cart::zzz, i);
            }
          }
          //------------------------------------------------------
        }

        // which ones do we want to store
        int cartoffset_alpha = shell_alpha->getCartesianOffset();
        int cartoffset_beta = shell_beta->getCartesianOffset();

        int cartnumFunc_alpha = shell_alpha->getCartesianNumFunc();
        int cartnumFunc_beta = shell_beta->getCartesianNumFunc();

        const Eigen::MatrixXd trafo_beta = AOTransform::getTrafo(gaussian_beta);
        const Eigen::MatrixXd trafo_alpha =
            AOTransform::getTrafo(gaussian_alpha);

        for (int i_alpha = 0; i_alpha < shell_alpha->getNumFunc(); i_alpha++) {
          for (int i_beta = 0; i_beta < shell_beta->getNumFunc(); i_beta++) {
            for (int i_gamma = 0; i_gamma < shell_gamma->getNumFunc();
                 i_gamma++) {
              for (int i_beta_t = 0; i_beta_t < cartnumFunc_beta; i_beta_t++) {
                for (int i_alpha_t = 0; i_alpha_t < cartnumFunc_alpha;
                     i_alpha_t++) {
                  double coeff = R(i_alpha_t + cartoffset_alpha,
                                   i_beta_t + cartoffset_beta, i_gamma) *
                                 trafo_alpha(i_alpha_t, i_alpha) *
                                 trafo_beta(i_beta_t, i_beta);
                  if (alphabetaswitch) {
                    threec_block(i_gamma, i_beta, i_alpha) += coeff;
                  } else {
                    threec_block(i_gamma, i_alpha, i_beta) += coeff;
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
}

}  // namespace xtp
}  // namespace votca
