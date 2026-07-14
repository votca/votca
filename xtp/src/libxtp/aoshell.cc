/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/aoshell.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/checkpointtable.h"

namespace votca {
namespace xtp {

AOGaussianPrimitive::AOGaussianPrimitive(const GaussianPrimitive& gaussian)
    : decay_(gaussian.decay()), contraction_(gaussian.contraction()) {
  powfactor_ = CalcPowFactor(decay_);
}

void AOGaussianPrimitive::SetupCptTable(CptTable& table) {
  table.addCol<Index>("atomidx", HOFFSET(data, atomid));
  table.addCol<Index>("L", HOFFSET(data, l));
  table.addCol<Index>("startidx", HOFFSET(data, startindex));
  table.addCol<double>("decay", HOFFSET(data, decay));
  table.addCol<double>("contr", HOFFSET(data, contraction));
  table.addCol<double>("pos.x", HOFFSET(data, x));
  table.addCol<double>("pos.y", HOFFSET(data, y));
  table.addCol<double>("pos.z", HOFFSET(data, z));
  table.addCol<double>("scale", HOFFSET(data, scale));
}

void AOGaussianPrimitive::WriteData(data& d, const AOShell& s) const {
  d.atomid = s.getAtomIndex();
  d.l = static_cast<Index>(s.getL());
  d.startindex = s.getStartIndex();
  d.decay = getDecay();
  d.contraction = getContraction();
  d.x = s.getPos().x();
  d.y = s.getPos().y();
  d.z = s.getPos().z();
}

AOShell::AOShell(const Shell& shell, const QMAtom& atom, Index startIndex)
    : l_(shell.getL()),
      startIndex_(startIndex),
      pos_(atom.getPos()),
      atomindex_(atom.getId()) {
  ;
}

libint2::Shell AOShell::LibintShell() const {
  libint2::svector<libint2::Shell::real_t> decays;
  libint2::svector<libint2::Shell::Contraction> contractions;
  const Eigen::Vector3d& pos = getPos();
  libint2::Shell::Contraction contr;
  contr.l = static_cast<int>(getL());
  contr.pure = true;
  for (const auto& primitive : gaussians_) {
    decays.push_back(primitive.getDecay());
    contr.coeff.push_back(primitive.getContraction());
  }
  contractions.push_back(contr);
  std::array<libint2::Shell::real_t, 3> libintpos = {pos[0], pos[1], pos[2]};
  return libint2::Shell(decays, contractions, libintpos);
}

void AOShell::normalizeContraction() {
  AOOverlap overlap;
  Eigen::MatrixXd block = overlap.singleShellOverlap(*this);
  double norm = std::sqrt(block(0, 0));
  for (auto& gaussian : gaussians_) {
    gaussian.contraction_ /= norm;
  }
  return;
}

AOShell::AOValues AOShell::EvalAOspace(const Eigen::Vector3d& grid_pos) const {

  // need position of shell
  const Eigen::Vector3d center = (grid_pos - pos_);
  const double distsq = center.squaredNorm();
  AOShell::AOValues AO(getNumFunc());
  Eigen::VectorXd& AOvalues = AO.values;
  Eigen::MatrixX3d& gradAOvalues = AO.derivatives;

  // iterate over Gaussians in this shell
  for (const AOGaussianPrimitive& gaussian : gaussians_) {

    const double alpha = gaussian.getDecay();
    const double contraction = gaussian.getContraction();

    const double expofactor =
        gaussian.getPowfactor() * std::exp(-alpha * distsq);
    const Eigen::Vector3d second_term = -2.0 * alpha * center;

    switch (l_) {
      case L::S: {
        double AOvalue = contraction * expofactor;
        AOvalues(0) += AOvalue;                        // s-function
        gradAOvalues.row(0) += second_term * AOvalue;  // gradient of s-function
      } break;
      case L::P: {
        const double factor = 2. * sqrt(alpha) * contraction * expofactor;

        double AOvalue = factor * center.y();  // Y 1,-1
        AOvalues(0) += AOvalue;
        gradAOvalues.row(0) += second_term * AOvalue;
        gradAOvalues(0, 1) += factor;

        AOvalue = factor * center.z();  // Y 1,0
        AOvalues(1) += AOvalue;
        gradAOvalues.row(1) += second_term * AOvalue;
        gradAOvalues(1, 2) += factor;

        AOvalue = factor * center.x();  // Y 1,1
        AOvalues(2) += AOvalue;
        gradAOvalues(2, 0) += factor;
        gradAOvalues.row(2) += second_term * AOvalue;  // y gradient
      } break;
      case L::D: {
        const double factor = 2. * alpha * contraction * expofactor;
        const double factor_1 = factor / sqrt(3.);

        double AOvalue = 2. * factor * (center.x() * center.y());  // Y 2,-2
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {2 * center.y(), 2 * center.x(), 0};
        gradAOvalues.row(0) += factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = 2. * factor * (center.y() * center.z());  // Y 2,-1
        AOvalues(1) += AOvalue;
        coeff = {0, 2 * center.z(), 2 * center.y()};
        gradAOvalues.row(1) += factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_1 * (3. * center.z() * center.z() - distsq);  // Y 2,0
        AOvalues(2) += AOvalue;
        coeff = {-2, -2, 4};
        gradAOvalues.row(2) += (factor_1 * coeff * center.array()).matrix() +
                               second_term * AOvalue;

        AOvalue = 2. * factor * (center.x() * center.z());  // Y 2,1
        AOvalues(3) += AOvalue;
        coeff = {2 * center.z(), 0, 2 * center.x()};
        gradAOvalues.row(3) += factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor *
                  (center.x() * center.x() - center.y() * center.y());  // Y 2,2
        AOvalues(4) += AOvalue;
        coeff = {2 * center.x(), -2 * center.y(), 0};
        gradAOvalues.row(4) += factor * coeff.matrix() + second_term * AOvalue;
      } break;
      case L::F: {
        const double factor = 2. * pow(alpha, 1.5) * contraction * expofactor;
        const double factor_1 = factor * 2. / sqrt(15.);
        const double factor_2 = factor * sqrt(2.) / sqrt(5.);
        const double factor_3 = factor * sqrt(2.) / sqrt(3.);
        AxA c(center);

        double AOvalue =
            factor_3 * center.y() * (3. * c.xx() - c.yy());  // Y 3,-3
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {6. * c.xy(), 3. * (c.xx() - c.yy()), 0};
        gradAOvalues.row(0) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;

        AOvalue = 4. * factor * center.x() * center.y() * center.z();  // Y 3,-2
        AOvalues(1) += AOvalue;
        coeff = {c.yz(), c.xz(), c.xy()};
        gradAOvalues.row(1) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_2 * center.y() * (5. * c.zz() - distsq);  // Y 3,-1
        AOvalues(2) += AOvalue;
        coeff = {-2. * c.xy(), 4. * c.zz() - c.xx() - 3. * c.yy(), 8. * c.yz()};
        gradAOvalues.row(2) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_1 * center.z() * (5. * c.zz() - 3. * distsq);  // Y 3,0
        AOvalues(3) += AOvalue;
        coeff = {-6. * c.xz(), -6. * c.yz(), 3. * (3. * c.zz() - distsq)};
        gradAOvalues.row(3) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_2 * center.x() * (5. * c.zz() - distsq);  // Y 3,1
        AOvalues(4) += AOvalue;
        coeff = {4. * c.zz() - c.yy() - 3. * c.xx(), -2. * c.xy(), 8. * c.xz()};
        gradAOvalues.row(4) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        AOvalue = 2. * factor * center.z() * (c.xx() - c.yy());  // Y 3,2
        AOvalues(5) += AOvalue;
        coeff = {2. * c.xz(), -2. * c.yz(), c.xx() - c.yy()};
        gradAOvalues.row(5) +=
            2 * factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_3 * center.x() * (c.xx() - 3. * c.yy());  // Y 3,3
        AOvalues(6) += AOvalue;
        coeff = {3. * (c.xx() - c.yy()), -6. * c.xy(), 0};
        gradAOvalues.row(6) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;
      } break;
      case L::G: {
        const double factor =
            2. / sqrt(3.) * alpha * alpha * contraction * expofactor;
        const double factor_1 = factor / sqrt(35.);
        const double factor_2 = factor * 4. / sqrt(14.);
        const double factor_3 = factor * 2. / sqrt(7.);
        const double factor_4 = factor * 2. * sqrt(2.);
        AxA c(center);

        double AOvalue = 4. * factor * c.xy() * (c.xx() - c.yy());  // Y 4,-4
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {center.y() * (3. * c.xx() - c.yy()),
                                center.x() * (c.xx() - 3. * c.yy()), 0};
        gradAOvalues.row(0) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_4 * c.yz() * (3. * c.xx() - c.yy());  // Y 4,-3
        AOvalues(1) += AOvalue;
        coeff = {6. * center.x() * c.yz(), 3. * center.z() * (c.xx() - c.yy()),
                 center.y() * (3. * c.xx() - c.yy())};
        gradAOvalues.row(1) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;

        AOvalue = 2. * factor_3 * c.xy() * (7. * c.zz() - distsq);  // Y 4,-2
        AOvalues(2) += AOvalue;
        coeff = {center.y() * (6. * c.zz() - 3. * c.xx() - c.yy()),
                 center.x() * (6. * c.zz() - c.xx() - 3. * c.yy()),
                 12. * center.z() * c.xy()};
        gradAOvalues.row(2) +=
            2 * factor_3 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_2 * c.yz() * (7. * c.zz() - 3. * distsq);  // Y 4,-1
        AOvalues(3) += AOvalue;
        coeff = {(-6. * center.x() * c.yz()),
                 center.z() * (4. * c.zz() - 3. * c.xx() - 9. * c.yy()),
                 3. * center.y() * (5. * c.zz() - distsq)};
        gradAOvalues.row(3) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_1 * (35. * c.zz() * c.zz() - 30. * c.zz() * distsq +
                              3. * distsq * distsq);  // Y 4,0
        AOvalues(4) += AOvalue;
        coeff = {12. * center.x() * (distsq - 5. * c.zz()),
                 12. * center.y() * (distsq - 5. * c.zz()),
                 16. * center.z() * (5. * c.zz() - 3. * distsq)};

        gradAOvalues.row(4) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_2 * c.xz() * (7. * c.zz() - 3. * distsq);  // Y 4,1
        AOvalues(5) += AOvalue;
        coeff = {center.z() * (4. * c.zz() - 9. * c.xx() - 3. * c.yy()),
                 (-6. * center.y() * c.xz()),
                 3. * center.x() * (5. * c.zz() - distsq)};
        gradAOvalues.row(5) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        AOvalue =
            factor_3 * (c.xx() - c.yy()) * (7. * c.zz() - distsq);  // Y 4,2
        AOvalues(6) += AOvalue;
        coeff = {4. * center.x() * (3. * c.zz() - c.xx()),
                 4. * center.y() * (c.yy() - 3. * c.zz()),
                 12. * center.z() * (c.xx() - c.yy())};
        gradAOvalues.row(6) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor_4 * c.xz() * (c.xx() - 3. * c.yy());  // Y 4,3
        AOvalues(7) += AOvalue;
        coeff = {3. * center.z() * (c.xx() - c.yy()),
                 (-6. * center.y() * c.xz()),
                 center.x() * (c.xx() - 3. * c.yy())};
        gradAOvalues.row(7) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;

        AOvalue = factor * (c.xx() * c.xx() - 6. * c.xx() * c.yy() +
                            c.yy() * c.yy());  // Y 4,4
        AOvalues(8) += AOvalue;
        coeff = {center.x() * (c.xx() - 3. * c.yy()),
                 center.y() * (c.yy() - 3. * c.xx()), 0};
        gradAOvalues.row(8) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;
      } break;
      default:
        throw std::runtime_error("Shell type:" + EnumToString(l_) +
                                 " not known");
        break;
    }
  }  // contractions
  return AO;
}

AOShell::AOValuesHessian AOShell::EvalAOspaceHessian(
    const Eigen::Vector3d& grid_pos) const {

  const Eigen::Vector3d center = (grid_pos - pos_);
  const double distsq = center.squaredNorm();
  AOShell::AOValuesHessian AO(getNumFunc());
  Eigen::VectorXd& AOvalues = AO.values;
  Eigen::MatrixX3d& gradAOvalues = AO.derivatives;
  std::vector<Eigen::Matrix3d>& hessians = AO.hessians;

  // Accumulates the Hessian contribution from ONE Gaussian primitive to
  // function index k, given that primitive's own prefactor, alpha, and
  // the shell-function's polynomial P (via its value P_val, gradient
  // dP_vec, and Hessian d2P_mat at this point) -- NOT the running,
  // multi-primitive-accumulated AOvalue/gradient, since a contracted
  // basis function sums primitives with DIFFERENT alpha, and this
  // formula's alpha-dependent terms must use each primitive's own
  // alpha, not some already-summed value. See the STATUS comment on
  // AOValuesHessian in aoshell.h for the formula and its derivation.
  auto addHessianContribution =
      [&](Index k, double prefactor, double P_val,
          const Eigen::Vector3d& dP_vec, const Eigen::Matrix3d& d2P_mat,
          double alpha) {
        double AOvalue_local = prefactor * P_val;
        Eigen::Vector3d second_term_local = -2.0 * alpha * center;
        Eigen::Vector3d grad_local =
            prefactor * dP_vec + second_term_local * AOvalue_local;
        Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
        for (Index i = 0; i < 3; ++i) {
          for (Index j = 0; j < 3; ++j) {
            double delta_ij = (i == j) ? 1.0 : 0.0;
            H(i, j) = prefactor * d2P_mat(i, j) -
                      2.0 * alpha * delta_ij * AOvalue_local -
                      2.0 * alpha * center(i) * grad_local(j) -
                      2.0 * alpha * center(j) * grad_local(i) -
                      4.0 * alpha * alpha * center(i) * center(j) *
                          AOvalue_local;
          }
        }
        hessians[k] += H;
      };

  for (const AOGaussianPrimitive& gaussian : gaussians_) {

    const double alpha = gaussian.getDecay();
    const double contraction = gaussian.getContraction();

    const double expofactor =
        gaussian.getPowfactor() * std::exp(-alpha * distsq);
    const Eigen::Vector3d second_term = -2.0 * alpha * center;

    switch (l_) {
      case L::S: {
        double AOvalue = contraction * expofactor;
        AOvalues(0) += AOvalue;
        gradAOvalues.row(0) += second_term * AOvalue;
        // P = 1, dP = 0, d2P = 0
        addHessianContribution(0, contraction * expofactor, 1.0,
                               Eigen::Vector3d::Zero(),
                               Eigen::Matrix3d::Zero(), alpha);
      } break;
      case L::P: {
        const double factor = 2. * sqrt(alpha) * contraction * expofactor;
        Eigen::Matrix3d zero3 = Eigen::Matrix3d::Zero();

        double AOvalue = factor * center.y();  // Y 1,-1
        AOvalues(0) += AOvalue;
        gradAOvalues.row(0) += second_term * AOvalue;
        gradAOvalues(0, 1) += factor;
        addHessianContribution(0, factor, center.y(),
                               Eigen::Vector3d(0, 1, 0), zero3, alpha);

        AOvalue = factor * center.z();  // Y 1,0
        AOvalues(1) += AOvalue;
        gradAOvalues.row(1) += second_term * AOvalue;
        gradAOvalues(1, 2) += factor;
        addHessianContribution(1, factor, center.z(),
                               Eigen::Vector3d(0, 0, 1), zero3, alpha);

        AOvalue = factor * center.x();  // Y 1,1
        AOvalues(2) += AOvalue;
        gradAOvalues(2, 0) += factor;
        gradAOvalues.row(2) += second_term * AOvalue;
        addHessianContribution(2, factor, center.x(),
                               Eigen::Vector3d(1, 0, 0), zero3, alpha);
      } break;
      case L::D: {
        const double factor = 2. * alpha * contraction * expofactor;
        const double factor_1 = factor / sqrt(3.);
        Eigen::Matrix3d d2P;

        double AOvalue = 2. * factor * (center.x() * center.y());  // Y 2,-2
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {2 * center.y(), 2 * center.x(), 0};
        gradAOvalues.row(0) += factor * coeff.matrix() + second_term * AOvalue;
        d2P << 0, 2, 0,  2, 0, 0,  0, 0, 0;
        addHessianContribution(0, factor, 2. * center.x() * center.y(),
                               coeff.matrix(), d2P, alpha);

        AOvalue = 2. * factor * (center.y() * center.z());  // Y 2,-1
        AOvalues(1) += AOvalue;
        coeff = {0, 2 * center.z(), 2 * center.y()};
        gradAOvalues.row(1) += factor * coeff.matrix() + second_term * AOvalue;
        d2P << 0, 0, 0,  0, 0, 2,  0, 2, 0;
        addHessianContribution(1, factor, 2. * center.y() * center.z(),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_1 * (3. * center.z() * center.z() - distsq);  // Y 2,0
        AOvalues(2) += AOvalue;
        coeff = {-2, -2, 4};
        gradAOvalues.row(2) += (factor_1 * coeff * center.array()).matrix() +
                               second_term * AOvalue;
        d2P << -2, 0, 0,  0, -2, 0,  0, 0, 4;
        addHessianContribution(
            2, factor_1, 3. * center.z() * center.z() - distsq,
            (coeff * center.array()).matrix(), d2P, alpha);

        AOvalue = 2. * factor * (center.x() * center.z());  // Y 2,1
        AOvalues(3) += AOvalue;
        coeff = {2 * center.z(), 0, 2 * center.x()};
        gradAOvalues.row(3) += factor * coeff.matrix() + second_term * AOvalue;
        d2P << 0, 0, 2,  0, 0, 0,  2, 0, 0;
        addHessianContribution(3, factor, 2. * center.x() * center.z(),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor *
                  (center.x() * center.x() - center.y() * center.y());  // Y 2,2
        AOvalues(4) += AOvalue;
        coeff = {2 * center.x(), -2 * center.y(), 0};
        gradAOvalues.row(4) += factor * coeff.matrix() + second_term * AOvalue;
        d2P << 2, 0, 0,  0, -2, 0,  0, 0, 0;
        addHessianContribution(
            4, factor, center.x() * center.x() - center.y() * center.y(),
            coeff.matrix(), d2P, alpha);
      } break;
      case L::F: {
        const double factor = 2. * pow(alpha, 1.5) * contraction * expofactor;
        const double factor_1 = factor * 2. / sqrt(15.);
        const double factor_2 = factor * sqrt(2.) / sqrt(5.);
        const double factor_3 = factor * sqrt(2.) / sqrt(3.);
        AxA c(center);
        Eigen::Matrix3d d2P;
        double x = center.x(), y = center.y(), z = center.z();

        double AOvalue =
            factor_3 * center.y() * (3. * c.xx() - c.yy());  // Y 3,-3
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {6. * c.xy(), 3. * (c.xx() - c.yy()), 0};
        gradAOvalues.row(0) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;
        d2P << 6*y, 6*x, 0,  6*x, -6*y, 0,  0, 0, 0;
        addHessianContribution(0, factor_3, y * (3. * c.xx() - c.yy()),
                               coeff.matrix(), d2P, alpha);

        AOvalue = 4. * factor * center.x() * center.y() * center.z();  // Y 3,-2
        AOvalues(1) += AOvalue;
        coeff = {c.yz(), c.xz(), c.xy()};
        gradAOvalues.row(1) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;
        d2P << 0, z, y,  z, 0, x,  y, x, 0;
        addHessianContribution(1, 4. * factor, x * y * z, coeff.matrix(),
                               d2P, alpha);

        AOvalue = factor_2 * center.y() * (5. * c.zz() - distsq);  // Y 3,-1
        AOvalues(2) += AOvalue;
        coeff = {-2. * c.xy(), 4. * c.zz() - c.xx() - 3. * c.yy(), 8. * c.yz()};
        gradAOvalues.row(2) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;
        d2P << -2*y, -2*x, 0,  -2*x, -6*y, 8*z,  0, 8*z, 8*y;
        addHessianContribution(2, factor_2, y * (5. * c.zz() - distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_1 * center.z() * (5. * c.zz() - 3. * distsq);  // Y 3,0
        AOvalues(3) += AOvalue;
        coeff = {-6. * c.xz(), -6. * c.yz(), 3. * (3. * c.zz() - distsq)};
        gradAOvalues.row(3) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;
        d2P << -6*z, 0, -6*x,  0, -6*z, -6*y,  -6*x, -6*y, 12*z;
        addHessianContribution(3, factor_1,
                               z * (5. * c.zz() - 3. * distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_2 * center.x() * (5. * c.zz() - distsq);  // Y 3,1
        AOvalues(4) += AOvalue;
        coeff = {4. * c.zz() - c.yy() - 3. * c.xx(), -2. * c.xy(), 8. * c.xz()};
        gradAOvalues.row(4) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;
        d2P << -6*x, -2*y, 8*z,  -2*y, -2*x, 0,  8*z, 0, 8*x;
        addHessianContribution(4, factor_2, x * (5. * c.zz() - distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue = 2. * factor * center.z() * (c.xx() - c.yy());  // Y 3,2
        AOvalues(5) += AOvalue;
        coeff = {2. * c.xz(), -2. * c.yz(), c.xx() - c.yy()};
        gradAOvalues.row(5) +=
            2 * factor * coeff.matrix() + second_term * AOvalue;
        d2P << 2*z, 0, 2*x,  0, -2*z, -2*y,  2*x, -2*y, 0;
        addHessianContribution(5, 2. * factor, z * (c.xx() - c.yy()),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_3 * center.x() * (c.xx() - 3. * c.yy());  // Y 3,3
        AOvalues(6) += AOvalue;
        coeff = {3. * (c.xx() - c.yy()), -6. * c.xy(), 0};
        gradAOvalues.row(6) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;
        d2P << 6*x, -6*y, 0,  -6*y, -6*x, 0,  0, 0, 0;
        addHessianContribution(6, factor_3, x * (c.xx() - 3. * c.yy()),
                               coeff.matrix(), d2P, alpha);
      } break;
      case L::G: {
        const double factor =
            2. / sqrt(3.) * alpha * alpha * contraction * expofactor;
        const double factor_1 = factor / sqrt(35.);
        const double factor_2 = factor * 4. / sqrt(14.);
        const double factor_3 = factor * 2. / sqrt(7.);
        const double factor_4 = factor * 2. * sqrt(2.);
        AxA c(center);
        Eigen::Matrix3d d2P;
        double x = center.x(), y = center.y(), z = center.z();

        double AOvalue = 4. * factor * c.xy() * (c.xx() - c.yy());  // Y 4,-4
        AOvalues(0) += AOvalue;
        Eigen::Array3d coeff = {center.y() * (3. * c.xx() - c.yy()),
                                center.x() * (c.xx() - 3. * c.yy()), 0};
        gradAOvalues.row(0) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;
        d2P << 6*x*y, 3*x*x-3*y*y, 0,  3*x*x-3*y*y, -6*x*y, 0,  0, 0, 0;
        addHessianContribution(0, 4. * factor, x * y * (c.xx() - c.yy()),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_4 * c.yz() * (3. * c.xx() - c.yy());  // Y 4,-3
        AOvalues(1) += AOvalue;
        coeff = {6. * center.x() * c.yz(), 3. * center.z() * (c.xx() - c.yy()),
                 center.y() * (3. * c.xx() - c.yy())};
        gradAOvalues.row(1) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;
        d2P << 6*y*z, 6*x*z, 6*x*y,  6*x*z, -6*y*z, 3*x*x-3*y*y,
               6*x*y, 3*x*x-3*y*y, 0;
        addHessianContribution(1, factor_4, y * z * (3. * c.xx() - c.yy()),
                               coeff.matrix(), d2P, alpha);

        AOvalue = 2. * factor_3 * c.xy() * (7. * c.zz() - distsq);  // Y 4,-2
        AOvalues(2) += AOvalue;
        coeff = {center.y() * (6. * c.zz() - 3. * c.xx() - c.yy()),
                 center.x() * (6. * c.zz() - c.xx() - 3. * c.yy()),
                 12. * center.z() * c.xy()};
        gradAOvalues.row(2) +=
            2 * factor_3 * coeff.matrix() + second_term * AOvalue;
        d2P << -6*x*y, -3*x*x-3*y*y+6*z*z, 12*y*z,
               -3*x*x-3*y*y+6*z*z, -6*x*y, 12*x*z,
               12*y*z, 12*x*z, 12*x*y;
        addHessianContribution(2, 2. * factor_3,
                               x * y * (7. * c.zz() - distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_2 * c.yz() * (7. * c.zz() - 3. * distsq);  // Y 4,-1
        AOvalues(3) += AOvalue;
        coeff = {(-6. * center.x() * c.yz()),
                 center.z() * (4. * c.zz() - 3. * c.xx() - 9. * c.yy()),
                 3. * center.y() * (5. * c.zz() - distsq)};
        gradAOvalues.row(3) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;
        d2P << -6*y*z, -6*x*z, -6*x*y,
               -6*x*z, -18*y*z, -3*x*x-9*y*y+12*z*z,
               -6*x*y, -3*x*x-9*y*y+12*z*z, 24*y*z;
        addHessianContribution(3, factor_2,
                               y * z * (7. * c.zz() - 3. * distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor_1 * (35. * c.zz() * c.zz() - 30. * c.zz() * distsq +
                              3. * distsq * distsq);  // Y 4,0
        AOvalues(4) += AOvalue;
        coeff = {12. * center.x() * (distsq - 5. * c.zz()),
                 12. * center.y() * (distsq - 5. * c.zz()),
                 16. * center.z() * (5. * c.zz() - 3. * distsq)};
        gradAOvalues.row(4) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;
        d2P << 36*x*x+12*y*y-48*z*z, 24*x*y, -96*x*z,
               24*x*y, 12*x*x+36*y*y-48*z*z, -96*y*z,
               -96*x*z, -96*y*z, -48*x*x-48*y*y+96*z*z;
        addHessianContribution(
            4, factor_1,
            35. * c.zz() * c.zz() - 30. * c.zz() * distsq +
                3. * distsq * distsq,
            coeff.matrix(), d2P, alpha);

        AOvalue = factor_2 * c.xz() * (7. * c.zz() - 3. * distsq);  // Y 4,1
        AOvalues(5) += AOvalue;
        coeff = {center.z() * (4. * c.zz() - 9. * c.xx() - 3. * c.yy()),
                 (-6. * center.y() * c.xz()),
                 3. * center.x() * (5. * c.zz() - distsq)};
        gradAOvalues.row(5) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;
        d2P << -18*x*z, -6*y*z, -9*x*x-3*y*y+12*z*z,
               -6*y*z, -6*x*z, -6*x*y,
               -9*x*x-3*y*y+12*z*z, -6*x*y, 24*x*z;
        addHessianContribution(5, factor_2,
                               x * z * (7. * c.zz() - 3. * distsq),
                               coeff.matrix(), d2P, alpha);

        AOvalue =
            factor_3 * (c.xx() - c.yy()) * (7. * c.zz() - distsq);  // Y 4,2
        AOvalues(6) += AOvalue;
        coeff = {4. * center.x() * (3. * c.zz() - c.xx()),
                 4. * center.y() * (c.yy() - 3. * c.zz()),
                 12. * center.z() * (c.xx() - c.yy())};
        gradAOvalues.row(6) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;
        d2P << -12*x*x+12*z*z, 0, 24*x*z,
               0, 12*y*y-12*z*z, -24*y*z,
               24*x*z, -24*y*z, 12*x*x-12*y*y;
        addHessianContribution(
            6, factor_3, (c.xx() - c.yy()) * (7. * c.zz() - distsq),
            coeff.matrix(), d2P, alpha);

        AOvalue = factor_4 * c.xz() * (c.xx() - 3. * c.yy());  // Y 4,3
        AOvalues(7) += AOvalue;
        coeff = {3. * center.z() * (c.xx() - c.yy()),
                 (-6. * center.y() * c.xz()),
                 center.x() * (c.xx() - 3. * c.yy())};
        gradAOvalues.row(7) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;
        d2P << 6*x*z, -6*y*z, 3*x*x-3*y*y,
               -6*y*z, -6*x*z, -6*x*y,
               3*x*x-3*y*y, -6*x*y, 0;
        addHessianContribution(7, factor_4, x * z * (c.xx() - 3. * c.yy()),
                               coeff.matrix(), d2P, alpha);

        AOvalue = factor * (c.xx() * c.xx() - 6. * c.xx() * c.yy() +
                            c.yy() * c.yy());  // Y 4,4
        AOvalues(8) += AOvalue;
        coeff = {center.x() * (c.xx() - 3. * c.yy()),
                 center.y() * (c.yy() - 3. * c.xx()), 0};
        gradAOvalues.row(8) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;
        d2P << 12*x*x-12*y*y, -24*x*y, 0,  -24*x*y, -12*x*x+12*y*y, 0,
               0, 0, 0;
        // NOTE: dP_vec scaled by 4 here specifically, NOT prefactor --
        // this function's own gradient update uses coeff at a "quarter
        // scale" (4*factor*coeff.matrix(), not factor*coeff.matrix()
        // like every other function in this shell), confirmed by
        // cross-checking every AOvalue/gradient multiplier pair in this
        // file against each other (the only mismatch found across all
        // of D/F/G). d2P itself is unaffected (computed directly from
        // the exact polynomial, independent of coeff's scaling
        // convention) -- only dP_vec needs the correction, so that
        // addHessianContribution's internal grad_local reconstruction
        // (prefactor*dP_vec + second_term*AOvalue_local) matches the
        // real gradient (factor*(4*coeff) + second_term*AOvalue =
        // 4*factor*coeff + second_term*AOvalue, exactly what
        // gradAOvalues.row(8) above actually accumulates).
        addHessianContribution(
            8, factor,
            c.xx() * c.xx() - 6. * c.xx() * c.yy() + c.yy() * c.yy(),
            4.0 * coeff.matrix(), d2P, alpha);
      } break;
      default:
        throw std::runtime_error("Shell type:" + EnumToString(l_) +
                                 " not known (Hessian evaluation)");
        break;
    }
  }  // contractions
  return AO;
}

std::ostream& operator<<(std::ostream& out, const AOShell& shell) {
  out << "AtomIndex:" << shell.getAtomIndex();
  out << " Shelltype:" << EnumToString(shell.getL())
      << " StartIndex:" << shell.getStartIndex()
      << " MinDecay:" << shell.getMinDecay() << "\n";
  for (const auto& gaussian : shell) {
    out << " Gaussian Decay: " << gaussian.getDecay();
    out << " Contraction: " << gaussian.getContraction();
    out << "\n";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
