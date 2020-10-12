/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#include <libint2/shell.h>
namespace votca {
namespace xtp {

AOGaussianPrimitive::AOGaussianPrimitive(const GaussianPrimitive& gaussian,
                                         const AOShell& aoshell)
    : _decay(gaussian.decay()),
      _contraction(gaussian.contraction()),
      _aoshell(aoshell) {
  _powfactor = CalcPowFactor(_decay);
}

AOGaussianPrimitive::AOGaussianPrimitive(const AOGaussianPrimitive& gaussian,
                                         const AOShell& aoshell)
    : _decay(gaussian._decay),
      _contraction(gaussian._contraction),
      _aoshell(aoshell),
      _powfactor(gaussian._powfactor) {
  ;
}

void AOGaussianPrimitive::SetupCptTable(CptTable& table) const {
  table.addCol(getShell().getAtomIndex(), "atomidx", HOFFSET(data, atomid));
  table.addCol(static_cast<Index>(getShell().getL()), "L", HOFFSET(data, l));
  table.addCol(getShell().getStartIndex(), "startidx",
               HOFFSET(data, startindex));
  table.addCol(getDecay(), "decay", HOFFSET(data, decay));
  table.addCol(getContraction(), "contr", HOFFSET(data, contraction));
  table.addCol(getShell().getPos().x(), "pos.x", HOFFSET(data, x));
  table.addCol(getShell().getPos().y(), "pos.y", HOFFSET(data, y));
  table.addCol(getShell().getPos().z(), "pos.z", HOFFSET(data, z));
  table.addCol(getShell().getScale(), "scale", HOFFSET(data, scale));
}

void AOGaussianPrimitive::WriteData(data& d) const {
  d.atomid = getShell().getAtomIndex();
  d.l = static_cast<Index>(getShell().getL());
  d.startindex = getShell().getStartIndex();
  d.decay = getDecay();
  d.contraction = getContraction();
  d.x = getShell().getPos().x();
  d.y = getShell().getPos().y();
  d.z = getShell().getPos().z();
  d.scale = getShell().getScale();
}

AOShell::AOShell(const Shell& shell, const QMAtom& atom, Index startIndex)
    : _l(shell.getL()),
      _scale(shell.getScale()),
      _startIndex(startIndex),
      _pos(atom.getPos()),
      _atomindex(atom.getId()) {
  ;
}

AOShell::AOShell(const AOShell& shell) {
  _l = shell._l;
  _scale = shell._scale;
  _mindecay = shell._mindecay;
  _startIndex = shell._startIndex;
  _pos = shell._pos;
  _atomindex = shell._atomindex;
  _gaussians.reserve(shell._gaussians.size());
  for (const auto& gaus : shell._gaussians) {
    _gaussians.push_back(AOGaussianPrimitive(gaus, *this));
  }
}

libint2::Shell AOShell::LibintShell() const {
  libint2::svector<libint2::Shell::real_t> decays;
  libint2::svector<libint2::Shell::Contraction> contractions;
  const Eigen::Vector3d& pos = getPos();
  libint2::Shell::Contraction contr;
  contr.l = static_cast<int>(getL());
  contr.pure = true;
  for (const auto& primitive : _gaussians) {
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
  for (auto& gaussian : _gaussians) {
    gaussian._contraction /= norm;
  }
  return;
}

void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                          Eigen::Block<Eigen::MatrixX3d>& gradAOvalues,
                          const Eigen::Vector3d& grid_pos) const {

  // need position of shell
  const Eigen::Vector3d center = (grid_pos - _pos);
  const double distsq = center.squaredNorm();

  // iterate over Gaussians in this shell
  for (const AOGaussianPrimitive& gaussian : _gaussians) {

    const double alpha = gaussian.getDecay();
    const double contraction = gaussian.getContraction();

    const double expofactor =
        gaussian.getPowfactor() * std::exp(-alpha * distsq);
    const Eigen::Vector3d second_term = -2.0 * alpha * center;

    switch (_l) {
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
        throw std::runtime_error("Shell type:" + EnumToString(_l) +
                                 " not known");
        break;
    }
  }  // contractions
  return;
}  // namespace xtp

void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                          const Eigen::Vector3d& grid_pos) const {

  Eigen::MatrixX3d temp = Eigen::MatrixX3d::Zero(AOvalues.size(), 3);
  Eigen::Block<Eigen::MatrixX3d> temp2 =
      temp.block(0, 0, temp.rows(), temp.cols());
  EvalAOspace(AOvalues, temp2, grid_pos);
}

std::ostream& operator<<(std::ostream& out, const AOShell& shell) {
  out << "AtomIndex:" << shell.getAtomIndex();
  out << " Shelltype:" << EnumToString(shell.getL())
      << " StartIndex:" << shell.getStartIndex()
      << " Scale:" << shell.getScale() << " MinDecay:" << shell.getMinDecay()
      << "\n";
  for (const auto& gaussian : shell) {
    out << " Gaussian Decay: " << gaussian.getDecay();
    out << " Contraction: " << gaussian.getContraction();
    out << "\n";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
