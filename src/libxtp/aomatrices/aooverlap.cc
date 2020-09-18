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
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aotransform.h"

namespace votca {
namespace xtp {

void AOOverlap::Fill(const AOBasis& aobasis) {
  libint2::initialize();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  libint2::Operator obtype = libint2::Operator::overlap;
  libint2::Engine engine(obtype, aobasis.getMaxNprim(), aobasis.getMaxL(), 0);

  _aomatrix.resize(aobasis.AOBasisSize(), aobasis.AOBasisSize());

  auto shell2bf = aobasis.getMapToBasisFunctions();

  const auto& buf = engine.results();

  for (auto s1 = 0; s1 != shells.size(); ++s1) {

    auto bf1 = shell2bf[s1];  // first basis function in this shell
    auto n1 = shells[s1].size();

    for (auto s2 = 0; s2 <= s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair
      engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
      _aomatrix.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                     // block, note the transpose!
        _aomatrix.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }
  libint2::finalize();
}

// HERE SO OLD STUFF STILL COMPILES, WILL BE DELETED SOON
Eigen::MatrixXd AOOverlap::Primitive_Overlap(const AOGaussianPrimitive& g_row,
                                             const AOGaussianPrimitive& g_col,
                                             Index l_offset) const {
  Eigen::MatrixXd ol = Eigen::MatrixXd(5, 5);
  return ol;
}

// HERE SO OLD STUFF STILL COMPILES, WILL BE DELETED SOON
void AOOverlap::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                          const AOShell& shell_row,
                          const AOShell& shell_col) const {

  ;
}

// Needed to compute the normalization of the contractions
Eigen::MatrixXd AOOverlap::FillShell(const AOShell& shell) const {

  libint2::initialize();
  libint2::Operator obtype = libint2::Operator::overlap;
  libint2::Engine engine(obtype, shell.getNumFunc() + 3,
                         static_cast<Index>(shell.getL()), 0);

  const auto& buf = engine.results();

  engine.compute(toLibintShell(shell), toLibintShell(shell));

  Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], shell.getSize(),
                                            shell.getSize());

  libint2::finalize();
  return buf_mat;
}

// This is ugly here I will probably write it in the aoshell class, can then
// also be used in the aobasis
libint2::Shell AOOverlap::toLibintShell(const AOShell& shell) const {
  libint2::svector<libint2::Shell::real_t> decays;
  libint2::svector<libint2::Shell::Contraction> contractions;
  const Eigen::Vector3d& pos = shell.getPos();
  libint2::Shell::Contraction contr;
  contr.l = static_cast<int>(shell.getL());
  contr.pure = true;
  for (const auto& primitive : shell) {
    decays.push_back(primitive.getDecay());
    contr.coeff.push_back(primitive.getContraction());
  }
  contractions.push_back(contr);
  std::array<libint2::Shell::real_t, 3> libintpos = {pos[0], pos[1], pos[2]};
  libint2::Shell libintshell(decays, contractions, libintpos);
  return libintshell;
}

// Still need this bit
Eigen::MatrixXd AOOverlap::Pseudo_InvSqrt(double etol) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.eigenvalues().size());
  removedfunctions = 0;
  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.eigenvalues()(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.eigenvalues()(i));
    }
  }

  return es.eigenvectors() * diagonal.asDiagonal() *
         es.eigenvectors().transpose();
}

// And this as well
Eigen::MatrixXd AOOverlap::Sqrt() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  return es.operatorSqrt();
}

}  // namespace xtp
}  // namespace votca
