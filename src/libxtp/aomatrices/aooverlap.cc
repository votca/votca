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
  _aomatrix = computeOneBodyIntegrals<libint2::Operator::overlap>(aobasis)[0];
  libint2::finalize();
}

Eigen::MatrixXd AOOverlap::singleShellOverlap(const AOShell& shell) const {
  libint2::Shell::do_enforce_unit_normalization(false);
  libint2::initialize();
  libint2::Operator obtype = libint2::Operator::overlap;
  libint2::Engine engine(obtype, shell.getSize(),
                         static_cast<int>(shell.getL()), 0);

  const libint2::Engine::target_ptr_vec& buf = engine.results();

  libint2::Shell s = shell.LibintShell();
  engine.compute(s, s);

  Eigen::Map<const MatrixLibInt> buf_mat(buf[0], shell.getNumFunc(),
                                         shell.getNumFunc());

  libint2::finalize();
  return buf_mat;
}

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

Eigen::MatrixXd AOOverlap::Sqrt() {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
  smallestEigenvalue = es.eigenvalues()(0);
  return es.operatorSqrt();
}

}  // namespace xtp
}  // namespace votca
