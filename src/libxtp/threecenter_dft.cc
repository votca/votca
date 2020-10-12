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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/symmetric_matrix.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

void TCMatrix_dft::Fill(const AOBasis& auxbasis, const AOBasis& dftbasis) {

  AOCoulomb auxAOcoulomb;
  auxAOcoulomb.Fill(auxbasis);
  _inv_sqrt = auxAOcoulomb.Pseudo_InvSqrt(1e-8);
  _removedfunctions = auxAOcoulomb.Removedfunctions();

  _matrix = std::vector<Symmetric_Matrix>(
      auxbasis.AOBasisSize(), Symmetric_Matrix(dftbasis.AOBasisSize()));
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

#pragma omp parallel for schedule(dynamic)
  for (Index is = dftbasis.getNumofShells() - 1; is >= 0; is--) {

    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    const libint2::Shell& dftshell = dftshells[is];
    Index start = shell2bf[is];
    std::vector<Eigen::MatrixXd> block(dftshell.size());
    for (Index i = 0; i < Index(dftshell.size()); i++) {
      Index size = start + i + 1;
      block[i] = Eigen::MatrixXd::Zero(auxbasis.AOBasisSize(), size);
    }

    for (Index aux = 0; aux < auxbasis.getNumofShells(); aux++) {
      const libint2::Shell& auxshell = auxshells[aux];
      Index aux_start = auxshell2bf[aux];

      for (Index dis = 0; dis <= is; dis++) {

        const libint2::Shell& shell_col = dftshells[dis];
        Index col_start = shell2bf[dis];
        engine.compute(auxshell, dftshell, shell_col);

        if (buf[0] == nullptr) {
          continue;
        }
        Eigen::TensorMap<Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
            result(buf[0], auxshell.size(), dftshell.size(), shell_col.size());

        for (size_t left = 0; left < dftshell.size(); left++) {
          for (size_t auxf = 0; auxf < auxshell.size(); auxf++) {
            for (size_t col = 0; col < shell_col.size(); col++) {
              // symmetry
              if ((col_start + col) > (start + left)) {
                break;
              }
              block[left](aux_start + auxf, col_start + col) =
                  result(auxf, left, col);
            }
          }
        }
      }
    }

    for (Index i = 0; i < Index(block.size()); ++i) {
      Eigen::MatrixXd temp = _inv_sqrt * block[i];
      for (Index mu = 0; mu < temp.rows(); ++mu) {
        for (Index j = 0; j < temp.cols(); ++j) {
          _matrix[mu](i + start, j) = temp(mu, j);
        }
      }
    }
  }

  return;
}

}  // namespace xtp
}  // namespace votca
