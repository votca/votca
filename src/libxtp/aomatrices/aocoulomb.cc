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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aotransform.h"

namespace votca {
namespace xtp {

void AOCoulomb::Fill(const AOBasis& aobasis) {
  libint2::initialize();
  computeCoulombIntegrals(aobasis);
  libint2::finalize();
}

void AOCoulomb::computeCoulombIntegrals(const AOBasis& aobasis) {
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  const Index n = aobasis.AOBasisSize();
  const Index nshells = aobasis.getNumofShells();

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);

  // build engines for each thread
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb, aobasis.getMaxNprim(),
                      static_cast<int>(aobasis.getMaxL()), 0);
  engines[0].set(libint2::BraKet::xs_xs);
  for (Index i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = aobasis.getMapToBasisFunctions();
  auto unitshell = libint2::Shell::unit();

  auto compute = [&](Index thread_id) {
    const auto& buf = engines[thread_id].results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = shells[s1].size();

      for (auto s2 = 0; s2 <= s1; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = shells[s2].size();

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(shells[s1], shells[s2]);
        if (buf[0] == nullptr)
          continue;  // if all integrals screened out, skip to next shell set

        // "map" buffer to a const Eigen Matrix, and copy it to the
        // corresponding blocks of the result
        Eigen::Map<const MatrixLibInt> buf_mat(buf[0], n1, n2);
        result.block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                       // block, note the transpose!
          result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  };

  parallel_do(compute);

  _aomatrix = result;
}

// This converts V into ((S-1/2 V S-1/2)-1/2 S-1/2)T, which is needed to
// construct 4c integrals,
Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                                double etol) {

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(auxoverlap.Matrix());
  removedfunctions = 0;
  Eigen::VectorXd diagonal_overlap =
      Eigen::VectorXd::Zero(eo.eigenvalues().size());
  for (Index i = 0; i < diagonal_overlap.size(); ++i) {
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

  for (Index i = 0; i < diagonal.size(); ++i) {
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

}  // namespace xtp
}  // namespace votca
