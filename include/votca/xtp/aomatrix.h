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

#pragma once
#ifndef VOTCA_XTP_AOMATRIX_H
#define VOTCA_XTP_AOMATRIX_H

// Local VOTCA includes
#include "aobasis.h"
#include <unordered_map>

#include <libint2.hpp>

namespace votca {
namespace xtp {

class AOMatrix {
 public:
  virtual void Fill(const AOBasis& aobasis) = 0;
  virtual Index Dimension() = 0;

 protected:
  // libint uses rowmajor storage for the computed integrals
  using MatrixLibInt =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  template <libint2::Operator obtype,
            typename OperatorParams =
                typename libint2::operator_traits<obtype>::oper_params_type>
  std::array<MatrixLibInt, libint2::operator_traits<obtype>::nopers>
      computeOneBodyIntegrals(const AOBasis& basis,
                              OperatorParams oparams = OperatorParams());
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final {
    libint2::initialize();
    _aomatrix = computeOneBodyIntegrals<libint2::Operator::kinetic>(aobasis)[0];
    libint2::finalize();
  }
  Index Dimension() final { return _aomatrix.rows(); }
  const Eigen::MatrixXd& Matrix() const { return _aomatrix; }

 private:
  Eigen::MatrixXd _aomatrix;
};

// derived class for atomic orbital overlap
class AOOverlap : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return _aomatrix.rows(); }
  const Eigen::MatrixXd& Matrix() const { return _aomatrix; }

  Eigen::MatrixXd singleShellOverlap(const AOShell& shell) const;
  Index Removedfunctions() const { return removedfunctions; }
  double SmallestEigenValue() const { return smallestEigenvalue; }

  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Eigen::MatrixXd Sqrt();

 private:
  Index removedfunctions;
  double smallestEigenvalue;
  Eigen::MatrixXd _aomatrix;
};

// derived class for atomic orbital Coulomb interaction
class AOCoulomb : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return _aomatrix.rows(); }
  const Eigen::MatrixXd& Matrix() const { return _aomatrix; }

  Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                       double etol);
  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Index Removedfunctions() const { return removedfunctions; }

 private:
  void computeCoulombIntegrals(const AOBasis& aobasis);
  Index removedfunctions;
  Eigen::MatrixXd _aomatrix;
};

/* derived class for atomic orbital electrical dipole matrices, required for
 * electrical transition dipoles
 */
class AODipole : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return _aomatrix[0].rows(); }
  const std::array<Eigen::MatrixXd, 3>& Matrix() const { return _aomatrix; }

  void setCenter(const Eigen::Vector3d& r) {
    for (Index i = 0; i < 3; i++) {
      _r[0] = r[0];
    }
  }  // definition of a center around which the moment should be calculated

 private:
  std::array<Eigen::MatrixXd, 3> _aomatrix;
  std::array<libint2::Shell::real_t, 3> _r = {0, 0, 0};
};

/**********************************
 *
 * TEMPLATE DEFINITIONS
 *
 * *******************************/
template <libint2::Operator obtype, typename OperatorParams>
std::array<AOMatrix::MatrixLibInt, libint2::operator_traits<obtype>::nopers>
    AOMatrix::computeOneBodyIntegrals(const AOBasis& aobasis,
                                      OperatorParams oparams) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();

  std::vector<std::vector<Index>> shellpair_list = aobasis.ComputeShellPairs();

  Index nopers = static_cast<Index>(libint2::operator_traits<obtype>::nopers);
  std::array<MatrixLibInt, libint2::operator_traits<obtype>::nopers> result;
  for (MatrixLibInt& r : result) {
    r = MatrixLibInt::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  }

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(obtype, aobasis.getMaxNprim(),
                               static_cast<int>(aobasis.getMaxL()), 0);
  engines[0].set_params(oparams);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0l; s1 < aobasis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();

    for (Index s2 : shellpair_list[s1]) {
      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();

      engine.compute(shells[s1], shells[s2]);

      for (unsigned int op = 0; op != nopers; ++op) {
        Eigen::Map<const MatrixLibInt> buf_mat(buf[op], n1, n2);
        result[op].block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding
                       // {s2,s1} block, note the transpose!
          result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  }

  return result;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
