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
  Index Dimension() { return _aomatrix.rows(); }
  const Eigen::MatrixXd& Matrix() const { return _aomatrix; }
  void Fill(const AOBasis& aobasis);

 protected:
 template <typename Lambda>
  void parallel_do(Lambda& lambda);

  template <libint2::Operator obtype>
  void computeOneBodyIntegrals(const AOBasis& basis);

  std::unordered_map<Index, std::vector<Index>> compute_shellpairs(
      const AOBasis& bs1, const double threshold = 1e-20);

  virtual void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                         const AOShell& shell_row,
                         const AOShell& shell_col) const = 0;
  Eigen::MatrixXd _aomatrix;
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix {
  public:
  void Fill(const AOBasis& aobasis) {
    libint2::initialize();
    computeOneBodyIntegrals<libint2::Operator::kinetic>(aobasis);
    libint2::finalize();
  }
 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;
};

// derived class for atomic orbital overlap
class AOOverlap : public AOMatrix {
 public:
  Eigen::MatrixXd FillShell(const AOShell& shell) const;
  Index Removedfunctions() const { return removedfunctions; }
  double SmallestEigenValue() const { return smallestEigenvalue; }
  void Fill(const AOBasis& aobasis);

  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Eigen::MatrixXd Sqrt();

  Eigen::MatrixXd Primitive_Overlap(const AOGaussianPrimitive& g_row,
                                    const AOGaussianPrimitive& g_col,
                                    Index l_offset = 0) const;

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;

 private:
  Index removedfunctions;
  double smallestEigenvalue;
};

// derived class for atomic orbital Coulomb interaction
class AOCoulomb : public AOMatrix {
 public:
  Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                       double etol);
  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Index Removedfunctions() const { return removedfunctions; }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;

 private:
  Index removedfunctions;
};


/**********************************
 * 
 * TEMPLATE DEFINITIONS
 * 
 * *******************************/

template <typename Lambda>
void AOMatrix::parallel_do(Lambda& lambda) {
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
}

template <libint2::Operator obtype>
void AOMatrix::computeOneBodyIntegrals(const AOBasis& aobasis) {

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  const Index n = aobasis.AOBasisSize();
  const Index nshells = aobasis.getNumofShells();

  auto shellpair_list = compute_shellpairs(aobasis);

  Index nopers = static_cast<Index>(libint2::operator_traits<obtype>::nopers);
  std::array<Eigen::MatrixXd, libint2::operator_traits<obtype>::nopers> result;
  for (auto& r : result) {
    r = Eigen::MatrixXd::Zero(n, n);
  }

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(obtype, aobasis.getMaxNprim(),
                               static_cast<int>(aobasis.getMaxL()), 0);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = aobasis.getMapToBasisFunctions();

  auto compute = [&](int thread_id) {
    const auto& buf = engines[thread_id].results();

    for (auto s1 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = shells[s1].size();

      auto s1_offset = s1 * (s1 + 1) / 2;
      for (auto s2 : shellpair_list[s1]) {
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = shells[s2].size();

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(shells[s1], shells[s2]);

        for (unsigned int op = 0; op != nopers; ++op) {
          // "map" buffer to a const Eigen Matrix, and copy it to the
          // corresponding blocks of the result
          Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> buf_mat(buf[op], n1, n2);
          result[op].block(bf1, bf2, n1, n2) = buf_mat;
          if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding
                         // {s2,s1} block, note the transpose!
            result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
      }
    }
  };
  parallel_do(compute);

  _aomatrix = result[0];
}


}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
