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

namespace votca {
namespace xtp {

class AOMatrix {
 public:
  virtual void Fill(const AOBasis& aobasis) = 0;
  virtual Index Dimension() = 0;
  using MatrixLibInt =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
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
      _r[i] = r[i];
    }
  }  // definition of a center around which the moment should be calculated

 private:
  std::array<Eigen::MatrixXd, 3> _aomatrix;
  std::array<libint2::Shell::real_t, 3> _r = {0, 0, 0};
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
