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

#pragma once
#ifndef VOTCA_XTP_AOMATRIX_H
#define VOTCA_XTP_AOMATRIX_H

#include <votca/xtp/aobasis.h>

namespace votca {
namespace xtp {

class AOMatrix {
 public:
  int Dimension() { return _aomatrix.rows(); }
  const Eigen::MatrixXd& Matrix() const { return _aomatrix; }
  void Fill(const AOBasis& aobasis);

 protected:
  virtual void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                         const AOShell& shell_row,
                         const AOShell& shell_col) const = 0;
  Eigen::MatrixXd _aomatrix;
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix {
 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;
};

// derived class for atomic orbital overlap
class AOOverlap : public AOMatrix {
 public:
  Eigen::MatrixXd FillShell(const AOShell& shell) const;
  int Removedfunctions() const { return removedfunctions; }
  double SmallestEigenValue() const { return smallestEigenvalue; }

  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Eigen::MatrixXd Sqrt();

  Eigen::MatrixXd Primitive_Overlap(const AOGaussianPrimitive& g_row,
                                    const AOGaussianPrimitive& g_col,
                                    int l_offset = 0) const;

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  int removedfunctions;
  double smallestEigenvalue;
};

// derived class for atomic orbital Coulomb interaction
class AOCoulomb : public AOMatrix {
 public:
  Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                       double etol);
  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  int Removedfunctions() const { return removedfunctions; }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  int removedfunctions;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
