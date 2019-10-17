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
#ifndef VOTCA_XTP_AOMATRIX3D_H
#define VOTCA_XTP_AOMATRIX3D_H

#include <votca/xtp/aobasis.h>

namespace votca {
namespace xtp {

class AOMatrix3D {
 public:
  const std::array<Eigen::MatrixXd, 3>& Matrix() const { return _aomatrix; }
  void Fill(const AOBasis& aobasis);
  int Dimension() { return _aomatrix[0].rows(); }

 protected:
  std::array<Eigen::MatrixXd, 3> _aomatrix;
  virtual void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                         const AOShell& shell_row,
                         const AOShell& shell_col) const = 0;
};

/* derived class for atomic orbital gradient matrices, required for
 * momentum transition dipoles
 */
class AOMomentum : public AOMatrix3D {
 protected:
  void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;
};

/* derived class for atomic orbital electrical dipole matrices, required for
 * electrical transition dipoles
 */
class AODipole : public AOMatrix3D {
 public:
  void setCenter(const Eigen::Vector3d& r) {
    _r = r;
  }  // definition of a center around which the moment should be calculated
 protected:
  void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;

 private:
  Eigen::Vector3d _r = Eigen::Vector3d::Zero();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX3D_H
