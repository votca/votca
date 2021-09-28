/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#ifndef VOTCA_XTP_AOPOTENTIAL_H
#define VOTCA_XTP_AOPOTENTIAL_H

// Local VOTCA includes
#include "aobasis.h"
#include "ecpaobasis.h"
#include "staticsite.h"

namespace votca {
namespace xtp {

class QMMolecule;

// base class for 1D atomic orbital matrix types (overlap, Coulomb, ESP)
template <class T>
class AOPotential {
 public:
  Index Dimension() { return _aopotential.rows(); }
  const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Matrix() const {
    return _aopotential;
  }

 protected:
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Fill(
      const AOBasis& aobasis) const;
  virtual void FillBlock(
      Eigen::Block<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrix,
      const AOShell& shell_row, const AOShell& shell_col) const = 0;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _aopotential;
};

// derived class for Effective Core Potentials
class AOECP : public AOPotential<double> {
 public:
  void FillPotential(const AOBasis& aobasis, const ECPAOBasis& ecp);

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>&, const AOShell&,
                 const AOShell&) const final{};

 private:
};

class AOMultipole : public AOPotential<double> {
 public:
  void FillPotential(const AOBasis& aobasis, const QMMolecule& atoms);
  void FillPotential(const AOBasis& aobasis, const Eigen::Vector3d& r);
  void FillPotential(
      const AOBasis& aobasis,
      const std::vector<std::unique_ptr<StaticSite>>& externalsites);

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;

 private:
  void setSite(const StaticSite* site) { _site = site; };

  const StaticSite* _site;
};

class AOPlanewave : public AOPotential<std::complex<double>> {
 public:
  void FillPotential(const AOBasis& aobasis,
                     const std::vector<Eigen::Vector3d>& kpoints);

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXcd>& matrix,
                 const AOShell& shell_row,
                 const AOShell& shell_col) const override;

 private:
  void setkVector(const Eigen::Vector3d& k) { _k = k; };
  Eigen::Vector3d _k;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOPOTENTIAL_H
