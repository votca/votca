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
#ifndef VOTCA_XTP_AOMATRIX_H
#define VOTCA_XTP_AOMATRIX_H

// Local VOTCA includes
#include "aobasis.h"
#include <votca/xtp/ewald/ewaldcontainer.h>

namespace votca {
namespace xtp {

class AOMatrix {
 public:
  virtual void Fill(const AOBasis& aobasis) = 0;
  virtual Index Dimension() = 0;
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

 private:
  Eigen::MatrixXd aomatrix_;
};

// derived class for atomic orbital overlap
class AOOverlap : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  Eigen::MatrixXd singleShellOverlap(const AOShell& shell) const;
  Index Removedfunctions() const { return removedfunctions; }
  double SmallestEigenValue() const { return smallestEigenvalue; }

  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Eigen::MatrixXd Sqrt();

 private:
  Index removedfunctions;
  double smallestEigenvalue;
  Eigen::MatrixXd aomatrix_;
};

// derived class for atomic orbital Coulomb interaction
class AOCoulomb : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                       double etol);
  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Index Removedfunctions() const { return removedfunctions; }

 private:
  void computeCoulombIntegrals(const AOBasis& aobasis);
  Index removedfunctions;
  Eigen::MatrixXd aomatrix_;
};

/* derived class for atomic orbital electrical dipole matrices, required for
 * electrical transition dipoles
 */
class AODipole : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_[0].rows(); }
  const std::array<Eigen::MatrixXd, 3>& Matrix() const { return aomatrix_; }

  void setCenter(const Eigen::Vector3d& r) {
    for (Index i = 0; i < 3; i++) {
      r_[i] = r[i];
    }
  }  // definition of a center around which the moment should be calculated

 private:
  std::array<Eigen::MatrixXd, 3> aomatrix_;
  std::array<libint2::Shell::real_t, 3> r_ = {0, 0, 0};
};

/* derived class for atomic orbital overlap and electrical dipole matrices, required for
 * Ewald Shape Correction
 */
class AOEwaldShapeCorrection : public AOMatrix {
 public:
  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_[0].rows(); }
  const std::array<Eigen::MatrixXd, 4>& Matrix() const { return aomatrix_; }

  void setCenter(const Eigen::Vector3d& r) {
    for (Index i = 0; i < 3; i++) {
      r_[i] = r[i];
    }
  }  // definition of a center around which the moment should be calculated

 private:
  std::array<Eigen::MatrixXd, 4> aomatrix_;
  std::array<libint2::Shell::real_t, 3> r_ = {0, 0, 0};
};

// derived class for atomic orbital real-space Ewald potential from point charges
class AOEwaldRealSpaceCharges : public AOMatrix {
 public:

  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  void setEta(double eta) { eta_ = eta; }
  double getEta() const { return eta_; }

  void setCharges(const std::vector<ewaldcontainer::PointCharge>& charges) { charges_ = charges; }
  void addCharge(double charge, const Eigen::Vector3d& pos) {
    charges_.push_back({charge, pos});
  }
  const std::vector<ewaldcontainer::PointCharge>& Charges() const { return charges_; }

 private:
  double eta_ = 0.0;
  Eigen::MatrixXd aomatrix_;
  std::vector<ewaldcontainer::PointCharge> charges_;
};

// derived class for atomic orbital real-space Ewald Foreground Correction from point charges
class AOEwaldForegroundCharges : public AOMatrix {
 public:

  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  void setEta(double eta) { eta_ = eta; }
  double getEta() const { return eta_; }

  void setCharges(const std::vector<ewaldcontainer::PointCharge>& charges) { charges_ = charges; }
  void addCharge(double charge, const Eigen::Vector3d& pos) {
    charges_.push_back({charge, pos});
  }
  const std::vector<ewaldcontainer::PointCharge>& Charges() const { return charges_; }

 private:
  double eta_ = 0.0;
  Eigen::MatrixXd aomatrix_;
  std::vector<ewaldcontainer::PointCharge> charges_;
};

// derived class for atomic orbital real-space Ewald potential from point dipoles
class AOEwaldRealSpaceDipoles : public AOMatrix {
 public:
  struct PointDipole {
    Eigen::Vector3d dipole = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  };

  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  void setEta(double eta) { eta_ = eta; }
  void setDipoles(const std::vector<PointDipole>& dipoles) { dipoles_ = dipoles; }
  void addDipole(const Eigen::Vector3d& dipole, const Eigen::Vector3d& pos) {
    dipoles_.push_back({dipole, pos});
  }

 private:
  double eta_ = 0.0;
  Eigen::MatrixXd aomatrix_;
  std::vector<PointDipole> dipoles_;
};

/*
class AOEwaldRealSpaceDipoles : public AOMatrix {
 public:
  struct PointDipole {
    Eigen::Vector3d dipole = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  };

  void Fill(const AOBasis& aobasis) final;
  Index Dimension() final { return aomatrix_.rows(); }
  const Eigen::MatrixXd& Matrix() const { return aomatrix_; }

  void setEta(double eta) { eta_ = eta; }
  double getEta() const { return eta_; }

  // central finite-difference step for d/dR of erfc(eta r)/r
  void setDelta(double delta) { delta_ = delta; }
  double getDelta() const { return delta_; }

  void setDipoles(const std::vector<PointDipole>& dipoles) { dipoles_ = dipoles; }
  void addDipole(const Eigen::Vector3d& dipole, const Eigen::Vector3d& pos) {
    dipoles_.push_back({dipole, pos});
  }
  const std::vector<PointDipole>& Dipoles() const { return dipoles_; }

 private:
  double eta_ = 0.0;
  double delta_ = 1e-5;
  Eigen::MatrixXd aomatrix_;
  std::vector<PointDipole> dipoles_;
};*/

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
