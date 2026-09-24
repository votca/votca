
/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#include <Eigen/Core>
#include <complex>
#include <vector>
#include <stdexcept>
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {
namespace ewaldcontainer {

struct PointCharge {
  double charge = 0.0;
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
};

struct PointDipole {
  Eigen::Vector3d dipole = Eigen::Vector3d::Zero();
  Eigen::Vector3d position = Eigen::Vector3d::Zero();
};

struct ReciprocalTerm {
  Eigen::Vector3d G = Eigen::Vector3d::Zero();
  std::complex<double> coefficient = {0.0, 0.0};
};

class PotentialData {
 public:
  PotentialData() = default;

  explicit PotentialData(double eta) : eta_(eta) {}

  /* ---------- eta ---------- */

  void setEta(double eta) {
    if (eta <= 0.0) {
      throw std::runtime_error("Ewald eta must be positive");
    }
    // eta in Ewald is in nm^-1; QM need bohr^-1
    eta_ = eta * nm_inv_to_bohr_inv;
  }

  double eta() const { return eta_; }

  /* ---------- charges ---------- */

  void addCharge(double q, const Eigen::Vector3d& pos) {
    charges_.push_back({q, pos*tools::conv::nm2bohr});
  }

  const std::vector<PointCharge>& charges() const { return charges_; }
  std::vector<PointCharge>& charges() { return charges_; }

  size_t numCharges() const { return charges_.size(); }

  /* ---------- dipoles ---------- */

  void addDipole(const Eigen::Vector3d& mu, const Eigen::Vector3d& pos) {
    dipoles_.push_back({mu*tools::conv::nm2bohr, pos*tools::conv::nm2bohr});
  }

  const std::vector<PointDipole>& dipoles() const { return dipoles_; }
  std::vector<PointDipole>& dipoles() { return dipoles_; }

  size_t numDipoles() const { return dipoles_.size(); }

  /* ---------- reciprocal terms ---------- */

  void addReciprocalTerm(const Eigen::Vector3d& G,
                         const std::complex<double>& coeff) {
    reciprocal_terms_.push_back({G* nm_inv_to_bohr_inv, coeff * nm_inv_to_bohr_inv});
  }

  const std::vector<ReciprocalTerm>& reciprocalTerms() const {
    return reciprocal_terms_;
  }

  std::vector<ReciprocalTerm>& reciprocalTerms() {
    return reciprocal_terms_;
  }

  size_t numReciprocalTerms() const { return reciprocal_terms_.size(); }

  /* ---------- utilities ---------- */

  bool empty() const {
    return charges_.empty() &&
           dipoles_.empty() &&
           reciprocal_terms_.empty();
  }

  void clear() {
    charges_.clear();
    dipoles_.clear();
    reciprocal_terms_.clear();
    shape_factors_.setZero();
  }

  /* total charge useful for sanity checks */
  double totalCharge() const {
    double q = 0.0;
    for (const auto& c : charges_) {
      q += c.charge;
    }
    return q;
  }

  void setShapeFactors(Eigen::Vector4d vector){
    shape_factors_ = vector *tools::conv::bohr2nm ;
  } 

  Eigen::Vector4d& shapeFactors(){
    return shape_factors_;
  } 

  const Eigen::Vector4d& shapeFactors() const {
    return shape_factors_;
  } 

 private:
  double eta_ = 0.0;
  const double nm_inv_to_bohr_inv = tools::conv::bohr2nm;


  std::vector<PointCharge> charges_;
  std::vector<PointDipole> dipoles_;
  std::vector<ReciprocalTerm> reciprocal_terms_;
  Eigen::Vector4d shape_factors_ = Eigen::VectorXd::Zero(4);


};

}  // namespace ewaldcontainer
}  // namespace xtp
}  // namespace votca