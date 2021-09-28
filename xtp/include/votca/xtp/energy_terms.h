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
#ifndef VOTCA_XTP_ENERGY_TERMS_H
#define VOTCA_XTP_ENERGY_TERMS_H

// Local VOTCA includes
#include "eeinteractor.h"
#include "eigen.h"

/**
 * \brief Small container for the individual energy terms in a polar region
 *
 *
 */

namespace votca {
namespace xtp {

class Energy_terms {
 public:
  Energy_terms& operator+=(const Energy_terms& right) {
    this->data_ += right.data_;
    return *this;
  }

  Energy_terms operator+(Energy_terms right) const {
    right.data_ += this->data_;
    return right;
  }

  Energy_terms operator-(Energy_terms right) const {
    right.data_ *= (-1);
    right.data_ += this->data_;
    return right;
  }

  void addInternalPolarContrib(const eeInteractor::E_terms& induction_terms) {
    data_.segment<3>(0) += induction_terms.data();
  }

  double Etotal() const { return data_.sum(); }  // total energy
  double Epolar() const {
    return data_.segment<4>(0).sum();
  }  // all polar inside region and from outside contributions
     // dQ-dQ,Q-dQ,E_internal
  double Estatic() const {
    return data_.segment<2>(4).sum();
  }  // all static contributions Q-Q inside region and from outside
  double Eextern() const {
    return data_.segment<2>(3).sum();
  }  // all external contributions
  double Eintern() const {
    return Etotal() - Eextern();
  }  // all internal contributions

  double& E_indu_indu() { return data_[0]; }  // dQ-dQ inside region
  double& E_indu_stat() { return data_[1]; }  // dQ-Q inside region
  double& E_internal() { return data_[2]; }   // e_internal
  double& E_polar_ext() {
    return data_[3];
  }  // dQ-Q and dQ-dQ from outside regions
  double& E_static_ext() { return data_[4]; }     // Q-Q from outside regions
  double& E_static_static() { return data_[5]; }  // Q-Q inside region

  const Eigen::Matrix<double, 6, 1>& data() const { return data_; }
  Eigen::Matrix<double, 6, 1>& data() { return data_; }

 private:
  Eigen::Matrix<double, 6, 1> data_ = Eigen::Matrix<double, 6, 1>::Zero();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ENERGY_TERMS_H
