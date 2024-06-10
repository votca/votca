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
#ifndef VOTCA_XTP_EEINTERACTOR_H
#define VOTCA_XTP_EEINTERACTOR_H

// Local VOTCA includes
#include "classicalsegment.h"
#include "eigen.h"

namespace votca {
namespace xtp {

enum Estatic : bool {
  noE_V = true,
  V = false,
};

/**
 * \brief Mediates interaction between polar and static sites
 */
class eeInteractor {
 public:
  explicit eeInteractor() = default;
  explicit eeInteractor(double expdamping) : expdamping_(expdamping) {};

  Eigen::Matrix3d FillTholeInteraction(const PolarSite& site1,
                                       const PolarSite& site2) const;

  Eigen::VectorXd Cholesky_IntraSegment(const PolarSegment& seg) const;

  template <class T, enum Estatic>
  double ApplyStaticField(const T& segment1, PolarSegment& segment2) const;

  template <enum Estatic>
  double ApplyInducedField(const PolarSegment& segment1,
                           PolarSegment& segment2) const;

  template <class S1, class S2>
  double CalcStaticEnergy(const S1& segment1, const S2& segment2) const;

  class E_terms {
   public:
    E_terms& operator+=(const E_terms& right) {
      this->data_ += right.data_;
      return *this;
    }

    E_terms operator+(E_terms right) const {
      right.data_ += this->data_;
      return right;
    }

    double sum() { return data_.sum(); }

    double& E_indu_indu() { return data_.x(); }
    double& E_indu_stat() { return data_.y(); }
    double& E_internal() { return data_.z(); }

    const Eigen::Vector3d& data() const { return data_; }

   private:
    Eigen::Vector3d data_ = Eigen::Vector3d::Zero();
  };

  template <class S1, class S2>
  E_terms CalcPolarEnergy(const S1& segment1, const S2& segment2) const;

  template <class S>
  double CalcStaticEnergy_IntraSegment(const S& seg) const;

  double CalcPolarEnergy_IntraSegment(const PolarSegment& seg) const;

  double CalcStaticEnergy_site(const StaticSite& site1,
                               const StaticSite& site2) const;

 private:
  template <int N>
  Eigen::Matrix<double, N, 1> VSiteA(const StaticSite& site1,
                                     const StaticSite& site2) const;
  template <enum Estatic>
  double ApplyInducedField_site(const PolarSite& site1, PolarSite& site2) const;
  template <enum Estatic>
  double ApplyStaticField_site(const StaticSite& site1, PolarSite& site2) const;

  double CalcPolar_stat_Energy_site(const PolarSite& site1,
                                    const StaticSite& site2) const;

  E_terms CalcPolarEnergy_site(const PolarSite& site1,
                               const StaticSite& site2) const;

  E_terms CalcPolarEnergy_site(const PolarSite& site1,
                               const PolarSite& site2) const;

  double expdamping_ = 0.39;  // dimensionless
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EEINTERACTOR_H
