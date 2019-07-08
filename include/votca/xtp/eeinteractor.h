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
#ifndef VOTCA_XTP_EEINTERACTOR_H
#define VOTCA_XTP_EEINTERACTOR_H

#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/eigen.h>
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
  explicit eeInteractor(){};
  explicit eeInteractor(double expdamping) : _expdamping(expdamping){};

  Eigen::Matrix3d FillTholeInteraction(const PolarSite& site1,
                                       const PolarSite& site2) const;

  Eigen::Vector3d VThole(const PolarSite& site1, const PolarSite& site2,
                         const Eigen::Vector3d& dQ) const;

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
      this->_data += right._data;
      return *this;
    }

    E_terms operator+(E_terms right) const {
      right._data += this->_data;
      return right;
    }

    double sum() { return _data.sum(); }

    double& E_indu_indu() { return _data.x(); }
    double& E_indu_stat() { return _data.y(); }
    double& E_internal() { return _data.z(); }

    const Eigen::Vector3d& data() const { return _data; }

   private:
    Eigen::Vector3d _data = Eigen::Vector3d::Zero();
  };

  template <class S1, class S2>
  E_terms CalcPolarEnergy(const S1& segment1, const S2& segment2) const;

  template <class S>
  double CalcStaticEnergy_IntraSegment(const S& seg) const;

  double CalcPolarEnergy_IntraSegment(const PolarSegment& seg) const;

 private:
  class AxA {
   public:
    AxA(const Eigen::Vector3d& a) {
      _data.segment<3>(0) = a.x() * a;
      _data.segment<2>(3) = a.y() * a.segment<2>(1);
      _data[5] = a.z() * a.z();
    }
    inline const double& xx() const { return _data[0]; }
    inline const double& xy() const { return _data[1]; }
    inline const double& xz() const { return _data[2]; }
    inline const double& yy() const { return _data[3]; }
    inline const double& yz() const { return _data[4]; }
    inline const double& zz() const { return _data[5]; }

   private:
    Eigen::Matrix<double, 6, 1> _data;
  };

  template <int N>
  Eigen::Matrix<double, N, 1> VSiteA(const StaticSite& site1,
                                     const StaticSite& site2) const;
  template <enum Estatic>
  double ApplyInducedField_site(const PolarSite& site1, PolarSite& site2) const;
  template <enum Estatic>
  double ApplyStaticField_site(const StaticSite& site1, PolarSite& site2) const;

  double CalcStaticEnergy_site(const StaticSite& site1,
                               const StaticSite& site2) const;

  double CalcPolar_stat_Energy_site(const PolarSite& site1,
                                    const StaticSite& site2) const;

  E_terms CalcPolarEnergy_site(const PolarSite& site1,
                               const StaticSite& site2) const;

  E_terms CalcPolarEnergy_site(const PolarSite& site1,
                               const PolarSite& site2) const;

  double _expdamping = 0.39;  // dimensionless
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EEINTERACTOR_H
