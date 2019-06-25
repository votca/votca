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

/**
 * \brief Mediates interaction between polar and static sites
 */
class eeInteractor {
 public:
  explicit eeInteractor(){};
  explicit eeInteractor(double expdamping) : _expdamping(expdamping){};

  Eigen::Matrix3d FillTholeInteraction(const PolarSite& site1,
                                       const PolarSite& site2) const;

  Eigen::VectorXd Cholesky_IntraSegment(const PolarSegment& seg) const;

  template <class T>
  void ApplyStaticField(const T& segment1, PolarSegment& segment2) const;
  void ApplyStaticField_IntraSegment(PolarSegment& seg) const;

  void ApplyInducedField(const PolarSegment& segment1,
                         PolarSegment& segment2) const;

  template <class S1, class S2>
  double CalcStaticEnergy(const S1& segment1, const S2& segment2) const;

  template <class S1, class S2>
  double CalcPolarEnergy(const S1& segment1, const S2& segment2) const;

  template <class S>
  double CalcStaticEnergy_IntraSegment(const S& seg) const;

  double CalcPolarEnergy_IntraSegment(const PolarSegment& seg) const;

 private:
  template <int N, int M>
  Eigen::Matrix<double, N, M> FillInteraction(const StaticSite& site1,
                                              const StaticSite& site2) const;

  void ApplyInducedField_site(const PolarSite& site1, PolarSite& site2) const;
  void ApplyStaticField_site(const StaticSite& site1, PolarSite& site2) const;

  double CalcStaticEnergy_site(const StaticSite& site1,
                               const StaticSite& site2) const;

  double CalcPolar_stat_Energy_site(const PolarSite& site1,
                                    const StaticSite& site2) const;

  double CalcPolarEnergy_site(const PolarSite& site1,
                              const StaticSite& site2) const;

  double CalcPolarEnergy_site(const PolarSite& site1,
                              const PolarSite& site2) const;

  double _expdamping = 0.39;  // dimensionless
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EEINTERACTOR_H
