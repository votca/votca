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
  explicit eeInteractor(double expdamping = 0.39) : _expdamping(expdamping){};

  template <class T1, class T2>
  double InteractStatic(T1& seg1, T2& seg2) const;

  template <class T>
  double InteractStatic_IntraSegment(T& seg) const;

  double InteractPolar_IntraSegment(PolarSegment& seg1) const;

  double InteractPolar(PolarSegment& seg1, PolarSegment& seg2) const;

  double InteractPolar(const PolarSegment& seg1, PolarSegment& seg2) const;

  Eigen::Matrix3d FillTholeInteraction_diponly(const PolarSite& site1,
                                               const PolarSite& site2) const;

  Eigen::VectorXd Cholesky_IntraSegment(const PolarSegment& seg) const;

 private:
  double InteractStatic_site(StaticSite& site1, StaticSite& site2) const;
  double InteractStatic_site(const StaticSite& site1, StaticSite& site2) const;

  double InteractPolar_site(PolarSite& site1, PolarSite& site2) const;
  void InteractPolar_site_small(PolarSite& site1, PolarSite& site2) const;
  double InteractPolar_site(const PolarSite& site1, PolarSite& site2) const;
  Matrix9d FillInteraction(const StaticSite& site1,
                           const StaticSite& site2) const;
  Matrix9d FillTholeInteraction(const PolarSite& site1,
                                const PolarSite& site2) const;

  Eigen::Matrix4d FillInteraction_noQuadrupoles(const StaticSite& site1,
                                                const StaticSite& site2) const;

  Eigen::Matrix4d FillTholeInteraction_noQuadrupoles(
      const PolarSite& site1, const PolarSite& site2) const;

  double _expdamping = 0.39;  // dimensionless
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EEINTERACTOR_H
