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
#ifndef _VOTCA_XTP_SIGMA_PPM_H
#define _VOTCA_XTP_SIGMA_PPM_H
#include <votca/xtp/ppm.h>
#include <votca/xtp/sigma_base.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_PPM : public Sigma_base {
 public:
  Sigma_PPM(TCMatrix_gwbse& Mmn, RPA& rpa) : Sigma_base(Mmn, rpa){};

  // Sets up the screening parametrisation
  void PrepareScreening() override;
  // Calculates Sigma_c diagonal elements
  virtual double CalcCorrelation(Index gw_level,
                                 double frequency) const override;
  // Calculates Sigma_c off-diagonal elements
  virtual double CalcCorrelation(Index gw_level1, Index gw_level2,
                                 double frequency1,
                                 double frequency2) const override;

 private:
  inline void Stabilize(Eigen::ArrayXd& denom) const;
  PPM _ppm;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_SIGMA_H */
