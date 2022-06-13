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
#ifndef VOTCA_XTP_ERDIABATIZATION_H
#define VOTCA_XTP_ERDIABATIZATION_H

#include "votca/xtp/ERIs.h"
//#include "votca/xtp/logger.h"
#include "logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmtool.h"
#include <cstdio>
#include <votca/tools/property.h>
#include <votca/tools/types.h>
#include <votca/xtp/aobasis.h>
//#include <votca/xtp/aomatrix.h>

namespace votca {
namespace xtp {

class ERDiabatization {
 public:
  ERDiabatization();

  ERDiabatization(Orbitals& orbitals1, Orbitals& orbitals2, Logger* log)
      : orbitals1_(orbitals1), orbitals2_(orbitals2), pLog_(log){};

  // Write a function to set up all the matrices we need
  void setUpMatrices();

  // Option for ERdiabatization
  struct options_erdiabatization {

    Index state_idx_1;
    Index state_idx_2;
    std::string qmtype;
    std::string outfile;
  };

  void configure(const options_erdiabatization& opt);

  double Calculate_angle() const;
  Eigen::MatrixXd Calculate_diabatic_H(const double angle) const;

 private:
  Orbitals& orbitals1_;
  Orbitals& orbitals2_;
  QMStateType qmtype_;
  Logger* pLog_;

  ERIs eris_;
  AOBasis dftbasis_;
  AOBasis auxbasis_;

  Index bse_cmax_;
  Index bse_cmin_;
  Index bse_vmax_;
  Index bse_vmin_;
  Index bse_vtotal_;
  Index bse_ctotal_;
  Index basissize_;

  Eigen::MatrixXd occlevels1_;
  Eigen::MatrixXd virtlevels1_;

  Eigen::MatrixXd occlevels2_;
  Eigen::MatrixXd virtlevels2_;

  bool useRI_;
  options_erdiabatization opt_;

  double E1_;
  double E2_;

  template <bool AR>
  Eigen::MatrixXd CalculateD(Index stateindex1, Index stateindex2) const;

  Eigen::MatrixXd CalculateD_R(Index stateindex1, Index stateindex2) const {
    return CalculateD<false>(stateindex1, stateindex2);
  }
  Eigen::MatrixXd CalculateD_AR(Index stateindex1, Index stateindex2) const {
    return CalculateD<true>(stateindex1, stateindex2);
  }

  double CalculateR(const Eigen::MatrixXd& D_JK,
                    const Eigen::MatrixXd& D_LM) const;
  Eigen::MatrixXd CalculateU(const double phi) const;
  Eigen::Tensor<double, 4> CalculateRtensor() const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERDIABATIZATION_H
