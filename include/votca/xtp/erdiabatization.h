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

namespace votca {
namespace xtp {

class ERDiabatization {
 public:
  ERDiabatization();

  ERDiabatization(Orbitals& orbitals, Logger* log)
      : _orbitals(orbitals), _pLog(log){};

  // Write a function to set up all the matrices we need
  void setUpMatrices();

  // Option for ERdiabatization
  struct options_erdiabatization {

    Index state_idx_1;
    Index state_idx_2;
    std::string qmtype;
    std::string outfile;
  };

  // Edit options
  void configure(const options_erdiabatization& opt);
  // Set up options

Eigen::VectorXd CalculateER(const Orbitals& orb, QMStateType type) const;
 Eigen::MatrixXd Calculate_diabatic_H(const double E1, const double E2, const double angle) const;
 private:
  Logger* _pLog;

  Orbitals& _orbitals;
  ERIs _eris;
  AOBasis _dftbasis;
  AOBasis _auxbasis;

  std::string _type;
  Index _nostates;
  Index _bse_cmax;
  Index _bse_cmin;
  Index _bse_vmax;
  Index _bse_vmin;
  Index _bse_vtotal;
  Index _bse_ctotal;
  Index _basis;
  Index _bse_size_ao;
  Eigen::VectorXd _occlevels;
  Eigen::VectorXd _virtlevels;

  options_erdiabatization _opt;
  Eigen::MatrixXd CalculateD(const Orbitals& orb, QMStateType type,
                             Index stateindex1, Index stateindex2) const;
  
  double CalculateR(const Eigen::MatrixXd& D_JK,
                    const Eigen::MatrixXd& D_LM) const;
  Eigen::MatrixXd CalculateU(const double phi) const;
  Eigen::Tensor<double, 4> CalculateRtensor(const Orbitals& orb,
                                            QMStateType type) const;
  
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ERDIABATIZATION_H
