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
#ifndef VOTCA_XTP_SIGMA_BASE_UKS_H
#define VOTCA_XTP_SIGMA_BASE_UKS_H

#include "eigen.h"
#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

class Sigma_base_UKS {
 public:
  Sigma_base_UKS(TCMatrix_gwbse_spin& Mmn, const RPA_UKS& rpa,
                 TCMatrix::SpinChannel spin)
      : Mmn_spin_(Mmn), Mmn_(Mmn[spin]), rpa_(rpa), spin_(spin) {}

  virtual ~Sigma_base_UKS() = default;

  struct options {
    Index homo;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double eta;
    std::string quadrature_scheme;
    Index order;
    double alpha;
  };

  void configure(options opt) {
    opt_ = opt;
    qptotal_ = opt.qpmax - opt.qpmin + 1;
    rpatotal_ = opt.rpamax - opt.rpamin + 1;
  }

  Eigen::MatrixXd CalcExchangeMatrix() const;
  Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const;
  Eigen::MatrixXd CalcCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const;

  virtual void PrepareScreening() = 0;
  virtual double CalcCorrelationDiagElementDerivative(
      Index gw_level, double frequency) const = 0;
  virtual double CalcCorrelationDiagElement(Index gw_level,
                                            double frequency) const = 0;
  virtual double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                               double frequency1,
                                               double frequency2) const = 0;

 protected:
  const Eigen::VectorXd& getSpinRPAInputEnergies() const {
    return (spin_ == TCMatrix::SpinChannel::Alpha)
               ? rpa_.getRPAInputEnergiesAlpha()
               : rpa_.getRPAInputEnergiesBeta();
  }

  options opt_;
  TCMatrix_gwbse_spin& Mmn_spin_;
  TCMatrix_gwbse& Mmn_;
  const RPA_UKS& rpa_;
  TCMatrix::SpinChannel spin_;

  Index qptotal_ = 0;
  Index rpatotal_ = 0;
};

}  // namespace xtp
}  // namespace votca

#endif