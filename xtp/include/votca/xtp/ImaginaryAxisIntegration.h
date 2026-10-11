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

#ifndef VOTCA_XTP_IMAGINARYAXISINTEGRATION_H
#define VOTCA_XTP_IMAGINARYAXISINTEGRATION_H

#include "eigen.h"
#include "quadrature_factory.h"
#include "rpa.h"
#include "rpa_uks.h"
#include <memory>

// Computes the contribution from the Gauss-Laguerre quadrature to the
// self-energy expectation matrix for given RPA and frequencies
namespace votca {
namespace xtp {

class ImaginaryAxisIntegration {

 public:
  struct options {
    Index order;
    Index qptotal;
    Index qpmin;
    Index homo;
    Index rpamin;
    Index rpamax;
    std::string quadrature_scheme;
    double alpha;
  };

  ImaginaryAxisIntegration(const Eigen::VectorXd& energies,
                           const TCMatrix_gwbse& Mmn);

  void configure(options opt, const RPA& rpa,
                 const Eigen::MatrixXd& kDielMxInv_zero);
  void configure(options opt, const RPA_UKS& rpa,
                 const Eigen::MatrixXd& kDielMxInv_zero);

  double SigmaGQDiag(double frequency, Index gw_level, double eta) const;

  /// q(i) = M.row(i) K M.row(i)^T for all rows of M.
  static Eigen::VectorXd RowQuadraticForms(const Eigen::MatrixXd& M,
                                           const Eigen::MatrixXd& K);

 private:
  options opt_;

  std::unique_ptr<GaussianQuadratureBase> gq_ = nullptr;

  // This function calculates and stores inverses of the microscopic dielectric
  // matrix in a matrix vector
  template <class RPAType>
  void CalcDielInvVector(const RPAType& rpa,
                         const Eigen::MatrixXd& kDielMxInv_zero);
  // Imx W_j Imx^T enters the integrand only through its diagonal, so it is
  // taken once per screening: node_terms_[gw_level](i, j) for the
  // quadrature node j. The node matrices W_j are not kept.
  void CalcNodeTerms(const std::vector<Eigen::MatrixXd>& dielinv_matrices);

  const Eigen::VectorXd& energies_;
  std::vector<Eigen::MatrixXd> node_terms_;
  const TCMatrix_gwbse& Mmn_;
};
}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IMAGINARYAXISINTEGRATION_H
