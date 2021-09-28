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

  double SigmaGQDiag(double frequency, Index gw_level, double eta) const;

 private:
  options opt_;

  std::unique_ptr<GaussianQuadratureBase> gq_ = nullptr;

  // This function calculates and stores inverses of the microscopic dielectric
  // matrix in a matrix vector
  void CalcDielInvVector(const RPA& rpa,
                         const Eigen::MatrixXd& kDielMxInv_zero);
  const Eigen::VectorXd& energies_;
  std::vector<Eigen::MatrixXd> dielinv_matrices_r_;
  const TCMatrix_gwbse& Mmn_;
};
}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IMAGINARYAXISINTEGRATION_H
