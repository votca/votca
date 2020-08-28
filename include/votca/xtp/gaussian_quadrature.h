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

#ifndef __XTP_GAUSSIAN_QUADRATURE__H
#define __XTP_GAUSSIAN_QUADRATURE__H

#include <votca/xtp/eigen.h>
#include <votca/xtp/rpa.h>

// Computes the contribution from the Gauss-Laguerre quadrature to the
// self-energy expectation matrix for given RPA and frequencies
namespace votca {
namespace xtp {

class GaussianQuadrature {

 public:
  struct options {
    Index order = 12;
    Index qptotal;
    Index qpmin;
    Index homo;
    Index rpamin;
    Index rpamax;
    double alpha = 0.001;
    std::string quadrature_scheme;
  };

  GaussianQuadrature(const Eigen::VectorXd& energies,
                     const TCMatrix_gwbse& Mmn);

  void configure(options opt, const RPA& rpa);

  double SigmaGQDiag(double frequency, Index gw_level, double eta) const;
 
 private:
  options _opt;

  // This function calculates and stores inverses of the microscopic dielectric
  // matrix in a matrix vector
  void CalcDielInvVector(const RPA& rpa);
  const Eigen::VectorXd& _energies;

  std::vector<Eigen::MatrixXcd> _dielinv_matrices_r;
  const TCMatrix_gwbse& _Mmn;
  Eigen::VectorXd _quadpoints;
  Eigen::VectorXd _quadadaptedweights;
};
}  // namespace xtp
}  // namespace votca
#endif /* GAUSSIAN_QUADRATURE_H */
