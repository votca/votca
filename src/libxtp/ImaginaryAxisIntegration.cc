
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

#include <votca/tools/constants.h>

#include "votca/xtp/ImaginaryAxisIntegration.h"
#include "votca/xtp/quadrature_factory.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

ImaginaryAxisIntegration::ImaginaryAxisIntegration(
    const Eigen::VectorXd& energies, const TCMatrix_gwbse& Mmn)
    : _energies(energies), _Mmn(Mmn) {}

void ImaginaryAxisIntegration::configure(
    options opt, const RPA& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  _opt = opt;
  QuadratureFactory factory;
  _gq = std::unique_ptr<GaussianQuadratureBase>(
      factory.Create(_opt.quadrature_scheme));
  _gq->configure(_opt.order);

  CalcDielInvVector(rpa, kDielMxInv_zero);
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
void ImaginaryAxisIntegration::CalcDielInvVector(
    const RPA& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  _dielinv_matrices_r.resize(_gq->Order());

  for (Index j = 0; j < _gq->Order(); j++) {
    double newpoint = _gq->ScaledPoint(j);
    Eigen::MatrixXd eps_inv_j = rpa.calculate_epsilon_i(newpoint).inverse();
    eps_inv_j.diagonal().array() -= 1.0;
    _dielinv_matrices_r[j] =
        -eps_inv_j +
        kDielMxInv_zero * std::exp(-std::pow(_opt.alpha * newpoint, 2));
  }
}

class FunctionEvaluation {
 public:
  FunctionEvaluation(const Eigen::MatrixXd& Imx, const Eigen::ArrayXcd& DeltaE,
                     const std::vector<Eigen::MatrixXd>& dielinv_matrices_r)
      : _Imx(Imx), _DeltaE(DeltaE), _dielinv_matrices_r(dielinv_matrices_r){};

  double operator()(Index j, double point, bool symmetry) const {
    Eigen::VectorXcd denominator;
    const std::complex<double> cpoint(0.0, point);
    if (symmetry) {
      denominator =
          (_DeltaE + cpoint).cwiseInverse() + (_DeltaE - cpoint).cwiseInverse();
    } else {
      denominator = (_DeltaE + cpoint).cwiseInverse();
    }
    return 0.5 / tools::conv::Pi *
           ((_Imx * (_dielinv_matrices_r[j].conjugate()))
                .cwiseProduct(denominator.asDiagonal() * _Imx))
               .sum()
               .real();
  }

 private:
  const Eigen::MatrixXd& _Imx;
  const Eigen::ArrayXcd& _DeltaE;
  const std::vector<Eigen::MatrixXd>& _dielinv_matrices_r;
};

double ImaginaryAxisIntegration::SigmaGQDiag(double frequency, Index gw_level,
                                             double eta) const {
  Index lumo = _opt.homo + 1;
  const Index occ = lumo - _opt.rpamin;
  const Index unocc = _opt.rpamax - _opt.homo;
  Index gw_level_offset = gw_level + _opt.qpmin - _opt.rpamin;
  const Eigen::MatrixXd& Imx = _Mmn[gw_level_offset];
  Eigen::ArrayXcd DeltaE = frequency - _energies.array();
  DeltaE.imag().head(occ) = eta;
  DeltaE.imag().tail(unocc) = -eta;
  FunctionEvaluation f(Imx, DeltaE, _dielinv_matrices_r);
  return _gq->Integrate(f);
}

}  // namespace xtp

}  // namespace votca
