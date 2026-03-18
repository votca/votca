
/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <set>
#include <sstream>

#include "votca/xtp/ImaginaryAxisIntegration.h"
#include "votca/xtp/quadrature_factory.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

namespace {

bool CdaDebugEnabled() {
  const char* env = std::getenv("VOTCA_XTP_CDA_DEBUG");
  return env != nullptr && std::string(env) == "1";
}

bool CdaDebugLevelMatch(Index level) {
  const char* env = std::getenv("VOTCA_XTP_CDA_DEBUG_LEVELS");
  if (env == nullptr) {
    return true;
  }

  std::set<int> allowed;
  std::stringstream ss(env);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (!token.empty()) {
      allowed.insert(std::stoi(token));
    }
  }
  return allowed.count(static_cast<int>(level)) > 0;
}

bool CdaDebugMatch(Index level) {
  return CdaDebugEnabled() && CdaDebugLevelMatch(level);
}

}  // namespace

ImaginaryAxisIntegration::ImaginaryAxisIntegration(
    const Eigen::VectorXd& energies, const TCMatrix_gwbse& Mmn)
    : energies_(energies), Mmn_(Mmn) {}

void ImaginaryAxisIntegration::configure(
    options opt, const RPA& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  opt_ = opt;
  gq_ = QuadratureFactory().Create(opt_.quadrature_scheme);
  gq_->configure(opt_.order);

  CalcDielInvVector(rpa, kDielMxInv_zero);
}

void ImaginaryAxisIntegration::configure(
    options opt, const RPA_UKS& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  opt_ = opt;
  gq_ = QuadratureFactory().Create(opt_.quadrature_scheme);
  gq_->configure(opt_.order);

  CalcDielInvVector(rpa, kDielMxInv_zero);
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
void ImaginaryAxisIntegration::CalcDielInvVector(
    const RPA& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  dielinv_matrices_r_.resize(gq_->Order());

  for (Index j = 0; j < gq_->Order(); j++) {
    double newpoint = gq_->ScaledPoint(j);
    Eigen::MatrixXd eps_inv_j = rpa.calculate_epsilon_i(newpoint).inverse();
    eps_inv_j.diagonal().array() -= 1.0;
    dielinv_matrices_r_[j] =
        -eps_inv_j +
        kDielMxInv_zero * std::exp(-std::pow(opt_.alpha * newpoint, 2));
  }
}

void ImaginaryAxisIntegration::CalcDielInvVector(
    const RPA_UKS& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  dielinv_matrices_r_.resize(gq_->Order());

  for (Index j = 0; j < gq_->Order(); j++) {
    double newpoint = gq_->ScaledPoint(j);
    Eigen::MatrixXd eps_inv_j = rpa.calculate_epsilon_i(newpoint).inverse();
    eps_inv_j.diagonal().array() -= 1.0;
    dielinv_matrices_r_[j] =
        -eps_inv_j +
        kDielMxInv_zero * std::exp(-std::pow(opt_.alpha * newpoint, 2));
  }
}

class FunctionEvaluation {
 public:
  FunctionEvaluation(const Eigen::MatrixXd& Imx, const Eigen::ArrayXcd& DeltaE,
                     const std::vector<Eigen::MatrixXd>& dielinv_matrices_r)
      : Imx_(Imx), DeltaE_(DeltaE), dielinv_matrices_r_(dielinv_matrices_r) {};

  double operator()(Index j, double point, bool symmetry) const {
    Eigen::VectorXcd denominator;
    const std::complex<double> cpoint(0.0, point);
    if (symmetry) {
      denominator =
          (DeltaE_ + cpoint).cwiseInverse() + (DeltaE_ - cpoint).cwiseInverse();
    } else {
      denominator = (DeltaE_ + cpoint).cwiseInverse();
    }
    return 0.5 / tools::conv::Pi *
           ((Imx_ * (dielinv_matrices_r_[j].conjugate()))
                .cwiseProduct(denominator.asDiagonal() * Imx_))
               .sum()
               .real();
  }

 private:
  const Eigen::MatrixXd& Imx_;
  const Eigen::ArrayXcd& DeltaE_;
  const std::vector<Eigen::MatrixXd>& dielinv_matrices_r_;
};

double ImaginaryAxisIntegration::SigmaGQDiag(double frequency, Index gw_level,
                                             double eta) const {
  Index lumo = opt_.homo + 1;
  const Index occ = lumo - opt_.rpamin;
  const Index unocc = opt_.rpamax - opt_.homo;
  Index gw_level_offset = gw_level + opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd& Imx = Mmn_[gw_level_offset];
  Eigen::ArrayXcd DeltaE = frequency - energies_.array();
  DeltaE.imag().head(occ) = eta;
  DeltaE.imag().tail(unocc) = -eta;

  if (CdaDebugMatch(gw_level)) {
    double min_abs = std::numeric_limits<double>::max();
    for (Index i = 0; i < DeltaE.size(); ++i) {
      min_abs = std::min(min_abs, std::abs(DeltaE(i)));
    }

    std::cout << "\n [CDA-DEBUG] ImagAxis level=" << gw_level
              << " freq=" << frequency << " eta=" << eta
              << " min|DeltaE|=" << min_abs << " ||Imx||=" << Imx.norm();

    for (Index i = 0; i < DeltaE.size(); ++i) {
      std::cout << "\n [CDA-DEBUG] ImagAxis i=" << i
                << " occ=" << (i < occ ? 1 : 0)
                << " DeltaE=" << DeltaE(i);
    }
    std::cout << std::flush;
  }

  FunctionEvaluation f(Imx, DeltaE, dielinv_matrices_r_);
  return gq_->Integrate(f);
}

}  // namespace xtp
}  // namespace votca