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

#include <boost/progress.hpp>

#include <votca/tools/constants.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>
#include <votca/xtp/gauss_laguerre_quadrature_constants.h>
#include <votca/xtp/gauss_legendre_quadrature_constants.h>
#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

// Constructor
GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& energies,
                                       const TCMatrix_gwbse& Mmn)
    : _energies(energies), _Mmn(Mmn) {}

void GaussianQuadrature::configure(options opt, const RPA& rpa,
                                   const Eigen::MatrixXd& kDielMxInv_zero) {
  _opt = opt;
  if (_opt.quadrature_scheme == "laguerre") {
    Gauss_Laguerre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa, kDielMxInv_zero);
  } else if (_opt.quadrature_scheme == "legendre") {
    Gauss_Legendre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa, kDielMxInv_zero);
  } else if (_opt.quadrature_scheme == "modified_legendre") {
    Gauss_Legendre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa, kDielMxInv_zero);
  } else if (_opt.quadrature_scheme == "hermite") {
    Gauss_Hermite_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa, kDielMxInv_zero);
  } else {
    std::cout << "There no such a thing as the integration scheme you asked"
              << std::endl;
  }
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
void GaussianQuadrature::CalcDielInvVector(
    const RPA& rpa, const Eigen::MatrixXd& kDielMxInv_zero) {
  _dielinv_matrices_r.resize(_opt.order);
  Eigen::MatrixXd eps_inv_j;
  // the matrix here
  Eigen::MatrixXd eps_inv_j_tail;
  double halfpi = 0.5 * votca::tools::conv::Pi;
  double newpoint = 0.0;

  std::cout << "\n... ... Preparing RPA for Gaussian quadrature along "
               "imaginary axis with "
            << _opt.order << " points" << std::endl;
  boost::progress_display progress(_opt.order, std::cout, "... ... ",
                                   "... ... ", "... ... ");
  for (Index j = 0; j < _opt.order; j++) {
    if (_opt.quadrature_scheme == "legendre") {
      newpoint = std::tan(halfpi * _quadpoints(j));
    } else if (_opt.quadrature_scheme == "modified_legendre") {
      double exponent = (1.0 + _quadpoints(j)) / (1.0 - _quadpoints(j));
      newpoint = std::pow(0.5, exponent);
    } else {
      newpoint = _quadpoints(j);
    }
    eps_inv_j = rpa.calculate_epsilon_i(newpoint).inverse();
    eps_inv_j.diagonal().array() -= 1.0;

    _dielinv_matrices_r[j] =
        -eps_inv_j +
        kDielMxInv_zero * std::exp(-std::pow(_opt.alpha * newpoint, 2));
    ++progress;
  }
}

double GaussianQuadrature::SigmaGQDiag(double frequency, Index gw_level,
                                       double eta) const {
  Index lumo = _opt.homo + 1;
  const Index occ = lumo - _opt.rpamin;
  const Index unocc = _opt.rpamax - _opt.homo;
  Index gw_level_offset = gw_level + _opt.qpmin - _opt.rpamin;
  const Eigen::MatrixXd& Imx = _Mmn[gw_level_offset];
  Eigen::ArrayXcd DeltaE = frequency - _energies.array();
  DeltaE.head(occ).imag() = eta;
  DeltaE.tail(unocc).imag() = -eta;

  std::complex<double> result(0.0, 0.0);
  if (_opt.quadrature_scheme == "laguerre") {
    // The laguerre method is suitable for integration limits a = 0 b =
    // +infty. Here we have value1 and value2 because we split the original
    // integral from -infty to +infty in two parts
    for (Index j = 0; j < _opt.order; ++j) {
      std::complex<double> newpoint(0.0, _quadpoints(j));
      Eigen::VectorXcd denominator = (DeltaE + newpoint).cwiseInverse() +
                                     (DeltaE - newpoint).cwiseInverse();
      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j].conjugate()))
               .cwiseProduct(denominator.asDiagonal() * Imx))
              .sum();
      result += _quadadaptedweights(j) * value;
    }

  } else if (_opt.quadrature_scheme == "modified_legendre") {
    // The modified legendre method is suitable for integration limits a = 0 b =
    // +infty. Here we have value1 and value2 because we split the original
    // integral from -infty to +infty in two parts. Original legendre quadrature
    // is meant for integratal with integration limits of -1 and 1. To overcome
    // this we use the transformation x' = 0.5 ^ (1+x/1-x)
    for (Index j = 0; j < _opt.order; ++j) {
      double exponent = (1.0 + _quadpoints(j)) / (1.0 - _quadpoints(j));
      std::complex<double> newpoint(0.0, std::pow(0.5, exponent));
      double den =
          (1.0 - _quadadaptedweights(j)) * (1.0 - _quadadaptedweights(j));
      double newweight = (2.0 * _quadadaptedweights(j) * 0.5) / den;
      Eigen::VectorXcd denominator = (DeltaE + newpoint).cwiseInverse() +
                                     (DeltaE - newpoint).cwiseInverse();
      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j].conjugate()))
               .cwiseProduct(denominator.asDiagonal() * Imx))
              .sum();
      result += newweight * value;
    }

  } else if (_opt.quadrature_scheme == "hermite") {
    // The hermite quadrature method is suitable for integration limits a =
    // -infty b = +infty. Here we don't do any modification to the original
    // hermite quadrature method. Points and weights are not transformed
    for (Index j = 0; j < _opt.order; ++j) {
      std::complex<double> newpoint(0.0, _quadpoints(j));
      Eigen::VectorXcd denominator = (DeltaE + newpoint).cwiseInverse();
      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator.asDiagonal() * Imx))
              .sum();
      result += _quadadaptedweights(j) * value;
    }

  } else if (_opt.quadrature_scheme == "legendre") {
    // This particular legendre quadrature is suitable for integration limits a
    // = -infty b = +infty. Original Legendre quadrature is meant for
    // integration limits -1 and +1. The change of variables is x' = tan (pi/2 *
    // x). When x=-1 we have x'=-infty. When x=1 we have x'=+infty
    double halfpi = 0.5 * tools::conv::Pi;
    for (Index j = 0; j < _opt.order; ++j) {
      std::complex<double> newpoint(0.0, std::tan(halfpi * _quadpoints(j)));
      double den =
          std::cos(halfpi * _quadpoints(j)) * std::cos(halfpi * _quadpoints(j));
      double num = halfpi / den;
      Eigen::VectorXcd denominator = (DeltaE + newpoint).cwiseInverse();
      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator.asDiagonal() * Imx))
              .sum();
      result += num * _quadadaptedweights(j) * value;
    }
  }
  return 0.5 * result.real() / (tools::conv::Pi);
}

}  // namespace xtp

}  // namespace votca
