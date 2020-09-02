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

void GaussianQuadrature::configure(options opt, const RPA& rpa) {
  _opt = opt;
  if (_opt.quadrature_scheme == "laguerre") {
    Gauss_Laguerre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa);
  } else if (_opt.quadrature_scheme == "legendre") {
    Gauss_Legendre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa);
  } else if (_opt.quadrature_scheme == "modified_legendre") {
    Gauss_Legendre_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa);
  } else if (_opt.quadrature_scheme == "hermite") {
    Gauss_Hermite_Quadrature_Constants glqc;
    _quadpoints = glqc.getPoints(_opt.order);
    _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
    CalcDielInvVector(rpa);
  } else {
    std::cout << "There no such a thing as the integration scheme you asked"
              << std::endl;
  }
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
void GaussianQuadrature::CalcDielInvVector(const RPA& rpa) {
  _dielinv_matrices_r.resize(_opt.order);

  double halfpi = 0.5 * votca::tools::conv::Pi;
#pragma openmp parallel schedule(guided)
  for (Index j = 0; j < _opt.order; j++) {
    if (_opt.quadrature_scheme == "legendre") {
      double newpoint = std::tan(halfpi * _quadpoints(j));
      Eigen::MatrixXcd eps_inv_j =
          rpa.calculate_epsilon_complex(0.0, newpoint).inverse();
      eps_inv_j.diagonal().array() -= 1.0;
      _dielinv_matrices_r[j] = -eps_inv_j;
    } else if (_opt.quadrature_scheme == "modified_legendre") {
      double exponent = (1.0 + _quadpoints(j)) / (1.0 - _quadpoints(j));
      double newpoint = std::pow(0.5, exponent);
      Eigen::MatrixXcd eps_inv_j =
          rpa.calculate_epsilon_complex(0.0, newpoint).inverse();
      eps_inv_j.diagonal().array() -= 1.0;
      _dielinv_matrices_r[j] = -eps_inv_j;
    } else {
      double newpoint = _quadpoints(j);
      Eigen::MatrixXcd eps_inv_j =
          rpa.calculate_epsilon_complex(0.0, newpoint).inverse();
      eps_inv_j.diagonal().array() -= 1.0;
      _dielinv_matrices_r[j] = -eps_inv_j;
    }
  }
}

double GaussianQuadrature::SigmaGQDiag(double frequency, Index gw_level,
                                       double eta) const {
  Index lumo = _opt.homo + 1;
  const Index occ = lumo - _opt.rpamin;
  const Index unocc = _opt.rpamax - _opt.homo;
  
  const Eigen::MatrixXd& Imx = _Mmn[gw_level];
  Eigen::ArrayXcd DeltaE = frequency - _energies.array();
  std::complex<double> result = std::complex<double>(0.0, 0.0);
  if (_opt.quadrature_scheme == "laguerre") {
    // The laguerre quadrature is suitable for integration limits a = 0 b =
    // +infty

    DeltaE.head(occ).imag() = eta;
    DeltaE.tail(unocc).imag() = -eta;

    for (Index j = 0; j < _opt.order; ++j) {
      double newpoint = _quadpoints(j);
      double newweight = _quadadaptedweights(j);
      Eigen::VectorXcd denominator1 =
          (DeltaE + std::complex<double>(0.0, newpoint)).cwiseInverse();
      // Eigen::MatrixXcd Amx = denominator1.asDiagonal() * Imx;
      // Eigen::MatrixXcd Cmx = Imx * (_dielinv_matrices_r[j]);
      std::complex<double> value1 =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator1.asDiagonal() * Imx))
              .sum();
      Eigen::VectorXcd denominator2 =
          (DeltaE + std::complex<double>(0.0, -newpoint)).cwiseInverse();
      // Eigen::MatrixXcd Dmx = denominator2.asDiagonal() * Imx;
      // Eigen::MatrixXcd Emx = Imx * (_dielinv_matrices_r[j].conjugate());
      std::complex<double> value2 =
          ((Imx * (_dielinv_matrices_r[j].conjugate()))
               .cwiseProduct(denominator2.asDiagonal() * Imx))
              .sum();

      result += newweight * (value1 + value2);
    }

    result *= 0.5;
  } else if (_opt.quadrature_scheme == "modified_legendre") {
    // The modified legendre method is suitable for integration limits a = 0 b =
    // +infty

    DeltaE.head(occ).imag() = eta;
    DeltaE.tail(unocc).imag() = -eta;

    for (Index j = 0; j < _opt.order; ++j) {
      double exponent = (1.0 + _quadpoints(j)) / (1.0 - _quadpoints(j));
      double newpoint = std::pow(0.5, exponent);
      double den =
          (1.0 - _quadadaptedweights(j)) * (1.0 - _quadadaptedweights(j));
      double newweight = (2.0 * _quadadaptedweights(j) * 0.5) / den;
      Eigen::VectorXcd denominator1 =
          (DeltaE + std::complex<double>(0.0, newpoint)).cwiseInverse();
      std::complex<double> value1 =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator1.asDiagonal() * Imx))
              .sum();
      Eigen::VectorXcd denominator2 =
          (DeltaE + std::complex<double>(0.0, -newpoint)).cwiseInverse();
      std::complex<double> value2 =
          ((Imx * (_dielinv_matrices_r[j].conjugate()))
               .cwiseProduct(denominator2.asDiagonal() * Imx))
              .sum();

      result += newweight * (value1 + value2);
    }

    result *= 0.5;

  } else if (_opt.quadrature_scheme == "hermite") {
    // The hermite quadrature method is suitable for integration limits a =
    // -infty b = +infty

    DeltaE.head(occ).imag() = eta;
    DeltaE.tail(unocc).imag() = -eta;

    for (Index j = 0; j < _opt.order; ++j) {

      Eigen::VectorXcd denominator =
          (DeltaE + std::complex<double>(0.0, _quadpoints(j))).cwiseInverse();
      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator.asDiagonal() * Imx))
              .sum();
      result += _quadadaptedweights(j) * value;
    }

    result *= 0.5;
  } else if (_opt.quadrature_scheme == "laguerre") {
    // The laguerre quadrature is suitable for integration limits a = 0 b =
    // +infty
    for (Index j = 0; j < _opt.order; ++j) {
      Eigen::VectorXcd coeffs1 =
          (DeltaE) *
          (DeltaE.square() + std::pow(_quadpoints(j), 2)).cwiseInverse();

      std::complex<double> value =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(coeffs1.asDiagonal() * Imx))
              .sum();
      result += _quadadaptedweights(j) * value;
    }
  } else if (_opt.quadrature_scheme == "legendre") {
    // This particular legendre quadrature is suitable for integration limits a
    // = -infty b = +infty

    DeltaE.head(occ).imag() = eta;
    DeltaE.tail(unocc).imag() = -eta;

    double halfpi = std::acos(0.0);
    for (Index j = 0; j < _opt.order; ++j) {

      double newpoint = std::tan(halfpi * _quadpoints(j));
      double den =
          std::cos(halfpi * _quadpoints(j)) * std::cos(halfpi * _quadpoints(j));
      double num = halfpi / den;
      double newweight = _quadadaptedweights(j);
      Eigen::VectorXcd denominator1 =
          (num) * (DeltaE + std::complex<double>(0.0, newpoint)).cwiseInverse();
      // Eigen::MatrixXcd Amx = denominator1.asDiagonal() * Imx;
      // Eigen::MatrixXcd Cmx = Imx * (_dielinv_matrices_r[j]);
      std::complex<double> value1 =
          ((Imx * (_dielinv_matrices_r[j]))
               .cwiseProduct(denominator1.asDiagonal() * Imx))
              .sum();
      result += newweight * (value1);
    }

    result *= 0.5;
  }

  return result.real() / (tools::conv::Pi);
}

}  // namespace xtp

}  // namespace votca
