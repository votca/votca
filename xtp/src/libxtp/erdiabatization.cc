/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// VOTCA includes
#include "votca/xtp/erdiabatization.h"
#include "votca/xtp/ERIs.h"
#include "votca/xtp/transition_densities.h"
#include <votca/tools/constants.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void ERDiabatization::setUpMatrices() {

  // Check on dftbasis
  if (orbitals1_.getDFTbasisName() != orbitals2_.getDFTbasisName()) {
    throw std::runtime_error("Different DFT basis for the two input file.");
  } else {
    dftbasis_ = orbitals1_.getDftBasis();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Data was created with basis set "
        << orbitals1_.getDFTbasisName() << flush;
  }

  // Check on auxbasis
  if (orbitals1_.hasAuxbasisName() && orbitals2_.hasAuxbasisName()) {
    if (orbitals1_.getAuxbasisName() == orbitals2_.getAuxbasisName()) {
      auxbasis_ = orbitals1_.getAuxBasis();
      hasRI_ = true;
    } else {
      throw std::runtime_error(
          "Different DFT aux-basis for the two input file.");
    }
  } else {
    XTP_LOG(Log::error, *pLog_)
        << "No auxbasis for file1: This will affect perfomances " << flush;
    hasRI_ = false;
  }

  // check BSE indices
  if (orbitals1_.getBSEvmin() != orbitals2_.getBSEvmin()) {
    throw std::runtime_error("Different BSE vmin for the two input file.");
  }

  if (orbitals1_.getBSEvmax() != orbitals2_.getBSEvmax()) {
    throw std::runtime_error("Different BSE vmax for the two input file.");
  }

  if (orbitals1_.getBSEcmin() != orbitals2_.getBSEcmin()) {
    throw std::runtime_error("Different BSE cmin for the two input file.");
  }

  if (orbitals1_.getBSEcmax() != orbitals2_.getBSEcmax()) {
    throw std::runtime_error("Different BSE cmax for the two input file.");
  }

  // Use different RI initialization according to what is in the orb files.
  if (hasRI_ && useRI_) {
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Using RI " << flush;
    eris_.Initialize(dftbasis_, auxbasis_);
  } else {
    eris_.Initialize_4c(dftbasis_);
  }
}

void ERDiabatization::configure() {
  qmtype_.FromString(qmtype_str_);

  if (qmtype_ == QMStateType::Singlet) {
    E1_ = orbitals1_.BSESinglets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSESinglets().eigenvalues()[state_idx_2_ - 1];
  } else {
    E1_ = orbitals1_.BSETriplets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSETriplets().eigenvalues()[state_idx_2_ - 1];
  }
}

double ERDiabatization::CalculateR(const Eigen::MatrixXd& D_JK,
                                   const Eigen::MatrixXd& D_LM) const {

  Eigen::MatrixXd contracted;
  if (hasRI_ && useRI_) {
    contracted = eris_.CalculateERIs_3c(D_LM);
  } else {
    contracted = eris_.CalculateERIs_4c(D_LM, 1e-12);
  }

  return D_JK.cwiseProduct(contracted).sum();
}

Eigen::MatrixXd ERDiabatization::CalculateU(const double phi) const {
  Eigen::MatrixXd U(2, 2);
  U(0, 0) = std::cos(phi);
  U(0, 1) = -1.0 * std::sin(phi);
  U(1, 0) = std::sin(phi);
  U(1, 1) = std::cos(phi);
  return U;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd>
    ERDiabatization::Calculate_diabatic_H(const double angle) const {
  Eigen::VectorXd ad_energies(2);
  ad_energies << E1_, E2_;

  XTP_LOG(Log::debug, *pLog_)
      << format("Adiabatic energies: %1$+1.12f eV and %2$+1.12f eV") %
             (E1_ * votca::tools::conv::hrt2ev) %
             (E2_ * votca::tools::conv::hrt2ev)
      << flush;

  XTP_LOG(Log::debug, *pLog_)
      << TimeStamp() << "Rotation angle (degrees) " << angle * 57.2958 << flush;

  Eigen::MatrixXd U = CalculateU(angle);
  return std::pair<Eigen::VectorXd, Eigen::MatrixXd>(
      ad_energies, U.transpose() * ad_energies.asDiagonal() * U);
}

double ERDiabatization::Calculate_angle() const {
  Eigen::Tensor<double, 4> rtensor = CalculateRtensor();

  double A_12 =
      rtensor(0, 1, 0, 1) - 0.25 * (rtensor(0, 0, 0, 0) + rtensor(1, 1, 1, 1) -
                                    2. * rtensor(0, 0, 1, 1));
  double B_12 = rtensor(0, 0, 0, 1) - rtensor(1, 1, 0, 1);

  double cos4alpha = -A_12 / (std::sqrt(A_12 * A_12 + B_12 * B_12));

  // also calculate sin4alpha
  double sin4alpha = B_12 / (std::sqrt(A_12 * A_12 + B_12 * B_12));

  // In some cases acos may give the wrong solution, that das not maych
  // sin4alpha Instead, follow procedure in original ER paper
  double sqroot = std::sqrt(1.0 - 0.5 * (1.0 - cos4alpha));
  double x_squared_plus = 0.5 * (1.0 + sqroot);
  double x_squared_minus = 0.5 * (1.0 - sqroot);

  double x1 = sqrt(x_squared_plus);
  double y1 = sqrt(1.0 - x_squared_plus);

  double x2 = sqrt(x_squared_minus);
  double y2 = sqrt(1.0 - x_squared_minus);

  double test1 = 4.0 * x1 * y1 * (x1 * x1 - y1 * y1);
  double test2 = 4.0 * x2 * y2 * (x2 * x2 - y2 * y2);

  double cos_alpha0 = 0;
  double sin_alpha0 = 0;

  // check which (x,y) pair matches sin(4alpha)
  if (std::abs(test1 - sin4alpha) < 1e-4) {
    cos_alpha0 = x1;
    sin_alpha0 = y1;
  } else if (std::abs(test2 - sin4alpha) < 1e-4) {
    cos_alpha0 = x2;
    sin_alpha0 = y2;
  } else {
    XTP_LOG(Log::debug, *pLog_) << "Can't find correct angle " << flush;
  }

  XTP_LOG(Log::debug, *pLog_)
      << "Coupling element directly: "
      << cos_alpha0 * sin_alpha0 * (E2_ - E1_) * votca::tools::conv::hrt2ev
      << flush;

  double angle = std::acos(cos_alpha0);

  XTP_LOG(Log::debug, *pLog_) << "B12 " << B_12 << flush;
  XTP_LOG(Log::debug, *pLog_) << "A12 " << A_12 << flush;

  XTP_LOG(Log::debug, *pLog_)
      << "angle MAX (degrees) " << angle * 57.2958 << flush;

  return angle;
}

Eigen::Tensor<double, 4> ERDiabatization::CalculateRtensor() const {
  Eigen::Tensor<double, 4> r_tensor(2, 2, 2, 2);

  // get a transition density calculator object
  TransitionDensities tdmat(orbitals1_, orbitals2_, pLog_);
  tdmat.configure();

  // define QM states
  QMState state1 = QMState(qmtype_, state_idx_1_ - 1, false);
  QMState state2 = QMState(qmtype_, state_idx_2_ - 1, false);
  std::vector<QMState> states;
  states.push_back(state1);
  states.push_back(state2);

  // determine tensor elements
  for (Index J = 0; J < 2; J++) {
    for (Index K = 0; K < 2; K++) {
      Eigen::MatrixXd D_JK = tdmat.Matrix(states[J], states[K]);
      if (J == 0 && K == 0) {
        D_JK += orbitals1_.DensityMatrixGroundState();
      }
      if (J == 1 && K == 1) {
        D_JK += orbitals2_.DensityMatrixGroundState();
      }
      for (Index L = 0; L < 2; L++) {
        for (Index M = 0; M < 2; M++) {
          Eigen::MatrixXd D_LM = tdmat.Matrix(states[L], states[M]);
          if (L == 0 && M == 0) {
            D_LM += orbitals1_.DensityMatrixGroundState();
          }
          if (L == 1 && M == 1) {
            D_LM += orbitals2_.DensityMatrixGroundState();
          }
          r_tensor(J, K, L, M) = CalculateR(D_JK, D_LM);
        }
      }
    }
  }
  return r_tensor;
}

}  // namespace xtp
}  // namespace votca
