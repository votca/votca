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
#include <votca/tools/constants.h>

using std::flush;

namespace votca {
namespace xtp {

void ERDiabatization::setUpMatrices() {

  // Check on dftbasis
  if (_orbitals1.getDFTbasisName() != _orbitals2.getDFTbasisName()) {
    throw std::runtime_error("Different DFT basis for the two input file.");
  } else {

    XTP_LOG(Log::error, *_pLog) << "Data was created with basis set "
                                << _orbitals1.getDFTbasisName() << flush;
  }

  // Check on auxbasis
  if (_orbitals1.hasAuxbasisName()) {
    // this->_auxbasis1 = _orbitals1.SetupAuxBasis();
    _auxbasis1 = _orbitals1.getAuxBasis();
    _useRI = true;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << "No auxbasis for file1: This will affect perfomances " << flush;
    _useRI = false;
  }

  if (_orbitals2.hasAuxbasisName()) {
    _auxbasis2 = _orbitals2.getAuxBasis();
    _useRI = true;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << "No auxbasis for file2: This will affect perfomances " << flush;
    _useRI = false;
  }

  if (_orbitals1.getAuxbasisName() != _orbitals2.getAuxbasisName()) {
    throw std::runtime_error("Different DFT aux-basis for the two input file.");
  } else {
    if (_useRI) {
      XTP_LOG(Log::error, *_pLog) << "RI was used with Auxbasis set "
                                  << _orbitals1.getAuxbasisName() << flush;
    }
  }

  XTP_LOG(Log::debug, *_pLog) << "Setting up basis" << flush;

  _dftbasis1 = _orbitals1.getDftBasis();
  _bse_cmax1 = _orbitals1.getBSEcmax();
  _bse_cmin1 = _orbitals1.getBSEcmin();
  _bse_vmax1 = _orbitals1.getBSEvmax();
  _bse_vmin1 = _orbitals1.getBSEvmin();
  _bse_vtotal1 = _bse_vmax1 - _bse_vmin1 + 1;
  _bse_ctotal1 = _bse_cmax1 - _bse_cmin1 + 1;
  _basissize1 = _orbitals1.getBasisSetSize();
  //_bse_size_ao1 = _basis1 * _basis1;

  _occlevels1 = _orbitals1.MOs().eigenvectors().block(
      0, _bse_vmin1, _basissize1, _bse_vtotal1);
  _virtlevels1 = _orbitals1.MOs().eigenvectors().block(
      0, _bse_cmin1, _basissize1, _bse_ctotal1);

  // does this make any sense? Should be completely the same as in mol1
  _dftbasis2 = _orbitals2.getDftBasis();
  _bse_cmax2 = _orbitals2.getBSEcmax();
  _bse_cmin2 = _orbitals2.getBSEcmin();
  _bse_vmax2 = _orbitals2.getBSEvmax();
  _bse_vmin2 = _orbitals2.getBSEvmin();
  _bse_vtotal2 = _bse_vmax2 - _bse_vmin2 + 1;
  _bse_ctotal2 = _bse_cmax2 - _bse_cmin2 + 1;
  _basissize2 = _orbitals2.getBasisSetSize();
  _occlevels2 = _orbitals2.MOs().eigenvectors().block(
      0, _bse_vmin2, _basissize2, _bse_vtotal2);
  _virtlevels2 = _orbitals2.MOs().eigenvectors().block(
      0, _bse_cmin2, _basissize2, _bse_ctotal2);

  // Use different RI initialization according to what is in the orb files.
  if (_useRI) {
    XTP_LOG(Log::error, *_pLog) << "Using RI " << flush;
    _eris.Initialize(_dftbasis1, _auxbasis1);
  } else {
    XTP_LOG(Log::error, *_pLog) << "Not using RI. It can be slow. " << flush;
    _eris.Initialize_4c(_dftbasis1);
  }
}

void ERDiabatization::configure(const options_erdiabatization& opt) {
  _opt = opt;
}

double ERDiabatization::CalculateR(const Eigen::MatrixXd& D_JK,
                                   const Eigen::MatrixXd& D_LM) const {

  // Here I want to do \sum_{kl} (ij|kl) D^{LM}_{jk}. Is it right?
  Eigen::MatrixXd contracted;
  // I still have to figure how to try 3c and if it fails go to 4c
  if (_useRI) {
    contracted = _eris.CalculateERIs_3c(D_LM);
  } else {
    contracted = _eris.CalculateERIs_4c(D_LM, 1e-12);
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

Eigen::MatrixXd ERDiabatization::Calculate_diabatic_H(
    const double E1, const double E2, const double angle) const {
  Eigen::VectorXd ad_energies(2);
  ad_energies << E1, E2;

  // To uncomment
  XTP_LOG(Log::error, *_pLog)
      << "Adiabatic energies in eV "
      << "E1: " << E1 * votca::tools::conv::hrt2ev
      << " E2: " << E2 * votca::tools::conv::hrt2ev << flush;
  XTP_LOG(Log::error, *_pLog)
      << "Rotation angle (degrees) " << angle * 57.2958 << flush;

  Eigen::MatrixXd U = CalculateU(angle);
  return U.transpose() * ad_energies.asDiagonal() * U;
}

double ERDiabatization::Calculate_angle(const Orbitals& orb1,
                                        const Orbitals& orb2,
                                        QMStateType type) const {
  Eigen::Tensor<double, 4> rtensor = CalculateRtensor(orb1, orb2, type);

  double A_12 =
      rtensor(0, 1, 0, 1) - 0.25 * (rtensor(0, 0, 0, 0) + rtensor(1, 1, 1, 1) -
                                    2. * rtensor(0, 0, 1, 1));
  double B_12 = rtensor(0, 0, 0, 1) - rtensor(1, 1, 0, 1);

  // I keep these definitions here for the sake of debugging
  // double tan4alpha = -1.0*(B_12/A_12);
  // double sin4alpha = B_12 /(std::sqrt(A_12*A_12 + B_12*B_12));

  double cos4alpha = -A_12 / (std::sqrt(A_12 * A_12 + B_12 * B_12));

  double angle = 0.25 * std::acos(cos4alpha);

  XTP_LOG(Log::error, *_pLog) << "B12 " << B_12 << flush;
  XTP_LOG(Log::error, *_pLog) << "A12 " << A_12 << flush;

  XTP_LOG(Log::error, *_pLog)
      << "angle MAX (degrees) " << angle * 57.2958 << flush;

  return angle;
}

void ERDiabatization::Print_ERfunction(Eigen::VectorXd results) const {
  // // TO DO: This loop should be printed on a file
  std::cout << "\n" << std::endl;
  double Pi = votca::tools::conv::Pi;
  // Initial mixing angle
  double phi_in = 0.;
  // Final mixing angle
  double phi_fin = .5 * Pi;
  // We divide the interval into equal bits
  double step = (phi_fin - phi_in) / double(results.size());

  for (Index n = 0; n < results.size(); n++) {
    std::cout << (57.2958) * (phi_in + double(n + 1) * step) << " "
              << results(n) << std::endl;
  }

  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculation done. Selecting maximum " << flush;

  // Get all the ingredients we need for evaluating the diabatic Hamiltonian
  // We need the angle that maximise the ER functional
  Index pos_min;
  Index pos_max;

  double maxval = results.maxCoeff(&pos_max);
  double minval = results.minCoeff(&pos_min);

  XTP_LOG(Log::error, *_pLog)
      << "Maximum EF is: " << maxval << " at position " << pos_max << flush;
  XTP_LOG(Log::error, *_pLog)
      << "Minimum EF is: " << minval << " at position " << pos_min << flush;

  double angle = phi_in + double(pos_max + 1) * step;
  double angle_min = phi_in + double(pos_min + 1) * step;
  XTP_LOG(Log::error, *_pLog)
      << "From plot: "
      << "MAX " << angle * 57.2958 << " MIN " << angle_min * 57.2958 << flush;
}

Eigen::Tensor<double, 4> ERDiabatization::CalculateRtensor(
    const Orbitals& orb1, const Orbitals& orb2, QMStateType type) const {
  XTP_LOG(Log::error, *_pLog) << "Computing R tensor" << flush;
  Eigen::Tensor<double, 4> r_tensor(2, 2, 2, 2);
  for (Index J = 0; J < 2; J++) {
    for (Index K = 0; K < 2; K++) {
      Eigen::MatrixXd D_JK = CalculateD_R(orb1, orb2, type, J, K);
      if (!orb1.getTDAApprox() and !orb2.getTDAApprox()) {
        D_JK -= CalculateD_AR(orb1, orb2, type, J, K);
      }
      for (Index L = 0; L < 2; L++) {
        for (Index M = 0; M < 2; M++) {
          Eigen::MatrixXd D_LM = CalculateD_R(orb1, orb2, type, L, M);
          if (!orb1.getTDAApprox() and !orb2.getTDAApprox()) {
            D_LM -= CalculateD_AR(orb1, orb2, type, L, M);
          }
          r_tensor(J, K, L, M) = CalculateR(D_JK, D_LM);
        }
      }
    }
  }
  return r_tensor;
}

Eigen::VectorXd ERDiabatization::CalculateER(const Orbitals& orb1,
                                             const Orbitals& orb2,
                                             QMStateType type) const {

  Eigen::Tensor<double, 4> R_JKLM = CalculateRtensor(orb1, orb2, type);
  const double pi = votca::tools::conv::Pi;
  // Scanning through angles
  Eigen::VectorXd results = Eigen::VectorXd::Zero(360);
  // Initial mixing angle
  double phi_in = 0.;
  // Final mixing angle
  double phi_fin = 0.5 * pi;
  // We divide the interval into equal bits
  double step = (phi_fin - phi_in) / double(results.size());
  // Define angle we are iterating
  double phi;
  for (Index n = 0; n < results.size(); n++) {
    // Update angle
    phi = phi_in + double(n + 1) * step;
    Eigen::MatrixXd U = CalculateU(phi);
    // Complicated loop to handle. Can we make it a bit better?
    double result = 0.;
    for (Index I = 0; I < 2; I++) {
      for (Index J = 0; J < 2; J++) {
        for (Index K = 0; K < 2; K++) {
          for (Index L = 0; L < 2; L++) {
            for (Index M = 0; M < 2; M++) {
              result +=
                  U(I, J) * U(I, K) * U(I, L) * U(I, M) * R_JKLM(J, K, L, M);
            }
          }
        }
      }
    }
    results(n) = result;
  }
  return results;
}

template <bool AR>
Eigen::MatrixXd ERDiabatization::CalculateD(const Orbitals& orb1,
                                            const Orbitals& orb2,
                                            QMStateType type, Index stateindex1,
                                            Index stateindex2) const {

  // D matrix depends on 2 indeces. These can be either 0 or 1.
  // Index=0 means "take the first excited state" as Index=1 means "take the
  // second excitate state"
  // This is the reason for this
  Index index1;
  Index index2;

  if (stateindex1 == 0) {
    index1 = _opt.state_idx_1;
  } else if (stateindex1 == 1) {
    index1 = _opt.state_idx_2;
  } else {
    throw std::runtime_error("Invalid state index specified.");
  }
  if (stateindex2 == 0) {
    index2 = _opt.state_idx_1;
  } else if (stateindex2 == 1) {
    index2 = _opt.state_idx_2;
  } else {
    throw std::runtime_error("Invalid state index specified.");
  }

  Eigen::VectorXd exciton1;
  Eigen::VectorXd exciton2;

  // If AR==True we want Bs from BSE. If AR==False we want As from BSE.
  if (AR) {
    if (type == QMStateType::Singlet) {
      exciton1 = orb1.BSESinglets().eigenvectors2().col(index1 - 1);
      exciton2 = orb2.BSESinglets().eigenvectors2().col(index2 - 1);
    } else {
      exciton1 = orb1.BSETriplets().eigenvectors2().col(index1 - 1);
      exciton2 = orb2.BSETriplets().eigenvectors2().col(index2 - 1);
    }
  } else {
    if (type == QMStateType::Singlet) {
      exciton1 = orb1.BSESinglets().eigenvectors().col(index1 - 1);
      exciton2 = orb2.BSESinglets().eigenvectors().col(index2 - 1);
    } else {
      exciton1 = orb1.BSETriplets().eigenvectors().col(index1 - 1);
      exciton2 = orb2.BSETriplets().eigenvectors().col(index2 - 1);
    }
  }

  Eigen::Map<const Eigen::MatrixXd> mat1(exciton1.data(), _bse_ctotal1,
                                         _bse_vtotal1);
  Eigen::Map<const Eigen::MatrixXd> mat2(exciton2.data(), _bse_ctotal2,
                                         _bse_vtotal2);

  Eigen::MatrixXd AuxMat_vv = mat1.transpose() * mat2;

  Eigen::MatrixXd AuxMat_cc = mat1 * mat2.transpose();

  Eigen::MatrixXd results =
      _virtlevels1 * AuxMat_cc * _virtlevels2.transpose() -
      _occlevels1 * AuxMat_vv * _occlevels2.transpose();

  // This adds the GS density
  if (stateindex1 == stateindex2) {
    if (stateindex1 == 0) {
      results += orb1.DensityMatrixGroundState();
    }
    if (stateindex1 == 1) {
      results += orb2.DensityMatrixGroundState();
    }
  }

  return results;
}

}  // namespace xtp
}  // namespace votca
