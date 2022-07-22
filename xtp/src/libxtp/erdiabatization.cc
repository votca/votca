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
    basissize_ = orbitals1_.getBasisSetSize();
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
  } else {
    bse_vmin_ = orbitals1_.getBSEvmin();
  }

  if (orbitals1_.getBSEvmax() != orbitals2_.getBSEvmax()) {
    throw std::runtime_error("Different BSE vmax for the two input file.");
  } else {
    bse_vmax_ = orbitals1_.getBSEvmax();
  }

  if (orbitals1_.getBSEcmin() != orbitals2_.getBSEcmin()) {
    throw std::runtime_error("Different BSE cmin for the two input file.");
  } else {
    bse_cmin_ = orbitals1_.getBSEcmin();
  }

  if (orbitals1_.getBSEcmax() != orbitals2_.getBSEcmax()) {
    throw std::runtime_error("Different BSE cmax for the two input file.");
  } else {
    bse_cmax_ = orbitals1_.getBSEcmax();
  }

  bse_vtotal_ = bse_vmax_ - bse_vmin_ + 1;
  bse_ctotal_ = bse_cmax_ - bse_cmin_ + 1;

  occlevels1_ = orbitals1_.MOs().eigenvectors().block(0, bse_vmin_, basissize_,
                                                      bse_vtotal_);
  virtlevels1_ = orbitals1_.MOs().eigenvectors().block(0, bse_cmin_, basissize_,
                                                       bse_ctotal_);

  occlevels2_ = orbitals2_.MOs().eigenvectors().block(0, bse_vmin_, basissize_,
                                                      bse_vtotal_);
  virtlevels2_ = orbitals2_.MOs().eigenvectors().block(0, bse_cmin_, basissize_,
                                                       bse_ctotal_);

  // Use different RI initialization according to what is in the orb files.
  if (hasRI_ && opt_.use_RI) {
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Using RI " << flush;
    eris_.Initialize(dftbasis_, auxbasis_);
  } else {
    eris_.Initialize_4c(dftbasis_);
  }
}

void ERDiabatization::configure(const options_erdiabatization& opt) {
  opt_ = opt;
  qmtype_.FromString(opt.qmtype);

  if (qmtype_ == QMStateType::Singlet) {
    E1_ = orbitals1_.BSESinglets().eigenvalues()[opt_.state_idx_1 - 1];
    E2_ = orbitals2_.BSESinglets().eigenvalues()[opt_.state_idx_2 - 1];
  } else {
    E1_ = orbitals1_.BSETriplets().eigenvalues()[opt_.state_idx_1 - 1];
    E2_ = orbitals2_.BSETriplets().eigenvalues()[opt_.state_idx_2 - 1];
  }
}

double ERDiabatization::CalculateR(const Eigen::MatrixXd& D_JK,
                                   const Eigen::MatrixXd& D_LM) const {

  Eigen::MatrixXd contracted;
  if (hasRI_ && opt_.use_RI) {
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

Eigen::MatrixXd ERDiabatization::Calculate_diabatic_H(
    const double angle) const {
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
  return U.transpose() * ad_energies.asDiagonal() * U;
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
  

  double angle = 0.25 * std::acos(cos4alpha);
  double sin_angle = 0.25 * std::asin(sin4alpha);

 XTP_LOG(Log::debug, *pLog_) << "cos angle " << angle *  57.2958 << flush;
 XTP_LOG(Log::debug, *pLog_) << "sin angle " <<   sin_angle  *  57.2958 << flush;


 XTP_LOG(Log::debug, *pLog_) << "cos(4 alpha) " <<  cos4alpha   << flush;
 XTP_LOG(Log::debug, *pLog_) << "cos(4*sin angle) " <<   std::cos(4. * sin_angle)   << flush;

// testing strange procedure according to original ER paper
double sqroot = std::sqrt( 1.0 - 0.5 * ( 1.0 - cos4alpha ) );
double x_squared_plus =  0.5 * ( 1.0 + sqroot );
double x_squared_minus =  0.5 * ( 1.0 - sqroot );

double x1 = sqrt(x_squared_plus);
double y1 = sqrt(1.0 - x_squared_plus);

double x2 = sqrt(x_squared_minus);
double y2 = sqrt(1.0 - x_squared_minus);

double test1 = 4.0 * x1*y1*(x1*x1-y1*y1);
double test2 = 4.0 * x2*y2*(x2*x2-y2*y2);

double cos_alpha0 = 0;
double sin_alpha0 = 0;

// check which (x,y) pair matches sin(4alpha)
if ( std::abs( test1 - sin4alpha ) < 1e-4 ) {

  cos_alpha0 = x1;
  sin_alpha0 = y1; 

} else if ( std::abs( test2 - sin4alpha ) < 1e-4 ) {

  cos_alpha0 = x2;
  sin_alpha0 = y2; 

} else{

  XTP_LOG(Log::debug, *pLog_) << "Can't find correct angle " << flush;
}


XTP_LOG(Log::debug, *pLog_) << "Coupling element directly with cos_alpha0: "
                              << cos_alpha0 * sin_alpha0 *
                                     (E2_ - E1_) * votca::tools::conv::hrt2ev
                              << flush;




XTP_LOG(Log::debug, *pLog_) << "test1 " << test1 << " test2 " << test2 << " " << " sin(4a) " << sin4alpha  << flush;



  XTP_LOG(Log::debug, *pLog_) << "B12 " << B_12 << flush;
  XTP_LOG(Log::debug, *pLog_) << "A12 " << A_12 << flush;

  XTP_LOG(Log::debug, *pLog_)
      << "angle MAX (degrees) " << angle * 57.2958 << flush;

  XTP_LOG(Log::debug, *pLog_) << "Coupling element directly: "
                              << 0.5 * std::sqrt(0.5 - 0.5 * cos4alpha) *
                                     (E2_ - E1_) * votca::tools::conv::hrt2ev
                              << flush;





  return angle;
}

Eigen::Tensor<double, 4> ERDiabatization::CalculateRtensor() const {
  Eigen::Tensor<double, 4> r_tensor(2, 2, 2, 2);
  for (Index J = 0; J < 2; J++) {
    for (Index K = 0; K < 2; K++) {
      Eigen::MatrixXd D_JK = CalculateD_R(J, K);
      if (!orbitals1_.getTDAApprox() and !orbitals2_.getTDAApprox()) {
        D_JK -= CalculateD_AR(J, K);
      }
      for (Index L = 0; L < 2; L++) {
        for (Index M = 0; M < 2; M++) {
          Eigen::MatrixXd D_LM = CalculateD_R(L, M);
          if (!orbitals1_.getTDAApprox() and !orbitals2_.getTDAApprox()) {
            D_LM -= CalculateD_AR(L, M);
          }
          r_tensor(J, K, L, M) = CalculateR(D_JK, D_LM);
          std::cout << J << " " << K << " " << L << " " << M << " " << r_tensor(J, K, L, M) << "\n" << std::endl;
        }
      }
    }
  }
  return r_tensor;
}

template <bool AR>
Eigen::MatrixXd ERDiabatization::CalculateD(Index stateindex1,
                                            Index stateindex2) const {

  // D matrix depends on 2 indeces. These can be either 0 or 1.
  // Index=0 means "take the first excited state" as Index=1 means "take the
  // second excitate state"
  // This is the reason for this
  Index index1;
  Index index2;

  if (stateindex1 == 0) {
    index1 = opt_.state_idx_1;
  } else if (stateindex1 == 1) {
    index1 = opt_.state_idx_2;
  } else {
    throw std::runtime_error("Invalid state index specified.");
  }
  if (stateindex2 == 0) {
    index2 = opt_.state_idx_1;
  } else if (stateindex2 == 1) {
    index2 = opt_.state_idx_2;
  } else {
    throw std::runtime_error("Invalid state index specified.");
  }

  Eigen::VectorXd exciton1;
  Eigen::VectorXd exciton2;

  // If AR==True we want Bs from BSE. If AR==False we want As from BSE.
  if (AR) {
    if (qmtype_ == QMStateType::Singlet) {
      exciton1 = orbitals1_.BSESinglets().eigenvectors2().col(index1 - 1);
      exciton2 = orbitals2_.BSESinglets().eigenvectors2().col(index2 - 1);
    } else {
      exciton1 = orbitals1_.BSETriplets().eigenvectors2().col(index1 - 1);
      exciton2 = orbitals2_.BSETriplets().eigenvectors2().col(index2 - 1);
    }
  } else {
    if (qmtype_ == QMStateType::Singlet) {
      exciton1 = orbitals1_.BSESinglets().eigenvectors().col(index1 - 1);
      exciton2 = orbitals2_.BSESinglets().eigenvectors().col(index2 - 1);
    } else {
      exciton1 = orbitals1_.BSETriplets().eigenvectors().col(index1 - 1);
      exciton2 = orbitals2_.BSETriplets().eigenvectors().col(index2 - 1);
    }
  }

  // single BSEcoefficient vectors are treated as a    x     matrix
  Eigen::Map<const Eigen::MatrixXd> mat1(exciton1.data(), bse_ctotal_,
                                         bse_vtotal_);

  //std::cout << "mat1 has " << mat1.rows() << " rows and " << mat1.cols() << " columns " << std::endl;  


  Eigen::Map<const Eigen::MatrixXd> mat2(exciton2.data(), bse_ctotal_,
                                         bse_vtotal_);

  Eigen::MatrixXd AuxMat_vv = mat1.transpose() * mat2;

  //std::cout << "AuxMax_vv has " << AuxMat_vv.rows() << " rows and " << AuxMat_vv.cols() << " columns " << std::endl;  


  Eigen::MatrixXd AuxMat_cc = mat1 * mat2.transpose();

  //std::cout << "AuxMax_cc has " << AuxMat_cc.rows() << " rows and " << AuxMat_cc.cols() << " columns " << std::endl;  


  Eigen::MatrixXd results =
      virtlevels1_ * AuxMat_cc * virtlevels2_.transpose() -
      occlevels1_ * AuxMat_vv * occlevels2_.transpose();

  // This adds the GS density
  /*if (stateindex1 == stateindex2) {
    if (stateindex1 == 0) {
      results += orbitals1_.DensityMatrixGroundState();
    }
    if (stateindex1 == 1) {
      results += orbitals2_.DensityMatrixGroundState();
    }
  }*/

  //std::cout << "result has " << results.rows() << " rows and " << results.cols() << " columns " << std::endl;  


  return results;
}

}  // namespace xtp
}  // namespace votca
