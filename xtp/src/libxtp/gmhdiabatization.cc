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
#include "votca/xtp/gmhdiabatization.h"
#include "votca/xtp/aomatrix.h"
#include <votca/tools/constants.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void GMHDiabatization::configure() {

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

  qmtype_.FromString(qmstate_str_);

  if (qmtype_ == QMStateType::Singlet) {
    E1_ = orbitals1_.BSESinglets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSESinglets().eigenvalues()[state_idx_2_ - 1];
  } else {
    E1_ = orbitals1_.BSETriplets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSETriplets().eigenvalues()[state_idx_2_ - 1];
  }
}

std::pair<double, double> GMHDiabatization::calculate_coupling() {

  // calculate electric dipole moment in state 1
  QMState state1 = QMState(qmtype_, state_idx_1_ - 1, false);
  Eigen::Vector3d state1_dipole = orbitals1_.CalcElDipole(state1);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Dipole Moment excited state 1 " << state1_dipole[0]
      << "\t" << state1_dipole[1] << "\t" << state1_dipole[2] << flush;

  // calculate electric dipole moment in state 2
  QMState state2 = QMState(qmtype_, state_idx_2_ - 1, false);
  Eigen::Vector3d state2_dipole = orbitals2_.CalcElDipole(state2);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Dipole Moment excited state 2 " << state2_dipole[0]
      << "\t" << state2_dipole[1] << "\t" << state2_dipole[2] << flush;

  // calculate transition dipole moment from 1 to 2
  Eigen::Vector3d transition_dip = transition_dipole(state1, state2);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Transition Dipole Moment " << transition_dip[0]
      << "\t" << transition_dip[1] << "\t" << transition_dip[2] << flush;

  // Generalized Mulliken-Hush coupling
  double abs_transition_dip = transition_dip.norm();
  double dipole_diff_norm = (state1_dipole - state2_dipole).norm();
  double coupling = (abs_transition_dip * (E2_ - E1_)) /
                    std::sqrt(std::pow(dipole_diff_norm, 2) +
                              4.0 * std::pow(abs_transition_dip, 2));

  // with the "projection" from Stanford people
  Eigen::Vector3d CT_direction =
      (state1_dipole - state2_dipole) / (state1_dipole - state2_dipole).norm();
  double proj_transition_dip = transition_dip.dot(CT_direction);
  double coupling_proj = (proj_transition_dip * (E2_ - E1_)) /
                         std::sqrt(std::pow(dipole_diff_norm, 2) +
                                   4.0 * std::pow(proj_transition_dip, 2));

  return std::pair<double, double>(coupling, coupling_proj);
}

Eigen::Vector3d GMHDiabatization::transition_dipole(QMState state1,
                                                    QMState state2) {

  // get the resonant part of exciton1
  Eigen::VectorXd BSEcoeffs1_res =
      orbitals1_.BSESinglets().eigenvectors().col(state1.StateIdx());

  // get the resonant part of exciton2
  Eigen::VectorXd BSEcoeffs2_res =
      orbitals2_.BSESinglets().eigenvectors().col(state2.StateIdx());

  // view BSEcoeffs as matrix
  Eigen::Map<const Eigen::MatrixXd> exciton1_res(BSEcoeffs1_res.data(),
                                                 bse_ctotal_, bse_vtotal_);

  Eigen::Map<const Eigen::MatrixXd> exciton2_res(BSEcoeffs2_res.data(),
                                                 bse_ctotal_, bse_vtotal_);

  Eigen::MatrixXd AuxMat_vv = exciton1_res.transpose() * exciton2_res;
  Eigen::MatrixXd AuxMat_cc = exciton1_res * exciton2_res.transpose();

  // subtract the antiresonant part, if applicable
  if (!orbitals1_.getTDAApprox()) {
    Eigen::VectorXd BSEcoeffs1_antires =
        orbitals1_.BSESinglets().eigenvectors2().col(state1.StateIdx());

    Eigen::Map<const Eigen::MatrixXd> exciton1_antires(
        BSEcoeffs1_antires.data(), bse_ctotal_, bse_vtotal_);

    Eigen::VectorXd BSEcoeffs2_antires =
        orbitals2_.BSESinglets().eigenvectors2().col(state1.StateIdx());

    Eigen::Map<const Eigen::MatrixXd> exciton2_antires(
        BSEcoeffs2_antires.data(), bse_ctotal_, bse_vtotal_);

    AuxMat_vv -= exciton1_antires.transpose() * exciton2_antires;
    AuxMat_cc -= exciton1_antires * exciton2_antires.transpose();
  }

  AOBasis basis = orbitals1_.getDftBasis();
  AODipole dipole;
  dipole.setCenter(orbitals1_.QMAtoms().getPos());
  dipole.Fill(basis);
  Eigen::MatrixXd transition_dmat =
      virtlevels1_ * AuxMat_cc * virtlevels2_.transpose() -
      occlevels1_ * AuxMat_vv * occlevels2_.transpose();

  AOOverlap overlap;
  overlap.Fill(basis);
  double nelec = transition_dmat.cwiseProduct(overlap.Matrix()).sum();

  Eigen::Vector3d electronic_dip;
  for (Index i = 0; i < 3; ++i) {
    electronic_dip(i) = transition_dmat.cwiseProduct(dipole.Matrix()[i]).sum();
  }

  return electronic_dip;
}

}  // namespace xtp
}  // namespace votca
