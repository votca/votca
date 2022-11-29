/*
 *            Copyright 2009-2022 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/transition_densities.h"

namespace votca {
namespace xtp {

void TransitionDensities::configure() {

  // Check on dftbasis
  if (orbitals1_.getDFTbasisName() != orbitals2_.getDFTbasisName()) {
    throw std::runtime_error("Different DFT basis for the two input file.");
  } else {
    dftbasis_ = orbitals1_.getDftBasis();
    basissize_ = orbitals1_.getBasisSetSize();
    XTP_LOG(Log::error, *log_)
        << TimeStamp() << " Data was created with basis set "
        << orbitals1_.getDFTbasisName() << std::flush;
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
}

Eigen::MatrixXd TransitionDensities::Matrix(QMState state1, QMState state2) {

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

  Eigen::MatrixXd transition_dmat =
      virtlevels1_ * AuxMat_cc * virtlevels2_.transpose() -
      occlevels1_ * AuxMat_vv * occlevels2_.transpose();

  return transition_dmat;
}

}  // namespace xtp
}  // namespace votca
