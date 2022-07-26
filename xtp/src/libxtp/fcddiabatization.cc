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
#include "votca/xtp/fcddiabatization.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/populationanalysis.h"
#include "votca/xtp/transition_densities.h"
#include <votca/tools/constants.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void FCDDiabatization::configure() {

  // Check on dftbasis
  if (orbitals1_.getDFTbasisName() != orbitals2_.getDFTbasisName()) {
    throw std::runtime_error("Different DFT basis for the two input file.");
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

  qmtype_.FromString(qmstate_str_);

  if (qmtype_ == QMStateType::Singlet) {
    E1_ = orbitals1_.BSESinglets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSESinglets().eigenvalues()[state_idx_2_ - 1];
  } else {
    E1_ = orbitals1_.BSETriplets().eigenvalues()[state_idx_1_ - 1];
    E2_ = orbitals2_.BSETriplets().eigenvalues()[state_idx_2_ - 1];
  }
}

double FCDDiabatization::calculate_coupling() {

  // calculate fragment charge differences in state 1
  QMState state1 = QMState(qmtype_, state_idx_1_ - 1, false);
  Lowdin pops;
  pops.CalcChargeperFragment(fragments_, orbitals1_, state1.Type());

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Fragment 1 ground state "
                              << fragments_[0].value().Gs << flush;
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Fragment 2 ground state "
                              << fragments_[1].value().Gs << flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 1 poulation electron "
      << fragments_[0].value().E(state_idx_1_ - 1) << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 1 poulation hole "
      << fragments_[0].value().H(state_idx_1_ - 1) << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 2 poulation electron "
      << fragments_[1].value().E(state_idx_1_ - 1) << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 2 poulation hole "
      << fragments_[1].value().H(state_idx_1_ - 1) << flush;

  double dQ11 = fragments_[0].value().H(state_idx_1_ - 1) +
                fragments_[0].value().E(state_idx_1_ - 1) -
                fragments_[1].value().H(state_idx_1_ - 1) -
                fragments_[1].value().E(state_idx_1_ - 1);

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " dQ11 " << dQ11 << flush;

  // calculate fragment charge differences in state 2
  QMState state2 = QMState(qmtype_, state_idx_2_ - 1, false);
  pops.CalcChargeperFragment(fragments_, orbitals2_, state2.Type());
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 1 poulation electron "
      << fragments_[0].value().E(state_idx_2_ - 1) << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Fragment 1 poulation hole "
      << fragments_[0].value().H(state_idx_2_ - 1) << flush;

  double dQ22 = fragments_[0].value().H(state_idx_2_ - 1) +
                fragments_[0].value().E(state_idx_2_ - 1) -
                fragments_[1].value().H(state_idx_2_ - 1) -
                fragments_[1].value().E(state_idx_2_ - 1);

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " dQ22 " << dQ22 << flush;

  // calculate fragment charge difference from 1 to 2 transition density
  TransitionDensities tdmat(orbitals1_, orbitals2_, pLog_);
  tdmat.configure();
  Eigen::MatrixXd tmat = tdmat.Matrix(state1, state2);
  pops.CalcChargeperFragmentTransition(fragments_, orbitals1_, tmat);
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Fragment 1 transition "
                              << fragments_[0].value().Gs << flush;
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Fragment 2 transition "
                              << fragments_[1].value().Gs << flush;

  double dQ12 = fragments_[0].value().Gs - fragments_[1].value().Gs;

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " dQ12 " << dQ12 << flush;

  // Fragment Charge Difference coupling
  double coupling =
      (std::abs(dQ12) * (E2_ - E1_)) /
      std::sqrt(std::pow(dQ11 - dQ22, 2) + 4.0 * std::pow(dQ12, 2));

  return coupling;
}

}  // namespace xtp
}  // namespace votca
