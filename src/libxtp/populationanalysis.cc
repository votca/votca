/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/xtp/populationanalysis.h>

namespace votca {
namespace xtp {

template <bool T>
StaticSegment Populationanalysis<T>::CalcChargeperAtom(
    const Orbitals& orbitals, const QMState& state) const {

  AOBasis basis = orbitals.SetupDftBasis();
  Eigen::MatrixXd dmat = orbitals.DensityMatrixFull(state);
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::VectorXd charges = -CalcElecChargeperAtom(dmat, overlap, basis);
  if (!state.isTransition()) {
    charges += CalcNucChargeperAtom(orbitals.QMAtoms());
  }

  StaticSegment seg =
      StaticSegment(orbitals.QMAtoms().getType(), orbitals.QMAtoms().getId());
  for (int i = 0; i < orbitals.QMAtoms().size(); ++i) {
    seg.push_back(StaticSite(orbitals.QMAtoms()[i], charges(i)));
  }
  return seg;
}

template <bool T>
void Populationanalysis<T>::CalcChargeperFragment(
    std::vector<QMFragment<BSE_Population> >& frags, const Orbitals& orbitals,
    QMStateType type) const {
  if (!type.isExciton()) {
    throw std::runtime_error(
        "CalcChargeperFragment: QmStateType must be an exciton");
  }
  AOBasis basis = orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::VectorXd nuccharges = CalcNucChargeperAtom(orbitals.QMAtoms());
  Eigen::MatrixXd dmatgs = orbitals.DensityMatrixGroundState();
  Eigen::VectorXd gscharges =
      nuccharges - CalcElecChargeperAtom(dmatgs, overlap, basis);
  int numofstates = orbitals.NumberofStates(type);
  for (auto& frag : frags) {
    frag.value().Initialize(numofstates);
    frag.value().Gs = frag.ExtractFromVector(gscharges);
  }
  for (int i_state = 0; i_state < numofstates; i_state++) {
    QMState state(type, i_state, false);
    std::array<Eigen::MatrixXd, 2> dmat_ex =
        orbitals.DensityMatrixExcitedState(state);
    Eigen::VectorXd atom_h = CalcElecChargeperAtom(dmat_ex[0], overlap, basis);
    Eigen::VectorXd atom_e = -CalcElecChargeperAtom(dmat_ex[1], overlap, basis);
    for (auto& frag : frags) {
      frag.value().E(i_state) = frag.ExtractFromVector(atom_e);
      frag.value().H(i_state) = frag.ExtractFromVector(atom_h);
    }
  }
}

template <bool T>
Eigen::VectorXd Populationanalysis<T>::CalcNucChargeperAtom(
    const QMMolecule& mol) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(mol.size());
  for (long i = 0; i < mol.size(); i++) {
    result(i) = mol[i].getNuccharge();
  }
  return result;
}

template <bool T>
Eigen::VectorXd Populationanalysis<T>::CalcElecChargeperAtom(
    const Eigen::MatrixXd& dmat, AOOverlap& overlap,
    const AOBasis& basis) const {
  Eigen::MatrixXd prodmat;
  if (T) {
    Eigen::MatrixXd Smsqrt = overlap.Sqrt();
    prodmat = Smsqrt * dmat * Smsqrt;
  } else {
    prodmat = dmat * overlap.Matrix();
  }
  std::vector<int> funcperatom = basis.getFuncPerAtom();
  long noofatoms = long(funcperatom.size());
  Eigen::VectorXd charges = Eigen::VectorXd::Zero(noofatoms);
  int start = 0;
  for (long i = 0; i < charges.size(); ++i) {
    int nofunc = funcperatom[i];
    charges(i) = prodmat.diagonal().segment(start, nofunc).sum();
    start += nofunc;
  }
  return charges;
}

template class Populationanalysis<true>;
template class Populationanalysis<false>;

}  // namespace xtp
}  // namespace votca
