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

#pragma once
#ifndef VOTCA_XTP_POPULATIONANALYSIS_H
#define VOTCA_XTP_POPULATIONANALYSIS_H

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmfragment.h>

/**
 * \brief Takes a list of atoms, and the corresponding density matrix and puts
 * out a table of wavefunction partial charges
 *
 *
 *
 */

namespace votca {
namespace xtp {

struct BSE_Population {
  Eigen::VectorXd H;
  Eigen::VectorXd E;
  double Gs = 0;

  void Initialize(int size) {
    H = Eigen::VectorXd::Zero(size);
    E = Eigen::VectorXd::Zero(size);
    Gs = 0;
  }
  friend std::ostream& operator<<(std::ostream& out,
                                  const BSE_Population& pop) {
    if (pop.H.size() < 1) {
      return out;
    }
    Eigen::VectorXd diff = pop.H - pop.E;
    out << "GroundstateCharge:" << pop.Gs << std::endl;
    out << "Index hole electron dQ Qeff" << std::endl;
    for (int i = 0; i < pop.H.size(); ++i) {
      out << i << " " << pop.H(i) << " " << pop.E(i) << " " << diff(i) << " "
          << diff(i) + pop.Gs << std::endl;
    }
    return out;
  }
};

template <bool T>
class Populationanalysis {
 public:
  void CalcChargeperAtom(Orbitals& orbitals, const QMState& state) const {

    AOBasis basis = orbitals.SetupDftBasis();
    Eigen::MatrixXd dmat = orbitals.DensityMatrixFull(state);
    AOOverlap overlap;
    overlap.Fill(basis);
    Eigen::VectorXd charges = CalcElecChargeperAtom(dmat, overlap, basis);
    if (!state.isTransition()) {
      charges += CalcNucChargeperAtom(orbitals.QMAtoms());
    }

    StaticSegment seg =
        StaticSegment(orbitals.QMAtoms().getName(), orbitals.QMAtoms().getId());
    for (int i = 0; i < orbitals.QMAtoms().size(); ++i) {
      seg.push_back(StaticSite(orbitals.QMAtoms()[i], charges(i)));
    }
    orbitals.Multipoles() = seg;
    return;
  }

  void CalcChargeperFragment(std::vector<QMFragment<BSE_Population> >& frags,
                             const Orbitals& orbitals, QMStateType type) {
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
      std::vector<Eigen::MatrixXd> dmat_ex =
          orbitals.DensityMatrixExcitedState(state);
      Eigen::VectorXd atom_h =
          CalcElecChargeperAtom(dmat_ex[0], overlap, basis);
      Eigen::VectorXd atom_e =
          CalcElecChargeperAtom(dmat_ex[1], overlap, basis);
      for (auto& frag : frags) {
        frag.value().E(i_state) = frag.ExtractFromVector(atom_e);
        frag.value().H(i_state) = frag.ExtractFromVector(atom_h);
      }
    }
  }

 private:
  Eigen::VectorXd CalcNucChargeperAtom(const QMMolecule& mol) const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(mol.size());
    for (int i = 0; i < mol.size(); i++) {
      result(i) = mol[i].getNuccharge();
    }
    return result;
  }

  Eigen::VectorXd CalcElecChargeperAtom(const Eigen::MatrixXd& dmat,
                                        AOOverlap& overlap,
                                        const AOBasis& basis) const {
    Eigen::MatrixXd prodmat;
    if (T) {
      Eigen::MatrixXd Smsqrt = overlap.Sqrt();
      prodmat = Smsqrt * dmat * Smsqrt;
    } else {
      prodmat = dmat * overlap.Matrix();
    }
    int noofatoms = basis.getFuncPerAtom().size();
    Eigen::VectorXd charges = Eigen::VectorXd::Zero(noofatoms);
    int start = 0;
    for (int i = 0; i < charges.size(); ++i) {
      int nofunc = basis.getFuncPerAtom()[i];
      charges(i) = prodmat.diagonal().segment(start, nofunc).sum();
      start += nofunc;
    }
    return charges;
  }
};

typedef Populationanalysis<false> Mulliken;
typedef Populationanalysis<true> Lowdin;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_POPULATIONANALYSIS_H
