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

#pragma once
#ifndef VOTCA_XTP_DFTCOUPLING_H
#define VOTCA_XTP_DFTCOUPLING_H

// Local VOTCA includes
#include "couplingbase.h"

namespace votca {
namespace xtp {

/**
 * \brief Evaluates electronic coupling elements
 *
 * B. Baumeier, J. Kirkpatrick, D. Andrienko,
 * Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
 *
 */

class DFTcoupling : public CouplingBase {
 public:
  std::string Identify() const { return "dftcoupling"; }

  void CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                          const Orbitals& orbitalsAB) override;

  void Initialize(tools::Property&) override;

  void Addoutput(tools::Property& type_summary, const Orbitals& orbitalsA,
                 const Orbitals& orbitalsB) const override;

 private:
  void WriteToProperty(tools::Property& type_summary, const Orbitals& orbitalsA,
                       const Orbitals& orbitalsB, Index a, Index b) const;

  double getCouplingElement(Index levelA, Index levelB,
                            const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) const;

  std::pair<Index, Index> DetermineRangeOfStates(const Orbitals& orbital,
                                                 Index numberofstates) const;

  /**
   * \brief Write an Eigen matrix as rows of space-separated values into an
   *        XML property node. Each row becomes an attribute "row_N".
   *
   * @param prop        parent Property node to attach the matrix node to
   * @param name        name of the new child node
   * @param mat         matrix to write
   * @param conversion  optional unit conversion factor (default 1.0)
   */
  static void WriteMatrixToProperty(tools::Property& prop,
                                    const std::string& name,
                                    const Eigen::MatrixXd& mat,
                                    double conversion = 1.0);

  // --- existing: Lowdin-orthogonalized effective coupling matrix [Hrt] ---
  // JAB(i,j) where i in [0, levelsA) indexes fragment A MOs and
  // j in [levelsA, levelsA+levelsB) indexes fragment B MOs.
  // Both hole and electron blocks share this matrix, partitioned by range.
  Eigen::MatrixXd JAB;

  // --- new: raw (pre-Lowdin) TB matrices [Hrt for H, dimensionless for S] ---
  // Stored separately for holes and electrons since they are computed from
  // different MO ranges and assembled into the TB Hamiltonian independently.
  //
  // Block structure of each matrix (hole or electron):
  //   H = [ H_AA  H_AB ]    S = [ S_AA  S_AB ]
  //       [ H_BA  H_BB ]        [ S_BA  S_BB ]
  // where A = fragment A MOs, B = fragment B MOs for that carrier type.
  // H is in Hartree; conversion to eV happens in Addoutput.
  Eigen::MatrixXd JAB_dimer_hole_;  // pre-Lowdin H, hole block [Hrt]
  Eigen::MatrixXd S_AxB_hole_;      // overlap, hole block
  Eigen::MatrixXd JAB_dimer_elec_;  // pre-Lowdin H, electron block [Hrt]
  Eigen::MatrixXd S_AxB_elec_;      // overlap, electron block

  // --- new: monomer MO energies [eV] ---
  // Isolated monomer quasiparticle energies for the relevant MO ranges.
  // These provide pairwise-consistent TB site energies: the same monomer
  // calculation is used regardless of which dimer partner is processed,
  // unlike the dimer Hamiltonian diagonal which depends on the partner.
  // Environmental corrections (e.g. from polarisation embedding) should
  // be added on top of these values externally.
  //
  // KS: Kohn-Sham DFT eigenvalues, always available.
  // QP: perturbative GW quasiparticle energies, present only when a GW
  //     calculation was run on the monomer. Stored as empty vector otherwise.
  //     Indexed identically to the KS vectors — same MO range, same ordering.
  Eigen::VectorXd moEnergiesA_hole_KS_;  // fragment A hole MOs, KS [eV]
  Eigen::VectorXd moEnergiesB_hole_KS_;
  Eigen::VectorXd moEnergiesA_elec_KS_;  // fragment A electron MOs, KS [eV]
  Eigen::VectorXd moEnergiesB_elec_KS_;

  Eigen::VectorXd moEnergiesA_hole_QP_;  // fragment A hole MOs, QPpert [eV]
  Eigen::VectorXd moEnergiesB_hole_QP_;  // empty if GW not available
  Eigen::VectorXd moEnergiesA_elec_QP_;
  Eigen::VectorXd moEnergiesB_elec_QP_;

  // --- new: basis diagnostics ---
  // Minimum eigenvalue of the overlap matrix S_AxB for each carrier type.
  // Small values indicate near-linear dependence in the MO projection basis,
  // which would make the Lowdin orthogonalization ill-conditioned.
  double min_S_eigenvalue_hole_ = 1.0;
  double min_S_eigenvalue_elec_ = 1.0;

  double degeneracy_ = 0.0;
  Index numberofstatesA_ = 1;
  Index numberofstatesB_ = 1;
  // When true, write monomer MO energies, raw H/S TB matrices, and
  // diagnostics to the XML output. Set false (default) for KMC/rate
  // workflows that only need scalar effective couplings.
  bool output_tb_ = false;

  std::pair<Index, Index> Range_orbA;
  std::pair<Index, Index> Range_orbB;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTCOUPLING_H
