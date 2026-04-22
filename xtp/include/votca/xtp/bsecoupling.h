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
#ifndef VOTCA_XTP_BSECOUPLING_H
#define VOTCA_XTP_BSECOUPLING_H

// Local VOTCA includes
#include "couplingbase.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

/**
 * \brief Evaluates electronic coupling elements
 *
 * J. Wehner,B. Baumeier,
 * JCTC DOI: 10.1021/acs.jctc.6b00935
 *
 */

class BSECoupling : public CouplingBase {
 public:
  void Initialize(tools::Property& options) override;
  std::string Identify() const { return "bsecoupling"; }

  void Addoutput(tools::Property& type_summary, const Orbitals& orbitalsA,
                 const Orbitals& orbitalsB) const override;

  /**
   * \brief evaluates electronic couplings
   *
   * @param  orbitalsA molecular orbitals of molecule A
   * @param  orbitalsB molecular orbitals of molecule B
   * @param  orbitalsAB molecular orbitals of the dimer AB
   */
  void CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                          const Orbitals& orbitalsAB) override;

 private:
  /**
   * \brief Diagnostic quantities for assessing whether CT downfolding is safe.
   *
   * xi:                  max over FE-CT pairs of |H_FE_CT| / |E_FE - E_CT|.
   *                      Small xi (< ~0.3) means perturbative downfolding is
   *                      reliable.
   * pt_rm_discrepancy:   max |J_pert - J_diag| over all FE pairs [Hrt].
   *                      Large discrepancy flags near-resonant CT states.
   * downfolding_safe:    true if both xi and pt_rm_discrepancy are below
   *                      their respective thresholds.
   */
  struct Diagnostics {
    double xi = 0.0;
    double pt_rm_discrepancy = 0.0;
    bool downfolding_safe = true;
  };

  void WriteToProperty(tools::Property& summary, const QMState& stateA,
                       const QMState& stateB) const;

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

  double getSingletCouplingElement(Index levelA, Index levelB,
                                   Index methodindex) const;

  double getTripletCouplingElement(Index levelA, Index levelB,
                                   Index methodindex) const;

  Eigen::MatrixXd SetupCTStates(Index bseA_vtotal, Index bseB_vtotal,
                                Index bseAB_vtotal, Index bseAB_ctotal,
                                const Eigen::MatrixXd& A_AB,
                                const Eigen::MatrixXd& B_AB) const;

  Eigen::MatrixXd ProjectFrenkelExcitons(const Eigen::MatrixXd& BSE_Coeffs,
                                         const Eigen::MatrixXd& X_AB,
                                         Index bseX_vtotal, Index bseX_ctotal,
                                         Index bseAB_vtotal,
                                         Index bseAB_ctotal) const;

  /**
   * \brief Compute J_pert and J_diag, and store raw J_dimer/S_dimer and
   *        diagnostics for later output.
   *
   * @param FE_AB       projected FE states in dimer basis
   * @param CTStates    CT states in dimer basis
   * @param H           BSE Hamiltonian operator
   * @param J_dimer_out raw (pre-Lowdin) Hamiltonian matrix [out]
   * @param S_dimer_out raw overlap matrix [out]
   * @param diag_out    diagnostics struct [out]
   * @return            array of {J_pert, J_diag} matrices
   */
  template <class BSE_OPERATOR>
  std::array<Eigen::MatrixXd, 2> ProjectExcitons(
      Eigen::MatrixXd& FE_AB, Eigen::MatrixXd& CTStates, BSE_OPERATOR H,
      Eigen::MatrixXd& J_dimer_out, Eigen::MatrixXd& S_dimer_out,
      Diagnostics& diag_out) const;

  /**
   * \brief Form J_dimer and S_dimer from the projection, then Lowdin
   *        orthogonalize to produce J_ortho.
   *
   * @param H           BSE Hamiltonian operator
   * @param projection  merged FE+CT projection matrix (consumed)
   * @param J_dimer_out raw Hamiltonian matrix before Lowdin [out]
   * @param S_dimer_out raw overlap matrix before Lowdin [out]
   * @return            Lowdin-orthogonalized Hamiltonian J_ortho
   */
  template <class BSE_OPERATOR>
  Eigen::MatrixXd CalcJ_dimer(BSE_OPERATOR& H, Eigen::MatrixXd& projection,
                              Eigen::MatrixXd& J_dimer_out,
                              Eigen::MatrixXd& S_dimer_out) const;

  /**
   * \brief Merge FE and CT projection vectors into a single projection matrix
   *        without pre-orthogonalization.
   *
   * The non-orthogonality between FE and CT states is handled correctly
   * downstream by the joint Lowdin in CalcJ_dimer (for the reduction method)
   * or by the generalized eigenvalue problem / SRG (for TB use of J_dimer
   * and S_dimer directly).
   *
   * Pre-orthogonalizing CT states against FE states (previous behaviour) used
   * the incorrect projector P = F*F^T, which is only exact when FE_AB columns
   * are orthonormal. This introduced a systematic error that grows with the
   * number of CT states included. The function name is retained for interface
   * compatibility but the orthogonalization step has been removed.
   */
  Eigen::MatrixXd OrthogonalizeCTs(Eigen::MatrixXd& FE_AB,
                                   Eigen::MatrixXd& CTStates) const;

  Eigen::MatrixXd Fulldiag(const Eigen::MatrixXd& J_dimer) const;

  Eigen::MatrixXd Perturbation(const Eigen::MatrixXd& J_dimer) const;

  /**
   * \brief Compute diagnostics from raw J_dimer and the two effective
   *        coupling matrices.
   *
   * @param J_dimer  raw (pre-Lowdin) Hamiltonian in FE+CT basis
   * @param J_pert   effective couplings from perturbation theory
   * @param J_diag   effective couplings from reduction method
   * @return         populated Diagnostics struct
   */
  Diagnostics ComputeDiagnostics(const Eigen::MatrixXd& J_dimer,
                                 const Eigen::MatrixXd& J_pert,
                                 const Eigen::MatrixXd& J_diag) const;

  // --- effective coupling results (existing) ---
  std::array<Eigen::MatrixXd, 2> JAB_singlet;  // [J_pert, J_diag]
  std::array<Eigen::MatrixXd, 2> JAB_triplet;

  // --- raw TB matrices (new) ---
  // These are the pre-Lowdin J and S in the full FE+CT basis,
  // suitable for assembling a TB Hamiltonian or feeding into SRG.
  Eigen::MatrixXd J_dimer_singlet_;
  Eigen::MatrixXd S_dimer_singlet_;
  Eigen::MatrixXd J_dimer_triplet_;
  Eigen::MatrixXd S_dimer_triplet_;

  // --- diagnostics (new) ---
  Diagnostics diag_singlet_;
  Diagnostics diag_triplet_;

  // --- monomer FE energies (new) ---
  // Isolated monomer excitation energies in eV, stored per spin channel.
  // These provide pairwise-consistent site energies for TB assembly:
  // the same monomer calculation is used regardless of which dimer pair
  // is being processed, unlike the dimer diagonal H_FE_FE[i,i] which
  // varies with the partner. Environmental corrections (electrostatic
  // embedding) should be added separately on top of these values.
  Eigen::VectorXd monomerA_energies_singlet_;
  Eigen::VectorXd monomerB_energies_singlet_;
  Eigen::VectorXd monomerA_energies_triplet_;
  Eigen::VectorXd monomerB_energies_triplet_;

  bool doTriplets_;
  bool doSinglets_;
  bool output_perturbation_;
  Index levA_;
  Index levB_;
  Index occA_;
  Index unoccA_;
  Index occB_;
  Index unoccB_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSECOUPLING_H
