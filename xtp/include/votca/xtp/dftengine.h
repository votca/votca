/*
 *            Copyright 2009-2023 The VOTCA Development Team
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
#ifndef VOTCA_XTP_DFTENGINE_H
#define VOTCA_XTP_DFTENGINE_H

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "ERIs.h"
#include "convergenceacc.h"
#include "uks_convergenceacc.h"

#include "ecpaobasis.h"
#include "extended_hueckel.h"
#include "logger.h"
#include "qmmolecule.h"
#include "staticsite.h"
#include "vxc_grid.h"
#include "vxc_potential.h"
namespace votca {
namespace xtp {
class Orbitals;
class DFTEngineTestAccess;

/**
 * \brief Electronic ground-state via Density-Functional Theory.
 *
 * This class assembles the one- and two-electron matrix contributions needed
 * for self-consistent Kohn-Sham calculations in a Gaussian AO basis. The SCF
 * machinery supports restricted closed-shell calculations and unrestricted
 * spin-polarized calculations, while reusing the same integral and numerical
 * XC infrastructure whenever possible.
 */

class DFTEngine {
 public:
  /// Read DFT, grid, and SCF settings from the user options tree.
  void Initialize(tools::Property& options);

  /// Attach the logger used for SCF progress and diagnostics.
  void setLogger(Logger* pLog) { pLog_ = pLog; }

  /// Provide external static sites whose electrostatic potential enters the
  /// Hamiltonian.
  void setExternalcharges(
      std::vector<std::unique_ptr<StaticSite> >* externalsites) {
    externalsites_ = externalsites;
  }

  /// Run a full ground-state DFT calculation and store the results in the
  /// orbital container.
  bool Evaluate(Orbitals& orb);

  /// Run an embedded active-region DFT calculation for the supplied orbital
  /// container.
  bool EvaluateActiveRegion(Orbitals& orb);
  /// Run the truncated active-region workflow used for reduced embedded
  /// calculations.
  bool EvaluateTruncatedActiveRegion(Orbitals& trunc_orb);

  /// Return the configured AO basis-set name for the DFT calculation.
  std::string getDFTBasisName() const { return dftbasis_name_; };

  /// Spin-resolved density matrices used to pass alpha and beta quantities
  /// together through the UKS workflow.
  struct SpinDensity {
    Eigen::MatrixXd alpha;
    Eigen::MatrixXd beta;

    /// Return the total density P = P^alpha + P^beta.
    Eigen::MatrixXd total() const { return alpha + beta; }

    /// Return the spin density P^alpha - P^beta.
    Eigen::MatrixXd spin() const { return alpha - beta; }
  };

  /// Report whether the current electron counts define a spin-polarized
  /// reference.
  bool IsRestrictedOpenShell() const {
    return num_alpha_electrons_ != num_beta_electrons_;
  }

  /// Return the number of spatial orbitals occupied by at least one electron
  /// in a restricted open-shell reference.
  Index NumberOfRestrictedOccupiedOrbitals() const {
    return std::max(num_alpha_electrons_, num_beta_electrons_);
  }

  /// Construct alpha and beta density matrices from a shared MO coefficient
  /// matrix and the current occupation metadata.
  SpinDensity BuildSpinDensity(const tools::EigenSystem& MOs) const;

  /// Assemble SCF acceleration settings consistent with the current spin
  /// treatment and occupation model.
  ConvergenceAcc::options BuildConvergenceOptions() const;

 private:
  friend class DFTEngineTestAccess;

  /// Initialize basis sets, integral engines, and electron counts before
  /// entering the SCF loop.
  void Prepare(Orbitals& orb, Index numofelectrons = -1);

  /// Build the numerical exchange-correlation integration object for the
  /// current molecule.
  Vxc_Potential<Vxc_Grid> SetupVxc(const QMMolecule& mol);

  /// Orthonormalize an initial MO guess with respect to the AO overlap matrix.
  Eigen::MatrixXd OrthogonalizeGuess(const Eigen::MatrixXd& GuessMOs) const;
  /// Print a one-spin list of orbital energies and occupations to the logger.
  void PrintMOs(const Eigen::VectorXd& MOEnergies, Log::Level level);
  /// Print separate alpha and beta orbital energies for a UKS calculation.
  void PrintMOsUKS(const Eigen::VectorXd& alpha_energies,
                   const Eigen::VectorXd& beta_energies,
                   Log::Level level) const;
  /// Evaluate and print the electronic dipole moment from the final density.
  void CalcElDipole(const Orbitals& orb) const;

  /// Build Coulomb and exact-exchange matrix contributions from the current MO
  /// coefficients and density.
  std::array<Eigen::MatrixXd, 2> CalcERIs_EXX(const Eigen::MatrixXd& MOCoeff,
                                              const Eigen::MatrixXd& Dmat,
                                              double error) const;

  /// Build the Coulomb matrix contribution from the current density matrix.
  Eigen::MatrixXd CalcERIs(const Eigen::MatrixXd& Dmat, double error) const;

  /// Propagate basis-set, XC, and metadata settings into the orbital container.
  void ConfigOrbfile(Orbitals& orb);
  /// Precompute AO matrices that remain unchanged throughout the SCF procedure.
  void SetupInvariantMatrices();
  /// Apply McWeeny purification to improve the idempotency of a density-matrix
  /// guess.
  Eigen::MatrixXd McWeenyPurification(Eigen::MatrixXd& Dmat_in,
                                      AOOverlap& overlap);

  /// Assemble the one-electron core Hamiltonian for the current molecule.
  Mat_p_Energy SetupH0(const QMMolecule& mol) const;
  /// Integrate the electrostatic potential generated by external multipoles
  /// into the AO basis.
  Mat_p_Energy IntegrateExternalMultipoles(
      const QMMolecule& mol,
      const std::vector<std::unique_ptr<StaticSite> >& multipoles) const;
  /// Integrate an external electron density represented by another orbital
  /// container.
  Mat_p_Energy IntegrateExternalDensity(const QMMolecule& mol,
                                        const Orbitals& extdensity) const;

  /// Integrate a homogeneous external electric field into the AO basis.
  Eigen::MatrixXd IntegrateExternalField(const QMMolecule& mol) const;

  /// Generate an initial guess by diagonalizing the core Hamiltonian only.
  tools::EigenSystem IndependentElectronGuess(const Mat_p_Energy& H0) const;
  /// Generate an initial guess from a model potential including numerical XC
  /// contributions.
  tools::EigenSystem ModelPotentialGuess(
      const Mat_p_Energy& H0, const QMMolecule& mol,
      const Vxc_Potential<Vxc_Grid>& vxcpotential) const;

  /// Build an atomic-density based initial guess in the AO basis.
  Eigen::MatrixXd AtomicGuess(const QMMolecule& mol) const;

  /// Build orbital energies used in the extended-Hückel starting guess.
  Eigen::VectorXd BuildEHTOrbitalEnergies(const QMMolecule& mol) const;
  /// Build the extended-Hückel Hamiltonian for the current molecule.
  Eigen::MatrixXd BuildEHTHamiltonian(const QMMolecule& mol) const;
  /// Generate an initial guess by diagonalizing the extended-Hückel
  /// Hamiltonian.
  tools::EigenSystem ExtendedHuckelGuess(const QMMolecule& mol) const;
  /// Generate an extended-Hückel based guess refined with the one-electron DFT
  /// Hamiltonian.
  tools::EigenSystem ExtendedHuckelDFTGuess(
      const Mat_p_Energy& H0, const QMMolecule& mol,
      const Vxc_Potential<Vxc_Grid>& vxcpotential) const;
  /// Run an unrestricted atomic reference calculation used in open-shell atomic
  /// guesses.
  Eigen::MatrixXd RunAtomicDFT_unrestricted(const QMAtom& uniqueAtom) const;

  /// Compute the classical nucleus-nucleus repulsion energy.
  double NuclearRepulsion(const QMMolecule& mol) const;
  /// Compute the classical interaction energy between nuclei and external
  /// multipoles.
  double ExternalRepulsion(
      const QMMolecule& mol,
      const std::vector<std::unique_ptr<StaticSite> >& multipoles) const;
  /// Average density-matrix elements over functions belonging to the same
  /// atomic shell.
  Eigen::MatrixXd SphericalAverageShells(const Eigen::MatrixXd& dmat,
                                         const AOBasis& dftbasis) const;

  /// Project the full-system Hamiltonian and densities onto the selected
  /// truncated active basis.
  void TruncateBasis(Orbitals& orb, std::vector<Index>& activeatoms,
                     Mat_p_Energy& H0,
                     Eigen::MatrixXd InitialActiveDensityMatrix,
                     Eigen::MatrixXd v_embedding,
                     Eigen::MatrixXd InitialInactiveMOs);

  /// Expand truncated active-region orbitals back into the full AO basis with
  /// zero padding.
  void TruncMOsFullBasis(Orbitals& orb, std::vector<Index> activeatoms,
                         std::vector<Index> numfuncpatom);
  /// Insert zero columns into an MO coefficient matrix at the requested
  /// position.
  Eigen::MatrixXd InsertZeroCols(Eigen::MatrixXd MOsMatrix, Index startidx,
                                 Index numofzerocols);
  /// Insert zero rows into an MO coefficient matrix at the requested position.
  Eigen::MatrixXd InsertZeroRows(Eigen::MatrixXd MOsMatrix, Index startidx,
                                 Index numofzerorows);

  /// Run the restricted closed-shell SCF loop and store the converged result.
  bool EvaluateClosedShell(Orbitals& orb, const Mat_p_Energy& H0,
                           const Vxc_Potential<Vxc_Grid>& vxcpotential);

  /// Run the unrestricted Kohn-Sham SCF loop and store alpha and beta orbitals
  /// separately.
  bool EvaluateUKS(Orbitals& orb, const Mat_p_Energy& H0,
                   const Vxc_Potential<Vxc_Grid>& vxcpotential);

  Logger* pLog_;

  // basis sets
  std::string auxbasis_name_;
  std::string dftbasis_name_;
  std::string ecp_name_;
  AOBasis dftbasis_;
  AOBasis auxbasis_;
  ECPAOBasis ecp_;

  Index fock_matrix_reset_;
  // Pre-screening
  double screening_eps_;

  // numerical integration Vxc
  std::string grid_name_;

  // AO Matrices
  AOOverlap dftAOoverlap_;

  std::string initial_guess_;

  // Convergence
  Index numofelectrons_ = 0;
  Index max_iter_;
  ConvergenceAcc::options conv_opt_;
  // DIIS variables
  ConvergenceAcc conv_accelerator_;
  // Electron repulsion integrals
  ERIs ERIs_;

  // external charges
  std::vector<std::unique_ptr<StaticSite> >* externalsites_ = nullptr;

  // exchange and correlation
  double ScaHFX_;
  std::string xc_functional_name_;

  bool integrate_ext_density_ = false;
  // integrate external density
  std::string orbfilename_;
  std::string gridquality_;
  std::string state_;

  QMMolecule activemol_ =
      QMMolecule("molecule made of atoms participating in Active region", 1);

  Eigen::Vector3d extfield_ = Eigen::Vector3d::Zero();
  bool integrate_ext_field_ = false;

  std::string active_atoms_as_string_;
  double active_threshold_;
  double levelshift_;

  // truncation
  Eigen::MatrixXd H0_trunc_;
  Eigen::MatrixXd InitialActiveDmat_trunc_;
  Eigen::MatrixXd v_embedding_trunc_;
  bool truncate_;
  Index active_electrons_;
  double Total_E_full_;
  double E_nuc_;
  double truncation_threshold_;
  std::vector<Index> active_and_border_atoms_;
  std::vector<Index> numfuncpatom_;

  // Spin-DFT Extension
  Index num_alpha_electrons_ = 0;
  Index num_beta_electrons_ = 0;
  Index num_docc_ = 0;
  Index num_socc_alpha_ = 0;
  Index spin_ = 1;
  Index charge_ = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTENGINE_H
