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

// Standard includes
#include <map>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "ERIs.h"
#include "convergenceacc.h"
#include "dfttimings.h"
#include "hirshfeldpartition.h"
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
/**
 * \brief Parts of the DFT setup that stay the same between runs on the same
 * molecule with the same basis sets, such as the outer iterations of QM/MM:
 * the RI three-centre integrals and the AO matrix of an external Ewald
 * potential. Owned by the caller (XTPDFT), lent to one DFTEngine per run.
 * Entries are reused only if their key matches; otherwise they are rebuilt.
 */
struct DFTSetupCache {
  std::string eris_key;
  bool has_eris = false;
  ERIs eris;
  std::string ewald_key;
  Eigen::MatrixXd ewald_matrix;
  // AO potential of the permanent moments of the external sites; the induced
  // dipoles change between QM/MM iterations and are integrated every run.
  std::string multipole_key;
  Eigen::MatrixXd multipole_sites;  // one row per site: pos, rank, Q
  Eigen::MatrixXd multipole_matrix;
};

class DFTEngineTestAccess;

/// True if libint2 supports first-derivative integrals for every category
/// the forces need (one-body, two-centre Coulomb metric, three-centre RI).
/// Compile-time check (LIBINT2_MAX_DERIV_ORDER and the LIBINT_INCLUDE_*
/// guards in libint2_derivative_calls.cc), not a probe of the library linked
/// at runtime. Lets callers, e.g. tests, skip force calculations before
/// running an SCF.
bool HasLibint2DerivativeSupport();

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

  void setEwaldgrid(const Vxc_Grid& ewaldgrid) {
    external_ewaldgrid_ = ewaldgrid;
    has_ewaldgrid_ = true;
  }

  // The nuclei's share of the external Ewald potential, sum_A Z_A
  // phi(R_A). Added to E0 alongside the grid's contribution to H0 -- the
  // grid integrates against the DENSITY and so reaches the electrons
  // only, exactly as IntegrateExternalMultipoles pairs its AO matrix with
  // ExternalRepulsion for the multipole route.
  void setEwaldNuclearEnergy(double energy) { ewald_nuclear_energy_ = energy; }

  /// Use the MOs already in the Orbitals passed to Evaluate as the SCF guess
  /// when they fit this calculation (same basis set and size, same electron
  /// counts); otherwise the configured initial_guess is used. Set by QM/MM
  /// for iterations after the first, where only the environment changed.
  void setWarmStart(bool warm_start) { warm_start_ = warm_start; }

  /// Reuse setup from, and keep it in, this cache (see DFTSetupCache). The
  /// RI integrals are taken out of the cache while the engine runs; call
  /// ReturnSetupCache() afterwards to put them back.
  void setSetupCache(DFTSetupCache* cache) { setup_cache_ = cache; }

  /// Raise the SCF convergence thresholds (energy in Hartree, DIIS error) to
  /// at least these values; call after Initialize.
  void setSCFToleranceFloor(double energy, double error);
  void ReturnSetupCache();

  /// Run a full ground-state DFT calculation and store the results in the
  /// orbital container.
  bool Evaluate(Orbitals& orb);

  /// Run a single, charge-constrained DFT (CDFT) calculation: finds the
  /// Lagrange multiplier lambda such that
  /// Tr[(P_alpha + P_beta) * constraint.weight_matrix] equals
  /// constraint.target_population, then converges the SCF at that
  /// lambda (Wu-Van Voorhis outer loop). Each trial SCF is warm-started
  /// from the previous trial's converged density via the "orbfile"
  /// initial guess, as in CP2K. Always uses EvaluateUKS, since a localized
  /// extra charge is usually an open-shell situation.
  ///
  /// constraint.lambda is the initial guess for the bisection (0.0 is a
  /// reasonable default) and holds the converged value on return. Returns
  /// false if no root can be bracketed (constraints_ cleared, lambda
  /// unchanged) or if the bisection exhausts max_cdft_iterations_ without
  /// reaching cdft_population_tolerance_ (lambda set to the last bracket
  /// midpoint). Throws std::runtime_error if an inner SCF does not converge,
  /// to distinguish that from running out of outer iterations.
  bool RunCDFT(Orbitals& orb, HirshfeldPartition::Constraint& constraint);

  /// Run an embedded active-region DFT calculation for the supplied orbital
  /// container.
  bool EvaluateActiveRegion(Orbitals& orb);
  /// Run the truncated active-region workflow used for reduced embedded
  /// calculations.
  bool EvaluateTruncatedActiveRegion(Orbitals& trunc_orb);

  /// Return the configured AO basis-set name for the DFT calculation.
  std::string getDFTBasisName() const { return dftbasis_name_; };

  /// Static linear response of the converged closed-shell ground state in
  /// orb (same basis, functional and grid as this engine) to one-electron
  /// perturbations h_k (AO matrices of an electron's potential energy, per
  /// unit strength): coupled-perturbed Kohn-Sham,
  ///   (e_a - e_i) U_ai + [C_v^T G[dP] C_o]_ai = -[C_v^T h C_o]_ai,
  ///   dP = 2 (C_v U C_o^T + C_o U^T C_v^T),
  /// G = J + 0.5 a_x K + f_xc (LDA/GGA, global hybrids), solved by
  /// preconditioned conjugate gradients for all h_k together until
  /// |residual| < tolerance |rhs|. Returns the first-order density matrices
  /// dP_k (total). Throws for open-shell states or if CG does not converge.
  std::vector<Eigen::MatrixXd> StaticResponse(
      const Orbitals& orb, const std::vector<Eigen::MatrixXd>& perturbations,
      double tolerance = 1e-7, Index max_iterations = 100);

  /// Name of the integration grid of the SCF (xtpdft.integration_grid).
  const std::string& getGridName() const { return grid_name_; }

  /// A static one-electron perturbation (AO matrix, Hartree) added to the
  /// core Hamiltonian of every following Evaluate: the potential energy
  /// operator of an electron, e.g. lambda * W_k for a Hirshfeld-localized
  /// potential. Its size must match the DFT basis; an empty matrix removes
  /// it. Used for finite-difference response calculations
  /// (ChargeResponseFit).
  void setExternalPotentialMatrix(const Eigen::MatrixXd& v) {
    ext_potential_matrix_ = v;
  }

  /// Free-atom reference densities of the Hirshfeld partition of mol (the
  /// ones CDFT uses), as atom-centred references for HirshfeldPartition.
  std::vector<HirshfeldPartition::AtomicReference> HirshfeldReferences(
      const QMMolecule& mol) const;
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

  // development path RKS vs UKS
  void setForceUKSPath(bool force) { force_uks_path_ = force; }

 private:
  friend class DFTEngineTestAccess;

  /// Initialize basis sets, integral engines, and electron counts before
  /// entering the SCF loop.
  void Prepare(Orbitals& orb, Index numofelectrons = -1);

  /// Build the numerical exchange-correlation integration object for the
  /// current molecule.
  Vxc_Potential<Vxc_Grid> SetupVxc(const QMMolecule& mol);

  /// C_v^T G[dP_k] C_o for dP_k = 2 (C_v U_k C_o^T + C_o U_k^T C_v^T), with
  /// G the closed-shell response of J, 0.5 a_x K and V_xc (f_xc at the
  /// density P0): the coupling part of the orbital Hessian, shared by
  /// StaticResponse and the second-order SCF.
  UKSConvergenceAcc::SpinFock BuildUKSFock(
      const Eigen::MatrixXd& H0, const Vxc_Potential<Vxc_Grid>& vxcpotential,
      const Eigen::MatrixXd& alpha_new, const Eigen::MatrixXd& beta_new);
  Eigen::MatrixXd RKSOrbitalHessian(const Vxc_Potential<Vxc_Grid>& vxc,
                                    const tools::EigenSystem& mos,
                                    const Eigen::MatrixXd& F,
                                    const Eigen::MatrixXd& X);
  Eigen::MatrixXd UKSOrbitalHessian(const Vxc_Potential<Vxc_Grid>& vxc,
                                    const tools::EigenSystem& mos_alpha,
                                    const tools::EigenSystem& mos_beta,
                                    const Eigen::MatrixXd& F_alpha,
                                    const Eigen::MatrixXd& F_beta,
                                    const Eigen::MatrixXd& X);

  std::vector<Eigen::MatrixXd> CoupledResponse(
      const Vxc_Potential<Vxc_Grid>& vxc, const Eigen::MatrixXd& P0,
      const Eigen::MatrixXd& C_occ, const Eigen::MatrixXd& C_virt,
      const std::vector<Eigen::MatrixXd>& U);

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

  /// Static polarizability of the converged closed-shell state in orb
  /// (StaticResponse with the dipole operators), stored in orb.
  void ComputePolarizability(Orbitals& orb);

  /// Evaluate() without the timing report.
  bool EvaluateAndTime(Orbitals& orb);
  /// Basis dimensions and the memory of the large intermediates.
  void ReportDimensionsAndMemory() const;

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

  /// Build a dimer guess (separate alpha and beta MOs) from two converged
  /// monomer .orb files (dimer_guess_orbA_name_/dimer_guess_orbB_name_),
  /// combined via Orbitals::PrepareDimerGuessMixedSpin.
  ///
  /// Checks each monomer against dimer_mol: (1) its element sequence must
  /// match the first N_A (monomer A) or the remaining (monomer B) atoms;
  /// (2) all intra-monomer interatomic distances must match those of the
  /// corresponding dimer atoms to a tight tolerance. Absolute positions are
  /// not compared, since a monomer may be translated/rotated in the dimer.
  Orbitals BuildDimerGuessFromMonomerFiles(const QMMolecule& dimer_mol) const;

  /// Run an unrestricted atomic reference calculation used in open-shell atomic
  /// guesses.
  ///
  /// use_hunds_rule_occupation: if true, take the alpha/beta electron counts
  /// from a Hund's-rule ground-state table for main-group (s/p-block)
  /// elements instead of the default parity split (odd Z -> one extra alpha
  /// electron, even Z -> alpha == beta). The parity split is wrong for e.g.
  /// carbon (triplet, alpha=4, beta=2), which is irrelevant for the SAD guess
  /// (AtomicGuess) but not for Hirshfeld reference densities in CDFT.
  /// Elements not in the table (d/f-block, where the ground state is
  /// functional-dependent) fall back to the parity split with a warning; see
  /// HundsRuleAlphaBetaElectrons in dftengine.cc.
  Eigen::MatrixXd RunAtomicDFT_unrestricted(
      const QMAtom& uniqueAtom, bool use_hunds_rule_occupation = false) const;

  /// Isolated-atom reference density (in the atom's own basis) for each
  /// unique element in mol, keyed by element symbol: the promolecular
  /// densities for Hirshfeld-based CDFT. Same structure as AtomicGuess, but
  /// with use_hunds_rule_occupation=true, and the per-element densities are
  /// returned instead of being assembled into a molecular AO matrix, since
  /// Hirshfeld evaluates them in real space on each atom's position.
  std::map<std::string, Eigen::MatrixXd> ComputeHirshfeldReferenceDensities(
      const QMMolecule& mol) const;

  /// CDFT settings parsed from the \<cdft\> options block in Initialize():
  /// atom indices and target charge only. The full
  /// HirshfeldPartition::Constraint needs the molecule and is built later by
  /// BuildCDFTConstraint.
  struct CDFTConstraintSpec {
    std::vector<Index> atom_indices;
    // Charge relative to the neutral fragment (sum of its nuclear charges);
    // +1.0 means one electron removed. BuildCDFTConstraint converts it to an
    // absolute electron count (CP2K's TARGET convention).
    double target_charge = 0.0;
    double initial_lambda = 0.0;
    // "warmstart" (default) or "fresh", see dftpackage.xml. Kept as a string
    // so that invalid values are rejected by the XML choices validation.
    std::string guess_strategy = "warmstart";
  };

  /// Build the HirshfeldPartition::Constraint for spec on mol. The weight
  /// matrix is the sum of BuildWeightMatrix over spec.atom_indices (Hirshfeld
  /// weights are additive, w_fragment(r) = sum_{i in fragment} w_i(r)). The
  /// absolute target population is sum_{i in fragment} Z_i -
  /// spec.target_charge.
  HirshfeldPartition::Constraint BuildCDFTConstraint(
      const QMMolecule& mol, const CDFTConstraintSpec& spec) const;

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

  /// Assemble the total ground-state nuclear gradient (one-electron
  /// [kinetic + nuclear attraction] + overlap "Pulay force" + nuclear
  /// repulsion + RI-J Coulomb + RI-K exact exchange [hybrid functionals]
  /// + XC, LDA or GGA) from a converged density matrix and store it in
  /// the orbital container via Orbitals::setForces(), as the physical
  /// force (-dE/dR, matching the convention external tools such as ASE
  /// expect from a Calculator's getForces()).
  ///
  /// The overlap Pulay term, -tr(W dS/dR) with the energy-weighted density
  /// W, arises because C^T S C = I depends on the geometry through the basis
  /// functions; it is distinct from DFTGradient::PulayGradient, the
  /// basis-function term of the XC integral.
  ///
  /// The exchange gradient uses the asymmetric V^-1 RI fit, which gives the
  /// same exchange energy as the symmetric V^-1/2 fit in
  /// ERIs::CalculateEXX_mos (exact identity for positive-definite V), so no
  /// derivative of V^-1/2 is needed.
  ///
  /// Requires RI (auxbasis_name_ non-empty); RIJGradient/RIKGradient have no
  /// 4-centre implementation. Only called if compute_forces_ is set
  /// (\<xtpdft\>\<compute_forces\>true\</compute_forces\>\</xtpdft\>,
  /// default false), since it adds non-trivial cost to every SCF.
  void ComputeAndStoreForces(Orbitals& orb, const Eigen::MatrixXd& Dmat,
                             const Vxc_Potential<Vxc_Grid>& vxcpotential) const;

  /// UKS overlap Pulay gradient with W = W_alpha + W_beta (no factor of 2,
  /// unlike RKS). Separate method because it cannot be checked against a
  /// fixed-C finite difference: it accounts for the implicit R-dependence of
  /// C through orthonormality and is only valid at an SCF stationary point.
  /// test_dftengine_private.cc checks it against the RKS formula for
  /// alpha == beta.
  Eigen::MatrixXd ComputeOverlapPulayGradientUKS(
      const QMMolecule& mol, const tools::EigenSystem& MOs_alpha,
      const tools::EigenSystem& MOs_beta) const;

  /// Non-XC part of the UKS gradient (nuclear repulsion, one-electron,
  /// overlap Pulay, RI-J, RI-K), see ComputeAndStoreForcesUKS. Returns
  /// dE/dR (natoms x 3), not the force; the sign flip is done on storage.
  Eigen::MatrixXd ComputeNonXCGradientUKS(
      const QMMolecule& mol, const UKSConvergenceAcc::SpinDensity& Dspin,
      const tools::EigenSystem& MOs_alpha,
      const tools::EigenSystem& MOs_beta) const;

  /// UKS (open-shell) analog of ComputeAndStoreForces; stores the force via
  /// Orbitals::setForces().
  ///
  /// Terms: nuclear repulsion; one-electron [kinetic + nuclear attraction]
  /// and RI-J with D_total = D_alpha + D_beta; overlap Pulay with
  /// W = W_alpha + W_beta (spin densities are not pre-doubled); RI-K for
  /// hybrids as 0.5 * ScaHFX_ * [RIKGradient(C_alpha_occ) +
  /// RIKGradient(C_beta_occ)], the 0.5 because UKS exchange uses
  /// ERIs::CalculateEXX_dmat(P) = 0.5 * ERIs::CalculateEXX_mos(C) for
  /// P = C C^T; and XC (PulayGradientUKS + GridWeightGradientUKS, LDA and
  /// GGA, including the sigma_aa/sigma_ab/sigma_bb terms), checked against
  /// finite differences of IntegrateVXCSpin.
  void ComputeAndStoreForcesUKS(
      Orbitals& orb, const UKSConvergenceAcc::SpinDensity& Dspin,
      const tools::EigenSystem& MOs_alpha, const tools::EigenSystem& MOs_beta,
      const Vxc_Potential<Vxc_Grid>& vxcpotential) const;

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
  // Eigenvalues of the AO overlap below this are removed from S^-1/2
  // (xtpdft.overlap_tolerance).
  double overlap_tolerance_ = 1e-8;
  // Shell pairs whose RI three-centre integrals are bounded by less than this
  // are not stored (xtpdft.ri_pair_threshold).
  double ri_pair_threshold_ = 1e-10;

  // numerical integration Vxc
  std::string grid_name_;

  // AO Matrices
  AOOverlap dftAOoverlap_;

  std::string initial_guess_;

  // Only used if initial_guess_ == "dimer_guess", see
  // BuildDimerGuessFromMonomerFiles.
  std::string dimer_guess_orbA_name_;
  std::string dimer_guess_orbB_name_;

  // Convergence
  Index numofelectrons_ = 0;
  Index max_iter_;
  ConvergenceAcc::options conv_opt_;
  // DIIS variables
  ConvergenceAcc conv_accelerator_;
  // Electron repulsion integrals
  ERIs ERIs_;
  // Wall-clock accounting, reported at the end of Evaluate()
  mutable DFTTimings timings_;
  bool warm_start_ = false;
  bool warm_started_ = false;
  DFTSetupCache* setup_cache_ = nullptr;
  std::string eris_key_;  // key of the RI integrals in ERIs_
  /// Identifies what the RI integrals depend on: basis sets (shell types,
  /// primitives, centres) and the pair threshold.
  std::string RISetupKey() const;
  /// Whether orb's MOs (computed in previous_basis) can serve as the guess
  /// for this calculation; if not, why not.
  bool UsableAsWarmStart(const Orbitals& orb, const std::string& previous_basis,
                         std::string& reason) const;

  // external charges
  std::vector<std::unique_ptr<StaticSite> >* externalsites_ = nullptr;

  // exchange and correlation
  double ScaHFX_ = 0.0;
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
  // xtpdft.polarizability: static polarizability after the SCF
  bool compute_polarizability_ = false;
  Eigen::MatrixXd ext_potential_matrix_;

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

  // QMEwald
  Vxc_Grid external_ewaldgrid_;
  bool has_ewaldgrid_ = false;
  double ewald_nuclear_energy_ = 0.0;
  // Spin-DFT Extension
  Index num_alpha_electrons_ = 0;
  Index num_beta_electrons_ = 0;
  Index num_docc_ = 0;
  Index num_socc_alpha_ = 0;
  Index spin_ = 1;
  Index charge_ = 0;
  bool force_uks_path_ = false;
  // xtpdft.compute_forces: compute the nuclear gradient after the SCF.
  // Off by default because of its cost (one-electron and overlap
  // derivatives, RI-J, XC and, for hybrids, RI-K gradients).
  bool compute_forces_ = false;

  // CDFT constraints added to the UKS Fock matrix; empty (no-op) for
  // ordinary runs. Set by RunCDFT, which updates each lambda between
  // warm-started EvaluateUKS calls.
  std::vector<HirshfeldPartition::Constraint> constraints_;

  // CDFT outer-loop control, used only by RunCDFT. Read from
  // cdft.max_iterations/cdft.population_tolerance in Initialize() if
  // cdft.enabled; the defaults apply when RunCDFT is called directly.
  Index max_cdft_iterations_ = 50;
  double cdft_population_tolerance_ = 1.e-4;

  bool cdft_enabled_ = false;
  CDFTConstraintSpec cdft_constraint_spec_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTENGINE_H
