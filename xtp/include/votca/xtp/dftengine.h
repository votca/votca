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
#include "hirshfeldpartition.h"
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

  /// Run a single, charge-constrained DFT (CDFT) calculation: finds the
  /// Lagrange multiplier lambda such that
  /// Tr[(P_alpha + P_beta) * constraint.weight_matrix] equals
  /// constraint.target_population, then converges the SCF at that
  /// lambda -- the standard Wu-Van Voorhis outer loop, warm-started
  /// (matching CP2K's own documented approach: each new trial's SCF is
  /// restarted from the PREVIOUS trial's converged density, not a cold
  /// start) via the existing "orbfile" initial-guess mechanism, reusing
  /// it exactly as written rather than building new warm-start
  /// machinery. Only ever wires through EvaluateUKS (never
  /// EvaluateClosedShell) -- CDFT charge constraints are built on UKS
  /// from the start, per the design discussion that preceded this: a
  /// localized extra charge is almost always naturally an open-shell/
  /// radical situation regardless of whether spin constraints are ever
  /// added later.
  ///
  /// constraint.lambda is used as the initial guess for the bisection
  /// search (0.0 is a reasonable default for most systems) and is left
  /// holding the converged value on return. Returns false (with
  /// constraints_ left populated, holding the last-attempted lambda)
  /// if EITHER a root cannot be bracketed at all, OR the outer
  /// bisection loop exhausts max_cdft_iterations_ without reaching
  /// cdft_population_tolerance_. If any individual INNER SCF call
  /// itself fails to converge, this throws std::runtime_error instead
  /// (does not return false) -- an inner SCF failure means something
  /// more fundamental than "the outer loop needs more iterations" is
  /// wrong (e.g. a genuinely bad initial guess, or too tight an SCF
  /// convergence threshold for this system), and silently returning
  /// false would look identical to the ordinary "ran out of outer
  /// iterations" case, which it is not.
  bool RunCDFT(Orbitals& orb, HirshfeldPartition::Constraint& constraint);

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
  ///
  /// use_hunds_rule_occupation (default false, preserving all EXISTING
  /// callers' behavior exactly): when true, use a small, explicit
  /// Hund's-rule ground-state alpha/beta electron-count table for
  /// common main-group (s/p-block) elements, instead of the simpler
  /// parity-based split (odd nuclear charge -> one extra alpha electron;
  /// even -> alpha == beta) used by default. That default split is
  /// wrong for many real ground states -- e.g. carbon (true ground
  /// state alpha=4,beta=2, a triplet) gets alpha=beta=3 (an artificial
  /// singlet) -- but this does not matter for a SAD initial-guess
  /// starting DENSITY MATRIX SHAPE (AtomicGuess, this function's only
  /// existing caller), since the full molecule's own SCF reshapes the
  /// density regardless of the isolated reference atom's spin state.
  /// It DOES matter for promolecular reference densities used in
  /// Hirshfeld-based CDFT constraints, which is what this parameter
  /// exists for. Falls back to the default, parity-based split (with a
  /// logged warning) for any element not covered by the table --
  /// currently d/f-block only, where the ground-state configuration is
  /// genuinely ambiguous/functional-dependent rather than a simple,
  /// textbook Hund's-rule case; see HundsRuleAlphaBetaElectrons's own
  /// comment in dftengine.cc.
  Eigen::MatrixXd RunAtomicDFT_unrestricted(
      const QMAtom& uniqueAtom, bool use_hunds_rule_occupation = false) const;

  /// Build one isolated-atom reference density per unique element in
  /// mol, keyed by element symbol -- the promolecular densities
  /// Hirshfeld-based CDFT constraints need. Mirrors AtomicGuess's own
  /// "find unique elements, run RunAtomicDFT_unrestricted once each,
  /// cache by element" structure exactly, but (a) always passes
  /// use_hunds_rule_occupation=true (unlike AtomicGuess's own call,
  /// which never does), and (b) returns the per-element densities
  /// directly rather than assembling them into one combined,
  /// molecule-sized AO-basis matrix -- Hirshfeld only ever needs each
  /// reference density evaluated as a real-space scalar function,
  /// using that element's own (small, atom-only) basis re-centered on
  /// each real atom's actual position, never embedded into the full
  /// molecule's AO basis at all, so there is no molecule-sized object
  /// to assemble here in the first place.
  std::map<std::string, Eigen::MatrixXd> ComputeHirshfeldReferenceDensities(
      const QMMolecule& mol) const;

  /// Parsed directly from the <cdft> options block at Initialize()
  /// time -- atom indices and target charge only, NOT yet a full
  /// HirshfeldPartition::Constraint (which needs the reference
  /// densities/weight matrix, neither of which exist until Evaluate()
  /// actually has a real molecule to work with). BuildCDFTConstraint
  /// (just below) does that later conversion. Deliberately defined
  /// HERE, before BuildCDFTConstraint's own declaration -- a type must
  /// already be visible before it is used as a parameter type in a
  /// member function DECLARATION (unlike a function body, which can
  /// freely reference members declared later in the same class, since
  /// the whole class body is parsed first); this was a real compile
  /// error caught directly (previously defined near
  /// cdft_constraint_spec_ at the end of this class, well after this
  /// point).
  struct CDFTConstraintSpec {
    std::vector<Index> atom_indices;
    // Relative to the fragment's own neutral reference state (the sum
    // of its atoms' nuclear charges) -- e.g. +1.0 means one electron
    // REMOVED (a cation). Converted to an absolute target electron
    // count once, inside BuildCDFTConstraint, matching CP2K's own
    // internal (absolute) TARGET convention exactly -- only the
    // user-facing options syntax is charge-relative, per the earlier
    // design discussion on this.
    double target_charge = 0.0;
    double initial_lambda = 0.0;
  };

  /// Converts a parsed CDFTConstraintSpec (atom indices + relative
  /// target charge, from Initialize()'s own <cdft> options parsing)
  /// into a fully-built HirshfeldPartition::Constraint, given the real
  /// molecule this calculation is actually running on. Builds the
  /// weight matrix as the SUM of BuildWeightMatrix over every atom in
  /// spec.atom_indices -- Hirshfeld weights are additive across atoms
  /// in a fragment (w_fragment(r) = sum_{i in fragment} w_i(r)), so
  /// this generalizes correctly to a multi-atom region, not just a
  /// single atom. The absolute target_population is computed as
  /// (sum of the fragment atoms' own nuclear charges) -
  /// spec.target_charge, matching CP2K's own internal (absolute)
  /// TARGET convention -- only the OPTIONS-file syntax is
  /// charge-relative, per the earlier design discussion on this; the
  /// underlying Constraint/RunCDFT machinery itself was never changed
  /// and still only ever deals in absolute populations.
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
  /// The one-electron term and the overlap Pulay force were initially
  /// MISSING entirely -- discovered by the first genuine end-to-end
  /// SCF+forces test (test_dftengine_forces.cc), since every earlier
  /// gradient test in this branch validated individual terms against
  /// fixed density matrices without ever checking the complete gradient
  /// against a real total SCF energy. The overlap Pulay force (arising
  /// because the MO orthonormality constraint C^T S C = I depends on
  /// geometry through the moving basis functions, weighted by orbital
  /// energies rather than occupation) is a standard part of any
  /// Gaussian-basis SCF gradient, distinct from the "PulayGradient"
  /// naming used elsewhere in this codebase for the XC-integral
  /// basis-function term.
  ///
  /// RI-K/hybrid-functional support was added after DFTGradient::
  /// RIKGradient's energy convention was fixed (a self-introduced
  /// regression during that work: a plausible-looking hand-algebra
  /// "correction" to a half-transformed structure was wrong, caught by
  /// directly, numerically simulating ERIs::CalculateEXX_mos's real
  /// algorithm before committing to it -- the original fully-MO-
  /// transformed structure was correct, needing only a missing factor
  /// of 2) and then verified via a real C++ finite-difference test
  /// against ERIs::CalculateEXX_mos itself, not just a self-consistent
  /// formula. Also confirmed (numerically, to machine precision):
  /// CalculateEXX_mos's symmetric V^-1/2 RI fitting and RIKGradient's
  /// simpler asymmetric V^-1 fitting give IDENTICAL exchange energies
  /// (an exact algebraic identity for symmetric positive-definite V),
  /// so no matrix square root derivative was ever actually needed.
  ///
  /// SCOPE, explicitly checked and logged rather than silently producing
  /// a wrong result: only supported when RI is actually in use for the
  /// SCF (auxbasis_name_ non-empty -- DFTGradient::RIJGradient/RIKGradient
  /// only implement the RI path, not conventional 4-center ERIs).
  ///
  /// OPT-IN: only called at all if compute_forces_ is true (see its
  /// declaration below), settable via
  /// \<xtpdft\>\<compute_forces\>true\</compute_forces\>\</xtpdft\> in
  /// the options tree, defaulting to false.
  /// Computing forces adds real, non-trivial cost to every converged
  /// SCF, so this is deliberately not silently always-on -- added after
  /// this was pointed out as an unflagged side effect of the original,
  /// unconditional wiring.
  void ComputeAndStoreForces(Orbitals& orb, const Eigen::MatrixXd& Dmat,
                             const Vxc_Potential<Vxc_Grid>& vxcpotential) const;

  /// UKS overlap Pulay force -- W = W_alpha + W_beta, each WITHOUT the
  /// factor of 2 RKS uses. Split out as its own method (rather than
  /// inlined in ComputeNonXCGradientUKS) because it needs a genuinely
  /// DIFFERENT validation strategy than the other four terms: it is NOT
  /// checkable against a fixed-C finite difference (confirmed directly
  /// by a failed attempt to do exactly that -- see git history), since
  /// it specifically corrects for C's implicit R-dependence through the
  /// orthonormality constraint, valid only at a genuine SCF stationary
  /// point. See test_dftengine_private.cc for how this is actually
  /// validated instead (reduction to the already-validated RKS formula
  /// when alpha==beta).
  Eigen::MatrixXd ComputeOverlapPulayGradientUKS(
      const QMMolecule& mol, const tools::EigenSystem& MOs_alpha,
      const tools::EigenSystem& MOs_beta) const;

  /// The four non-XC-adjacent gradient terms that generalize cleanly to
  /// UKS -- see the detailed derivation on ComputeAndStoreForcesUKS
  /// below, which calls this and then decides whether/how to report the
  /// result (currently: never stores it, since XC is missing). Returns
  /// the (natoms x 3) dE/dR gradient directly (NOT negated to the
  /// physical force convention -- that flip, if/when this becomes part
  /// of a complete, storable UKS gradient, belongs at the point of
  /// storage, same as the RKS path).
  Eigen::MatrixXd ComputeNonXCGradientUKS(
      const QMMolecule& mol, const UKSConvergenceAcc::SpinDensity& Dspin,
      const tools::EigenSystem& MOs_alpha,
      const tools::EigenSystem& MOs_beta) const;

  /// UKS (open-shell) analog of ComputeAndStoreForces.
  ///
  /// STATUS: PARTIAL, deliberately. Four of the five non-XC-adjacent
  /// terms generalize cleanly to UKS and are implemented here: nuclear
  /// repulsion (unchanged), one-electron [kinetic + nuclear attraction]
  /// (uses D_total = Dspin.alpha + Dspin.beta, exactly the same
  /// convention RKS's Dmat already uses), the overlap Pulay force
  /// (W = W_alpha + W_beta, each WITHOUT the factor of 2 RKS uses, since
  /// UKS spin densities are not pre-doubled), and RI-K exact exchange
  /// for hybrids (0.5 * ScaHFX_ * [RIKGradient(C_alpha_occ,...) +
  /// RIKGradient(C_beta_occ,...)] -- the extra factor of 0.5 relative to
  /// the naive guess of ScaHFX_*(...) confirmed both algebraically and
  /// numerically: ERIs::CalculateEXX_dmat(P) == 0.5 *
  /// ERIs::CalculateEXX_mos(C) when P = C*C^T, checked directly rather
  /// than assumed, since UKS's exact exchange goes through
  /// CalculateEXX_dmat, a different code path than the one
  /// RIKGradient/CalculateEXX_mos were validated against).
  ///
  /// The XC gradient (PulayGradientUKS + GridWeightGradientUKS, LDA and
  /// GGA) is now included too -- initially deferred as new derivation
  /// work (spin-polarized rho_alpha/rho_beta, and for GGA a genuinely
  /// new sigma_alpha-alpha/alpha-beta/beta-beta cross-term structure
  /// with no analog in the spin-restricted case), then completed and
  /// validated (Python-verified formulas first, then a real C++ finite-
  /// difference test against IntegrateVXCSpin -- caught and fixed one
  /// real transcription bug, a missing factor of 2 in the GGA sigma
  /// term's Hessian contraction, found by careful line-by-line
  /// comparison against the verified Python once the first real test
  /// run showed a partial, non-catastrophic discrepancy). With XC now
  /// included, this function DOES call Orbitals::setForces(), same as
  /// the RKS ComputeAndStoreForces.
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
  bool force_uks_path_ = false;
  // Default false to preserve existing performance for callers not
  // using this feature -- computing forces adds real, non-trivial cost
  // (kinetic/nuclear-attraction/overlap derivatives, RI-J gradient, full
  // XC gradient, RI-K for hybrids) to every converged SCF, so this must
  // be explicit opt-in, not silently always-on. Settable via the
  // \<xtpdft\> options block, which flows through unmodified from
  // XTPDFT::RunDFT() (options_) straight into DFTEngine::Initialize --
  // confirmed directly by reading XTPDFT::ParseSpecificOptions, which
  // only extracts a single unrelated field (temporary_file) and does
  // not filter/transform anything else -- so this option is
  // automatically available through the full QMPackage/XTPDFT flow with
  // no changes needed there.
  bool compute_forces_ = false;

  // Empty by default -- the ONLY thing a standard, non-CDFT run needs
  // to know about this member is that it is empty, checked via a
  // single, cheap constraints_.empty() guard inside
  // EvaluateUKS's own Fock-matrix assembly (see that function's own
  // comment at the point the constraint potential term is added).
  // When empty, that guard means the added term is a complete no-op:
  // the Hamiltonian is built exactly as it always was, with no
  // measurable overhead and no change in behavior whatsoever for any
  // run that never touches this member. Populated only by the
  // (not yet implemented) outer Lagrange-multiplier optimization loop,
  // which is expected to modify each Constraint's own lambda field in
  // place between successive, warm-started calls into EvaluateUKS --
  // per the design discussion this grew out of (CP2K's own documented
  // approach: restart the inner SCF from the previous trial's
  // converged density at each new lambda, rather than a cold start).
  std::vector<HirshfeldPartition::Constraint> constraints_;

  // CDFT outer-loop (Lagrange-multiplier) control, used only by
  // RunCDFT below -- never read by the ordinary Evaluate/EvaluateUKS
  // path at all, so these have no bearing on any standard run either.
  // max_cdft_iterations_/cdft_population_tolerance_ are also settable
  // from options (see Initialize()'s own cdft.max_iterations/
  // cdft.population_tolerance parsing) when cdft.enabled=true; their
  // defaults here are what a directly-constructed RunCDFT call (e.g.
  // from a test, bypassing Initialize() entirely) gets instead.
  Index max_cdft_iterations_ = 50;
  double cdft_population_tolerance_ = 1.e-4;

  bool cdft_enabled_ = false;
  CDFTConstraintSpec cdft_constraint_spec_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTENGINE_H
