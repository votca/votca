/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

// Third party includes
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <map>
#include <optional>
#include <string>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/IncrementalFockBuilder.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/dftgradient.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/mmregion.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/pmlocalization.h"
#include "votca/xtp/uks_convergenceacc.h"

namespace votca {
namespace xtp {

namespace {

void CanonicalizeOrbitalPhases(Eigen::MatrixXd& coeffs) {
  constexpr double tol = 1e-14;

  for (Index col = 0; col < coeffs.cols(); ++col) {
    Eigen::Index pivot = 0;
    const double maxabs = coeffs.col(col).cwiseAbs().maxCoeff(&pivot);

    if (maxabs <= tol) {
      continue;
    }

    if (coeffs(pivot, col) < 0.0) {
      coeffs.col(col) *= -1.0;
    }
  }
}

void CanonicalizeOrbitalPhases(tools::EigenSystem& mos) {
  CanonicalizeOrbitalPhases(mos.eigenvectors());
}

}  // namespace

// Defined in libint2_derivative_calls.cc -- forward declared here
// (rather than only near ComputeAndStoreForces further down, where the
// other libint2_derivative_calls.cc forward declarations live) because
// Initialize() below needs it too, to warn early -- at options-parsing
// time, before any SCF work at all -- if compute_forces=true was
// requested on a build that cannot actually do it. See that file's own
// compile-time guard (around LIBINT2_MAX_DERIV_ORDER) for the full
// explanation of why this check exists.
bool HasLibint2DerivativeSupport();

/**
 * Self-consistent Kohn-Sham implementation.
 *
 * The SCF cycle solves F C = S C eps in a Gaussian AO basis. In the
 * restricted branch a single density matrix is iterated, whereas in the UKS
 * branch separate alpha and beta densities are propagated while sharing the
 * same one-electron Hamiltonian, Coulomb term, and AO overlap matrix.
 *
 * Relative to the earlier restricted implementation, the UKS extension keeps
 * the spin channels separate only where the equations require it: exchange,
 * spin-resolved XC potentials, occupations, and convergence acceleration.
 */

void DFTEngine::Initialize(tools::Property& options) {

  const std::string key_xtpdft = "xtpdft";
  dftbasis_name_ = options.get(".basisset").as<std::string>();

  if (options.exists(".auxbasisset")) {
    auxbasis_name_ = options.get(".auxbasisset").as<std::string>();
  }

  if (!auxbasis_name_.empty()) {
    screening_eps_ = options.get(key_xtpdft + ".screening_eps").as<double>();
    fock_matrix_reset_ =
        options.get(key_xtpdft + ".fock_matrix_reset").as<Index>();
  }
  if (options.exists(".ecp")) {
    ecp_name_ = options.get(".ecp").as<std::string>();
  }

  if (options.exists(key_xtpdft + ".force_uks_path")) {
    force_uks_path_ = options.get(key_xtpdft + ".force_uks_path").as<bool>();
  }

  if (options.exists(key_xtpdft + ".compute_forces")) {
    compute_forces_ = options.get(key_xtpdft + ".compute_forces").as<bool>();
  }
  if (compute_forces_ && !HasLibint2DerivativeSupport()) {
    // Fail fast, at options-parsing time, rather than only discovering
    // this after a full (potentially expensive) SCF has already
    // converged. Genuinely throws now (previously only printed a
    // std::cerr WARNING and let Initialize() return normally, so the
    // full SCF still ran to completion regardless, wasting real
    // compute on a calculation that could never produce forces) --
    // throwing std::runtime_error here for an invalid/impossible
    // options combination matches this file's own, already-established
    // convention (see e.g. "Spin multiplicity must be >= 1." and
    // several other throw std::runtime_error(...) calls elsewhere in
    // this same function), not a new pattern.
    throw std::runtime_error(
        "compute_forces=true was requested, but the libint2 this was "
        "built against does not support derivative integrals for one "
        "or more operator categories it needs (one-body, the two-center "
        "Coulomb metric, or three-center RI -- see "
        "libint2_derivative_calls.cc's own compile guards for exactly "
        "which). Many pre-packaged libint2 builds (Homebrew, Ubuntu "
        "apt, etc.) do not enable derivative-integral support for all "
        "of these by default; rebuild libint2 with "
        "--enable-1body/--enable-eri2/--enable-eri3 to use this "
        "feature.");
  }
  if (compute_forces_ && !ecp_name_.empty()) {
    // A real, previously-unguarded gap: the SCF's own Hamiltonian
    // genuinely includes the ECP contribution (H0 = T + V_nuc + V_ECP
    // + V_ext, see the comment on that further down in this file, and
    // dftAOECP.FillPotential(dftbasis_, ecp_) actually called during
    // the SCF itself) -- but ComputeAndStoreForces/ComputeAndStoreForcesUKS
    // have no d(V_ECP)/dR term at all (confirmed directly: neither
    // function references ecp_ or ecp_name_ anywhere). Computing
    // forces anyway in this case would not fail cleanly the way the
    // libint2-support case above does -- it would silently produce a
    // physically INCOMPLETE result (missing the ECP contribution to
    // the force entirely) that looks like a normal, valid force
    // output, which is worse than refusing outright. Refuse instead.
    throw std::runtime_error(
        "compute_forces=true was requested together with an ECP ('" +
        ecp_name_ +
        "'), but analytic nuclear forces do not yet include the ECP "
        "contribution to the force (d(V_ECP)/dR) -- computing forces "
        "in this configuration would silently omit that term rather "
        "than fail visibly. Either drop the ECP or do not request "
        "compute_forces until this is implemented.");
  }

  if (options.exists(key_xtpdft + ".cdft.enabled")) {
    cdft_enabled_ = options.get(key_xtpdft + ".cdft.enabled").as<bool>();
  }
  if (cdft_enabled_) {
    // Deliberately parsed into a CDFTConstraintSpec (atom indices +
    // charge, both directly from the options tree) here, at
    // Initialize() time, rather than building the actual
    // HirshfeldPartition::Constraint (which needs the reference
    // densities and weight matrix) right away -- those need the
    // molecule and basis, neither of which exist yet at this point;
    // BuildCDFTConstraint (elsewhere in this file) does that
    // conversion later, once Evaluate() actually has an Orbitals
    // object with real QMAtoms to work with.
    std::string indices_str =
        options.get(key_xtpdft + ".cdft.indices").as<std::string>();
    if (indices_str.empty()) {
      throw std::runtime_error(
          "cdft.enabled=true was requested, but cdft.indices is empty -- "
          "specify which atoms (0-based, e.g. '1 3 13:17', same syntax "
          "already used for diabatization.xml's own fragment indices) "
          "make up the constrained fragment.");
    }
    cdft_constraint_spec_.atom_indices =
        IndexParser().CreateIndexVector(indices_str);
    cdft_constraint_spec_.target_charge =
        options.get(key_xtpdft + ".cdft.charge").as<double>();
    cdft_constraint_spec_.initial_lambda =
        options.get(key_xtpdft + ".cdft.initial_lambda").as<double>();
    max_cdft_iterations_ =
        options.get(key_xtpdft + ".cdft.max_iterations").as<Index>();
    cdft_population_tolerance_ =
        options.get(key_xtpdft + ".cdft.population_tolerance").as<double>();
    // Note: CDFT itself needs no derivative integral at all (it only
    // ever builds ENERGY-level quantities -- reference densities,
    // weight matrices, Fock-matrix potentials -- never the
    // deriv_order=1 machinery compute_forces needs), so there is no
    // separate HasLibint2DerivativeSupport() check needed here; if
    // compute_forces were ALSO left enabled alongside cdft.enabled,
    // the earlier compute_forces-specific check above already covers
    // that combination.
  }

  initial_guess_ = options.get(".initial_guess").as<std::string>();

  grid_name_ = options.get(key_xtpdft + ".integration_grid").as<std::string>();
  xc_functional_name_ = options.get(".functional").as<std::string>();

  if (options.exists(key_xtpdft + ".externaldensity")) {
    integrate_ext_density_ = true;
    orbfilename_ =
        options.get(key_xtpdft + ".externaldensity.orbfile").as<std::string>();
    gridquality_ = options.get(key_xtpdft + ".externaldensity.gridquality")
                       .as<std::string>();
    state_ =
        options.get(key_xtpdft + ".externaldensity.state").as<std::string>();
  }

  if (options.exists(".externalfield")) {
    integrate_ext_field_ = true;
    extfield_ = options.get(".externalfield").as<Eigen::Vector3d>();
  }

  conv_opt_.Econverged =
      options.get(key_xtpdft + ".convergence.energy").as<double>();
  conv_opt_.error_converged =
      options.get(key_xtpdft + ".convergence.error").as<double>();
  max_iter_ =
      options.get(key_xtpdft + ".convergence.max_iterations").as<Index>();

  std::string method =
      options.get(key_xtpdft + ".convergence.method").as<std::string>();
  if (method == "DIIS") {
    conv_opt_.usediis = true;
  } else if (method == "mixing") {
    conv_opt_.usediis = false;
  }
  if (!conv_opt_.usediis) {
    conv_opt_.histlength = 1;
    conv_opt_.maxout = false;
  }
  conv_opt_.mixingparameter =
      options.get(key_xtpdft + ".convergence.mixing").as<double>();
  conv_opt_.levelshift =
      options.get(key_xtpdft + ".convergence.levelshift").as<double>();
  conv_opt_.levelshiftend =
      options.get(key_xtpdft + ".convergence.levelshift_end").as<double>();
  conv_opt_.maxout =
      options.get(key_xtpdft + ".convergence.DIIS_maxout").as<bool>();
  conv_opt_.histlength =
      options.get(key_xtpdft + ".convergence.DIIS_length").as<Index>();
  conv_opt_.diis_start =
      options.get(key_xtpdft + ".convergence.DIIS_start").as<double>();
  conv_opt_.adiis_start =
      options.get(key_xtpdft + ".convergence.ADIIS_start").as<double>();

  if (options.exists(key_xtpdft + ".dft_in_dft.activeatoms")) {
    active_atoms_as_string_ =
        options.get(key_xtpdft + ".dft_in_dft.activeatoms").as<std::string>();
    active_threshold_ =
        options.get(key_xtpdft + ".dft_in_dft.threshold").as<double>();
    levelshift_ =
        options.get(key_xtpdft + ".dft_in_dft.levelshift").as<double>();
    truncate_ =
        options.get(key_xtpdft + ".dft_in_dft.truncate_basis").as<bool>();
    if (truncate_) {
      truncation_threshold_ =
          options.get(key_xtpdft + ".dft_in_dft.truncation_threshold")
              .as<double>();
    }
  }
}

void DFTEngine::PrintMOs(const Eigen::VectorXd& MOEnergies, Log::Level level) {
  XTP_LOG(level, *pLog_) << "  Orbital energies: " << std::flush;
  XTP_LOG(level, *pLog_) << "  index occupation energy(Hartree) " << std::flush;

  for (Index i = 0; i < MOEnergies.size(); ++i) {
    Index occupancy = 0;
    if (i < num_docc_) {
      occupancy = 2;
    } else if (i < num_docc_ + num_socc_alpha_) {
      occupancy = 1;
    }

    XTP_LOG(level, *pLog_) << (boost::format(" %1$5d      %2$1d   %3$+1.10f") %
                               i % occupancy % MOEnergies(i))
                                  .str()
                           << std::flush;
  }
  return;
}

void DFTEngine::PrintMOsUKS(const Eigen::VectorXd& alpha_energies,
                            const Eigen::VectorXd& beta_energies,
                            Log::Level level) const {
  XTP_LOG(level, *pLog_) << "  UKS orbital energies:" << std::flush;
  XTP_LOG(level, *pLog_) << "  index   occ   eps_a(Ha)         eps_b(Ha)"
                         << std::flush;

  const Index nrows =
      std::max<Index>(alpha_energies.size(), beta_energies.size());

  for (Index i = 0; i < nrows; ++i) {
    const bool occ_a = (i < num_alpha_electrons_);
    const bool occ_b = (i < num_beta_electrons_);

    std::string occ = "0";
    if (occ_a && occ_b) {
      occ = "2";
    } else if (occ_a) {
      occ = "a";
    } else if (occ_b) {
      occ = "b";
    }

    std::string eps_a = "     -";
    std::string eps_b = "     -";

    if (i < alpha_energies.size()) {
      eps_a = (boost::format("%+1.10f") % alpha_energies(i)).str();
    }
    if (i < beta_energies.size()) {
      eps_b = (boost::format("%+1.10f") % beta_energies(i)).str();
    }

    XTP_LOG(level, *pLog_) << (boost::format(
                                   " %1$5d   %2$1s   %3$15s   %4$15s") %
                               i % occ % eps_a % eps_b)
                                  .str()
                           << std::flush;
  }

  if (num_alpha_electrons_ > 0 &&
      num_alpha_electrons_ < alpha_energies.size()) {
    XTP_LOG(level, *pLog_) << (boost::format(
                                   "  alpha HOMO-LUMO gap: %+1.10f Ha") %
                               (alpha_energies(num_alpha_electrons_) -
                                alpha_energies(num_alpha_electrons_ - 1)))
                                  .str()
                           << std::flush;
  }

  if (num_beta_electrons_ > 0 && num_beta_electrons_ < beta_energies.size()) {
    XTP_LOG(level, *pLog_) << (boost::format(
                                   "  beta  HOMO-LUMO gap: %+1.10f Ha") %
                               (beta_energies(num_beta_electrons_) -
                                beta_energies(num_beta_electrons_ - 1)))
                                  .str()
                           << std::flush;
  }
}

void DFTEngine::CalcElDipole(const Orbitals& orb) const {
  QMState state = QMState("n");
  Eigen::Vector3d result = orb.CalcElDipole(state);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Electric Dipole is[e*bohr]:\n\t\t dx=" << result[0]
      << "\n\t\t dy=" << result[1] << "\n\t\t dz=" << result[2] << std::flush;
  return;
}

// Assembles the total ground-state gradient (nuclear repulsion + RI-J
// Coulomb + XC, LDA or GGA) from the converged density matrix, negates it
// to the physical force convention, and stores it via Orbitals::setForces().
// See the detailed SCOPE note on the declaration in dftengine.h for exactly
// which cases this does and does not support, and why.
//
// Every individual term here (NuclearRepulsionDerivative, RIJGradient,
// PulayGradient, GridWeightGradient) was separately derived and validated
// via finite-difference tests earlier in this branch (see
// test_dftgradient.cc and test_xcgradient.cc) -- this function's own new
// content is just the SUMMATION and the sign convention, not any new
// derivative math.
// Defined in libint2_derivative_calls.cc, not yet in any header (same
// STATUS noted throughout that file) -- forward declared here. Unlike
// DFTGradient::RIJGradient/PulayGradient/etc., these return RAW AO-matrix
// derivatives (d(matrix_munu)/dR), not already-contracted energy
// gradients -- the contraction with Dmat is done explicitly below.
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis& aobasis, const QMMolecule& mol);
// HasLibint2DerivativeSupport() (used below in both ComputeAndStoreForces
// and ComputeAndStoreForcesUKS) is already forward declared earlier in
// this file, near Initialize() -- see that declaration's own comment
// for why it needed to be that early.

void DFTEngine::ComputeAndStoreForces(
    Orbitals& orb, const Eigen::MatrixXd& Dmat,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {
  if (auxbasis_name_.empty()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping force calculation: RI-J gradient (DFTGradient::"
           "RIJGradient) only implements the RI path, but this SCF ran "
           "without an auxiliary basis (conventional 4-center ERIs)."
        << std::flush;
    return;
  }

  if (!HasLibint2DerivativeSupport()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping force calculation: the libint2 this was built "
           "against does not support derivative integrals for one or "
           "more operator categories it needs. Many pre-packaged "
           "libint2 builds (Homebrew, Ubuntu apt, etc.) do not enable "
           "this by default -- rebuild libint2 with "
           "--enable-1body/--enable-eri2/--enable-eri3 to use analytic "
           "forces."
        << std::flush;
    return;
  }

  if (!ecp_name_.empty()) {
    // Same reasoning as Initialize()'s own, earlier check (which
    // should already have caught this before any SCF work even
    // started) -- this is a defense-in-depth repeat, matching the
    // existing HasLibint2DerivativeSupport() re-check just above,
    // in case compute_forces_/ecp_name_ were ever set some other way
    // than through Initialize()'s own options parsing. Skips cleanly
    // (log + return) rather than throwing here, matching this
    // function's own existing style for the libint2-support case
    // above -- by the time SCF has already converged this far,
    // throwing would be a less graceful failure than simply not
    // storing forces, though Initialize()'s own check is the
    // preferred, much earlier place for this to actually be caught.
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping force calculation: an ECP ('" << ecp_name_
        << "') was used for this SCF, but analytic nuclear forces do "
           "not yet include the ECP contribution to the force "
           "(d(V_ECP)/dR) -- computing forces in this configuration "
           "would silently omit that term rather than fail visibly."
        << std::flush;
    return;
  }

  const QMMolecule& mol = orb.QMAtoms();
  Index natoms = mol.size();

  // One-electron (kinetic + nuclear attraction) contribution --
  // dEone/dR_A = Tr[Dmat . d(T+V_ne)/dR_A]. This was the piece
  // discovered MISSING from the total gradient by the first genuine
  // end-to-end SCF+forces test (test_dftengine_forces.cc): kinetic
  // derivatives were validated at the very start of this whole branch
  // and then never actually wired into any gradient assembly, and
  // nuclear attraction derivatives were never implemented at all until
  // that gap was found. See ComputeNuclearAttractionDerivatives in
  // libint2_derivative_calls.cc for the detailed derivation (sign
  // convention checked directly against AOMultipole's own,
  // already-validated energy-level code, not assumed).
  std::vector<AOMatrixDerivative> dT = ComputeKineticDerivatives(dftbasis_);
  std::vector<AOMatrixDerivative> dVne =
      ComputeNuclearAttractionDerivatives(dftbasis_, mol);
  Eigen::MatrixXd eone_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      eone_grad(a, xyz) = Dmat.cwiseProduct(dT[a][xyz] + dVne[a][xyz]).sum();
    }
  }

  // Overlap "Pulay force" -- a SECOND, genuinely distinct missing term,
  // found after the kinetic+nuclear-attraction fix improved but did not
  // fully resolve the discrepancy against the end-to-end finite-difference
  // test (magnitude dropped ~7x in the right direction, but still wrong
  // by roughly the size of a real missing term, not noise).
  //
  // Distinct from the earlier "PulayGradient" naming (which is about
  // basis functions inside the XC integral) -- this is the CLASSICAL
  // SCF Pulay/overlap force, present in essentially any Gaussian-basis
  // HF/DFT gradient: the MO coefficients C are only implicitly
  // R-independent because they satisfy the orthonormality constraint
  // C^T S C = I, and S itself depends on R (basis functions move). At
  // the SCF stationary point, the Lagrange multipliers for this
  // constraint are exactly the orbital energies (canonical MOs), giving
  // an extra term dE/dR_A|_overlap = -Tr[W . dS/dR_A], where
  // W = 2 * C_occ * diag(eps_occ) * C_occ^T (the "energy-weighted
  // density matrix", factor of 2 matching the same doubled convention
  // Dmat already uses for closed-shell restricted). Confirmed as a
  // standard, expected term by libint2's own reference SCF-gradient
  // example (compute_1body_ints_deriv<Operator::overlap> combined with
  // exactly this W construction, in
  // libint2/include/libint2/lcao/1body.h) -- not a novel derivation.
  //
  // ComputeOverlapDerivatives itself was validated (finite-difference
  // tested) at the very start of this whole branch and then never
  // actually used in any gradient assembly until now, same as kinetic.
  Index n_occ = num_docc_ + num_socc_alpha_;
  Eigen::MatrixXd C_occ = orb.MOs().eigenvectors().leftCols(n_occ);
  Eigen::VectorXd eps_occ = orb.MOs().eigenvalues().head(n_occ);
  Eigen::MatrixXd W = 2.0 * C_occ * eps_occ.asDiagonal() * C_occ.transpose();

  std::vector<AOMatrixDerivative> dS = ComputeOverlapDerivatives(dftbasis_);
  Eigen::MatrixXd overlap_pulay_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      overlap_pulay_grad(a, xyz) = -W.cwiseProduct(dS[a][xyz]).sum();
    }
  }

  Eigen::MatrixXd rij_term = DFTGradient::RIJGradient(Dmat, auxbasis_, dftbasis_);
  Eigen::MatrixXd pulay_term = vxcpotential.PulayGradient(Dmat, dftbasis_);
  Eigen::MatrixXd weight_term = vxcpotential.GridWeightGradient(Dmat, mol);
  Eigen::MatrixXd nucrep_term = DFTGradient::NuclearRepulsionDerivative(mol);

  Eigen::MatrixXd grad =
      nucrep_term +
      eone_grad +
      overlap_pulay_grad +
      rij_term +
      pulay_term +
      weight_term;

  // Exact-exchange (RI-K) gradient -- hybrid functionals only. Skipped
  // entirely (not just multiplied by a zero ScaHFX_) when not needed,
  // since RIKGradient is genuinely expensive (O(nocc^2 * naux) linear
  // solves) unlike the GGA sigma terms, which are cheap enough to
  // compute unconditionally.
  //
  // RIKGradient's own energy convention (E_K = -sum_ij c_ij.d_ij) was
  // confirmed, via direct numerical simulation of
  // ERIs::CalculateEXX_mos's real algorithm and then a real C++
  // finite-difference test against that same production function
  // (test_dftgradient.cc), to equal EXACTLY 0.25*Dmat.cwiseProduct(K).sum()
  // at ScaHFX_=1 -- so for general ScaHFX_, the contribution is
  // ScaHFX_ * RIKGradient(...), a direct scaling, matching exactly how
  // the real SCF energy scales its own exx term
  // (exx = 0.25*ScaHFX_*Dmat.cwiseProduct(K).sum()).
  //
  // This removes what was previously an explicit, logged SCOPE
  // limitation (hybrid functionals skipped entirely) -- see git history
  // for the full derivation/verification that led to this.
  if (ScaHFX_ > 0.0) {
    grad += ScaHFX_ * DFTGradient::RIKGradient(C_occ, auxbasis_, dftbasis_);
  }

  // Sanity check independent of the finite-difference tests already done
  // per-term: translational invariance means the TOTAL gradient must sum
  // to zero across all atoms. Logged rather than asserted/thrown --
  // deliberately not blocking a real SCF run over a force-only sanity
  // check, but worth knowing about if it ever fires.
  Eigen::Vector3d sum = grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " WARNING: computed forces do not sum to zero across atoms "
           "(translational invariance check failed, max component="
        << sum.cwiseAbs().maxCoeff() << ") -- treat these forces with "
        "caution."
        << std::flush;
  }

  // Physical force = -dE/dR, matching the convention external tools
  // (e.g. ASE's Calculator.get_forces()) expect -- NuclearRepulsionDerivative/
  // RIJGradient/PulayGradient/GridWeightGradient all return dE/dR directly
  // (the gradient, not the force), consistent with each other throughout
  // this branch; negating once here, at the point of storage, rather than
  // in each individual term, keeps that internal convention consistent
  // and puts the physical-force sign flip in exactly one place.
  Eigen::MatrixXd force = -grad;
  orb.setForces(force);

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Computed and stored ground-state nuclear forces."
      << std::flush;
  // Atomic units (Hartree/Bohr) -- deliberately NOT converted here, same
  // convention as what gets stored via setForces()/WriteToCpt above and
  // what the rest of this file's own log output uses for energies
  // (Hartree throughout; only the earlier "Molecule Coordinates" section
  // converts to Angstrom, for readability, and that conversion is
  // unrelated to this).
  XTP_LOG(Log::error, *pLog_) << " Forces [Ha/Bohr]" << std::flush;
  for (Index a = 0; a < natoms; ++a) {
    std::string output = (boost::format(" %1$s"
                                        "   %2$+1.6f %3$+1.6f %4$+1.6f") %
                          mol[a].getElement() % force(a, 0) % force(a, 1) %
                          force(a, 2))
                             .str();
    XTP_LOG(Log::error, *pLog_) << output << std::flush;
  }
}

Eigen::MatrixXd DFTEngine::ComputeOverlapPulayGradientUKS(
    const QMMolecule& mol, const tools::EigenSystem& MOs_alpha,
    const tools::EigenSystem& MOs_beta) const {
  // W = W_alpha + W_beta, each WITHOUT the factor of 2 RKS uses -- UKS
  // spin densities/MO occupations are not pre-doubled (each spin
  // channel already corresponds to its own electron count).
  //
  // NOTE ON VALIDATION: this term is NOT checkable against a
  // fixed-C finite difference the way the other UKS gradient terms are
  // -- confirmed directly by a failed attempt to do exactly that (see
  // git history). The overlap Pulay force specifically corrects for C's
  // IMPLICIT R-dependence through the orthonormality constraint
  // C^T S(R) C = I, valid only at a genuine variational stationary
  // point (the Lagrange-multiplier argument requires C to actually be a
  // converged SCF solution) -- a fixed, arbitrary C held constant across
  // displaced geometries never satisfies that constraint
  // self-consistently, so there is no fixed-C energy this term is
  // supposed to match. This mirrors exactly why the RKS version of this
  // term was only ever validated by a genuine, self-consistent
  // end-to-end SCF test (test_dftengine_forces.cc), never a
  // fixed-density-matrix unit test. See
  // compute_non_xc_gradient_uks_finite_difference and
  // overlap_pulay_gradient_uks_reduces_to_rks in
  // test_dftengine_private.cc for how this piece is actually checked
  // instead: the other four terms against a fixed-C finite difference,
  // and this term separately against the already-validated RKS formula
  // in the alpha==beta limit.
  Index n_occ_alpha = num_alpha_electrons_;
  Index n_occ_beta = num_beta_electrons_;
  Eigen::MatrixXd C_alpha_occ = MOs_alpha.eigenvectors().leftCols(n_occ_alpha);
  Eigen::MatrixXd C_beta_occ = MOs_beta.eigenvectors().leftCols(n_occ_beta);
  Eigen::VectorXd eps_alpha_occ = MOs_alpha.eigenvalues().head(n_occ_alpha);
  Eigen::VectorXd eps_beta_occ = MOs_beta.eigenvalues().head(n_occ_beta);
  Eigen::MatrixXd W = C_alpha_occ * eps_alpha_occ.asDiagonal() *
                          C_alpha_occ.transpose() +
                      C_beta_occ * eps_beta_occ.asDiagonal() *
                          C_beta_occ.transpose();

  Index natoms = mol.size();
  std::vector<AOMatrixDerivative> dS = ComputeOverlapDerivatives(dftbasis_);
  Eigen::MatrixXd overlap_pulay_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      overlap_pulay_grad(a, xyz) = -W.cwiseProduct(dS[a][xyz]).sum();
    }
  }
  return overlap_pulay_grad;
}

Eigen::MatrixXd DFTEngine::ComputeNonXCGradientUKS(
    const QMMolecule& mol, const UKSConvergenceAcc::SpinDensity& Dspin,
    const tools::EigenSystem& MOs_alpha,
    const tools::EigenSystem& MOs_beta) const {
  Index natoms = mol.size();
  const Eigen::MatrixXd D_total = Dspin.total();

  // One-electron and RI-J: identical formulas/conventions to the RKS
  // case, just built from D_total = Dspin.alpha + Dspin.beta -- matches
  // exactly how RKS's own Dmat is already alpha+beta (E_one and E_coul
  // in EvaluateUKS use D_total the same way EvaluateClosedShell's Eone/
  // Etwo use Dmat), confirmed directly by reading EvaluateUKS rather
  // than assumed.
  std::vector<AOMatrixDerivative> dT = ComputeKineticDerivatives(dftbasis_);
  std::vector<AOMatrixDerivative> dVne =
      ComputeNuclearAttractionDerivatives(dftbasis_, mol);
  Eigen::MatrixXd eone_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      eone_grad(a, xyz) =
          D_total.cwiseProduct(dT[a][xyz] + dVne[a][xyz]).sum();
    }
  }

  Eigen::MatrixXd overlap_pulay_grad =
      ComputeOverlapPulayGradientUKS(mol, MOs_alpha, MOs_beta);

  Eigen::MatrixXd grad = DFTGradient::NuclearRepulsionDerivative(mol) +
                         eone_grad + overlap_pulay_grad +
                         DFTGradient::RIJGradient(D_total, auxbasis_, dftbasis_);

  // Exact exchange (RI-K), hybrids only. Factor of 0.5*ScaHFX_ (not
  // ScaHFX_ alone) -- confirmed both algebraically and numerically
  // (Python, to ~1e-14) that ERIs::CalculateEXX_dmat(P) ==
  // 0.5*ERIs::CalculateEXX_mos(C) when P=C*C^T, and UKS's own exact
  // exchange goes through CalculateEXX_dmat (a DIFFERENT code path than
  // RIKGradient was validated against, which uses CalculateEXX_mos
  // directly) -- tracing that factor of 0.5 through both spin channels'
  // energy expressions gives dE_exx/dR =
  // 0.5*ScaHFX_*[RIKGradient(C_alpha_occ)+RIKGradient(C_beta_occ)], not
  // the naive ScaHFX_*(...) that would be a factor-of-2 error.
  if (ScaHFX_ > 0.0) {
    Eigen::MatrixXd C_alpha_occ =
        MOs_alpha.eigenvectors().leftCols(num_alpha_electrons_);
    Eigen::MatrixXd C_beta_occ =
        MOs_beta.eigenvectors().leftCols(num_beta_electrons_);
    grad += 0.5 * ScaHFX_ *
            (DFTGradient::RIKGradient(C_alpha_occ, auxbasis_, dftbasis_) +
             DFTGradient::RIKGradient(C_beta_occ, auxbasis_, dftbasis_));
  }
  return grad;
}

void DFTEngine::ComputeAndStoreForcesUKS(
    Orbitals& orb, const UKSConvergenceAcc::SpinDensity& Dspin,
    const tools::EigenSystem& MOs_alpha, const tools::EigenSystem& MOs_beta,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {
  if (auxbasis_name_.empty()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping UKS force calculation: RI-J gradient only "
           "implements the RI path, but this SCF ran without an "
           "auxiliary basis."
        << std::flush;
    return;
  }

  if (!HasLibint2DerivativeSupport()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping UKS force calculation: the libint2 this was "
           "built against does not support derivative integrals for "
           "one or more operator categories it needs. Many "
           "pre-packaged libint2 builds (Homebrew, Ubuntu apt, etc.) "
           "do not enable this by default -- rebuild libint2 with "
           "--enable-1body/--enable-eri2/--enable-eri3 to use "
           "analytic forces."
        << std::flush;
    return;
  }

  if (!ecp_name_.empty()) {
    // Same reasoning as the RKS ComputeAndStoreForces' own, identical
    // check just above (and Initialize()'s own, earlier, preferred
    // check) -- ECP forces are not implemented in either spin
    // channel's gradient assembly.
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Skipping UKS force calculation: an ECP ('" << ecp_name_
        << "') was used for this SCF, but analytic nuclear forces do "
           "not yet include the ECP contribution to the force "
           "(d(V_ECP)/dR) -- computing forces in this configuration "
           "would silently omit that term rather than fail visibly."
        << std::flush;
    return;
  }

  Eigen::MatrixXd grad =
      ComputeNonXCGradientUKS(orb.QMAtoms(), Dspin, MOs_alpha, MOs_beta);

  // XC gradient (LDA and GGA both supported -- see the detailed
  // derivation/validation history on this function's declaration in
  // dftengine.h and on PulayGradientUKS/GridWeightGradientUKS in
  // vxc_potential.h).
  grad += vxcpotential.PulayGradientUKS(Dspin.alpha, Dspin.beta, dftbasis_);
  grad += vxcpotential.GridWeightGradientUKS(Dspin.alpha, Dspin.beta,
                                             orb.QMAtoms());

  // Sanity check independent of the finite-difference tests already
  // done per-term: translational invariance means the TOTAL gradient
  // must sum to zero across all atoms. Logged rather than asserted/
  // thrown, same as the RKS path -- deliberately not blocking a real
  // SCF run over a force-only sanity check.
  Eigen::Vector3d sum = grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " WARNING: computed UKS forces do not sum to zero across "
           "atoms (translational invariance check failed, max "
           "component="
        << sum.cwiseAbs().maxCoeff() << ") -- treat these forces with "
        "caution."
        << std::flush;
  }

  // Physical force = -dE/dR, matching the RKS ComputeAndStoreForces
  // convention exactly -- all pieces above return dE/dR directly (the
  // gradient, not the force), negated once here at the point of
  // storage.
  Eigen::MatrixXd force = -grad;
  orb.setForces(force);

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " Computed and stored ground-state UKS nuclear forces."
      << std::flush;
  // Same convention as ComputeAndStoreForces (RKS): atomic units
  // (Hartree/Bohr), matching what gets stored via setForces() above.
  const QMMolecule& mol_for_print = orb.QMAtoms();
  XTP_LOG(Log::error, *pLog_) << " Forces [Ha/Bohr]" << std::flush;
  for (Index a = 0; a < force.rows(); ++a) {
    std::string output = (boost::format(" %1$s"
                                        "   %2$+1.6f %3$+1.6f %4$+1.6f") %
                          mol_for_print[a].getElement() % force(a, 0) %
                          force(a, 1) % force(a, 2))
                             .str();
    XTP_LOG(Log::error, *pLog_) << output << std::flush;
  }
}

// Build the Coulomb and exact-exchange contributions generated by the current
// density matrix. The returned pair is conventionally interpreted as
//
//   (J[P], -K[P]),
//
// so that the hybrid Fock update becomes F = H0 + J[P] + a_x (-K[P]) + V_xc.
// For RI/3c builds the occupied MO block is supplied when available to avoid an
// unnecessary reconstruction of exchange intermediates.
std::array<Eigen::MatrixXd, 2> DFTEngine::CalcERIs_EXX(
    const Eigen::MatrixXd& MOCoeff, const Eigen::MatrixXd& Dmat,
    double error) const {
  if (!auxbasis_name_.empty()) {
    if (conv_accelerator_.getUseMixing() || MOCoeff.rows() == 0) {
      return ERIs_.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), Dmat);
    } else {
      Eigen::MatrixXd occblock = MOCoeff.leftCols(num_docc_ + num_socc_alpha_);
      return ERIs_.CalculateERIs_EXX_3c(occblock, Dmat);
    }
  } else {
    return ERIs_.CalculateERIs_EXX_4c(Dmat, error);
  }
}

// Pure Coulomb contribution J[P] from the current AO density matrix. The code
// dispatches to either RI/3c or conventional 4-center integral evaluation.
Eigen::MatrixXd DFTEngine::CalcERIs(const Eigen::MatrixXd& Dmat,
                                    double error) const {
  if (!auxbasis_name_.empty()) {
    return ERIs_.CalculateERIs_3c(Dmat);
  } else {
    return ERIs_.CalculateERIs_4c(Dmat, error);
  }
}

tools::EigenSystem DFTEngine::IndependentElectronGuess(
    const Mat_p_Energy& H0) const {
  return conv_accelerator_.SolveFockmatrix(H0.matrix());
}

// Construct a self-consistent model-potential guess by starting from an
// atomic density P^(0), evaluating
//
//   F[P^(0)] = H0 + J[P^(0)] + a_x (-K[P^(0)]) + V_xc[P^(0)],
//
// and diagonalizing the resulting Fock matrix once.
tools::EigenSystem DFTEngine::ModelPotentialGuess(
    const Mat_p_Energy& H0, const QMMolecule& mol,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {
  Eigen::MatrixXd Dmat = AtomicGuess(mol);
  Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Vxc matrix " << std::flush;

  Eigen::MatrixXd H = H0.matrix() + e_vxc.matrix();

  if (ScaHFX_ > 0) {
    std::array<Eigen::MatrixXd, 2> both =
        CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), Dmat, 1e-12);
    H += both[0];
    H += ScaHFX_ * both[1];
  } else {
    H += CalcERIs(Dmat, 1e-12);
  }
  return conv_accelerator_.SolveFockmatrix(H);
}

bool DFTEngine::Evaluate(Orbitals& orb) {
  if (cdft_enabled_) {
    // Deliberately dispatched here, BEFORE any of the normal
    // Prepare/SetupH0/SetupVxc/ConfigOrbfile setup below -- RunCDFT
    // does that same setup internally itself (matching this
    // function's own structure exactly), so doing it here too would
    // just duplicate the work. BuildCDFTConstraint needs orb.QMAtoms()
    // to already be set (the same requirement Evaluate() itself has,
    // via SetupH0(orb.QMAtoms()) below), so this is not adding any new
    // requirement on the caller.
    HirshfeldPartition::Constraint constraint =
        BuildCDFTConstraint(orb.QMAtoms(), cdft_constraint_spec_);
    bool converged = RunCDFT(orb, constraint);

    if (converged && orb.hasForces()) {
      // RunCDFT's own final EvaluateUKS call already computed and
      // stored the ordinary DFT force (via ComputeAndStoreForcesUKS,
      // triggered internally whenever compute_forces_ is also set) --
      // this adds the CDFT-specific correction on top of it. Done
      // HERE, once, after RunCDFT's outer loop has fully converged --
      // deliberately NOT inside ComputeAndStoreForcesUKS itself (which
      // would otherwise redo this work, wastefully and riskily, at
      // EVERY outer CDFT iteration, since RunCDFT calls EvaluateUKS
      // repeatedly) and deliberately NOT by changing RunCDFT's own
      // signature to pass through the original per-atom fragment
      // indices (constraint only carries the already-SUMMED
      // weight_matrix, not which atoms went into it -- rebuilding here
      // instead, via cdft_constraint_spec_'s own atom_indices, avoids
      // touching RunCDFT's own, already-validated signature/behavior
      // at all).
      //
      // Rebuilds the same reference densities/atomic references/
      // basis/grid BuildCDFTConstraint itself already built internally
      // -- a redundant but cheap recomputation (no SCF involved),
      // accepted deliberately for this reason.
      std::map<std::string, Eigen::MatrixXd> reference_densities =
          ComputeHirshfeldReferenceDensities(orb.QMAtoms());
      AOBasis full_dftbasis;
      {
        BasisSet basisset;
        basisset.Load(dftbasis_name_);
        full_dftbasis.Fill(basisset, orb.QMAtoms());
      }
      Vxc_Grid grid;
      grid.GridSetup(grid_name_, orb.QMAtoms(), full_dftbasis);
      std::vector<HirshfeldPartition::AtomicReference> atoms =
          HirshfeldPartition::BuildAtomicReferences(
              orb.QMAtoms(), dftbasis_name_, reference_densities);

      std::array<Eigen::MatrixXd, 2> Dspin =
          orb.DensityMatrixGroundStateSpinResolved();
      // Total (alpha+beta) density -- matches the charge constraint's
      // own spin_alpha_coefficient=spin_beta_coefficient=+1.0
      // convention exactly (Tr[(D_alpha+D_beta)*W] = the same
      // population EvaluateMismatch itself computes inside RunCDFT).
      Eigen::MatrixXd density_total = Dspin[0] + Dspin[1];

      Eigen::MatrixXd cdft_gradient_correction = Eigen::MatrixXd::Zero(
          static_cast<Index>(orb.QMAtoms().size()), 3);
      for (Index atom_index : cdft_constraint_spec_.atom_indices) {
        cdft_gradient_correction +=
            HirshfeldPartition::ComputeCDFTForceContribution(
                atoms, atom_index, density_total, orb.QMAtoms(),
                full_dftbasis, grid);
      }
      // Physical force = -dE/dR (ComputeAndStoreForcesUKS's own,
      // already-established convention): the CDFT correction to the
      // GRADIENT is +lambda*d(Tr[D*W_c])/dR (added directly, matching
      // ComputeCDFTForceContribution's own gradient-convention
      // return), so the correction to the FORCE is -lambda times this
      // same quantity.
      orb.setForces(orb.getForces() -
                    constraint.lambda * cdft_gradient_correction);
    }
    return converged;
  }

  Prepare(orb);
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  if (force_uks_path_ || num_alpha_electrons_ != num_beta_electrons_) {
    if (force_uks_path_ && num_alpha_electrons_ == num_beta_electrons_) {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp()
          << " Forcing closed-shell singlet through UKS development path."
          << std::flush;
    }
    return EvaluateUKS(orb, H0, vxcpotential);
  }
  return EvaluateClosedShell(orb, H0, vxcpotential);
}

bool DFTEngine::RunCDFT(Orbitals& orb,
                        HirshfeldPartition::Constraint& constraint) {
  Prepare(orb);
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  // Restored on every exit path (converged or not) -- RunCDFT
  // deliberately overrides this member's own value between outer
  // iterations (to force the warm-start "orbfile" guess from the
  // second iteration onward), so it must not leak whatever value the
  // caller's own options actually specified.
  std::string saved_initial_guess = initial_guess_;

  constraints_ = {constraint};

  // Bisection bracket for lambda -- deliberately not Newton's method:
  // bisection needs only that the population is monotonic in lambda
  // (true for a well-behaved CDFT problem: increasing lambda always
  // pushes more density toward -- or away from, depending on sign --
  // the constrained region), never an explicit dN/dlambda derivative,
  // making this the more robust choice for a first implementation.
  // Starts centered on the caller's own initial guess (constraint.lambda,
  // 0.0 by default) and expands outward, doubling each time, until the
  // mismatch changes sign across the bracket or a hard iteration limit
  // is hit -- rather than assuming any single fixed bracket width is
  // always wide enough for every system.
  double lambda_lo = constraint.lambda - 0.1;
  double lambda_hi = constraint.lambda + 0.1;

  auto EvaluateMismatch = [&](double lambda) -> double {
    constraints_[0].lambda = lambda;
    bool scf_converged = EvaluateUKS(orb, H0, vxcpotential);
    if (!scf_converged) {
      throw std::runtime_error(
          "RunCDFT: inner SCF did not converge at lambda=" +
          std::to_string(lambda));
    }
    initial_guess_ = "orbfile";  // warm start every subsequent call
    std::array<Eigen::MatrixXd, 2> Dspin =
        orb.DensityMatrixGroundStateSpinResolved();
    double population =
        constraint.spin_alpha_coefficient *
            Dspin[0].cwiseProduct(constraint.weight_matrix).sum() +
        constraint.spin_beta_coefficient *
            Dspin[1].cwiseProduct(constraint.weight_matrix).sum();
    return population - constraint.target_population;
  };

  double mismatch_lo;
  double mismatch_hi;
  try {
    mismatch_lo = EvaluateMismatch(lambda_lo);
    mismatch_hi = EvaluateMismatch(lambda_hi);

    Index bracket_attempts = 0;
    constexpr Index kMaxBracketAttempts = 10;
    while (mismatch_lo * mismatch_hi > 0.0 &&
          bracket_attempts < kMaxBracketAttempts) {
      double width = lambda_hi - lambda_lo;
      lambda_lo -= 0.5 * width;
      lambda_hi += 0.5 * width;
      mismatch_lo = EvaluateMismatch(lambda_lo);
      mismatch_hi = EvaluateMismatch(lambda_hi);
      ++bracket_attempts;
    }
    if (mismatch_lo * mismatch_hi > 0.0) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " RunCDFT: could not bracket a root for the population "
             "mismatch after "
          << kMaxBracketAttempts
          << " bracket-expansion attempts -- the target population may "
             "be unreachable for this system, or the initial "
             "lambda guess may be far from the actual root."
          << std::flush;
      initial_guess_ = saved_initial_guess;
      constraints_.clear();
      return false;
    }

    for (Index outer_iter = 0; outer_iter < max_cdft_iterations_;
        ++outer_iter) {
      double lambda_mid = 0.5 * (lambda_lo + lambda_hi);
      double mismatch_mid = EvaluateMismatch(lambda_mid);

      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " CDFT outer iteration " << outer_iter + 1
          << " of " << max_cdft_iterations_ << ": lambda=" << lambda_mid
          << " population mismatch=" << mismatch_mid << std::flush;

      if (std::abs(mismatch_mid) < cdft_population_tolerance_) {
        constraint.lambda = lambda_mid;
        initial_guess_ = saved_initial_guess;
        XTP_LOG(Log::error, *pLog_)
            << TimeStamp() << " CDFT converged after " << outer_iter + 1
            << " outer iterations, lambda=" << lambda_mid << std::flush;
        return true;
      }

      if (mismatch_mid * mismatch_lo < 0.0) {
        lambda_hi = lambda_mid;
        mismatch_hi = mismatch_mid;
      } else {
        lambda_lo = lambda_mid;
        mismatch_lo = mismatch_mid;
      }
    }
  } catch (const std::runtime_error&) {
    initial_guess_ = saved_initial_guess;
    constraints_.clear();
    throw;
  }

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " RunCDFT: outer bisection loop did not converge "
                       "within "
      << max_cdft_iterations_ << " iterations." << std::flush;
  constraint.lambda = 0.5 * (lambda_lo + lambda_hi);
  initial_guess_ = saved_initial_guess;
  return false;
}

// Restricted SCF loop. The total energy is assembled as
//
//   E = Tr[P H0] + E_nuc + E_coul + E_xc + E_exx,
//
// with P = 2 C_occ C_occ^T. DIIS or mixing updates the density until both
// the energy change and the commutator error are converged.
bool DFTEngine::EvaluateClosedShell(
    Orbitals& orb, const Mat_p_Energy& H0,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) {

  tools::EigenSystem MOs;
  MOs.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());

  if (initial_guess_ == "orbfile") {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Reading guess from orbitals object/file"
        << std::flush;
    MOs = orb.MOs();
    MOs.eigenvectors() = OrthogonalizeGuess(MOs.eigenvectors());
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Setup Initial Guess using: " << initial_guess_
        << std::flush;
    if (initial_guess_ == "independent") {
      MOs = IndependentElectronGuess(H0);
    } else if (initial_guess_ == "atom") {
      MOs = ModelPotentialGuess(H0, orb.QMAtoms(), vxcpotential);
    } else if (initial_guess_ == "huckel") {
      MOs = ExtendedHuckelGuess(orb.QMAtoms());
    } else if (initial_guess_ == "huckel_dft") {
      MOs = ExtendedHuckelDFTGuess(H0, orb.QMAtoms(), vxcpotential);
    } else {
      throw std::runtime_error("Initial guess method not known/implemented");
    }
  }

  ConvergenceAcc::SpinDensity spin_dmat =
      conv_accelerator_.DensityMatrixSpinResolved(MOs);
  Eigen::MatrixXd Dmat = spin_dmat.total();

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Guess Matrix gives N=" << std::setprecision(9)
      << Dmat.cwiseProduct(dftAOoverlap_.Matrix()).sum() << " electrons."
      << std::flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " STARTING SCF cycle" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << " ----------------------------------------------"
         "----------------------------"
      << std::flush;

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(Dmat.rows(), Dmat.cols());
  Eigen::MatrixXd K;
  if (ScaHFX_ > 0) {
    K = Eigen::MatrixXd::Zero(Dmat.rows(), Dmat.cols());
  }

  double start_incremental_F_threshold = 1e-4;
  if (!auxbasis_name_.empty()) {
    start_incremental_F_threshold = 0.0;  // Disable if RI is used
  }
  IncrementalFockBuilder incremental_fock(*pLog_, start_incremental_F_threshold,
                                          fock_matrix_reset_);
  incremental_fock.Configure(Dmat);

  for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << max_iter_ << std::flush;

    Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled DFT Vxc matrix " << std::flush;

    Eigen::MatrixXd H = H0.matrix() + e_vxc.matrix();
    double Eone = Dmat.cwiseProduct(H0.matrix()).sum();
    double Etwo = e_vxc.energy();
    double exx = 0.0;

    incremental_fock.Start(this_iter, conv_accelerator_.getDIIsError());
    incremental_fock.resetMatrices(J, K, Dmat);
    incremental_fock.UpdateCriteria(conv_accelerator_.getDIIsError(),
                                    this_iter);

    double integral_error =
        std::min(conv_accelerator_.getDIIsError() * 1e-5, 1e-5);

    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> both = CalcERIs_EXX(
          MOs.eigenvectors(), incremental_fock.getDmat_diff(), integral_error);
      J += both[0];
      H += J;
      Etwo += 0.5 * Dmat.cwiseProduct(J).sum();
      K += both[1];
      H += 0.5 * ScaHFX_ * K;
      exx = 0.25 * ScaHFX_ * Dmat.cwiseProduct(K).sum();
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " Filled F+K matrix " << std::flush;
    } else {
      J += CalcERIs(incremental_fock.getDmat_diff(), integral_error);
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " Filled F matrix " << std::flush;
      H += J;
      Etwo += 0.5 * Dmat.cwiseProduct(J).sum();
    }

    Etwo += exx;
    double totenergy = Eone + H0.energy() + Etwo;

    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Single particle energy "
                               << std::setprecision(12) << Eone << std::flush;
    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Two particle energy "
                               << std::setprecision(12) << Etwo << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << std::setprecision(12) << " Local Exc contribution "
        << e_vxc.energy() << std::flush;
    if (ScaHFX_ > 0) {
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << std::setprecision(12)
          << " Non local Ex contribution " << exx << std::flush;
    }
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Total Energy " << std::setprecision(12) << totenergy
        << std::flush;

    Dmat = conv_accelerator_.Iterate(Dmat, H, MOs, totenergy);
    incremental_fock.UpdateDmats(Dmat, conv_accelerator_.getDIIsError(),
                                 this_iter);

    PrintMOs(MOs.eigenvalues(), Log::info);

    if (num_docc_ + num_socc_alpha_ > 0 &&
        num_docc_ + num_socc_alpha_ < MOs.eigenvalues().size()) {
      XTP_LOG(Log::info, *pLog_)
          << "\t\tGAP "
          << MOs.eigenvalues()(num_docc_ + num_socc_alpha_) -
                 MOs.eigenvalues()(num_docc_ + num_socc_alpha_ - 1)
          << std::flush;
    }

    if (conv_accelerator_.isConverged()) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Total Energy has converged to "
          << std::setprecision(9) << conv_accelerator_.getDeltaE()
          << "[Ha] after " << this_iter + 1
          << " iterations. DIIS error is converged up to "
          << conv_accelerator_.getDIIsError() << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Final Single Point Energy "
          << std::setprecision(12) << totenergy << " Ha" << std::flush;
      XTP_LOG(Log::error, *pLog_) << TimeStamp() << std::setprecision(12)
                                  << " Final Local Exc contribution "
                                  << e_vxc.energy() << " Ha" << std::flush;
      if (ScaHFX_ > 0) {
        XTP_LOG(Log::error, *pLog_) << TimeStamp() << std::setprecision(12)
                                    << " Final Non Local Ex contribution "
                                    << exx << " Ha" << std::flush;
      }

      PrintMOs(MOs.eigenvalues(), Log::error);

      Index nuclear_charge = 0;
      for (const QMAtom& atom : orb.QMAtoms()) {
        nuclear_charge += atom.getNuccharge();
      }

      orb.setQMEnergy(totenergy);
      orb.MOs() = MOs;
      orb.setNumberOfAlphaElectrons(num_alpha_electrons_);
      orb.setNumberOfBetaElectrons(num_beta_electrons_);
      orb.setNumberOfOccupiedLevels(num_docc_ + num_socc_alpha_);
      orb.setChargeAndSpin(
          nuclear_charge - numofelectrons_,
          std::abs(num_alpha_electrons_ - num_beta_electrons_) + 1);

      if (compute_forces_) {
        ComputeAndStoreForces(orb, Dmat, vxcpotential);
      }

      CalcElDipole(orb);
      return true;
    } else if (this_iter == max_iter_ - 1) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " DFT calculation has not converged after "
          << max_iter_
          << " iterations. Use more iterations or another convergence "
             "acceleration scheme."
          << std::flush;
      return false;
    }
  }

  return true;
}

// Unrestricted SCF loop. The alpha and beta channels are iterated through
// separate Fock matrices
//
//   F^alpha = H0 + J[P^alpha + P^beta] + V_xc^alpha + K^alpha
//   F^beta  = H0 + J[P^alpha + P^beta] + V_xc^beta  + K^beta,
//
// while the total energy uses the spin-summed one-electron and Coulomb terms
// together with spin-resolved XC and exact-exchange contributions.
bool DFTEngine::EvaluateUKS(Orbitals& orb, const Mat_p_Energy& H0,
                            const Vxc_Potential<Vxc_Grid>& vxcpotential) {
  tools::EigenSystem MOs_alpha;
  tools::EigenSystem MOs_beta;

  MOs_alpha.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs_alpha.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());
  MOs_beta.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs_beta.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());

  UKSConvergenceAcc conv_uks;

  ConvergenceAcc::options opt_alpha = conv_opt_;
  opt_alpha.mode = ConvergenceAcc::KSmode::open;
  opt_alpha.numberofelectrons = num_alpha_electrons_;

  ConvergenceAcc::options opt_beta = conv_opt_;
  opt_beta.mode = ConvergenceAcc::KSmode::open;
  opt_beta.numberofelectrons = num_beta_electrons_;

  conv_uks.Configure(opt_alpha, opt_beta);
  conv_uks.setLogger(pLog_);
  conv_uks.setOverlap(dftAOoverlap_, 1e-8);

  if (initial_guess_ == "orbfile") {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Reading UKS guess from orbitals object/file"
        << std::flush;

    MOs_alpha = orb.MOs();
    MOs_alpha.eigenvectors() = OrthogonalizeGuess(MOs_alpha.eigenvectors());

    if (orb.hasBetaMOs()) {
      MOs_beta = orb.MOs_beta();
      MOs_beta.eigenvectors() = OrthogonalizeGuess(MOs_beta.eigenvectors());
    } else {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp()
          << " Orbital file has no beta MOs, using alpha guess for beta."
          << std::flush;
      MOs_beta = MOs_alpha;
    }
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Setup UKS Initial Guess using: " << initial_guess_
        << std::flush;

    tools::EigenSystem guess;
    if (initial_guess_ == "independent") {
      guess = IndependentElectronGuess(H0);
    } else if (initial_guess_ == "atom") {
      guess = ModelPotentialGuess(H0, orb.QMAtoms(), vxcpotential);
    } else if (initial_guess_ == "huckel") {
      guess = ExtendedHuckelGuess(orb.QMAtoms());
    } else if (initial_guess_ == "huckel_dft") {
      guess = ExtendedHuckelDFTGuess(H0, orb.QMAtoms(), vxcpotential);
    } else {
      throw std::runtime_error("Initial guess method not known/implemented");
    }

    MOs_alpha = guess;
    MOs_beta = guess;
  }

  // Build the initial spin densities P^alpha and P^beta from the chosen
  // starting orbitals before entering the coupled UKS iterations.
  UKSConvergenceAcc::SpinDensity Dspin =
      conv_uks.DensityMatrix(MOs_alpha, MOs_beta);

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " UKS guess gives Nalpha="
      << Dspin.alpha.cwiseProduct(dftAOoverlap_.Matrix()).sum()
      << " Nbeta=" << Dspin.beta.cwiseProduct(dftAOoverlap_.Matrix()).sum()
      << " Ntot=" << Dspin.total().cwiseProduct(dftAOoverlap_.Matrix()).sum()
      << std::flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " STARTING UKS SCF cycle" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << " ------------------------------------------------------------"
      << std::flush;

  for (Index this_iter = 0; this_iter < max_iter_; ++this_iter) {
    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << max_iter_ << std::flush;

    Eigen::MatrixXd H_alpha = H0.matrix();
    Eigen::MatrixXd H_beta = H0.matrix();

    // The Coulomb contribution depends only on the total density
    // P = P^alpha + P^beta, while exchange and XC remain spin resolved.
    const Eigen::MatrixXd D_total = Dspin.total();

    double E_one = Dspin.alpha.cwiseProduct(H0.matrix()).sum() +
                   Dspin.beta.cwiseProduct(H0.matrix()).sum();

    double E_coul = 0.0;
    double E_xc = 0.0;
    double E_exx = 0.0;

    double integral_error = std::min(conv_uks.getDIIsError() * 1e-5, 1e-5);

    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> both_alpha = CalcERIs_EXX(
          Eigen::MatrixXd::Zero(0, 0), Dspin.alpha, integral_error);
      std::array<Eigen::MatrixXd, 2> both_beta =
          CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), Dspin.beta, integral_error);

      Eigen::MatrixXd J = both_alpha[0] + both_beta[0];
      Eigen::MatrixXd K_alpha = both_alpha[1];
      Eigen::MatrixXd K_beta = both_beta[1];

      H_alpha += J + ScaHFX_ * K_alpha;
      H_beta += J + ScaHFX_ * K_beta;

      E_coul = 0.5 * D_total.cwiseProduct(J).sum();
      E_exx = 0.5 * ScaHFX_ *
              (Dspin.alpha.cwiseProduct(K_alpha).sum() +
               Dspin.beta.cwiseProduct(K_beta).sum());
    } else {
      Eigen::MatrixXd J = CalcERIs(D_total, integral_error);
      H_alpha += J;
      H_beta += J;
      E_coul = 0.5 * D_total.cwiseProduct(J).sum();
    }

    auto vxc = vxcpotential.IntegrateVXCSpin(Dspin.alpha, Dspin.beta);
    H_alpha += vxc.vxc_alpha;
    H_beta += vxc.vxc_beta;
    E_xc = vxc.energy;

    double totenergy = H0.energy() + E_one + E_coul + E_xc + E_exx;

    // CDFT constraint potential -- deliberately the LAST term added to
    // either Fock matrix, and gated by a single, cheap .empty() check:
    // for any standard, non-CDFT run (constraints_ left at its default,
    // empty state), this entire block is skipped, and both H_alpha and
    // H_beta are built exactly as they always were -- no measurable
    // overhead, no behavior change whatsoever. Adds
    // lambda_c * spin_alpha/beta_coefficient * W_c to the respective
    // Fock matrix for every active constraint c (a charge constraint
    // uses +1/+1, adding the identical potential to both channels; a
    // future spin constraint would use +1/-1 -- see Constraint's own
    // comment in hirshfeldpartition.h for why these are stored
    // separately rather than this code assuming "charge" specifically),
    // and the corresponding correction term to the reported total
    // energy: E_CDFT = E_KS + sum_c lambda_c * (N_c^computed -
    // N_c^target), the standard Wu-Van Voorhis Lagrangian.
    if (!constraints_.empty()) {
      for (const HirshfeldPartition::Constraint& c : constraints_) {
        H_alpha += (c.lambda * c.spin_alpha_coefficient) * c.weight_matrix;
        H_beta += (c.lambda * c.spin_beta_coefficient) * c.weight_matrix;
        double population =
            c.spin_alpha_coefficient *
                Dspin.alpha.cwiseProduct(c.weight_matrix).sum() +
            c.spin_beta_coefficient *
                Dspin.beta.cwiseProduct(c.weight_matrix).sum();
        totenergy += c.lambda * (population - c.target_population);
      }
    }

    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " One particle energy "
                               << std::setprecision(12) << E_one << std::flush;
    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Coulomb contribution "
                               << std::setprecision(12) << E_coul << std::flush;
    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " XC contribution "
                               << std::setprecision(12) << E_xc << std::flush;
    if (ScaHFX_ > 0) {
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " EXX contribution " << std::setprecision(12)
          << E_exx << std::flush;
    }
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Total Energy " << std::setprecision(12) << totenergy
        << std::flush;

    UKSConvergenceAcc::SpinFock Hspin{H_alpha, H_beta};
    Dspin = conv_uks.Iterate(Dspin, Hspin, MOs_alpha, MOs_beta, totenergy);
    if (force_uks_path_ && num_alpha_electrons_ == num_beta_electrons_) {
      MOs_beta = MOs_alpha;
      Dspin.beta = Dspin.alpha;
    }

    XTP_LOG(Log::info, *pLog_)
        << TimeStamp()
        << " Nalpha=" << Dspin.alpha.cwiseProduct(dftAOoverlap_.Matrix()).sum()
        << " Nbeta=" << Dspin.beta.cwiseProduct(dftAOoverlap_.Matrix()).sum()
        << std::flush;

    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " <Sz> = "
        << 0.5 * double(num_alpha_electrons_ - num_beta_electrons_)
        << std::flush;

    PrintMOsUKS(MOs_alpha.eigenvalues(), MOs_beta.eigenvalues(), Log::info);

    if (conv_uks.isConverged()) {
      Index nuclear_charge = 0;
      for (const QMAtom& atom : orb.QMAtoms()) {
        nuclear_charge += atom.getNuccharge();
      }

      CanonicalizeOrbitalPhases(MOs_alpha);
      CanonicalizeOrbitalPhases(MOs_beta);

      orb.setQMEnergy(totenergy);
      orb.MOs() = MOs_alpha;
      orb.MOs_beta() = MOs_beta;
      orb.setNumberOfAlphaElectrons(num_alpha_electrons_);
      orb.setNumberOfBetaElectrons(num_beta_electrons_);
      orb.setNumberOfOccupiedLevels(num_alpha_electrons_);
      orb.setNumberOfOccupiedLevelsBeta(num_beta_electrons_);
      orb.setChargeAndSpin(
          nuclear_charge - numofelectrons_,
          std::abs(num_alpha_electrons_ - num_beta_electrons_) + 1);

      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " UKS converged after " << this_iter + 1
          << " iterations. Delta E=" << conv_uks.getDeltaE()
          << " DIIS error=" << conv_uks.getDIIsError() << std::flush;

      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Final Single Point Energy "
          << std::setprecision(12) << totenergy << " Ha" << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << std::setprecision(12) << " Final XC contribution "
          << E_xc << " Ha" << std::flush;
      if (ScaHFX_ > 0) {
        XTP_LOG(Log::error, *pLog_)
            << TimeStamp() << std::setprecision(12)
            << " Final EXX contribution " << E_exx << " Ha" << std::flush;
      }

      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " <Sz> = "
          << 0.5 * double(num_alpha_electrons_ - num_beta_electrons_)
          << std::flush;

      PrintMOsUKS(MOs_alpha.eigenvalues(), MOs_beta.eigenvalues(), Log::error);

      if (compute_forces_) {
        ComputeAndStoreForcesUKS(orb, Dspin, MOs_alpha, MOs_beta, vxcpotential);
      }

      CalcElDipole(orb);
      return true;
    }

    if (this_iter == max_iter_ - 1) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " UKS calculation has not converged after "
          << max_iter_ << " iterations." << std::flush;
      return false;
    }
  }

  return false;
}

// One-electron core Hamiltonian and its constant energy offset.
//
// The matrix part is
//
//   H0 = T + V_nuc + V_ECP + V_ext,
//
// while the scalar energy collects all nucleus-nucleus and nucleus-external
// interaction terms that do not depend on the electronic density.
Mat_p_Energy DFTEngine::SetupH0(const QMMolecule& mol) const {

  AOKinetic dftAOkinetic;

  dftAOkinetic.Fill(dftbasis_);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Kinetic energy matrix ." << std::flush;

  AOMultipole dftAOESP;
  dftAOESP.FillPotential(dftbasis_, mol);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT nuclear potential matrix." << std::flush;

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Constructed independent particle hamiltonian "
      << std::flush;
  double E0 = NuclearRepulsion(mol);
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Nuclear Repulsion Energy is "
                              << std::setprecision(9) << E0 << std::flush;

  if (!ecp_name_.empty()) {
    AOECP dftAOECP;
    dftAOECP.FillPotential(dftbasis_, ecp_);
    H0 += dftAOECP.Matrix();
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled DFT ECP matrix" << std::flush;
  }

  if (externalsites_ != nullptr) {
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " " << externalsites_->size()
                                << " External sites" << std::flush;
    bool has_quadrupoles = std::any_of(
        externalsites_->begin(), externalsites_->end(),
        [](const std::unique_ptr<StaticSite>& s) { return s->getRank() == 2; });
    std::string header =
        " Name      Coordinates[a0]     charge[e]         dipole[e*a0]    ";
    if (has_quadrupoles) {
      header += "              quadrupole[e*a0^2]";
    }
    XTP_LOG(Log::error, *pLog_) << header << std::flush;
    Index limit = 50;
    Index counter = 0;
    for (const std::unique_ptr<StaticSite>& site : *externalsites_) {
      if (counter == limit) {
        break;
      }
      std::string output =
          (boost::format("  %1$s"
                         "   %2$+1.4f %3$+1.4f %4$+1.4f"
                         "   %5$+1.4f") %
           site->getElement() % site->getPos()[0] % site->getPos()[1] %
           site->getPos()[2] % site->getCharge())
              .str();
      const Eigen::Vector3d& dipole = site->getDipole();
      output += (boost::format("   %1$+1.4f %2$+1.4f %3$+1.4f") % dipole[0] %
                 dipole[1] % dipole[2])
                    .str();
      if (site->getRank() > 1) {
        Eigen::VectorXd quadrupole = site->Q().tail<5>();
        output +=
            (boost::format("   %1$+1.4f %2$+1.4f %3$+1.4f %4$+1.4f %5$+1.4f") %
             quadrupole[0] % quadrupole[1] % quadrupole[2] % quadrupole[3] %
             quadrupole[4])
                .str();
      }
      XTP_LOG(Log::error, *pLog_) << output << std::flush;
      counter++;
    }
    if (counter == limit) {
      XTP_LOG(Log::error, *pLog_)
          << "              ... (" << externalsites_->size() - limit
          << " sites not displayed)\n"
          << std::flush;
    }

    Mat_p_Energy ext_multipoles =
        IntegrateExternalMultipoles(mol, *externalsites_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Nuclei-external site interaction energy "
        << std::setprecision(9) << ext_multipoles.energy() << std::flush;
    E0 += ext_multipoles.energy();
    H0 += ext_multipoles.matrix();
  }

  if (integrate_ext_density_) {
    Orbitals extdensity;
    extdensity.ReadFromCpt(orbfilename_);
    Mat_p_Energy extdensity_result = IntegrateExternalDensity(mol, extdensity);
    E0 += extdensity_result.energy();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Nuclei-external density interaction energy "
        << std::setprecision(9) << extdensity_result.energy() << std::flush;
    H0 += extdensity_result.matrix();
  }

  if (integrate_ext_field_) {

    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Integrating external electric field with F[Hrt]="
        << extfield_.transpose() << std::flush;
    H0 += IntegrateExternalField(mol);
  }

  return Mat_p_Energy(E0, H0);
}

// Precompute SCF-invariant matrices: overlap for the generalized eigenvalue
// problem and the RI/4c electron-repulsion backend that later yields J[P] and
// K[P].
void DFTEngine::SetupInvariantMatrices() {
  dftAOoverlap_.Fill(dftbasis_);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Overlap matrix." << std::flush;

  conv_opt_.numberofelectrons = numofelectrons_;
  conv_opt_.number_alpha_electrons = num_alpha_electrons_;
  conv_opt_.number_beta_electrons = num_beta_electrons_;
  conv_opt_.mode = (num_alpha_electrons_ == num_beta_electrons_)
                       ? ConvergenceAcc::KSmode::closed
                       : ConvergenceAcc::KSmode::restricted_open;
  conv_accelerator_.Configure(conv_opt_);
  conv_accelerator_.setLogger(pLog_);
  conv_accelerator_.setOverlap(dftAOoverlap_, 1e-8);
  conv_accelerator_.PrintConfigOptions();

  if (!auxbasis_name_.empty()) {
    // prepare invariant part of electron repulsion integrals
    ERIs_.Initialize(dftbasis_, auxbasis_);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Inverted AUX Coulomb matrix, removed "
        << ERIs_.Removedfunctions() << " functions from aux basis"
        << std::flush;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Setup invariant parts of Electron Repulsion integrals "
        << std::flush;
  } else {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculating 4c diagonals. " << std::flush;
    ERIs_.Initialize_4c(dftbasis_);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculated 4c diagonals. " << std::flush;
  }

  return;
}

namespace {
// Hund's-rule ground-state (alpha electrons, beta electrons) for the
// main-group (s/p-block) elements most relevant to organic systems --
// H through Kr, plus the heavier halogens (Br, I) via their own,
// separately-computed period-5 entries. Explicitly does NOT cover
// d-block (Sc-Zn, Y-Cd) or f-block elements: the d^n s^2 vs d^(n+1) s^1
// (and worse, f-block) ground-state competition is genuinely subtle
// and functional-dependent -- exactly why CP2K's own isolated-atom
// ("ATOM") program requires explicit, manual per-subshell occupation
// specification rather than trusting any automatic rule (confirmed
// directly: HORTON's own CP2K pro-atom documentation states "The ATOM
// program of CP2K does not simply follow the Aufbau rule to assign
// orbital occupations"). Returns std::nullopt for anything not
// explicitly covered, so callers can fall back to the existing,
// simpler parity-based logic with a clear warning rather than silently
// guessing.
//
// Method: standard Aufbau filling order (1s,2s,2p,3s,3p,4s,3d,4p,5s,
// 4d,5p) up to (but explicitly skipping) each d-block range, applying
// Hund's rule within any open p subshell (spread across all 3 p
// orbitals with parallel/majority spin first, only pairing once every
// orbital in that subshell already has one) -- for p^n, n<=3 gives n
// alpha/0 beta in that subshell; n>3 gives 3 alpha/(n-3) beta. Every
// entry below was computed by hand from this rule and can be checked
// against any standard table of atomic ground-state term symbols
// (all are unambiguous, textbook Hund's-rule cases for main-group
// atoms -- no functional-dependent ambiguity of the kind that affects
// d/f-block).
std::optional<std::pair<Index, Index>> HundsRuleAlphaBetaElectrons(
    Index nuclear_charge) {
  switch (nuclear_charge) {
    case 1:  return std::make_pair(1, 0);    // H:  1s1
    case 2:  return std::make_pair(1, 1);    // He: 1s2
    case 3:  return std::make_pair(2, 1);    // Li: [He] 2s1
    case 4:  return std::make_pair(2, 2);    // Be: 2s2
    case 5:  return std::make_pair(3, 2);    // B:  2p1
    case 6:  return std::make_pair(4, 2);    // C:  2p2 (2a)
    case 7:  return std::make_pair(5, 2);    // N:  2p3 (3a)
    case 8:  return std::make_pair(5, 3);    // O:  2p4 (3a+1b)
    case 9:  return std::make_pair(5, 4);    // F:  2p5 (3a+2b)
    case 10: return std::make_pair(5, 5);    // Ne: 2p6
    case 11: return std::make_pair(6, 5);    // Na: [Ne] 3s1
    case 12: return std::make_pair(6, 6);    // Mg: 3s2
    case 13: return std::make_pair(7, 6);    // Al: 3p1
    case 14: return std::make_pair(8, 6);    // Si: 3p2 (2a)
    case 15: return std::make_pair(9, 6);    // P:  3p3 (3a)
    case 16: return std::make_pair(9, 7);    // S:  3p4 (3a+1b)
    case 17: return std::make_pair(9, 8);    // Cl: 3p5 (3a+2b)
    case 18: return std::make_pair(9, 9);    // Ar: 3p6
    case 19: return std::make_pair(10, 9);   // K:  [Ar] 4s1
    case 20: return std::make_pair(10, 10);  // Ca: 4s2
    // 21-30 (Sc-Zn): 3d block -- deliberately NOT covered.
    case 31: return std::make_pair(16, 15);  // Ga: [Zn] 4p1
    case 32: return std::make_pair(17, 15);  // Ge: 4p2 (2a)
    case 33: return std::make_pair(18, 15);  // As: 4p3 (3a)
    case 34: return std::make_pair(18, 16);  // Se: 4p4 (3a+1b)
    case 35: return std::make_pair(18, 17);  // Br: 4p5 (3a+2b)
    case 36: return std::make_pair(18, 18);  // Kr: 4p6
    // 39-48 (Y-Cd): 4d block -- deliberately NOT covered.
    case 49: return std::make_pair(25, 24);  // In: [Cd] 5p1
    case 50: return std::make_pair(26, 24);  // Sn: 5p2 (2a)
    case 51: return std::make_pair(27, 24);  // Sb: 5p3 (3a)
    case 52: return std::make_pair(27, 25);  // Te: 5p4 (3a+1b)
    case 53: return std::make_pair(27, 26);  // I:  5p5 (3a+2b)
    case 54: return std::make_pair(27, 27);  // Xe: 5p6
    default: return std::nullopt;
  }
}
}  // namespace

Eigen::MatrixXd DFTEngine::RunAtomicDFT_unrestricted(
    const QMAtom& uniqueAtom, bool use_hunds_rule_occupation) const {
  bool with_ecp = !ecp_name_.empty();
  if (uniqueAtom.getElement() == "H" || uniqueAtom.getElement() == "He") {
    with_ecp = false;
  }

  QMMolecule atom = QMMolecule("individual_atom", 0);
  atom.push_back(uniqueAtom);

  BasisSet basisset;
  basisset.Load(dftbasis_name_);
  AOBasis dftbasis;
  dftbasis.Fill(basisset, atom);
  Vxc_Grid grid;
  grid.GridSetup(grid_name_, atom, dftbasis);
  Vxc_Potential<Vxc_Grid> gridIntegration(grid);
  gridIntegration.setXCfunctional(xc_functional_name_);

  ECPAOBasis ecp;
  if (with_ecp) {
    ECPBasisSet ecps;
    ecps.Load(ecp_name_);
    ecp.Fill(ecps, atom);
  }

  Index numofelectrons = uniqueAtom.getNuccharge();
  Index alpha_e = 0;
  Index beta_e = 0;

  // Deliberately opt-in, defaulting to false: this changes ONLY which
  // total alpha/beta split is used for the reference atom's own SCF,
  // not the SphericalAverageShells step below (kept unconditionally,
  // for both modes) -- the existing SAD-initial-guess caller
  // (AtomicGuess) is not changed at all by this parameter existing, and
  // continues to use the simpler, parity-based split exactly as
  // before. A physically correct free-atom ground state is not needed
  // for a good SCF starting guess (the molecule's own overall spin
  // state, and the SCF that follows, will reshape this regardless);
  // it matters for the promolecular reference densities Hirshfeld-based
  // CDFT will need instead, which is what this parameter exists for.
  if (use_hunds_rule_occupation) {
    auto hunds_rule = HundsRuleAlphaBetaElectrons(numofelectrons);
    if (hunds_rule.has_value()) {
      alpha_e = hunds_rule->first;
      beta_e = hunds_rule->second;
    } else {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp() << " No Hund's-rule ground-state occupation table "
                            "entry for nuclear charge "
          << numofelectrons
          << " (d/f-block elements are not covered -- see "
             "HundsRuleAlphaBetaElectrons's own comment for why) -- "
             "falling back to the simpler, parity-based alpha/beta split."
          << std::flush;
      use_hunds_rule_occupation = false;
    }
  }
  if (!use_hunds_rule_occupation) {
    if ((numofelectrons % 2) != 0) {
      alpha_e = numofelectrons / 2 + numofelectrons % 2;
      beta_e = numofelectrons / 2;
    } else {
      alpha_e = numofelectrons / 2;
      beta_e = alpha_e;
    }
  }

  AOOverlap dftAOoverlap;
  AOKinetic dftAOkinetic;
  AOMultipole dftAOESP;
  AOECP dftAOECP;
  ERIs ERIs_atom;

  dftAOoverlap.Fill(dftbasis);
  dftAOkinetic.Fill(dftbasis);

  dftAOESP.FillPotential(dftbasis, atom);
  ERIs_atom.Initialize_4c(dftbasis);

  UKSConvergenceAcc conv_uks;
  ConvergenceAcc::options opt_alpha = conv_opt_;
  opt_alpha.mode = ConvergenceAcc::KSmode::open;
  opt_alpha.histlength = 20;
  opt_alpha.levelshift = 0.1;
  opt_alpha.levelshiftend = 0.0;
  opt_alpha.usediis = true;
  // adiis_start/diis_start deliberately NOT overridden here (previously
  // both hardcoded to 0.0) -- confirmed via ConvergenceAcc::Iterate's
  // own gating logic (the "diiserror_ < opt_.adiis_start ||
  // diiserror_ < opt_.diis_start" check) that 0.0 makes this condition
  // permanently false, since diiserror_ is a norm and can never be
  // negative. That silently disabled BOTH DIIS and ADIIS for the
  // entire atomic SCF regardless of usediis=true just above -- an
  // internal inconsistency, not an intentional design choice -- and is
  // the likely root cause of this function's own, separately reported
  // slow convergence (falling back to plain, level-shift-damped linear
  // mixing every single iteration, with no (A)DIIS extrapolation ever
  // actually applied). Inheriting these two fields from conv_opt_
  // (already the base for opt_alpha above, via "= conv_opt_") matches
  // the main UKS SCF path's own behavior, which never overrides them
  // at all.
  opt_alpha.numberofelectrons = alpha_e;

  ConvergenceAcc::options opt_beta = opt_alpha;
  opt_beta.numberofelectrons = beta_e;

  Logger log;
  // Single, shared accelerator -- previously two fully independent
  // ConvergenceAcc instances (Convergence_alpha, Convergence_beta),
  // each with its own separate DIIS history/error tracking/coefficient
  // calculation, with no coupling between the two spin channels at
  // all. The main UKS SCF path (DFTEngine::EvaluateUKS) instead uses
  // exactly this UKSConvergenceAcc class, which builds ONE combined
  // error metric from both spin channels together and applies a
  // single, jointly-derived set of (A)DIIS coefficients to both --
  // confirmed directly from uks_convergenceacc.cc's own comment ("one
  // shared DIIS/ADIIS history length") and its Iterate()'s own
  // "diis_.Update(maxerrorindex_, err_alpha, err_beta)" call. This is
  // the standard, textbook-correct formulation of UKS DIIS; the
  // previous two-independent-accelerators approach was not wrong in
  // the sense of being internally inconsistent (unlike the
  // adiis_start/diis_start bug above), but it let each spin channel's
  // extrapolation disagree with the other's, which is not how UKS DIIS
  // is meant to work.
  conv_uks.Configure(opt_alpha, opt_beta);
  conv_uks.setLogger(&log);
  conv_uks.setOverlap(dftAOoverlap, 1e-8);

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  if (with_ecp) {
    dftAOECP.FillPotential(dftbasis, ecp);
    H0 += dftAOECP.Matrix();
  }

  tools::EigenSystem MOs_alpha = conv_uks.SolveFockmatrix(H0);

  if (uniqueAtom.getElement() == "H") {
    // H has no beta electrons at all (beta_e=0 above) -- nocclevels_beta_
    // will be 0 once Configure() runs, and
    // DensityMatrixGroundState_unres already returns a zero matrix
    // whenever nocclevels==0 regardless of what MOs are passed in, so
    // reusing MOs_alpha as a dummy beta argument here is safe: only
    // Dspin_H.alpha is ever used below.
    UKSConvergenceAcc::SpinDensity Dspin_H =
        conv_uks.DensityMatrix(MOs_alpha, MOs_alpha);
    return Dspin_H.alpha;
  }

  tools::EigenSystem MOs_beta = conv_uks.SolveFockmatrix(H0);
  UKSConvergenceAcc::SpinDensity Dspin = conv_uks.DensityMatrix(MOs_alpha, MOs_beta);

  Index maxiter = 80;
  for (Index this_iter = 0; this_iter < maxiter; this_iter++) {
    Eigen::MatrixXd H_alpha = H0;
    Eigen::MatrixXd H_beta = H0;

    double E_coul = 0.0;
    double E_exx = 0.0;
    double E_xc = 0.0;

    // Matches EvaluateUKS's own formula exactly (a single combined
    // DIIS error now, not an average of two independent ones).
    double integral_error = std::min(conv_uks.getDIIsError() * 1e-5, 1e-5);

    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> both_alpha =
          ERIs_atom.CalculateERIs_EXX_4c(Dspin.alpha, integral_error);
      std::array<Eigen::MatrixXd, 2> both_beta =
          ERIs_atom.CalculateERIs_EXX_4c(Dspin.beta, integral_error);

      Eigen::MatrixXd Hartree = both_alpha[0] + both_beta[0];
      H_alpha += Hartree + ScaHFX_ * both_alpha[1];
      H_beta += Hartree + ScaHFX_ * both_beta[1];

      E_coul = 0.5 * Dspin.total().cwiseProduct(Hartree).sum();
      E_exx = 0.5 * ScaHFX_ *
              (both_alpha[1].cwiseProduct(Dspin.alpha).sum() +
               both_beta[1].cwiseProduct(Dspin.beta).sum());
    } else {
      Eigen::MatrixXd Hartree =
          ERIs_atom.CalculateERIs_4c(Dspin.total(), integral_error);
      H_alpha += Hartree;
      H_beta += Hartree;
      E_coul = 0.5 * Dspin.total().cwiseProduct(Hartree).sum();
    }

    auto vxc = gridIntegration.IntegrateVXCSpin(Dspin.alpha, Dspin.beta);
    H_alpha += vxc.vxc_alpha;
    H_beta += vxc.vxc_beta;
    E_xc = vxc.energy;

    double E_one_alpha = Dspin.alpha.cwiseProduct(H0).sum();
    double E_one_beta = Dspin.beta.cwiseProduct(H0).sum();
    double totenergy = E_one_alpha + E_one_beta + E_coul + E_exx + E_xc;

    UKSConvergenceAcc::SpinFock Hspin{H_alpha, H_beta};
    Dspin = conv_uks.Iterate(Dspin, Hspin, MOs_alpha, MOs_beta, totenergy);

    XTP_LOG(Log::debug, *pLog_)
        << TimeStamp() << " Iter " << this_iter << " of " << maxiter << " Etot "
        << totenergy << " diise " << conv_uks.getDIIsError() << "\n\t\t a_gap "
        << MOs_alpha.eigenvalues()(alpha_e) -
               MOs_alpha.eigenvalues()(alpha_e - 1)
        << " b_gap "
        << MOs_beta.eigenvalues()(beta_e) - MOs_beta.eigenvalues()(beta_e - 1)
        << " Nalpha=" << dftAOoverlap.Matrix().cwiseProduct(Dspin.alpha).sum()
        << " Nbeta=" << dftAOoverlap.Matrix().cwiseProduct(Dspin.beta).sum()
        << std::flush;

    bool converged = conv_uks.isConverged();
    if (converged || this_iter == maxiter - 1) {
      if (converged) {
        XTP_LOG(Log::info, *pLog_)
            << TimeStamp() << " Converged after " << this_iter + 1
            << " iterations" << std::flush;
      } else {
        XTP_LOG(Log::info, *pLog_)
            << TimeStamp() << " Not converged after " << this_iter + 1
            << " iterations. Unconverged density.\n\t\t\t"
            << " DIIsError=" << conv_uks.getDIIsError() << std::flush;
      }
      break;
    }
  }

  Eigen::MatrixXd avgmatrix =
      SphericalAverageShells(Dspin.total(), dftbasis);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Atomic density Matrix for " << uniqueAtom.getElement()
      << " gives N=" << std::setprecision(9)
      << avgmatrix.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons."
      << std::flush;
  return avgmatrix;
}

Eigen::MatrixXd DFTEngine::AtomicGuess(const QMMolecule& mol) const {

  std::vector<std::string> elements = mol.FindUniqueElements();
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Scanning molecule of size " << mol.size()
      << " for unique elements" << std::flush;
  QMMolecule uniqueelements = QMMolecule("uniqueelements", 0);
  for (auto element : elements) {
    uniqueelements.push_back(QMAtom(0, element, Eigen::Vector3d::Zero()));
  }

  XTP_LOG(Log::info, *pLog_) << TimeStamp() << " " << uniqueelements.size()
                             << " unique elements found" << std::flush;
  std::vector<Eigen::MatrixXd> uniqueatom_guesses;
  for (QMAtom& unique_atom : uniqueelements) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Calculating atom density for "
        << unique_atom.getElement() << std::flush;
    Eigen::MatrixXd dmat_unrestricted = RunAtomicDFT_unrestricted(unique_atom);
    uniqueatom_guesses.push_back(dmat_unrestricted);
  }

  Eigen::MatrixXd guess =
      Eigen::MatrixXd::Zero(dftbasis_.AOBasisSize(), dftbasis_.AOBasisSize());
  Index start = 0;
  for (const QMAtom& atom : mol) {
    Index index = 0;
    for (; index < uniqueelements.size(); index++) {
      if (atom.getElement() == uniqueelements[index].getElement()) {
        break;
      }
    }
    Eigen::MatrixXd& dmat_unrestricted = uniqueatom_guesses[index];
    guess.block(start, start, dmat_unrestricted.rows(),
                dmat_unrestricted.cols()) = dmat_unrestricted;
    start += dmat_unrestricted.rows();
  }

  return guess;
}

std::map<std::string, Eigen::MatrixXd>
DFTEngine::ComputeHirshfeldReferenceDensities(const QMMolecule& mol) const {
  std::vector<std::string> elements = mol.FindUniqueElements();
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Scanning molecule of size " << mol.size()
      << " for unique elements (Hirshfeld reference densities)"
      << std::flush;

  std::map<std::string, Eigen::MatrixXd> reference_densities;
  for (const std::string& element : elements) {
    QMAtom unique_atom(0, element, Eigen::Vector3d::Zero());
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Calculating Hirshfeld reference density for "
        << element << std::flush;
    // use_hunds_rule_occupation=true unconditionally here -- this is
    // the one and only caller that should ever request it; AtomicGuess
    // just above, the pre-existing SAD-guess caller, never does.
    reference_densities[element] =
        RunAtomicDFT_unrestricted(unique_atom, /*use_hunds_rule_occupation=*/true);
  }
  return reference_densities;
}

HirshfeldPartition::Constraint DFTEngine::BuildCDFTConstraint(
    const QMMolecule& mol, const CDFTConstraintSpec& spec) const {
  std::map<std::string, Eigen::MatrixXd> reference_densities =
      ComputeHirshfeldReferenceDensities(mol);

  AOBasis full_dftbasis;
  {
    BasisSet basisset;
    basisset.Load(dftbasis_name_);
    full_dftbasis.Fill(basisset, mol);
  }

  Vxc_Grid grid;
  grid.GridSetup(grid_name_, mol, full_dftbasis);

  std::vector<HirshfeldPartition::AtomicReference> atoms =
      HirshfeldPartition::BuildAtomicReferences(mol, dftbasis_name_,
                                                reference_densities);

  HirshfeldPartition::Constraint constraint;
  constraint.weight_matrix = Eigen::MatrixXd::Zero(full_dftbasis.AOBasisSize(),
                                                    full_dftbasis.AOBasisSize());
  double neutral_reference_population = 0.0;
  for (Index atom_index : spec.atom_indices) {
    if (atom_index < 0 || atom_index >= static_cast<Index>(mol.size())) {
      throw std::runtime_error(
          "BuildCDFTConstraint: cdft.indices contains atom index " +
          std::to_string(atom_index) + ", but this molecule only has " +
          std::to_string(mol.size()) +
          " atoms (0-based indexing -- valid range is 0.." +
          std::to_string(mol.size() - 1) + ").");
    }
    // Hirshfeld weights are additive across atoms in a fragment --
    // w_fragment(r) = sum_{i in fragment} w_i(r) -- so the fragment's
    // own weight matrix is just the sum of each atom's own
    // BuildWeightMatrix result, and the neutral reference population
    // (needed to convert the options file's charge-relative target
    // into RunCDFT's own absolute-population convention) is just the
    // sum of the fragment atoms' own nuclear charges.
    constraint.weight_matrix += HirshfeldPartition::BuildWeightMatrix(
        atoms, atom_index, full_dftbasis, grid);
    neutral_reference_population +=
        static_cast<double>(mol[atom_index].getNuccharge());
  }

  constraint.target_population =
      neutral_reference_population - spec.target_charge;
  constraint.lambda = spec.initial_lambda;
  constraint.spin_alpha_coefficient = 1.0;
  constraint.spin_beta_coefficient = 1.0;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " CDFT constraint: " << spec.atom_indices.size()
      << " atom(s), neutral reference population="
      << neutral_reference_population << ", requested relative charge="
      << spec.target_charge
      << ", absolute target population=" << constraint.target_population
      << std::flush;

  return constraint;
}

void DFTEngine::ConfigOrbfile(Orbitals& orb) {
  if (initial_guess_ == "orbfile") {

    if (orb.hasDFTbasisName()) {
      if (orb.getDFTbasisName() != dftbasis_name_) {
        throw std::runtime_error(
            (boost::format("Basisset Name in guess orb file "
                           "and in dftengine option file differ %1% vs %2%") %
             orb.getDFTbasisName() % dftbasis_name_)
                .str());
      }
    } else {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " WARNING: "
             "Orbital file has no basisset information,"
             "using it as a guess might work or not for calculation with "
          << dftbasis_name_ << std::flush;
    }
  }

  const Index target_charge = orb.getCharge();
  const Index multiplicity = orb.getSpin();

  orb.setChargeAndSpin(target_charge, multiplicity);
  orb.setNumberOfAlphaElectrons(num_alpha_electrons_);
  orb.setNumberOfBetaElectrons(num_beta_electrons_);

  orb.setXCFunctionalName(xc_functional_name_);
  orb.setXCGrid(grid_name_);
  orb.setScaHFX(ScaHFX_);
  if (!ecp_name_.empty()) {
    orb.setECPName(ecp_name_);
  }
  if (!auxbasis_name_.empty()) {
    orb.SetupAuxBasis(auxbasis_name_);
  }

  if (initial_guess_ == "orbfile") {
    if (orb.hasECPName() || !ecp_name_.empty()) {
      if (orb.getECPName() != ecp_name_) {
        throw std::runtime_error(
            (boost::format("ECPs in orb file: %1% and options %2% differ") %
             orb.getECPName() % ecp_name_)
                .str());
      }
    }
    if (orb.getNumberOfAlphaElectrons() != num_alpha_electrons_ ||
        orb.getNumberOfBetaElectrons() != num_beta_electrons_) {
      throw std::runtime_error(
          (boost::format("Number of electrons in guess orb file "
                         "and in dftengine differ: "
                         "alpha %1% vs %2%, beta %3% vs %4%.") %
           orb.getNumberOfAlphaElectrons() % num_alpha_electrons_ %
           orb.getNumberOfBetaElectrons() % num_beta_electrons_)
              .str());
    }
    if (orb.getBasisSetSize() != dftbasis_.AOBasisSize()) {
      throw std::runtime_error(
          (boost::format("Number of levels in guess orb file: "
                         "%1% and in dftengine: %2% differ.") %
           orb.getBasisSetSize() % dftbasis_.AOBasisSize())
              .str());
    }
  } else {
    orb.setNumberOfOccupiedLevels(num_alpha_electrons_);
    orb.setNumberOfOccupiedLevelsBeta(num_beta_electrons_);
  }
  return;
}

void DFTEngine::Prepare(Orbitals& orb, Index numofelectrons) {
  QMMolecule& mol = orb.QMAtoms();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Using " << OPENMP::getMaxThreads() << " threads"
      << std::flush;

  if (XTP_HAS_MKL_OVERLOAD()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using MKL overload for Eigen " << std::flush;
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Using native Eigen implementation, no BLAS overload "
        << std::flush;
  }

  XTP_LOG(Log::error, *pLog_) << " Molecule Coordinates [A] " << std::flush;
  for (const QMAtom& atom : mol) {
    const Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    std::string output = (boost::format("  %1$s"
                                        "   %2$+1.4f %3$+1.4f %4$+1.4f") %
                          atom.getElement() % pos[0] % pos[1] % pos[2])
                             .str();

    XTP_LOG(Log::error, *pLog_) << output << std::flush;
  }

  orb.SetupDftBasis(dftbasis_name_);
  dftbasis_ = orb.getDftBasis();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Loaded DFT Basis Set " << dftbasis_name_ << " with "
      << dftbasis_.AOBasisSize() << " functions" << std::flush;

  if (!auxbasis_name_.empty()) {
    BasisSet auxbasisset;
    auxbasisset.Load(auxbasis_name_);
    auxbasis_.Fill(auxbasisset, mol);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Loaded AUX Basis Set " << auxbasis_name_ << " with "
        << auxbasis_.AOBasisSize() << " functions" << std::flush;
  }
  if (!ecp_name_.empty()) {
    ECPBasisSet ecpbasisset;
    ecpbasisset.Load(ecp_name_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Loaded ECP library " << ecp_name_ << std::flush;

    std::vector<std::string> results = ecp_.Fill(ecpbasisset, mol);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled ECP Basis" << std::flush;
    if (results.size() > 0) {
      std::string message = "";
      for (const std::string& element : results) {
        message += " " + element;
      }
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Found no ECPs for elements" << message
          << std::flush;
    }
  }

  numofelectrons_ = 0;
  num_alpha_electrons_ = 0;
  num_beta_electrons_ = 0;
  num_docc_ = 0;
  num_socc_alpha_ = 0;

  Index nuclear_charge = 0;
  for (const QMAtom& atom : mol) {
    nuclear_charge += atom.getNuccharge();
  }

  Index target_charge = orb.getCharge();
  Index multiplicity = orb.getSpin();

  if (multiplicity < 1) {
    throw std::runtime_error("Spin multiplicity must be >= 1.");
  }

  if (numofelectrons >= 0) {
    numofelectrons_ = numofelectrons;
  } else {
    numofelectrons_ = nuclear_charge - target_charge;
  }

  Index spin_excess = multiplicity - 1;

  if (numofelectrons_ < 0) {
    throw std::runtime_error("Computed a negative number of electrons.");
  }

  if (spin_excess > numofelectrons_) {
    throw std::runtime_error(
        "Spin multiplicity incompatible with total number of electrons.");
  }

  if (((numofelectrons_ + spin_excess) % 2) != 0) {
    throw std::runtime_error(
        "Charge and spin multiplicity imply non-integer alpha/beta "
        "occupations.");
  }

  num_alpha_electrons_ = (numofelectrons_ + spin_excess) / 2;
  num_beta_electrons_ = (numofelectrons_ - spin_excess) / 2;

  num_docc_ = std::min(num_alpha_electrons_, num_beta_electrons_);
  num_socc_alpha_ = num_alpha_electrons_ - num_docc_;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Total number of electrons: " << numofelectrons_
      << " (charge=" << target_charge << ", multiplicity=" << multiplicity
      << ", alpha=" << num_alpha_electrons_ << ", beta=" << num_beta_electrons_
      << ", docc=" << num_docc_ << ", socc=" << num_socc_alpha_ << ")"
      << std::flush;

  SetupInvariantMatrices();
  return;
}

Vxc_Potential<Vxc_Grid> DFTEngine::SetupVxc(const QMMolecule& mol) {
  ScaHFX_ = Vxc_Potential<Vxc_Grid>::getExactExchange(xc_functional_name_);
  if (ScaHFX_ > 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using hybrid functional with alpha=" << ScaHFX_
        << std::flush;
  }
  Vxc_Grid grid;
  grid.GridSetup(grid_name_, mol, dftbasis_);
  Vxc_Potential<Vxc_Grid> vxc(grid);
  vxc.setXCfunctional(xc_functional_name_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Setup numerical integration grid " << grid_name_
      << " for vxc functional " << xc_functional_name_ << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << "\t\t "
      << " with " << grid.getGridSize() << " points"
      << " divided into " << grid.getBoxesSize() << " boxes" << std::flush;
  return vxc;
}

double DFTEngine::NuclearRepulsion(const QMMolecule& mol) const {
  double E_nucnuc = 0.0;

  for (Index i = 0; i < mol.size(); i++) {
    const Eigen::Vector3d& r1 = mol[i].getPos();
    double charge1 = double(mol[i].getNuccharge());
    for (Index j = 0; j < i; j++) {
      const Eigen::Vector3d& r2 = mol[j].getPos();
      double charge2 = double(mol[j].getNuccharge());
      E_nucnuc += charge1 * charge2 / (r1 - r2).norm();
    }
  }
  return E_nucnuc;
}

// spherically average the density matrix belonging to two shells
Eigen::MatrixXd DFTEngine::SphericalAverageShells(
    const Eigen::MatrixXd& dmat, const AOBasis& dftbasis) const {
  Eigen::MatrixXd avdmat = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  for (const AOShell& shellrow : dftbasis) {
    Index size_row = shellrow.getNumFunc();
    Index start_row = shellrow.getStartIndex();
    for (const AOShell& shellcol : dftbasis) {
      Index size_col = shellcol.getNumFunc();
      Index start_col = shellcol.getStartIndex();
      Eigen::MatrixXd shelldmat =
          dmat.block(start_row, start_col, size_row, size_col);
      if (shellrow.getL() == shellcol.getL()) {
        double diagavg = shelldmat.diagonal().sum() / double(shelldmat.rows());
        Index offdiagelements =
            shelldmat.rows() * shelldmat.cols() - shelldmat.cols();
        double offdiagavg = (shelldmat.sum() - shelldmat.diagonal().sum()) /
                            double(offdiagelements);
        avdmat.block(start_row, start_col, size_row, size_col).array() =
            offdiagavg;
        avdmat.block(start_row, start_col, size_row, size_col)
            .diagonal()
            .array() = diagavg;
      } else {
        double avg = shelldmat.sum() / double(shelldmat.size());
        avdmat.block(start_row, start_col, size_row, size_col).array() = avg;
      }
    }
  }
  return avdmat;
}

double DFTEngine::ExternalRepulsion(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  if (multipoles.size() == 0) {
    return 0;
  }

  double E_ext = 0;
  eeInteractor interactor;
  for (const QMAtom& atom : mol) {
    StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
    for (const std::unique_ptr<StaticSite>& site : *externalsites_) {
      if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
        XTP_LOG(Log::error, *pLog_) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << std::flush;
        continue;
      }
      E_ext += interactor.CalcStaticEnergy_site(*site, nucleus);
    }
  }
  return E_ext;
}

Eigen::MatrixXd DFTEngine::IntegrateExternalField(const QMMolecule& mol) const {

  AODipole dipole;
  dipole.setCenter(mol.getPos());
  dipole.Fill(dftbasis_);
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(dipole.Dimension(), dipole.Dimension());
  for (Index i = 0; i < 3; i++) {
    result -= dipole.Matrix()[i] * extfield_[i];
  }
  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalMultipoles(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  Mat_p_Energy result(dftbasis_.AOBasisSize(), dftbasis_.AOBasisSize());
  AOMultipole dftAOESP;

  dftAOESP.FillPotential(dftbasis_, multipoles);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Filled DFT external multipole potential matrix"
      << std::flush;
  result.matrix() = dftAOESP.Matrix();
  result.energy() = ExternalRepulsion(mol, multipoles);

  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalDensity(
    const QMMolecule& mol, const Orbitals& extdensity) const {
  BasisSet basis;
  basis.Load(extdensity.getDFTbasisName());
  AOBasis aobasis;
  aobasis.Fill(basis, extdensity.QMAtoms());
  Vxc_Grid grid;
  grid.GridSetup(gridquality_, extdensity.QMAtoms(), aobasis);
  DensityIntegration<Vxc_Grid> numint(grid);
  Eigen::MatrixXd dmat = extdensity.DensityMatrixFull(state_);

  numint.IntegrateDensity(dmat);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated external density" << std::flush;
  Eigen::MatrixXd e_contrib = numint.IntegratePotential(dftbasis_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated potential from electron density"
      << std::flush;
  AOMultipole esp;
  esp.FillPotential(dftbasis_, extdensity.QMAtoms());

  double nuc_energy = 0.0;
  for (const QMAtom& atom : mol) {
    nuc_energy +=
        numint.IntegratePotential(atom.getPos()) * double(atom.getNuccharge());
    for (const QMAtom& extatom : extdensity.QMAtoms()) {
      const double dist = (atom.getPos() - extatom.getPos()).norm();
      nuc_energy +=
          double(atom.getNuccharge()) * double(extatom.getNuccharge()) / dist;
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated potential from nuclei" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Electrostatic: " << nuc_energy << std::flush;
  return Mat_p_Energy(nuc_energy, e_contrib + esp.Matrix());
}

Eigen::MatrixXd DFTEngine::OrthogonalizeGuess(
    const Eigen::MatrixXd& GuessMOs) const {
  Eigen::MatrixXd nonortho =
      GuessMOs.transpose() * dftAOoverlap_.Matrix() * GuessMOs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(nonortho);
  Eigen::MatrixXd result = GuessMOs * es.operatorInverseSqrt();
  return result;
}

/*************************************************************
 * Extended Hueckel Theory
 ************************************************************/
Eigen::VectorXd DFTEngine::BuildEHTOrbitalEnergies(
    const QMMolecule& mol) const {

  ExtendedHuckelParameters params;

  const Index nao = dftbasis_.AOBasisSize();
  Eigen::VectorXd eps = Eigen::VectorXd::Zero(nao);

  for (const AOShell& shell : dftbasis_) {

    int l = static_cast<int>(shell.getL());
    Index start = shell.getStartIndex();
    Index nfunc = shell.getNumFunc();

    const QMAtom& atom = mol[shell.getAtomIndex()];
    const std::string& element = atom.getElement();

    int used_l = l;
    double e = params.GetWithFallback(element, l, &used_l);

    for (Index i = 0; i < nfunc; ++i) {
      eps(start + i) = e;
    }
  }

  return eps;
}

Eigen::MatrixXd DFTEngine::BuildEHTHamiltonian(const QMMolecule& mol) const {

  const Eigen::MatrixXd& S = dftAOoverlap_.Matrix();
  const Index nao = S.rows();
  Eigen::VectorXd eps = BuildEHTOrbitalEnergies(mol);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nao, nao);
  constexpr double K = 1.75;

  for (Index mu = 0; mu < nao; ++mu) {
    H(mu, mu) = eps(mu);
    for (Index nu = 0; nu < mu; ++nu) {
      double hij = K * S(mu, nu) * 0.5 * (eps(mu) + eps(nu));
      H(mu, nu) = hij;
      H(nu, mu) = hij;
    }
  }

  return H;
}

tools::EigenSystem DFTEngine::ExtendedHuckelGuess(const QMMolecule& mol) const {

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Building Extended Huckel guess" << std::flush;

  Eigen::MatrixXd H = BuildEHTHamiltonian(mol);

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Solving EHT generalized eigenproblem" << std::flush;

  return conv_accelerator_.SolveFockmatrix(H);
}

tools::EigenSystem DFTEngine::ExtendedHuckelDFTGuess(
    const Mat_p_Energy& H0, const QMMolecule& mol,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {

  tools::EigenSystem eht = ExtendedHuckelGuess(mol);

  Eigen::MatrixXd Dmat = conv_accelerator_.DensityMatrix(eht);

  Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);

  Eigen::MatrixXd H = H0.matrix() + e_vxc.matrix();

  if (ScaHFX_ > 0) {
    std::array<Eigen::MatrixXd, 2> both =
        CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), Dmat, 1e-12);
    H += both[0];
    H += ScaHFX_ * both[1];
  } else {
    H += CalcERIs(Dmat, 1e-12);
  }

  return conv_accelerator_.SolveFockmatrix(H);
}

}  // namespace xtp
}  // namespace votca