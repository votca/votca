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
#include <boost/format.hpp>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/IncrementalFockBuilder.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/aogaussiansites.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/citations.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/dftgradient.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/ewald_potential.h"
#include "votca/xtp/gaussiankernel.h"
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

// Defined in libint2_derivative_calls.cc; used by Initialize() to reject
// compute_forces on builds without derivative integrals.
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

  if (options.exists(key_xtpdft + ".polarizability")) {
    compute_polarizability_ =
        options.get(key_xtpdft + ".polarizability").as<bool>();
  }

  if (options.exists(key_xtpdft + ".compute_forces")) {
    compute_forces_ = options.get(key_xtpdft + ".compute_forces").as<bool>();
  }
  if (compute_forces_ && !HasLibint2DerivativeSupport()) {
    // Fail before the SCF rather than after it.
    throw std::runtime_error(
        "compute_forces=true was requested, but libint2 was built without "
        "derivative integrals for one or more of the required operator "
        "categories (one-body, two-center Coulomb metric, three-center RI). "
        "Many packaged libint2 builds (Homebrew, Ubuntu apt, etc.) lack "
        "them; rebuild libint2 with --enable-1body/--enable-eri2/"
        "--enable-eri3 to use this feature.");
  }
  if (compute_forces_ && !ecp_name_.empty()) {
    // H0 contains V_ECP, but the force routines have no d(V_ECP)/dR term;
    // refuse rather than return incomplete forces.
    throw std::runtime_error(
        "compute_forces=true was requested together with an ECP ('" +
        ecp_name_ +
        "'), but the analytic nuclear forces do not include the ECP "
        "term d(V_ECP)/dR. Either drop the ECP or do not request "
        "compute_forces.");
  }

  if (options.exists(key_xtpdft + ".cdft.enabled")) {
    cdft_enabled_ = options.get(key_xtpdft + ".cdft.enabled").as<bool>();
  }
  if (cdft_enabled_) {
    // Only the spec is stored here; the Constraint needs the molecule and
    // is built in Evaluate() by BuildCDFTConstraint.
    std::string indices_str =
        options.get(key_xtpdft + ".cdft.indices").as<std::string>();
    if (indices_str.empty()) {
      throw std::runtime_error(
          "cdft.enabled=true was requested, but cdft.indices is empty; "
          "specify the atoms of the constrained fragment (0-based, e.g. "
          "'1 3 13:17', same syntax as the fragment indices in "
          "diabatization.xml).");
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
    cdft_constraint_spec_.guess_strategy =
        options.get(key_xtpdft + ".cdft.guess_strategy").as<std::string>();
    // CDFT itself needs no derivative integrals, so no
    // HasLibint2DerivativeSupport() check here.
  }

  initial_guess_ = options.get(".initial_guess").as<std::string>();

  if (initial_guess_ == "dimer_guess") {
    dimer_guess_orbA_name_ = options.get(".dimer_guess_orbA").as<std::string>();
    dimer_guess_orbB_name_ = options.get(".dimer_guess_orbB").as<std::string>();
    if (dimer_guess_orbA_name_.empty() || dimer_guess_orbB_name_.empty()) {
      throw std::runtime_error(
          "initial_guess=dimer_guess requires both dimer_guess_orbA and "
          "dimer_guess_orbB to be set to real monomer .orb file paths.");
    }
  }

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
  // Ceiling for adaptive damping, see convergenceacc.h.
  conv_opt_.mixingmax =
      options.get(key_xtpdft + ".convergence.mixing_max").as<double>();
  conv_opt_.levelshift =
      options.get(key_xtpdft + ".convergence.levelshift").as<double>();
  conv_opt_.levelshiftend =
      options.get(key_xtpdft + ".convergence.levelshift_end").as<double>();
  // Independent from adiis_start, see convergenceacc.h and
  // UKSConvergenceAcc::Iterate.
  conv_opt_.mixingend =
      options.get(key_xtpdft + ".convergence.mixing_end").as<double>();
  conv_opt_.maxout =
      options.get(key_xtpdft + ".convergence.DIIS_maxout").as<bool>();
  conv_opt_.histlength =
      options.get(key_xtpdft + ".convergence.DIIS_length").as<Index>();
  conv_opt_.diis_start =
      options.get(key_xtpdft + ".convergence.DIIS_start").as<double>();
  conv_opt_.adiis_start =
      options.get(key_xtpdft + ".convergence.ADIIS_start").as<double>();
  conv_opt_.davidson_max_iter =
      options.get(key_xtpdft + ".convergence.davidson_max_iter").as<Index>();
  conv_opt_.energy_reset = options.ifExistsReturnElseReturnDefault<double>(
      key_xtpdft + ".convergence.energy_reset", 1.0);
  conv_opt_.second_order = options.ifExistsReturnElseReturnDefault<std::string>(
      key_xtpdft + ".convergence.second_order", "auto");
  if (conv_opt_.second_order != "auto" && conv_opt_.second_order != "always" &&
      conv_opt_.second_order != "off") {
    throw std::runtime_error(
        "xtpdft.convergence.second_order must be auto, always or off, not " +
        conv_opt_.second_order);
  }
  overlap_tolerance_ = options.ifExistsReturnElseReturnDefault<double>(
      key_xtpdft + ".overlap_tolerance", 1e-8);
  ri_pair_threshold_ = options.ifExistsReturnElseReturnDefault<double>(
      key_xtpdft + ".ri_pair_threshold", 1e-10);

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

// Defined in libint2_derivative_calls.cc. Unlike DFTGradient::RIJGradient
// etc., these return AO-matrix derivatives d(M_munu)/dR, which are
// contracted with the density below.
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis& aobasis, const QMMolecule& mol);

// Total ground-state gradient (nuclear repulsion, one-electron, overlap
// Pulay, RI-J, XC and, for hybrids, RI-K) from the converged density matrix,
// stored as the force -dE/dR via Orbitals::setForces(). The individual terms
// are tested in test_dftgradient.cc and test_xcgradient.cc.
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
           "against lacks derivative integrals for one or more required "
           "operator categories; rebuild libint2 with "
           "--enable-1body/--enable-eri2/--enable-eri3 to use analytic "
           "forces."
        << std::flush;
    return;
  }

  if (!ecp_name_.empty()) {
    // Normally rejected in Initialize(); repeated here in case the members
    // were set otherwise. Skip instead of throwing after a converged SCF.
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Skipping force calculation: an ECP ('" << ecp_name_
        << "') was used for this SCF, but analytic nuclear forces do "
           "not include the ECP term d(V_ECP)/dR."
        << std::flush;
    return;
  }

  const QMMolecule& mol = orb.QMAtoms();
  Index natoms = mol.size();

  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Starting force calculation ("
                              << natoms << " atoms)" << std::flush;

  // One-electron (kinetic + nuclear attraction) contribution,
  // dE_one/dR_A = Tr[D d(T+V_ne)/dR_A]. See
  // ComputeNuclearAttractionDerivatives for the sign convention.
  XTP_LOG(Log::error, *pLog_) << TimeStamp()
                              << "   Computing one-electron (kinetic + nuclear "
                                 "attraction) derivatives"
                              << std::flush;
  std::vector<AOMatrixDerivative> dT = ComputeKineticDerivatives(dftbasis_);
  std::vector<AOMatrixDerivative> dVne =
      ComputeNuclearAttractionDerivatives(dftbasis_, mol);
  Eigen::MatrixXd eone_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      eone_grad(a, xyz) = Dmat.cwiseProduct(dT[a][xyz] + dVne[a][xyz]).sum();
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   One-electron derivatives done" << std::flush;

  // Overlap (SCF Pulay) term from the orthonormality constraint
  // C^T S(R) C = I. For canonical MOs at convergence the Lagrange
  // multipliers are the orbital energies, giving
  //   dE/dR_A|_overlap = -Tr[W dS/dR_A],  W = 2 C_occ diag(eps_occ) C_occ^T
  // (energy-weighted density, factor 2 for closed shells; cf. the libint2
  // SCF gradient example). Not to be confused with PulayGradient, the
  // basis-function term of the XC integral.
  Index n_occ = num_docc_ + num_socc_alpha_;
  Eigen::MatrixXd C_occ = orb.MOs().eigenvectors().leftCols(n_occ);
  Eigen::VectorXd eps_occ = orb.MOs().eigenvalues().head(n_occ);
  Eigen::MatrixXd W = 2.0 * C_occ * eps_occ.asDiagonal() * C_occ.transpose();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   Computing overlap (Pulay) derivatives"
      << std::flush;
  std::vector<AOMatrixDerivative> dS = ComputeOverlapDerivatives(dftbasis_);
  Eigen::MatrixXd overlap_pulay_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      overlap_pulay_grad(a, xyz) = -W.cwiseProduct(dS[a][xyz]).sum();
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   Overlap derivatives done" << std::flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   Computing RI-J (Coulomb) gradient" << std::flush;
  Eigen::MatrixXd rij_term =
      DFTGradient::RIJGradient(Dmat, auxbasis_, dftbasis_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   RI-J gradient done" << std::flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   Computing XC grid (Pulay + weight) gradient terms"
      << std::flush;
  Eigen::MatrixXd pulay_term = vxcpotential.PulayGradient(Dmat, dftbasis_);
  Eigen::MatrixXd weight_term = vxcpotential.GridWeightGradient(Dmat, mol);
  Eigen::MatrixXd nucrep_term = DFTGradient::NuclearRepulsionDerivative(mol);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   XC grid gradient terms done" << std::flush;

  Eigen::MatrixXd grad = nucrep_term + eone_grad + overlap_pulay_grad +
                         rij_term + pulay_term + weight_term;

  // Exact-exchange (RI-K) gradient, hybrids only; skipped otherwise since
  // RIKGradient costs O(nocc^2 * naux) linear solves. Its energy,
  // E_K = -sum_ij c_ij.d_ij, equals 0.25*Dmat.cwiseProduct(K).sum() at
  // ScaHFX_ = 1 (tested in test_dftgradient.cc), so the contribution is
  // ScaHFX_ * RIKGradient, as exx = 0.25*ScaHFX_*Dmat.cwiseProduct(K).sum().
  if (ScaHFX_ > 0.0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   Computing RI-K (exact exchange) gradient"
        << std::flush;
    grad += ScaHFX_ * DFTGradient::RIKGradient(C_occ, auxbasis_, dftbasis_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   RI-K gradient done" << std::flush;
  }

  // Translational invariance: the total gradient must sum to zero over
  // atoms. Only logged, not thrown.
  Eigen::Vector3d sum = grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " WARNING: computed forces do not sum to zero across atoms "
           "(translational invariance check failed, max component="
        << sum.cwiseAbs().maxCoeff()
        << ") -- treat these forces with "
           "caution."
        << std::flush;
  }

  // All terms above are dE/dR; store the force -dE/dR (the convention of
  // e.g. ASE's Calculator.get_forces()).
  Eigen::MatrixXd force = -grad;
  orb.setForces(force);

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Computed and stored ground-state nuclear forces."
      << std::flush;
  // Printed in atomic units (Hartree/Bohr), as stored.
  XTP_LOG(Log::error, *pLog_) << " Forces [Ha/Bohr]" << std::flush;
  for (Index a = 0; a < natoms; ++a) {
    std::string output =
        (boost::format(" %1$s"
                       "   %2$+1.6f %3$+1.6f %4$+1.6f") %
         mol[a].getElement() % force(a, 0) % force(a, 1) % force(a, 2))
            .str();
    XTP_LOG(Log::error, *pLog_) << output << std::flush;
  }
}

Eigen::MatrixXd DFTEngine::ComputeOverlapPulayGradientUKS(
    const QMMolecule& mol, const tools::EigenSystem& MOs_alpha,
    const tools::EigenSystem& MOs_beta) const {
  // W = W_alpha + W_beta without the closed-shell factor 2.
  //
  // This term accounts for the implicit R-dependence of C through
  // C^T S(R) C = I and is only valid at an SCF stationary point, so it has
  // no fixed-C finite-difference counterpart. test_dftengine_private.cc
  // checks it against the RKS formula for alpha == beta
  // (overlap_pulay_gradient_uks_reduces_to_rks) and the other non-XC terms
  // by fixed-C finite differences
  // (compute_non_xc_gradient_uks_finite_difference).
  Index n_occ_alpha = num_alpha_electrons_;
  Index n_occ_beta = num_beta_electrons_;
  Eigen::MatrixXd C_alpha_occ = MOs_alpha.eigenvectors().leftCols(n_occ_alpha);
  Eigen::MatrixXd C_beta_occ = MOs_beta.eigenvectors().leftCols(n_occ_beta);
  Eigen::VectorXd eps_alpha_occ = MOs_alpha.eigenvalues().head(n_occ_alpha);
  Eigen::VectorXd eps_beta_occ = MOs_beta.eigenvalues().head(n_occ_beta);
  Eigen::MatrixXd W =
      C_alpha_occ * eps_alpha_occ.asDiagonal() * C_alpha_occ.transpose() +
      C_beta_occ * eps_beta_occ.asDiagonal() * C_beta_occ.transpose();

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

  // One-electron and RI-J terms as for RKS, with D_total = D_alpha + D_beta
  // in place of the closed-shell Dmat.
  std::vector<AOMatrixDerivative> dT = ComputeKineticDerivatives(dftbasis_);
  std::vector<AOMatrixDerivative> dVne =
      ComputeNuclearAttractionDerivatives(dftbasis_, mol);
  Eigen::MatrixXd eone_grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      eone_grad(a, xyz) = D_total.cwiseProduct(dT[a][xyz] + dVne[a][xyz]).sum();
    }
  }

  Eigen::MatrixXd overlap_pulay_grad =
      ComputeOverlapPulayGradientUKS(mol, MOs_alpha, MOs_beta);

  Eigen::MatrixXd grad =
      DFTGradient::NuclearRepulsionDerivative(mol) + eone_grad +
      overlap_pulay_grad +
      DFTGradient::RIJGradient(D_total, auxbasis_, dftbasis_);

  // Exact exchange (RI-K), hybrids only. UKS exchange uses
  // ERIs::CalculateEXX_dmat(P) = 0.5*ERIs::CalculateEXX_mos(C) for
  // P = C C^T, while RIKGradient corresponds to CalculateEXX_mos, hence
  // dE_exx/dR = 0.5*ScaHFX_*[RIKGradient(C_alpha_occ)+RIKGradient(C_beta_occ)].
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
           "built against lacks derivative integrals for one or more "
           "required operator categories; rebuild libint2 with "
           "--enable-1body/--enable-eri2/--enable-eri3 to use "
           "analytic forces."
        << std::flush;
    return;
  }

  if (!ecp_name_.empty()) {
    // See ComputeAndStoreForces.
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Skipping UKS force calculation: an ECP ('"
        << ecp_name_
        << "') was used for this SCF, but analytic nuclear forces do "
           "not include the ECP term d(V_ECP)/dR."
        << std::flush;
    return;
  }

  Eigen::MatrixXd grad =
      ComputeNonXCGradientUKS(orb.QMAtoms(), Dspin, MOs_alpha, MOs_beta);

  // XC gradient (LDA and GGA), see PulayGradientUKS/GridWeightGradientUKS.
  grad += vxcpotential.PulayGradientUKS(Dspin.alpha, Dspin.beta, dftbasis_);
  grad += vxcpotential.GridWeightGradientUKS(Dspin.alpha, Dspin.beta,
                                             orb.QMAtoms());

  // Translational invariance check, logged only (as for RKS).
  Eigen::Vector3d sum = grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " WARNING: computed UKS forces do not sum to zero across "
           "atoms (translational invariance check failed, max "
           "component="
        << sum.cwiseAbs().maxCoeff()
        << ") -- treat these forces with "
           "caution."
        << std::flush;
  }

  // Store the force -dE/dR, as for RKS.
  Eigen::MatrixXd force = -grad;
  orb.setForces(force);

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Computed and stored ground-state UKS nuclear forces."
      << std::flush;
  // Printed in atomic units (Hartree/Bohr), as stored.
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
    std::array<Eigen::MatrixXd, 2> result;
    {
      auto t = timings_.Measure("J (RI)");
      result[0] = ERIs_.CalculateERIs_3c(Dmat);
    }
    if (conv_accelerator_.getUseMixing() || MOCoeff.rows() == 0) {
      auto t = timings_.Measure("K (RI, from density matrix)");
      result[1] = ERIs_.CalculateEXX_3c(Eigen::MatrixXd::Zero(0, 0), Dmat);
    } else {
      auto t = timings_.Measure("K (RI, from occupied MOs)");
      Eigen::MatrixXd occblock = MOCoeff.leftCols(num_docc_ + num_socc_alpha_);
      result[1] = ERIs_.CalculateEXX_3c(occblock, Dmat);
    }
    return result;
  } else {
    auto t = timings_.Measure("J+K (4c)");
    return ERIs_.CalculateERIs_EXX_4c(Dmat, error);
  }
}

// Pure Coulomb contribution J[P] from the current AO density matrix. The code
// dispatches to either RI/3c or conventional 4-center integral evaluation.
Eigen::MatrixXd DFTEngine::CalcERIs(const Eigen::MatrixXd& Dmat,
                                    double error) const {
  if (!auxbasis_name_.empty()) {
    auto t = timings_.Measure("J (RI)");
    return ERIs_.CalculateERIs_3c(Dmat);
  } else {
    auto t = timings_.Measure("J (4c)");
    return ERIs_.CalculateERIs_4c(Dmat, error);
  }
}

void DFTEngine::ReportDimensionsAndMemory() const {
  const double gb = 1024.0 * 1024.0 * 1024.0;
  const double n = double(dftbasis_.AOBasisSize());
  const double threads = double(OPENMP::getMaxThreads());
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " DFT dimensions: " << dftbasis_.AOBasisSize()
      << " basis functions, " << OPENMP::getMaxThreads() << " threads"
      << std::flush;
  if (!auxbasis_name_.empty()) {
    const double naux = double(ERIs_.AuxSize());
    // the stored basis-function pairs for every aux function
    const double tensor =
        naux * double(ERIs_.StoredPairs()) * sizeof(double) / gb;
    const double kept = double(ERIs_.StoredPairs()) /
                        double(std::max<Index>(1, ERIs_.AllPairs()));
    // per thread: the OpenMP reduction copy of J or K, and for K the
    // unpacked 3c slice plus the product temporaries
    const double scratch = threads * 4 * n * n * sizeof(double) / gb;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " RI: " << ERIs_.AuxSize() << " aux functions ("
        << ERIs_.Removedfunctions()
        << " removed from the metric); stored 3c tensor "
        << std::setprecision(3) << tensor << " GB (" << 100.0 * kept
        << "% of the pairs at ri_pair_threshold " << ri_pair_threshold_
        << "); J/K scratch up to about " << scratch << " GB" << std::flush;
  }
  double rss = DFTTimings::ResidentMemoryGB(false);
  if (rss >= 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Resident memory after setup: " << std::setprecision(3) << rss
        << " GB" << std::flush;
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
  Eigen::MatrixXd Dmat = [&]() {
    auto t = timings_.Measure("guess: atomic densities");
    return AtomicGuess(mol);
  }();
  Mat_p_Energy e_vxc = [&]() {
    auto t = timings_.Measure("Vxc");
    return vxcpotential.IntegrateVXC(Dmat);
  }();
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
  timings_.Reset();
  if (cdft_enabled_) {
    citations::Cite("Wu:2005");
  }
  // a tensor carried over from a guess would describe another state
  orb.clearPolarizability();
  bool success = EvaluateAndTime(orb);
  timings_.Report(*pLog_, Log::error);
  if (success && compute_polarizability_) {
    ComputePolarizability(orb);
  }
  double peak = DFTTimings::ResidentMemoryGB(true);
  if (peak >= 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Peak resident memory of this process so far: "
        << std::setprecision(3) << peak << " GB" << std::flush;
  }
  return success;
}

void DFTEngine::ComputePolarizability(Orbitals& orb) {
  if (num_alpha_electrons_ != num_beta_electrons_) {
    throw std::runtime_error(
        "xtpdft polarizability: closed-shell ground states only");
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " Static polarizability (coupled-perturbed KS, 3 uniform fields)"
      << std::flush;
  AODipole dipole;
  dipole.setCenter(orb.QMAtoms().getPos());
  dipole.Fill(dftbasis_);
  const std::vector<Eigen::MatrixXd> h(dipole.Matrix().begin(),
                                       dipole.Matrix().end());
  const std::vector<Eigen::MatrixXd> dP = StaticResponse(orb, h);
  Eigen::Matrix3d alpha;
  for (Index a = 0; a < 3; ++a) {
    for (Index c = 0; c < 3; ++c) {
      alpha(a, c) = -dipole.Matrix()[std::size_t(a)]
                         .cwiseProduct(dP[std::size_t(c)])
                         .sum();
    }
  }
  alpha = 0.5 * (alpha + alpha.transpose()).eval();
  orb.setPolarizability(alpha);
  const double A3 = std::pow(tools::conv::bohr2ang, 3);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Polarizability [a.u.]:\n"
      << alpha << "\n isotropic " << alpha.trace() / 3
      << " a.u. = " << alpha.trace() / 3 * A3 << " A^3" << std::flush;
}

bool DFTEngine::EvaluateAndTime(Orbitals& orb) {
  if (cdft_enabled_) {
    // Dispatched before the setup below, which RunCDFT does itself.
    HirshfeldPartition::Constraint constraint =
        BuildCDFTConstraint(orb.QMAtoms(), cdft_constraint_spec_);

    // No forces during the lambda trials of RunCDFT; only the converged
    // lambda's force is needed. Restored whether or not RunCDFT converges.
    bool original_compute_forces = compute_forces_;
    compute_forces_ = false;
    bool converged = RunCDFT(orb, constraint);
    compute_forces_ = original_compute_forces;

    if (converged && original_compute_forces) {
      // One more UKS evaluation at the converged lambda to compute the DFT
      // force. orb holds the converged MOs, so the SCF converges almost
      // immediately; only the setup is repeated, to keep RunCDFT's
      // interface unchanged.
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " CDFT converged -- computing the ordinary DFT force once, "
             "for the final, converged density only"
          << std::flush;
      Prepare(orb);
      Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
      Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
      ConfigOrbfile(orb);
      EvaluateUKS(orb, H0, vxcpotential);
    }

    if (converged && orb.hasForces()) {
      // Add the CDFT constraint term to the DFT force, once, at the
      // converged lambda. The constraint only carries the summed weight
      // matrix, so the per-atom contributions are rebuilt from
      // cdft_constraint_spec_.atom_indices; this repeats the (cheap, no SCF)
      // reference-density, basis and grid setup of BuildCDFTConstraint.
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
      // Total density, as the charge constraint has
      // spin_alpha_coefficient = spin_beta_coefficient = +1.
      Eigen::MatrixXd density_total = Dspin[0] + Dspin[1];

      Eigen::MatrixXd cdft_gradient_correction =
          Eigen::MatrixXd::Zero(static_cast<Index>(orb.QMAtoms().size()), 3);
      for (Index atom_index : cdft_constraint_spec_.atom_indices) {
        cdft_gradient_correction +=
            HirshfeldPartition::ComputeCDFTForceContribution(
                atoms, atom_index, density_total, orb.QMAtoms(), full_dftbasis,
                grid);
      }
      // ComputeCDFTForceContribution returns d(Tr[D W_c])/dR; the gradient
      // correction is +lambda times that, so the force gets -lambda times it.
      orb.setForces(orb.getForces() -
                    constraint.lambda * cdft_gradient_correction);
    }
    return converged;
  }

  // Prepare replaces orb's basis, so record which basis its MOs belong to.
  const std::string previous_basis =
      orb.hasDFTbasisName() ? orb.getDFTbasisName() : "";
  Prepare(orb);
  ReportDimensionsAndMemory();

  const std::string configured_guess = initial_guess_;
  warm_started_ = false;
  if (warm_start_ && initial_guess_ != "orbfile") {
    std::string reason;
    if (UsableAsWarmStart(orb, previous_basis, reason)) {
      initial_guess_ = "orbfile";
      warm_started_ = true;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " Starting from the orbitals of the previous QM/MM iteration"
          << std::flush;
    } else {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Previous orbitals not usable as guess (" << reason
          << "); using " << initial_guess_ << std::flush;
    }
  }

  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  Vxc_Potential<Vxc_Grid> vxcpotential = [&]() {
    auto t = timings_.Measure("setup: XC grid");
    return SetupVxc(orb.QMAtoms());
  }();
  ConfigOrbfile(orb);

  bool success = false;
  if (force_uks_path_ || num_alpha_electrons_ != num_beta_electrons_) {
    if (force_uks_path_ && num_alpha_electrons_ == num_beta_electrons_) {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp()
          << " Forcing closed-shell singlet through UKS development path."
          << std::flush;
    }
    success = EvaluateUKS(orb, H0, vxcpotential);
  } else {
    success = EvaluateClosedShell(orb, H0, vxcpotential);
  }
  initial_guess_ = configured_guess;
  warm_started_ = false;
  return success;
}

bool DFTEngine::UsableAsWarmStart(const Orbitals& orb,
                                  const std::string& previous_basis,
                                  std::string& reason) const {
  if (!orb.hasMOs()) {
    reason = "no MOs";
    return false;
  }
  const Index n = dftbasis_.AOBasisSize();
  if (orb.MOs().eigenvectors().rows() != n ||
      orb.MOs().eigenvectors().cols() != n) {
    reason = "basis size differs";
    return false;
  }
  if (!previous_basis.empty() && previous_basis != orb.getDFTbasisName()) {
    reason = "basis set differs";
    return false;
  }
  if (orb.getNumberOfAlphaElectrons() != num_alpha_electrons_ ||
      orb.getNumberOfBetaElectrons() != num_beta_electrons_) {
    reason = "electron count differs";
    return false;
  }
  return true;
}

bool DFTEngine::RunCDFT(Orbitals& orb,
                        HirshfeldPartition::Constraint& constraint) {
  Prepare(orb);
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  // initial_guess_ is switched to "orbfile" for warm starts; restored on
  // every exit path.
  std::string saved_initial_guess = initial_guess_;

  constraints_ = {constraint};

  // Bisection in lambda: needs only monotonicity of the population in
  // lambda, no dN/dlambda. The bracket is centred on constraint.lambda and
  // doubled in width until the mismatch changes sign or
  // kMaxBracketAttempts is reached.
  double lambda_lo = constraint.lambda - 0.1;
  double lambda_hi = constraint.lambda + 0.1;

  auto EvaluateMismatch = [&](double lambda) -> double {
    constraints_[0].lambda = lambda;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " CDFT: starting inner SCF at lambda=" << lambda
        << std::flush;
    bool scf_converged = EvaluateUKS(orb, H0, vxcpotential);
    if (!scf_converged) {
      throw std::runtime_error(
          "RunCDFT: inner SCF did not converge at lambda=" +
          std::to_string(lambda));
    }
    if (cdft_constraint_spec_.guess_strategy == "warmstart") {
      initial_guess_ = "orbfile";  // warm start every subsequent call
    }
    // "fresh": keep the configured guess for every trial, for cases where
    // consecutive lambdas give very different densities (dftpackage.xml).
    std::array<Eigen::MatrixXd, 2> Dspin =
        orb.DensityMatrixGroundStateSpinResolved();
    double population =
        constraint.spin_alpha_coefficient *
            Dspin[0].cwiseProduct(constraint.weight_matrix).sum() +
        constraint.spin_beta_coefficient *
            Dspin[1].cwiseProduct(constraint.weight_matrix).sum();
    return population - constraint.target_population;
  };

  try {
    double mismatch_lo = EvaluateMismatch(lambda_lo);
    double mismatch_hi = EvaluateMismatch(lambda_hi);

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
          << TimeStamp() << " CDFT outer iteration " << outer_iter + 1 << " of "
          << max_cdft_iterations_ << ": lambda=" << lambda_mid
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
      << TimeStamp()
      << " RunCDFT: outer bisection loop did not converge "
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
    } else if (initial_guess_ == "dimer_guess") {
      // Closed-shell dimer: the block-diagonal guess is restricted as long
      // as both monomers carry identical alpha and beta MOs (restricted,
      // closed-shell monomers). Open-shell monomers need the UKS path.
      Orbitals dimer_guess_orb = BuildDimerGuessFromMonomerFiles(orb.QMAtoms());
      if (dimer_guess_orb.getNumberOfAlphaElectrons() !=
              dimer_guess_orb.getNumberOfBetaElectrons() ||
          !(dimer_guess_orb.MOs().eigenvectors() ==
            dimer_guess_orb.MOs_beta().eigenvectors())) {
        throw std::runtime_error(
            "initial_guess=dimer_guess: this is a restricted (closed-shell) "
            "calculation, but at least one monomer .orb file is open-shell "
            "or unrestricted. Use force_uks_path to run the dimer "
            "unrestricted with this guess.");
      }
      if (dimer_guess_orb.getNumberOfAlphaElectrons() != num_alpha_electrons_) {
        throw std::runtime_error(
            "initial_guess=dimer_guess: the monomers have " +
            std::to_string(2 * dimer_guess_orb.getNumberOfAlphaElectrons()) +
            " electrons in total, but this calculation has " +
            std::to_string(2 * num_alpha_electrons_) +
            ". Check the monomer charges against the dimer charge.");
      }
      MOs = dimer_guess_orb.MOs();
      MOs.eigenvectors() = OrthogonalizeGuess(MOs.eigenvectors());
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

  // Closed-shell energy Hessian for second-order (TRAH) SCF steps:
  // d^2E/dx dx = 4 [F_vv X - X F_oo + C_v^T G[dP] C_o], with X(a, i) = x_ia,
  // dP = 2 (C_v X C_o^T + C_o X^T C_v^T) and G the response of J, K and
  // V_xc, as in StaticResponse (there in canonical orbitals).
  if (conv_opt_.mode == ConvergenceAcc::KSmode::closed) {
    conv_accelerator_.setOrbitalHessian(
        [this, &vxcpotential](const tools::EigenSystem& mos,
                              const Eigen::MatrixXd& F,
                              const Eigen::MatrixXd& X) {
          return RKSOrbitalHessian(vxcpotential, mos, F, X);
        });
  }

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
  conv_accelerator_.StartNewSCF();

  for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << max_iter_ << std::flush;

    Mat_p_Energy e_vxc = [&]() {
      auto t = timings_.Measure("Vxc");
      return vxcpotential.IntegrateVXC(Dmat);
    }();
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

    {
      auto t = timings_.Measure("DIIS/ADIIS + diagonalisation");
      Dmat = conv_accelerator_.Iterate(Dmat, H, MOs, totenergy);
    }
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
        auto t = timings_.Measure("forces");
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
  conv_uks.setOverlap(dftAOoverlap_, overlap_tolerance_);

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
  } else if (initial_guess_ == "dimer_guess") {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Building UKS guess from two monomer .orb files (dimer_guess)"
        << std::flush;
    Orbitals dimer_guess_orb = BuildDimerGuessFromMonomerFiles(orb.QMAtoms());
    MOs_alpha = dimer_guess_orb.MOs();
    MOs_alpha.eigenvectors() = OrthogonalizeGuess(MOs_alpha.eigenvectors());
    MOs_beta = dimer_guess_orb.MOs_beta();
    MOs_beta.eigenvectors() = OrthogonalizeGuess(MOs_beta.eigenvectors());
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

    auto vxc = [&]() {
      auto t = timings_.Measure("Vxc");
      return vxcpotential.IntegrateVXCSpin(Dspin.alpha, Dspin.beta);
    }();
    H_alpha += vxc.vxc_alpha;
    H_beta += vxc.vxc_beta;
    E_xc = vxc.energy;

    double totenergy = H0.energy() + E_one + E_coul + E_xc + E_exx;

    // CDFT constraint potential (no-op if constraints_ is empty): adds
    // lambda_c * spin_alpha/beta_coefficient * W_c to the alpha/beta Fock
    // matrices (charge constraint: +1/+1; a spin constraint would use
    // +1/-1) and uses the Wu-Van Voorhis Lagrangian
    //   E_CDFT = E_KS + sum_c lambda_c (N_c - N_c^target).
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

    // Coupled Fock builder: both spin densities enter the Coulomb/exchange
    // and XC terms together, so CoupledAugmentedHessianStep/
    // BuildCoupledSigmaVector see the alpha-beta coupling (shared Coulomb
    // potential, cross-spin XC kernel). Same H0 + J/K + V_xc sequence as
    // H_alpha/H_beta above, so the current density reproduces them.
    conv_uks.setCoupledFockBuilder(
        [this, &H0, &vxcpotential](const Eigen::MatrixXd& alpha_new,
                                   const Eigen::MatrixXd& beta_new) {
          return BuildUKSFock(H0.matrix(), vxcpotential, alpha_new, beta_new);
        });
    // Analytic Hessian products (spin-polarized f_xc); without second
    // derivatives in libxc UKSConvergenceAcc falls back to finite
    // differences of the coupled Fock builder.
    conv_uks.setOrbitalHessian(
        [this, &vxcpotential](
            const tools::EigenSystem& mos_alpha,
            const tools::EigenSystem& mos_beta, const Eigen::MatrixXd& F_alpha,
            const Eigen::MatrixXd& F_beta, const Eigen::MatrixXd& X) {
          return UKSOrbitalHessian(vxcpotential, mos_alpha, mos_beta, F_alpha,
                                   F_beta, X);
        });

    {
      auto t = timings_.Measure("DIIS/ADIIS + diagonalisation");
      Dspin = conv_uks.Iterate(Dspin, Hspin, MOs_alpha, MOs_beta, totenergy);
    }
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
        auto t = timings_.Measure("forces");
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
  auto h0_timer =
      std::make_unique<DFTTimings::Scope>(timings_, "setup: one-electron H0");

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
  h0_timer.reset();

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

    auto t = timings_.Measure("setup: external multipoles");
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

  if (ext_potential_matrix_.size() > 0) {
    if (ext_potential_matrix_.rows() != H0.rows() ||
        ext_potential_matrix_.cols() != H0.cols()) {
      throw std::runtime_error(
          "DFTEngine: the external potential matrix does not match the DFT "
          "basis");
    }
    H0 += ext_potential_matrix_;
  }

  if (integrate_ext_field_) {

    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Integrating external electric field with F[Hrt]="
        << extfield_.transpose() << std::flush;
    H0 += IntegrateExternalField(mol);
  }

  if (has_ewaldgrid_) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Integrating external Ewald Potential" << std::flush;
    auto t = timings_.Measure("setup: Ewald potential on grid");
    Vxc_Grid ewaldgrid;
    ewaldgrid.GridSetup(grid_name_, mol, dftbasis_);

    // The rebuild is required: external_ewaldgrid_ was built in QMRegion
    // against a local AOBasis, and a Vxc_Grid holds raw `const AOShell*`
    // into that basis, which now dangle. Only its coordinates and potential
    // values may be read (see QMRegion::ewaldgrid_).
    //
    // The grids agree only if both used the same grid name, molecule and
    // basis; checked, since a mismatch would silently give a wrong H0.
    if (ewaldgrid.getBoxesSize() != external_ewaldgrid_.getBoxesSize()) {
      throw std::runtime_error(
          "DFTEngine: the external Ewald potential grid has " +
          std::to_string(external_ewaldgrid_.getBoxesSize()) +
          " boxes but this molecule's own grid has " +
          std::to_string(ewaldgrid.getBoxesSize()) +
          ". The potential was evaluated on a different grid than the one "
          "being integrated over.");
    }
    // make sure the Potential values are copied from the external ewald grid to
    // this one
    for (Index i = 0; i < ewaldgrid.getBoxesSize(); ++i) {
      GridBox& box = ewaldgrid[i];
      const std::vector<double>& source =
          external_ewaldgrid_[i].getPotentialValues();
      if (Index(source.size()) != box.size()) {
        throw std::runtime_error("DFTEngine: box " + std::to_string(i) +
                                 " of the external Ewald "
                                 "potential grid holds " +
                                 std::to_string(source.size()) +
                                 " values for " + std::to_string(box.size()) +
                                 " grid points.");
      }
      std::vector<double>& values = box.getPotentialValues();
      values = source;
    }
    // The AO matrix depends on the basis, the grid and the potential values;
    // reuse it from the previous run if all three are the same.
    std::string ewald_key;
    if (setup_cache_ != nullptr) {
      double sum = 0.0;
      double sum2 = 0.0;
      Index points = 0;
      for (Index i = 0; i < ewaldgrid.getBoxesSize(); ++i) {
        for (double v : ewaldgrid[i].getPotentialValues()) {
          sum += v;
          sum2 += v * v;
          ++points;
        }
      }
      std::ostringstream key;
      key << std::setprecision(17) << RISetupKey() << "|" << grid_name_ << "|"
          << points << "|" << sum << "|" << sum2;
      ewald_key = key.str();
    }
    if (setup_cache_ != nullptr && setup_cache_->ewald_key == ewald_key &&
        setup_cache_->ewald_matrix.rows() == dftbasis_.AOBasisSize()) {
      H0 += setup_cache_->ewald_matrix;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " Reusing the Ewald potential matrix of the previous run"
          << std::flush;
    } else {
      Ewald_Potential<Vxc_Grid> EwaldIntegration(ewaldgrid);
      const Eigen::MatrixXd ewald_matrix =
          EwaldIntegration.IntegrateEwald(dftbasis_.AOBasisSize()).matrix();
      H0 += ewald_matrix;
      if (setup_cache_ != nullptr) {
        setup_cache_->ewald_key = ewald_key;
        setup_cache_->ewald_matrix = ewald_matrix;
      }
    }

    // The grid reaches the electron density only. The nuclei sit in the
    // same potential, and their share arrives as a scalar -- the same
    // pairing IntegrateExternalMultipoles has with ExternalRepulsion.
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Nuclei-external Ewald potential energy "
        << std::setprecision(9) << ewald_nuclear_energy_ << std::flush;
    E0 += ewald_nuclear_energy_;
  }

  return Mat_p_Energy(E0, H0);
}

std::string DFTEngine::RISetupKey() const {
  std::ostringstream key;
  key << std::setprecision(17) << dftbasis_name_ << "|" << auxbasis_name_ << "|"
      << ri_pair_threshold_;
  for (const AOBasis* basis : {&dftbasis_, &auxbasis_}) {
    key << "|";
    for (const AOShell& shell : *basis) {
      key << static_cast<int>(shell.getL()) << "," << shell.getSize() << ","
          << shell.getPos().x() << "," << shell.getPos().y() << ","
          << shell.getPos().z() << ";";
    }
  }
  return key.str();
}

void DFTEngine::setSCFToleranceFloor(double energy, double error) {
  if (energy <= conv_opt_.Econverged && error <= conv_opt_.error_converged) {
    return;
  }
  conv_opt_.Econverged = std::max(conv_opt_.Econverged, energy);
  conv_opt_.error_converged = std::max(conv_opt_.error_converged, error);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " SCF thresholds for this run: Delta E "
      << conv_opt_.Econverged << " Ha, DIIS error " << conv_opt_.error_converged
      << " (set by the caller)" << std::flush;
}

void DFTEngine::ReturnSetupCache() {
  if (setup_cache_ == nullptr || auxbasis_name_.empty() ||
      ERIs_.AuxSize() == 0) {
    return;
  }
  setup_cache_->eris = std::move(ERIs_);
  setup_cache_->eris_key = eris_key_;
  setup_cache_->has_eris = true;
}

// Precompute SCF-invariant matrices: overlap for the generalized eigenvalue
// problem and the RI/4c electron-repulsion backend that later yields J[P] and
// K[P].
void DFTEngine::SetupInvariantMatrices() {
  auto overlap_timer =
      std::make_unique<DFTTimings::Scope>(timings_, "setup: overlap, S^-1/2");
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
  conv_accelerator_.setOverlap(dftAOoverlap_, overlap_tolerance_);
  conv_accelerator_.PrintConfigOptions();
  overlap_timer.reset();

  if (!auxbasis_name_.empty() && setup_cache_ != nullptr &&
      setup_cache_->has_eris && setup_cache_->eris_key == RISetupKey()) {
    // same basis sets and geometry as the previous run (QM/MM iteration)
    ERIs_ = std::move(setup_cache_->eris);
    setup_cache_->has_eris = false;
    eris_key_ = setup_cache_->eris_key;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Reusing the RI integrals of the previous run (same basis sets "
           "and geometry)"
        << std::flush;
  } else if (!auxbasis_name_.empty()) {
    // prepare invariant part of electron repulsion integrals
    auto ri_start = DFTTimings::Clock::now();
    if (setup_cache_ != nullptr) {
      setup_cache_->has_eris = false;
      setup_cache_->eris = ERIs();  // free a stale tensor before building
      eris_key_ = RISetupKey();
    }
    ERIs_.Initialize(dftbasis_, auxbasis_, ri_pair_threshold_);
    double ri_seconds =
        std::chrono::duration<double>(DFTTimings::Clock::now() - ri_start)
            .count();
    timings_.Add("setup: RI metric V^-1/2", ERIs_.MetricSeconds());
    timings_.Add("setup: RI 3c integrals", ri_seconds - ERIs_.MetricSeconds());
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
    auto t = timings_.Measure("setup: 4c Schwarz screening");
    ERIs_.Initialize_4c(dftbasis_);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculated 4c diagonals. " << std::flush;
  }

  return;
}

namespace {
// Spherically averaged atom SCF used for the atomic guess and the
// Hirshfeld reference densities.
//
// Averaged over rotations, an operator on one atom keeps, between two
// shells of equal l, (trace / (2l+1)) times the identity and nothing
// between different l. The averaged problem therefore separates into one
// small "radial" problem per l over the shells of that l, whose levels are
// (2l+1)-fold degenerate. Electrons are filled into these levels by aufbau;
// a partly filled level holds them spread evenly over its 2l+1 components.
// The density is spherical by construction, so there is no orientation of
// a partly filled shell to choose and no symmetry breaking to converge
// through (with integer occupation of the individual orbitals, an open p
// or d shell made the SCF stall or oscillate).
class SphericalAtomSCF {
 public:
  explicit SphericalAtomSCF(const AOBasis& basis) : n_(basis.AOBasisSize()) {
    for (const AOShell& shell : basis) {
      const Index l = static_cast<Index>(shell.getL());
      Index g = 0;
      while (g < Index(l_.size()) && l_[g] != l) {
        ++g;
      }
      if (g == Index(l_.size())) {
        l_.push_back(l);
        starts_.emplace_back();
      }
      starts_[g].push_back(shell.getStartIndex());
    }
  }

  Index Groups() const { return Index(l_.size()); }

  /// Radial matrix of group g: trace of each shell block / (2l+1).
  Eigen::MatrixXd Reduce(const Eigen::MatrixXd& full, Index g) const {
    const std::vector<Index>& st = starts_[g];
    const Index size = 2 * l_[g] + 1;
    Eigen::MatrixXd red(st.size(), st.size());
    for (Index i = 0; i < Index(st.size()); ++i) {
      for (Index j = 0; j < Index(st.size()); ++j) {
        red(i, j) = full.block(st[i], st[j], size, size).trace() / double(size);
      }
    }
    return red;
  }

  /// Full matrix from radial matrices: red_g(i,j) times the identity on
  /// every shell block of group g, zero elsewhere.
  Eigen::MatrixXd Expand(const std::vector<Eigen::MatrixXd>& red) const {
    Eigen::MatrixXd full = Eigen::MatrixXd::Zero(n_, n_);
    for (Index g = 0; g < Groups(); ++g) {
      const std::vector<Index>& st = starts_[g];
      const Index size = 2 * l_[g] + 1;
      for (Index i = 0; i < Index(st.size()); ++i) {
        for (Index j = 0; j < Index(st.size()); ++j) {
          full.block(st[i], st[j], size, size)
              .diagonal()
              .setConstant(red[g](i, j));
        }
      }
    }
    return full;
  }

  /// Radial density matrices: occupation(k) electrons in the k-th lowest
  /// level of group g, spread evenly over its 2l+1 components.
  std::vector<Eigen::MatrixXd> Density(
      const std::vector<Eigen::MatrixXd>& fock,
      const std::vector<Eigen::MatrixXd>& overlap,
      const std::vector<Eigen::VectorXd>& occupation) const {
    std::vector<Eigen::MatrixXd> dens(Groups());
    for (Index g = 0; g < Groups(); ++g) {
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(fock[g],
                                                                   overlap[g]);
      const Index nocc = occupation[g].size();
      const Eigen::MatrixXd c = es.eigenvectors().leftCols(nocc);
      dens[g] = c * (occupation[g] / double(2 * l_[g] + 1)).asDiagonal() *
                c.transpose();
    }
    return dens;
  }

  /// Aufbau occupation of nelectrons over the levels of all groups, a level
  /// holding up to 2l+1 electrons; used if no configuration is available.
  std::vector<Eigen::VectorXd> AufbauOccupation(
      const std::vector<Eigen::MatrixXd>& fock,
      const std::vector<Eigen::MatrixXd>& overlap, double nelectrons) const {
    struct Level {
      double energy;
      Index group;
      Index index;
    };
    std::vector<Level> levels;
    std::vector<Eigen::VectorXd> occ(Groups());
    for (Index g = 0; g < Groups(); ++g) {
      Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(
          fock[g], overlap[g], Eigen::EigenvaluesOnly);
      occ[g] = Eigen::VectorXd::Zero(es.eigenvalues().size());
      for (Index k = 0; k < es.eigenvalues().size(); ++k) {
        levels.push_back({es.eigenvalues()(k), g, k});
      }
    }
    std::stable_sort(
        levels.begin(), levels.end(),
        [](const Level& a, const Level& b) { return a.energy < b.energy; });
    double left = nelectrons;
    for (const Level& lev : levels) {
      const double take = std::clamp(left, 0.0, double(2 * l_[lev.group] + 1));
      occ[lev.group](lev.index) = take;
      left -= take;
    }
    return occ;
  }

  /// Occupation per group and level from subshell occupations (l, k-th
  /// level of that l, electrons of one spin); false if the basis has no
  /// such level.
  bool OccupationFromSubshells(
      const std::vector<std::array<double, 3>>& subshells,
      const std::vector<Eigen::MatrixXd>& overlap,
      std::vector<Eigen::VectorXd>& occ) const {
    occ.assign(Groups(), Eigen::VectorXd());
    for (Index g = 0; g < Groups(); ++g) {
      occ[g] = Eigen::VectorXd::Zero(overlap[g].rows());
    }
    for (const auto& sub : subshells) {
      const Index l = Index(sub[0]);
      const Index k = Index(sub[1]);
      Index g = 0;
      while (g < Groups() && l_[g] != l) {
        ++g;
      }
      if (g == Groups() || k >= occ[g].size()) {
        return false;
      }
      occ[g](k) += sub[2];
    }
    return true;
  }

 private:
  Index n_;
  std::vector<Index> l_;
  std::vector<std::vector<Index>> starts_;
};

// Pulay DIIS on a flattened Fock vector with a flattened error vector.
class SimpleDIIS {
 public:
  explicit SimpleDIIS(Index maxhist) : maxhist_(maxhist) {}

  Eigen::VectorXd Extrapolate(const Eigen::VectorXd& fock,
                              const Eigen::VectorXd& error) {
    focks_.push_back(fock);
    errors_.push_back(error);
    if (Index(focks_.size()) > maxhist_) {
      focks_.erase(focks_.begin());
      errors_.erase(errors_.begin());
    }
    const Index m = Index(focks_.size());
    if (m < 2) {
      return fock;
    }
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(m + 1, m + 1);
    for (Index i = 0; i < m; ++i) {
      for (Index j = 0; j <= i; ++j) {
        b(i, j) = b(j, i) = errors_[i].dot(errors_[j]);
      }
      b(i, m) = b(m, i) = -1.0;
    }
    // Scale the error overlaps to order one: near convergence they are
    // ~1e-12 next to the -1 constraint entries, and a rank-revealing solver
    // would treat the whole block as zero (the extrapolation then stalls).
    // The scale only changes the Lagrange multiplier, not the coefficients.
    const double scale = b.topLeftCorner(m, m).diagonal().maxCoeff();
    if (scale > 0.0) {
      b.topLeftCorner(m, m) /= scale;
    }
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m + 1);
    rhs(m) = -1.0;
    const Eigen::VectorXd x = b.colPivHouseholderQr().solve(rhs);
    if (!x.allFinite()) {
      return fock;
    }
    Eigen::VectorXd result = Eigen::VectorXd::Zero(fock.size());
    for (Index i = 0; i < m; ++i) {
      result += x(i) * focks_[i];
    }
    return result;
  }

 private:
  Index maxhist_;
  std::vector<Eigen::VectorXd> focks_;
  std::vector<Eigen::VectorXd> errors_;
};

// Valence subshells of a neutral atom for the spherical atom SCF, as
// {l, k, alpha electrons, beta electrons}, k counting the subshells of that
// l from the lowest valence one. Subshells are filled by the Madelung (n+l,
// then n) rule; the ncore electrons an ECP replaces are taken from the
// lowest subshells in (n, l) order (def2: 28 = 1s-3d, 46, 60 = up to 4f,
// ...). Every subshell holds alpha and beta electrons equally except that
// alpha - beta = nalpha - nbeta is put into the open subshells, highest
// first. Empty if the ECP core does not end at a subshell boundary or the
// spin difference does not fit.
std::vector<std::array<double, 4>> ValenceConfiguration(Index z, Index ncore,
                                                        Index nalpha,
                                                        Index nbeta) {
  struct Sub {
    Index n;
    Index l;
    double electrons;
  };
  std::vector<Sub> subs;
  Index left = z;
  for (Index sum = 1; left > 0; ++sum) {
    for (Index l = (sum - 1) / 2; l >= 0 && left > 0; --l) {
      const Index take = std::min(left, 2 * (2 * l + 1));
      subs.push_back({sum - l, l, double(take)});
      left -= take;
    }
  }
  std::stable_sort(subs.begin(), subs.end(), [](const Sub& a, const Sub& b) {
    return a.n < b.n || (a.n == b.n && a.l < b.l);
  });
  Index removed = 0;
  std::size_t first = 0;
  while (first < subs.size() && removed < ncore) {
    removed += Index(subs[first].electrons);
    ++first;
  }
  if (removed != ncore) {
    return {};
  }
  std::vector<std::array<double, 4>> result;
  for (std::size_t i = first; i < subs.size(); ++i) {
    Index k = 0;
    for (std::size_t j = first; j < i; ++j) {
      k += (subs[j].l == subs[i].l) ? 1 : 0;
    }
    result.push_back({double(subs[i].l), double(k), 0.5 * subs[i].electrons,
                      0.5 * subs[i].electrons});
  }
  // spin polarisation into the open subshells, last filled first
  double excess = 0.5 * double(nalpha - nbeta);
  for (auto it = result.rbegin(); it != result.rend() && excess > 0.0; ++it) {
    const double capacity = double(2 * Index((*it)[0]) + 1);
    const double move = std::min(excess, capacity - (*it)[2]);
    if (move > 0.0 && (*it)[3] >= move) {
      (*it)[2] += move;
      (*it)[3] -= move;
      excess -= move;
    }
  }
  if (excess > 1e-12) {
    return {};
  }
  return result;
}

// Hund's-rule ground-state (alpha, beta) electron counts for the s/p-block
// elements H-Ca, Ga-Kr and In-Xe. d- and f-block elements are not covered:
// the d^n s^2 vs d^(n+1) s^1 competition is functional-dependent (CP2K's
// ATOM code likewise requires explicit occupations). Returns std::nullopt
// for those, and callers fall back to the parity split.
//
// Aufbau order 1s,2s,2p,3s,3p,4s,3d,4p,5s,4d,5p with Hund's rule in the
// open p shell: p^n gives n alpha for n <= 3, and 3 alpha + (n-3) beta for
// n > 3.
std::optional<std::pair<Index, Index>> HundsRuleAlphaBetaElectrons(
    Index nuclear_charge) {
  switch (nuclear_charge) {
    case 1:
      return std::make_pair(1, 0);  // H:  1s1
    case 2:
      return std::make_pair(1, 1);  // He: 1s2
    case 3:
      return std::make_pair(2, 1);  // Li: [He] 2s1
    case 4:
      return std::make_pair(2, 2);  // Be: 2s2
    case 5:
      return std::make_pair(3, 2);  // B:  2p1
    case 6:
      return std::make_pair(4, 2);  // C:  2p2 (2a)
    case 7:
      return std::make_pair(5, 2);  // N:  2p3 (3a)
    case 8:
      return std::make_pair(5, 3);  // O:  2p4 (3a+1b)
    case 9:
      return std::make_pair(5, 4);  // F:  2p5 (3a+2b)
    case 10:
      return std::make_pair(5, 5);  // Ne: 2p6
    case 11:
      return std::make_pair(6, 5);  // Na: [Ne] 3s1
    case 12:
      return std::make_pair(6, 6);  // Mg: 3s2
    case 13:
      return std::make_pair(7, 6);  // Al: 3p1
    case 14:
      return std::make_pair(8, 6);  // Si: 3p2 (2a)
    case 15:
      return std::make_pair(9, 6);  // P:  3p3 (3a)
    case 16:
      return std::make_pair(9, 7);  // S:  3p4 (3a+1b)
    case 17:
      return std::make_pair(9, 8);  // Cl: 3p5 (3a+2b)
    case 18:
      return std::make_pair(9, 9);  // Ar: 3p6
    case 19:
      return std::make_pair(10, 9);  // K:  [Ar] 4s1
    case 20:
      return std::make_pair(10, 10);  // Ca: 4s2
    // 21-30 (Sc-Zn): 3d block, not covered.
    case 31:
      return std::make_pair(16, 15);  // Ga: [Zn] 4p1
    case 32:
      return std::make_pair(17, 15);  // Ge: 4p2 (2a)
    case 33:
      return std::make_pair(18, 15);  // As: 4p3 (3a)
    case 34:
      return std::make_pair(18, 16);  // Se: 4p4 (3a+1b)
    case 35:
      return std::make_pair(18, 17);  // Br: 4p5 (3a+2b)
    case 36:
      return std::make_pair(18, 18);  // Kr: 4p6
    // 39-48 (Y-Cd): 4d block, not covered.
    case 49:
      return std::make_pair(25, 24);  // In: [Cd] 5p1
    case 50:
      return std::make_pair(26, 24);  // Sn: 5p2 (2a)
    case 51:
      return std::make_pair(27, 24);  // Sb: 5p3 (3a)
    case 52:
      return std::make_pair(27, 25);  // Te: 5p4 (3a+1b)
    case 53:
      return std::make_pair(27, 26);  // I:  5p5 (3a+2b)
    case 54:
      return std::make_pair(27, 27);  // Xe: 5p6
    default:
      return std::nullopt;
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

  // Electrons of the neutral atom that the basis describes: without the
  // core an ECP replaces (ecp.Fill sets it on the atom in `atom`, not on
  // uniqueAtom, which therefore must not be used here).
  const Index z = uniqueAtom.getElementNumber();
  const Index ncore = z - atom[0].getNuccharge();
  const Index numofelectrons = atom[0].getNuccharge();
  Index alpha_e = 0;
  Index beta_e = 0;

  // Total alpha/beta split. The SAD guess (AtomicGuess) uses the
  // parity-based split; the Hirshfeld reference densities of CDFT ask for
  // the Hund's-rule ground state (use_hunds_rule_occupation). Either way
  // the atom is spherical and its open subshells fractionally occupied.
  if (use_hunds_rule_occupation) {
    auto hunds_rule = HundsRuleAlphaBetaElectrons(z);
    if (hunds_rule.has_value()) {
      // the ECP core is closed-shell
      alpha_e = hunds_rule->first - ncore / 2;
      beta_e = hunds_rule->second - ncore / 2;
    } else {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp()
          << " No Hund's-rule ground-state occupation table "
             "entry for nuclear charge "
          << z
          << " (d/f-block elements are not covered); "
             "falling back to the parity-based alpha/beta split."
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

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  if (with_ecp) {
    dftAOECP.FillPotential(dftbasis, ecp);
    H0 += dftAOECP.Matrix();
  }

  if (uniqueAtom.getElement() == "H") {
    // One electron: the lowest orbital of H0, as before.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(
        H0, dftAOoverlap.Matrix());
    const Eigen::VectorXd c = es.eigenvectors().col(0);
    return c * c.transpose();
  }

  // Spherically averaged SCF with fractional occupation of open shells
  // (see SphericalAtomSCF).
  const SphericalAtomSCF sph(dftbasis);
  const Index ngroups = sph.Groups();
  using Radial = std::vector<Eigen::MatrixXd>;
  Radial S_red(ngroups);
  Radial H0_red(ngroups);
  // S^-1/2 of each group, to measure the error in an orthonormal basis
  Radial X_red(ngroups);
  for (Index g = 0; g < ngroups; ++g) {
    S_red[g] = sph.Reduce(dftAOoverlap.Matrix(), g);
    H0_red[g] = sph.Reduce(H0, g);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_red[g]);
    X_red[g] = es.eigenvectors() *
               es.eigenvalues().cwiseInverse().cwiseSqrt().asDiagonal() *
               es.eigenvectors().transpose();
  }

  // Energy and averaged Fock matrices of the spin densities (radial)
  struct State {
    Radial Da, Db, Fa, Fb;
    double energy = 0.0;
  };
  // from the functional, not ScaHFX_: that member is only set in SetupVxc,
  // and this function is also called without it (Hirshfeld references, tests)
  const double scahfx =
      Vxc_Potential<Vxc_Grid>::getExactExchange(xc_functional_name_);
  auto Build = [&](const Radial& Da, const Radial& Db) {
    const Eigen::MatrixXd D_alpha = sph.Expand(Da);
    const Eigen::MatrixXd D_beta = sph.Expand(Db);
    const Eigen::MatrixXd D_total = D_alpha + D_beta;
    Eigen::MatrixXd H_alpha = H0;
    Eigen::MatrixXd H_beta = H0;
    double e_two = 0.0;
    if (scahfx > 0) {
      std::array<Eigen::MatrixXd, 2> both_alpha =
          ERIs_atom.CalculateERIs_EXX_4c(D_alpha, 1e-12);
      std::array<Eigen::MatrixXd, 2> both_beta =
          ERIs_atom.CalculateERIs_EXX_4c(D_beta, 1e-12);
      const Eigen::MatrixXd hartree = both_alpha[0] + both_beta[0];
      H_alpha += hartree + scahfx * both_alpha[1];
      H_beta += hartree + scahfx * both_beta[1];
      e_two = 0.5 * D_total.cwiseProduct(hartree).sum() +
              0.5 * scahfx *
                  (both_alpha[1].cwiseProduct(D_alpha).sum() +
                   both_beta[1].cwiseProduct(D_beta).sum());
    } else {
      const Eigen::MatrixXd hartree =
          ERIs_atom.CalculateERIs_4c(D_total, 1e-12);
      H_alpha += hartree;
      H_beta += hartree;
      e_two = 0.5 * D_total.cwiseProduct(hartree).sum();
    }
    const auto vxc = gridIntegration.IntegrateVXCSpin(D_alpha, D_beta);
    H_alpha += vxc.vxc_alpha;
    H_beta += vxc.vxc_beta;
    State st;
    st.Da = Da;
    st.Db = Db;
    st.energy = D_total.cwiseProduct(H0).sum() + e_two + vxc.energy;
    st.Fa.resize(ngroups);
    st.Fb.resize(ngroups);
    for (Index g = 0; g < ngroups; ++g) {
      st.Fa[g] = sph.Reduce(H_alpha, g);
      st.Fb[g] = sph.Reduce(H_beta, g);
    }
    return st;
  };
  // Fixed occupation of the valence subshells (Madelung configuration); an
  // aufbau occupation that may change between iterations only if there is
  // none for this atom and basis.
  std::vector<Eigen::VectorXd> occ_alpha;
  std::vector<Eigen::VectorXd> occ_beta;
  bool fixed_occupation = false;
  {
    const std::vector<std::array<double, 4>> config =
        ValenceConfiguration(z, ncore, alpha_e, beta_e);
    std::vector<std::array<double, 3>> alpha_subshells;
    std::vector<std::array<double, 3>> beta_subshells;
    for (const auto& sub : config) {
      alpha_subshells.push_back({sub[0], sub[1], sub[2]});
      beta_subshells.push_back({sub[0], sub[1], sub[3]});
    }
    fixed_occupation =
        !config.empty() &&
        sph.OccupationFromSubshells(alpha_subshells, S_red, occ_alpha) &&
        sph.OccupationFromSubshells(beta_subshells, S_red, occ_beta);
    if (!fixed_occupation) {
      XTP_LOG(Log::warning, *pLog_)
          << TimeStamp() << " No subshell configuration for "
          << uniqueAtom.getElement()
          << " in this basis/ECP; atomic SCF uses aufbau occupation"
          << std::flush;
    }
  }
  auto Occupy = [&](const Radial& F, bool alpha) {
    if (fixed_occupation) {
      return sph.Density(F, S_red, alpha ? occ_alpha : occ_beta);
    }
    return sph.Density(
        F, S_red,
        sph.AufbauOccupation(F, S_red, double(alpha ? alpha_e : beta_e)));
  };
  // flattened Fock matrices and commutator errors F D S - S D F
  auto Flatten = [&](const State& st, Eigen::VectorXd& fock,
                     Eigen::VectorXd& error) {
    Index total = 0;
    for (Index g = 0; g < ngroups; ++g) {
      total += 2 * st.Fa[g].size();
    }
    fock.resize(total);
    error.resize(total);
    Index pos = 0;
    for (Index spin = 0; spin < 2; ++spin) {
      const Radial& F = (spin == 0) ? st.Fa : st.Fb;
      const Radial& D = (spin == 0) ? st.Da : st.Db;
      for (Index g = 0; g < ngroups; ++g) {
        const Eigen::MatrixXd fds = F[g] * D[g] * S_red[g];
        const Eigen::MatrixXd err =
            X_red[g] * (fds - fds.transpose()) * X_red[g];
        const Index size = F[g].size();
        fock.segment(pos, size) =
            Eigen::Map<const Eigen::VectorXd>(F[g].data(), size);
        error.segment(pos, size) =
            Eigen::Map<const Eigen::VectorXd>(err.data(), size);
        pos += size;
      }
    }
  };
  auto Unflatten = [&](const Eigen::VectorXd& fock, Radial& Fa, Radial& Fb) {
    Index pos = 0;
    for (Index spin = 0; spin < 2; ++spin) {
      Radial& F = (spin == 0) ? Fa : Fb;
      for (Index g = 0; g < ngroups; ++g) {
        const Index rows = S_red[g].rows();
        F[g] = Eigen::Map<const Eigen::MatrixXd>(fock.data() + pos, rows, rows);
        pos += rows * rows;
      }
    }
  };

  // Pulay DIIS on the averaged Fock matrices. With the occupation of each
  // subshell fixed, the density is a smooth function of the Fock matrix and
  // DIIS converges in about 5-15 iterations for H-Kr (including the 3d
  // metals, where aufbau occupation flips between 4s and 3d).
  State cur = Build(Occupy(H0_red, true), Occupy(H0_red, false));
  SimpleDIIS diis(8);
  const Index maxiter = 100;
  // The molecule's tolerances, but not tighter than the atom needs as a
  // starting density: 1e-7 is at the noise floor of the atomic integrals and
  // grid for the 3d metals, where DIIS then stalls.
  const double error_tolerance = std::max(conv_opt_.error_converged, 1e-6);
  const double energy_tolerance = std::max(conv_opt_.Econverged, 1e-8);
  double energy_old = cur.energy;
  bool converged = false;
  Index this_iter = 0;
  for (; this_iter < maxiter; this_iter++) {
    Eigen::VectorXd fock;
    Eigen::VectorXd error;
    Flatten(cur, fock, error);
    const double max_error = error.cwiseAbs().maxCoeff();
    XTP_LOG(Log::debug, *pLog_)
        << TimeStamp() << " Iter " << this_iter << " of " << maxiter << " Etot "
        << std::setprecision(12) << cur.energy << " error " << max_error
        << std::flush;
    if (this_iter > 0 && max_error < error_tolerance &&
        std::abs(cur.energy - energy_old) < energy_tolerance) {
      converged = true;
      break;
    }
    energy_old = cur.energy;
    Radial Fa(ngroups);
    Radial Fb(ngroups);
    Unflatten(diis.Extrapolate(fock, error), Fa, Fb);
    cur = Build(Occupy(Fa, true), Occupy(Fb, false));
  }
  if (converged) {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Converged after " << this_iter + 1
        << " iterations, Etot=" << std::setprecision(12) << cur.energy
        << std::flush;
  } else {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Not converged after " << maxiter
        << " iterations. Unconverged density." << std::flush;
  }
  const Radial& Da = cur.Da;
  const Radial& Db = cur.Db;

  const Eigen::MatrixXd density = sph.Expand(Da) + sph.Expand(Db);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Atomic density Matrix for " << uniqueAtom.getElement()
      << " gives N=" << std::setprecision(9)
      << density.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons."
      << std::flush;
  return density;
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
      << " for unique elements (Hirshfeld reference densities)" << std::flush;

  std::map<std::string, Eigen::MatrixXd> reference_densities;
  for (const std::string& element : elements) {
    QMAtom unique_atom(0, element, Eigen::Vector3d::Zero());
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Calculating Hirshfeld reference density for "
        << element << std::flush;
    // Hund's-rule ground state, unlike the SAD guess in AtomicGuess.
    reference_densities[element] = RunAtomicDFT_unrestricted(
        unique_atom, /*use_hunds_rule_occupation=*/true);
  }
  return reference_densities;
}

// Coupled UKS Fock matrices of the spin densities (alpha_new, beta_new):
// H0 + J + a_x K + V_xc, plus the CDFT constraint potential, exactly as in
// the EvaluateUKS loop. Used for the finite-difference Hessian products.
UKSConvergenceAcc::SpinFock DFTEngine::BuildUKSFock(
    const Eigen::MatrixXd& H0, const Vxc_Potential<Vxc_Grid>& vxcpotential,
    const Eigen::MatrixXd& alpha_new, const Eigen::MatrixXd& beta_new) {
  UKSConvergenceAcc::SpinFock H_new;
  H_new.alpha = H0;
  H_new.beta = H0;
  constexpr double kIntegralError = 1e-8;
  if (ScaHFX_ > 0) {
    std::array<Eigen::MatrixXd, 2> both_alpha_new =
        CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), alpha_new, kIntegralError);
    std::array<Eigen::MatrixXd, 2> both_beta_new =
        CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), beta_new, kIntegralError);
    Eigen::MatrixXd J_new = both_alpha_new[0] + both_beta_new[0];
    H_new.alpha += J_new + ScaHFX_ * both_alpha_new[1];
    H_new.beta += J_new + ScaHFX_ * both_beta_new[1];
  } else {
    Eigen::MatrixXd J_new = CalcERIs(alpha_new + beta_new, kIntegralError);
    H_new.alpha += J_new;
    H_new.beta += J_new;
  }
  auto vxc_new = [&]() {
    auto t = timings_.Measure("Vxc");
    return vxcpotential.IntegrateVXCSpin(alpha_new, beta_new);
  }();
  H_new.alpha += vxc_new.vxc_alpha;
  H_new.beta += vxc_new.vxc_beta;
  for (const HirshfeldPartition::Constraint& c : constraints_) {
    H_new.alpha += (c.lambda * c.spin_alpha_coefficient) * c.weight_matrix;
    H_new.beta += (c.lambda * c.spin_beta_coefficient) * c.weight_matrix;
  }
  return H_new;
}

// Closed-shell energy Hessian d^2E/dx dx (layout of trah.h) at the orbitals
// mos with Fock matrix F, applied to the columns of X. With X(a, i) = x_ia:
//   (H x) = 4 [F_vv X - X F_oo + C_v^T G[dP] C_o],
//   dP = 2 (C_v X C_o^T + C_o X^T C_v^T),
// G the response of J, 0.5 a_x K and V_xc (CoupledResponse).
Eigen::MatrixXd DFTEngine::RKSOrbitalHessian(const Vxc_Potential<Vxc_Grid>& vxc,
                                             const tools::EigenSystem& mos,
                                             const Eigen::MatrixXd& F,
                                             const Eigen::MatrixXd& X) {
  const Index nocc = num_docc_;
  const Eigen::MatrixXd& C = mos.eigenvectors();
  const Index nvirt = C.cols() - nocc;
  const Eigen::MatrixXd C_occ = C.leftCols(nocc);
  const Eigen::MatrixXd C_virt = C.rightCols(nvirt);
  const Eigen::MatrixXd F_oo = C_occ.transpose() * F * C_occ;
  const Eigen::MatrixXd F_vv = C_virt.transpose() * F * C_virt;
  const Eigen::MatrixXd P0 = 2.0 * C_occ * C_occ.transpose();
  using RowMajor =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::vector<Eigen::MatrixXd> U;
  for (Index k = 0; k < X.cols(); ++k) {
    const Eigen::VectorXd x = X.col(k);
    U.push_back(Eigen::Map<const RowMajor>(x.data(), nocc, nvirt).transpose());
  }
  const std::vector<Eigen::MatrixXd> G =
      CoupledResponse(vxc, P0, C_occ, C_virt, U);
  Eigen::MatrixXd HX(X.rows(), X.cols());
  for (Index k = 0; k < X.cols(); ++k) {
    const RowMajor R = 4.0 * (F_vv * U[std::size_t(k)] -
                              U[std::size_t(k)] * F_oo + G[std::size_t(k)])
                                 .transpose();
    HX.col(k) = Eigen::Map<const Eigen::VectorXd>(R.data(), R.size());
  }
  return HX;
}

// UKS energy Hessian d^2E/dx dx (layout of trah.h, alpha block first) at
// the orbitals mos_alpha/mos_beta with Fock matrices F_alpha/F_beta,
// applied to the columns of X. Per spin s, with X_s(a, i) = x_ia:
//   (H x)_s = 2 [F_vv X_s - X_s F_oo + C_v^T G_s C_o],
//   dP_s = C_v X_s C_o^T + C_o X_s^T C_v^T,
//   G_s = J[dP_a + dP_b] + a_x K[dP_s] + dV_xc,s[dP_a, dP_b].
// The constraint potential of CDFT is in F; it has no response.
Eigen::MatrixXd DFTEngine::UKSOrbitalHessian(
    const Vxc_Potential<Vxc_Grid>& vxc, const tools::EigenSystem& mos_alpha,
    const tools::EigenSystem& mos_beta, const Eigen::MatrixXd& F_alpha,
    const Eigen::MatrixXd& F_beta, const Eigen::MatrixXd& X) {
  using RowMajor =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  struct Channel {
    Index nocc;
    Index nvirt;
    Eigen::MatrixXd C_occ;
    Eigen::MatrixXd C_virt;
    Eigen::MatrixXd F_oo;
    Eigen::MatrixXd F_vv;
    Eigen::MatrixXd P;
  };
  auto setup = [](const tools::EigenSystem& mos, const Eigen::MatrixXd& F,
                  Index nocc) {
    Channel c;
    const Eigen::MatrixXd& C = mos.eigenvectors();
    c.nocc = nocc;
    c.nvirt = C.cols() - nocc;
    c.C_occ = C.leftCols(nocc);
    c.C_virt = C.rightCols(c.nvirt);
    c.F_oo = c.C_occ.transpose() * F * c.C_occ;
    c.F_vv = c.C_virt.transpose() * F * c.C_virt;
    c.P = c.C_occ * c.C_occ.transpose();
    return c;
  };
  const std::array<Channel, 2> ch = {
      setup(mos_alpha, F_alpha, num_alpha_electrons_),
      setup(mos_beta, F_beta, num_beta_electrons_)};
  const std::array<Index, 2> offset = {0, ch[0].nocc * ch[0].nvirt};

  const Index ncols = X.cols();
  std::vector<std::array<Eigen::MatrixXd, 2>> U(
      static_cast<std::size_t>(ncols));
  std::vector<std::array<Eigen::MatrixXd, 2>> dP(
      static_cast<std::size_t>(ncols));
  for (Index k = 0; k < ncols; ++k) {
    for (std::size_t s = 0; s < 2; ++s) {
      const Channel& c = ch[s];
      const Eigen::VectorXd x = X.col(k).segment(offset[s], c.nocc * c.nvirt);
      U[std::size_t(k)][s] =
          Eigen::Map<const RowMajor>(x.data(), c.nocc, c.nvirt).transpose();
      const Eigen::MatrixXd half =
          c.C_virt * U[std::size_t(k)][s] * c.C_occ.transpose();
      dP[std::size_t(k)][s] = half + half.transpose();
    }
  }
  std::vector<std::array<Eigen::MatrixXd, 2>> G;
  {
    auto t = timings_.Measure("second order: f_xc (spin)");
    G = vxc.IntegrateFXCSpin(ch[0].P, ch[1].P, dP);
  }
  Eigen::MatrixXd HX(X.rows(), ncols);
  for (Index k = 0; k < ncols; ++k) {
    std::array<Eigen::MatrixXd, 2>& g = G[std::size_t(k)];
    const std::array<Eigen::MatrixXd, 2>& d = dP[std::size_t(k)];
    if (!auxbasis_name_.empty()) {
      {
        auto t = timings_.Measure("second order: J (RI)");
        const Eigen::MatrixXd J = ERIs_.CalculateERIs_3c(d[0] + d[1]);
        g[0] += J;
        g[1] += J;
      }
      if (ScaHFX_ > 0) {
        auto t = timings_.Measure("second order: K (RI)");
        for (std::size_t s = 0; s < 2; ++s) {
          // dP_s = 1/2 [(A + C_o)(A + C_o)^T - (A - C_o)(A - C_o)^T]
          const Eigen::MatrixXd A = ch[s].C_virt * U[std::size_t(k)][s];
          g[s] += 0.5 * ScaHFX_ *
                  ERIs_.CalculateEXX_factors(A + ch[s].C_occ, A - ch[s].C_occ);
        }
      }
    } else {
      auto t = timings_.Measure("second order: J+K (4c)");
      const std::array<Eigen::MatrixXd, 2> jk_a =
          ERIs_.CalculateERIs_EXX_4c(d[0], 1e-12);
      const std::array<Eigen::MatrixXd, 2> jk_b =
          ERIs_.CalculateERIs_EXX_4c(d[1], 1e-12);
      g[0] += jk_a[0] + jk_b[0];
      g[1] += jk_a[0] + jk_b[0];
      if (ScaHFX_ > 0) {
        g[0] += ScaHFX_ * jk_a[1];
        g[1] += ScaHFX_ * jk_b[1];
      }
    }
    for (std::size_t s = 0; s < 2; ++s) {
      const Channel& c = ch[s];
      const Eigen::MatrixXd& u = U[std::size_t(k)][s];
      const RowMajor R = 2.0 * (c.F_vv * u - u * c.F_oo +
                                c.C_virt.transpose() * g[s] * c.C_occ)
                                   .transpose();
      HX.col(k).segment(offset[s], c.nocc * c.nvirt) =
          Eigen::Map<const Eigen::VectorXd>(R.data(), R.size());
    }
  }
  return HX;
}

std::vector<Eigen::MatrixXd> DFTEngine::CoupledResponse(
    const Vxc_Potential<Vxc_Grid>& vxc, const Eigen::MatrixXd& P0,
    const Eigen::MatrixXd& C_occ, const Eigen::MatrixXd& C_virt,
    const std::vector<Eigen::MatrixXd>& U) {
  auto density = [&](const Eigen::MatrixXd& u) {
    const Eigen::MatrixXd half = (C_virt * u) * C_occ.transpose();
    return Eigen::MatrixXd(2.0 * (half + half.transpose()));
  };
  std::vector<Eigen::MatrixXd> dP;
  for (const Eigen::MatrixXd& u : U) {
    dP.push_back(density(u));
  }
  std::vector<Eigen::MatrixXd> G;
  {
    auto t = timings_.Measure("response: f_xc");
    G = vxc.IntegrateFXC(P0, dP);
  }
  std::vector<Eigen::MatrixXd> result;
  for (std::size_t k = 0; k < U.size(); ++k) {
    if (!auxbasis_name_.empty()) {
      {
        auto t = timings_.Measure("response: J (RI)");
        G[k] += ERIs_.CalculateERIs_3c(dP[k]);
      }
      if (ScaHFX_ > 0) {
        auto t = timings_.Measure("response: K (RI)");
        const Eigen::MatrixXd A = C_virt * U[k];
        G[k] +=
            0.5 * ScaHFX_ * ERIs_.CalculateEXX_factors(A + C_occ, A - C_occ);
      }
    } else {
      auto t = timings_.Measure("response: J+K (4c)");
      const std::array<Eigen::MatrixXd, 2> jk =
          ERIs_.CalculateERIs_EXX_4c(dP[k], 1e-12);
      G[k] += jk[0];
      if (ScaHFX_ > 0) {
        G[k] += 0.5 * ScaHFX_ * jk[1];
      }
    }
    result.push_back(C_virt.transpose() * G[k] * C_occ);
  }
  return result;
}

std::vector<Eigen::MatrixXd> DFTEngine::StaticResponse(
    const Orbitals& orb, const std::vector<Eigen::MatrixXd>& perturbations,
    double tolerance, Index max_iterations) {
  timings_.Reset();
  Orbitals work = orb;
  Prepare(work);
  if (num_alpha_electrons_ != num_beta_electrons_) {
    throw std::runtime_error(
        "DFTEngine::StaticResponse: closed-shell ground states only");
  }
  const Index n_ao = dftbasis_.AOBasisSize();
  if (!orb.hasMOs() || orb.MOs().eigenvectors().rows() != n_ao) {
    throw std::runtime_error(
        "DFTEngine::StaticResponse: the orbitals do not belong to this basis");
  }
  const Vxc_Potential<Vxc_Grid> vxc = SetupVxc(work.QMAtoms());

  const Index n_occ = num_docc_;
  const Index n_mo = orb.MOs().eigenvectors().cols();
  const Index n_virt = n_mo - n_occ;
  const Eigen::MatrixXd C_occ = orb.MOs().eigenvectors().leftCols(n_occ);
  const Eigen::MatrixXd C_virt = orb.MOs().eigenvectors().rightCols(n_virt);
  const Eigen::VectorXd& e = orb.MOs().eigenvalues();
  // e_a - e_i
  Eigen::MatrixXd diag(n_virt, n_occ);
  for (Index i = 0; i < n_occ; ++i) {
    diag.col(i) = e.tail(n_virt).array() - e(i);
  }
  if (diag.minCoeff() <= 0.0) {
    throw std::runtime_error(
        "DFTEngine::StaticResponse: no gap between occupied and virtual "
        "orbitals");
  }
  const Eigen::MatrixXd P0 = 2.0 * C_occ * C_occ.transpose();

  const Index n_rhs = Index(perturbations.size());
  std::vector<Eigen::MatrixXd> b(static_cast<std::size_t>(n_rhs));
  for (Index k = 0; k < n_rhs; ++k) {
    const Eigen::MatrixXd& h = perturbations[std::size_t(k)];
    if (h.rows() != n_ao || h.cols() != n_ao) {
      throw std::runtime_error(
          "DFTEngine::StaticResponse: perturbation does not match the basis");
    }
    b[std::size_t(k)] = -(C_virt.transpose() * h * C_occ);
  }

  auto density = [&](const Eigen::MatrixXd& U) {
    const Eigen::MatrixXd half = (C_virt * U) * C_occ.transpose();
    return Eigen::MatrixXd(2.0 * (half + half.transpose()));
  };
  // the orbital Hessian (A + B) applied to the U of several right-hand sides
  auto apply = [&](const std::vector<Eigen::MatrixXd>& U) {
    std::vector<Eigen::MatrixXd> result =
        CoupledResponse(vxc, P0, C_occ, C_virt, U);
    for (std::size_t k = 0; k < U.size(); ++k) {
      result[k] += diag.cwiseProduct(U[k]);
    }
    return result;
  };

  // Galerkin projection on one subspace shared by all right-hand sides,
  // grown by the preconditioned residuals (block preconditioned CG for a
  // symmetric positive definite Hessian); the right-hand sides profit from
  // each other's directions.
  const double kDropNorm = 1e-10;
  std::vector<Eigen::MatrixXd> V;   // orthonormal (Frobenius)
  std::vector<Eigen::MatrixXd> AV;  // Hessian applied to V
  std::vector<Eigen::MatrixXd> x(static_cast<std::size_t>(n_rhs));
  std::vector<Eigen::MatrixXd> r(static_cast<std::size_t>(n_rhs));
  std::vector<double> bnorm(static_cast<std::size_t>(n_rhs), 0.0);
  std::vector<Index> active;
  for (Index k = 0; k < n_rhs; ++k) {
    const std::size_t kk = std::size_t(k);
    bnorm[kk] = b[kk].norm();
    x[kk] = Eigen::MatrixXd::Zero(n_virt, n_occ);
    r[kk] = b[kk];
    if (bnorm[kk] > 0.0) {
      active.push_back(k);
    }
  }
  Index iteration = 0;
  while (!active.empty()) {
    if (iteration == max_iterations) {
      throw std::runtime_error(
          "DFTEngine::StaticResponse: CPKS did not converge in " +
          std::to_string(max_iterations) + " iterations");
    }
    ++iteration;
    // new directions: preconditioned residuals, orthonormalized against the
    // subspace and each other (twice, for stability)
    std::vector<Eigen::MatrixXd> fresh;
    for (Index k : active) {
      Eigen::MatrixXd v = r[std::size_t(k)].cwiseQuotient(diag);
      const double norm0 = v.norm();
      for (Index pass = 0; pass < 2; ++pass) {
        for (const Eigen::MatrixXd& u : V) {
          v -= u.cwiseProduct(v).sum() * u;
        }
        for (const Eigen::MatrixXd& u : fresh) {
          v -= u.cwiseProduct(v).sum() * u;
        }
      }
      const double norm = v.norm();
      if (norm > kDropNorm * norm0) {
        fresh.push_back(v / norm);
      }
    }
    if (fresh.empty()) {
      throw std::runtime_error(
          "DFTEngine::StaticResponse: CPKS subspace stopped growing");
    }
    const std::vector<Eigen::MatrixXd> Afresh = apply(fresh);
    V.insert(V.end(), fresh.begin(), fresh.end());
    AV.insert(AV.end(), Afresh.begin(), Afresh.end());
    const Index m = Index(V.size());
    Eigen::MatrixXd M(m, m);
    for (Index i = 0; i < m; ++i) {
      for (Index j = 0; j <= i; ++j) {
        M(i, j) =
            0.5 * (V[std::size_t(i)].cwiseProduct(AV[std::size_t(j)]).sum() +
                   V[std::size_t(j)].cwiseProduct(AV[std::size_t(i)]).sum());
        M(j, i) = M(i, j);
      }
    }
    const Eigen::LDLT<Eigen::MatrixXd> solver(M);
    std::vector<Index> still;
    double worst = 0.0;
    for (Index k : active) {
      const std::size_t kk = std::size_t(k);
      Eigen::VectorXd proj(m);
      for (Index i = 0; i < m; ++i) {
        proj(i) = V[std::size_t(i)].cwiseProduct(b[kk]).sum();
      }
      const Eigen::VectorXd y = solver.solve(proj);
      x[kk].setZero();
      r[kk] = b[kk];
      for (Index i = 0; i < m; ++i) {
        x[kk] += y(i) * V[std::size_t(i)];
        r[kk] -= y(i) * AV[std::size_t(i)];
      }
      const double rel = r[kk].norm() / bnorm[kk];
      worst = std::max(worst, rel);
      if (rel >= tolerance) {
        still.push_back(k);
      }
    }
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " CPKS iteration " << iteration << ": "
        << fresh.size() << " new directions, subspace " << m
        << ", largest relative residual " << worst << std::flush;
    active = still;
  }

  std::vector<Eigen::MatrixXd> dP;
  for (const Eigen::MatrixXd& u : x) {
    dP.push_back(density(u));
  }
  timings_.Report(*pLog_, Log::error);
  return dP;
}

std::vector<HirshfeldPartition::AtomicReference> DFTEngine::HirshfeldReferences(
    const QMMolecule& mol) const {
  return HirshfeldPartition::BuildAtomicReferences(
      mol, dftbasis_name_, ComputeHirshfeldReferenceDensities(mol));
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
    // Hirshfeld weights are additive, w_fragment(r) = sum_{i in fragment}
    // w_i(r); the neutral reference population is the sum of the nuclear
    // charges.
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
      << neutral_reference_population
      << ", requested relative charge=" << spec.target_charge
      << ", absolute target population=" << constraint.target_population
      << std::flush;

  return constraint;
}

void DFTEngine::ConfigOrbfile(Orbitals& orb) {
  // A warm start was checked in UsableAsWarmStart, against the basis the MOs
  // were computed in.
  if (initial_guess_ == "orbfile" && !warm_started_) {

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

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Using " << EigenBackendDescription() << std::flush;

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
      << "\t\t " << " with " << grid.getGridSize() << " points"
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
// Average of an atom's density matrix over all rotations about the nucleus.
// In real spherical harmonics a rotation acts on each shell by an orthogonal
// matrix that depends only on l, so the average of the block between two
// shells is (trace / (2l+1)) times the identity if both have the same l, and
// zero otherwise. The result is independent of the orientation of the input
// (an open-shell atom's SCF converges to a symmetry-broken solution whose
// orientation is decided by round-off) and keeps the number of electrons:
// the overlap between shells of one atom is diagonal within equal l and zero
// between different l.
Eigen::MatrixXd DFTEngine::SphericalAverageShells(
    const Eigen::MatrixXd& dmat, const AOBasis& dftbasis) const {
  Eigen::MatrixXd avdmat = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  for (const AOShell& shellrow : dftbasis) {
    for (const AOShell& shellcol : dftbasis) {
      if (shellrow.getL() != shellcol.getL()) {
        continue;
      }
      const Index size = shellrow.getNumFunc();
      const double diagavg = dmat.block(shellrow.getStartIndex(),
                                        shellcol.getStartIndex(), size, size)
                                 .trace() /
                             double(size);
      avdmat
          .block(shellrow.getStartIndex(), shellcol.getStartIndex(), size, size)
          .diagonal()
          .setConstant(diagavg);
    }
  }
  return avdmat;
}

namespace {
// The charge-response sites among the external sites (PolarSite with
// isChargeResponse): their induced charge and dipole are Gaussians.
const PolarSite* AsChargeResponseSite(const StaticSite& site) {
  const PolarSite* p = dynamic_cast<const PolarSite*>(&site);
  return (p != nullptr && p->isChargeResponse()) ? p : nullptr;
}
}  // namespace

double DFTEngine::ExternalRepulsion(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite>>& multipoles) const {

  if (multipoles.size() == 0) {
    return 0;
  }

  double E_ext = 0;
  eeInteractor interactor;
  for (const QMAtom& atom : mol) {
    StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
    for (const std::unique_ptr<StaticSite>& site : multipoles) {
      if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
        XTP_LOG(Log::error, *pLog_) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << std::flush;
        continue;
      }
      // The nuclei must see the permanent multipoles plus the induced
      // dipole, as the electrons do (AOMultipole::FillBlock uses the virtual
      // getDipole(), which includes the induced dipole of a PolarSite).
      // eeInteractor::CalcStaticEnergy_site reads the non-virtual Q(), i.e.
      // permanent moments only, so the induced dipole is added here. Since
      // the electronic and nuclear terms nearly cancel for a neutral QM
      // region, omitting it gives large errors. VSiteA itself is correct:
      // the induction solver treats permanent and induced moments
      // separately, and changing it would double-count induced dipoles.
      Vector9d Q = site->Q();
      if (const PolarSite* p = AsChargeResponseSite(*site)) {
        // induced charge and dipole of a charge-response site: Gaussians,
        // as the electrons see them (IntegrateExternalMultipoles)
        const Eigen::Vector3d r = nucleus.getPos() - p->getPos();
        const GaussianKernel k(r.norm(), p->CRWidth());
        E_ext += nucleus.getCharge() * (p->Induced_Charge() * k.B0 +
                                        k.B1 * p->getInducedDipole().dot(r));
      } else {
        Q.segment<3>(1) += site->getInducedDipole();  // zero for a StaticSite
      }
      Index rank = site->getRank();
      if (rank < 1 && Q.segment<3>(1).norm() > 1e-12) {
        rank = 1;  // same promotion, same threshold, as AOMultipole
      }
      StaticSite effective(site->getId(), site->getElement(), site->getPos());
      effective.setMultipole(Q, rank);

      E_ext += interactor.CalcStaticEnergy_site(effective, nucleus);
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
    const std::vector<std::unique_ptr<StaticSite>>& multipoles) const {

  Mat_p_Energy result(dftbasis_.AOBasisSize(), dftbasis_.AOBasisSize());
  result.energy() = ExternalRepulsion(mol, multipoles);

  // Charge-response sites: permanent multipoles as points (below, with all
  // sites), induced charge and dipole as Gaussians (AOGaussianSites).
  std::vector<AOGaussianSites::Source> gaussian;
  for (const std::unique_ptr<StaticSite>& site : multipoles) {
    if (const PolarSite* p = AsChargeResponseSite(*site)) {
      gaussian.push_back({p->getPos(), p->CRWidth(), p->Induced_Charge(),
                          p->getInducedDipole()});
    }
  }
  if (!gaussian.empty()) {
    // the point sites' induced dipoles, without the charge-response sites
    std::vector<std::unique_ptr<StaticSite>> point_sites;
    for (const std::unique_ptr<StaticSite>& site : multipoles) {
      if (AsChargeResponseSite(*site) != nullptr) {
        continue;
      }
      const PolarSite* p = dynamic_cast<const PolarSite*>(site.get());
      if (p != nullptr) {
        point_sites.push_back(std::make_unique<PolarSite>(*p));
      } else {
        point_sites.push_back(std::make_unique<StaticSite>(*site));
      }
    }
    AOMultipole permanent;
    permanent.FillPotential(dftbasis_, multipoles,
                            AOMultipole::Moments::Permanent);
    AOMultipole induced;
    induced.FillPotential(dftbasis_, point_sites,
                          AOMultipole::Moments::Induced);
    result.matrix() = permanent.Matrix() + induced.Matrix() -
                      AOGaussianSites::Potential(dftbasis_, gaussian);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Filled DFT external multipole potential matrix ("
        << gaussian.size() << " charge-response sites as Gaussians)"
        << std::flush;
    return result;
  }

  if (setup_cache_ == nullptr) {
    AOMultipole dftAOESP;
    dftAOESP.FillPotential(dftbasis_, multipoles);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Filled DFT external multipole potential matrix"
        << std::flush;
    result.matrix() = dftAOESP.Matrix();
    return result;
  }

  // In QM/MM only the induced dipoles change between iterations. The
  // potential is linear in the moments, so the permanent part is kept in the
  // cache and reused as long as basis, positions and permanent moments are
  // exactly the same; the induced dipoles are integrated every run.
  Eigen::MatrixXd sites(Index(multipoles.size()), 13);
  for (Index i = 0; i < Index(multipoles.size()); ++i) {
    const StaticSite& site = *multipoles[i];
    sites.block<1, 3>(i, 0) = site.getPos().transpose();
    sites(i, 3) = double(site.getRank());
    sites.block<1, 9>(i, 4) = site.Q().transpose();
  }
  const std::string key = RISetupKey();
  if (setup_cache_->multipole_key == key &&
      setup_cache_->multipole_matrix.rows() == dftbasis_.AOBasisSize() &&
      setup_cache_->multipole_sites.rows() == sites.rows() &&
      setup_cache_->multipole_sites == sites) {
    result.matrix() = setup_cache_->multipole_matrix;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Reusing the permanent multipole potential matrix of the previous "
           "run"
        << std::flush;
  } else {
    AOMultipole permanent;
    permanent.FillPotential(dftbasis_, multipoles,
                            AOMultipole::Moments::Permanent);
    result.matrix() = permanent.Matrix();
    setup_cache_->multipole_key = key;
    setup_cache_->multipole_sites = sites;
    setup_cache_->multipole_matrix = permanent.Matrix();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Filled DFT permanent multipole potential matrix"
        << std::flush;
  }

  const bool has_induced =
      std::any_of(multipoles.begin(), multipoles.end(),
                  [](const std::unique_ptr<StaticSite>& site) {
                    return site->getInducedDipole().norm() > 1e-12;
                  });
  if (has_induced) {
    AOMultipole induced;
    induced.FillPotential(dftbasis_, multipoles, AOMultipole::Moments::Induced);
    result.matrix() += induced.Matrix();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Filled DFT induced dipole potential matrix"
        << std::flush;
  }

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

Orbitals DFTEngine::BuildDimerGuessFromMonomerFiles(
    const QMMolecule& dimer_mol) const {
  Orbitals monomerA;
  monomerA.ReadFromCpt(dimer_guess_orbA_name_);
  Orbitals monomerB;
  monomerB.ReadFromCpt(dimer_guess_orbB_name_);

  const QMMolecule& atomsA = monomerA.QMAtoms();
  const QMMolecule& atomsB = monomerB.QMAtoms();
  Index nA = atomsA.size();
  Index nB = atomsB.size();

  // Check 1: atom count and element sequence (before the geometry check,
  // since an ordering error gives a clearer message here).
  if (nA + nB != dimer_mol.size()) {
    throw std::runtime_error(
        "BuildDimerGuessFromMonomerFiles: monomer A (" + std::to_string(nA) +
        " atoms) + monomer B (" + std::to_string(nB) +
        " atoms) does not equal this calculation's molecule (" +
        std::to_string(dimer_mol.size()) +
        " atoms): wrong monomer file(s), or the molecule is not the "
        "concatenation of the two monomers.");
  }
  for (Index i = 0; i < nA; ++i) {
    if (atomsA[i].getElement() != dimer_mol[i].getElement()) {
      throw std::runtime_error(
          "BuildDimerGuessFromMonomerFiles: monomer A atom " +
          std::to_string(i) + " (" + atomsA[i].getElement() +
          ") does not match this calculation's atom " + std::to_string(i) +
          " (" + dimer_mol[i].getElement() +
          "); dimer_guess requires monomer A to be the first N_A atoms of "
          "the molecule, in the same order.");
    }
  }
  for (Index i = 0; i < nB; ++i) {
    if (atomsB[i].getElement() != dimer_mol[nA + i].getElement()) {
      throw std::runtime_error(
          "BuildDimerGuessFromMonomerFiles: monomer B atom " +
          std::to_string(i) + " (" + atomsB[i].getElement() +
          ") does not match this calculation's atom " + std::to_string(nA + i) +
          " (" + dimer_mol[nA + i].getElement() +
          "); dimer_guess requires monomer B to be the remaining atoms of "
          "the molecule (after monomer A), in the same order.");
    }
  }

  // Check 2: intra-monomer interatomic distances, which are invariant under
  // the rigid translation/rotation placing a monomer into the dimer.
  constexpr double kGeometryToleranceBohr = 1e-3;
  auto CheckInternalGeometry = [&](const QMMolecule& monomer_atoms,
                                   Index offset_in_dimer,
                                   const std::string& label) {
    Index n = monomer_atoms.size();
    for (Index i = 0; i < n; ++i) {
      for (Index j = i + 1; j < n; ++j) {
        double monomer_distance =
            (monomer_atoms[i].getPos() - monomer_atoms[j].getPos()).norm();
        double dimer_distance = (dimer_mol[offset_in_dimer + i].getPos() -
                                 dimer_mol[offset_in_dimer + j].getPos())
                                    .norm();
        double diff = std::abs(monomer_distance - dimer_distance);
        if (diff > kGeometryToleranceBohr) {
          throw std::runtime_error(
              "BuildDimerGuessFromMonomerFiles: " + label +
              "'s internal geometry does not match this calculation's "
              "molecule: distance between its atoms " +
              std::to_string(i) + " and " + std::to_string(j) + " is " +
              std::to_string(monomer_distance) +
              " Bohr in the monomer file, but " +
              std::to_string(dimer_distance) +
              " Bohr in this calculation's molecule (difference " +
              std::to_string(diff) + " Bohr, tolerance " +
              std::to_string(kGeometryToleranceBohr) + " Bohr).");
        }
      }
    }
  };
  CheckInternalGeometry(atomsA, 0, "Monomer A");
  CheckInternalGeometry(atomsB, nA, "Monomer B");

  // The MO coefficients are copied without rotating them, so the guess is
  // only exact if each monomer is translated, not rotated, into the dimer.
  auto MaxDeviationFromTranslation = [&](const QMMolecule& monomer_atoms,
                                         Index offset_in_dimer) {
    Eigen::Vector3d shift =
        dimer_mol[offset_in_dimer].getPos() - monomer_atoms[0].getPos();
    double max_dev = 0.0;
    for (Index i = 0; i < monomer_atoms.size(); ++i) {
      double dev = (dimer_mol[offset_in_dimer + i].getPos() -
                    monomer_atoms[i].getPos() - shift)
                       .norm();
      max_dev = std::max(max_dev, dev);
    }
    return max_dev;
  };
  auto WarnIfRotated = [&](const QMMolecule& monomer_atoms,
                           Index offset_in_dimer, const std::string& label) {
    double dev = MaxDeviationFromTranslation(monomer_atoms, offset_in_dimer);
    if (dev > kGeometryToleranceBohr) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " WARNING: " << label
          << " is rotated with respect to its .orb file (max deviation " << dev
          << " bohr after translation). Its MO coefficients are not "
             "rotated, so the dimer guess will be poor."
          << std::flush;
    }
  };
  WarnIfRotated(atomsA, 0, "Monomer A");
  WarnIfRotated(atomsB, nA, "Monomer B");

  Orbitals dimer_guess;
  // PrepareDimerGuessMixedSpin calls SetupDftBasis, which needs QMAtoms()
  // to be set (as in iqm.cc for PrepareDimerGuess).
  dimer_guess.QMAtoms() = dimer_mol;
  dimer_guess.PrepareDimerGuessMixedSpin(monomerA, monomerB);
  return dimer_guess;
}

}  // namespace xtp
}  // namespace votca