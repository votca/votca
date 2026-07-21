/*
 *            Copyright 2009-2024 The VOTCA Development Team
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

// ===========================================================================
// STATUS UPDATE (later, via a real pyxtp run, not a unit test): a second,
// independent, more serious bug was found AFTER everything below was
// believed validated. Every buffer-mapping test in this file (and every
// finite-difference test built on top of it, throughout the whole
// plain-DFT-forces branch) used H2 with a minimal, all-s-shell basis
// (3-21G). libint2 returns each shell-pair derivative buffer in
// ROW-MAJOR order, but every buffer mapping below used plain
// Eigen::Map<const Eigen::MatrixXd> (Eigen's default COLUMN-MAJOR) --
// identical to row-major for a 1x1 (pure s-shell) block, silently wrong
// for any block with more than one function per shell (p/d/f angular
// momentum). Invisible for H2/3-21G; a large (~2-4%), real rotational-
// covariance violation for a genuine CO/def2-tzvp calculation (both
// atoms have p and d shells) -- found by comparing the same molecule in
// two orientations term-by-term, not by any test in this codebase.
// Fixed by switching to MatrixLibInt (Eigen::RowMajor), matching the
// already-correct, already-validated convention libint2_calls.cc uses
// for the energy-level (deriv_order=0) integrals this file's functions
// are the derivatives of. See the MatrixLibInt definition below for the
// full explanation. NOT yet re-confirmed by a real rerun as of this
// comment -- the fix is reasoned from a clear, specific mechanism (not
// a guess), but treat the "strongly validated" claims elsewhere in this
// file's history below with real skepticism: they were true only for
// all-s-shell systems, which is a much narrower claim than it reads.
// ===========================================================================
//
// STATUS (ORIGINAL, PRE-DATING THE BUG ABOVE): derivative-integral code
// below is now strongly validated -- FOR H2/3-21G SPECIFICALLY, not in
// general; see the update above.
//
// Full history: build-level root cause (libint2 built without derivative
// support -- eventually traced to Homebrew resolving formulae from its
// hosted API rather than the locally-edited tap file, requiring
// HOMEBREW_NO_INSTALL_FROM_API=1 to actually take effect) is fixed.
// libint2::Engine confirmed (via Psi4's own integral-programming docs) to
// provide all centers' derivatives explicitly, not via translational
// invariance -- code reads all 6 buffers directly, matching the ORIGINAL
// (v1) assumption; the diagnostic print confirms this directly at runtime
// (buf.size()=6, all six non-null).
//
// The finite-difference test then showed a discrepancy that turned out to
// be a bug in the TEST, not this file: BuildH2() displaces atoms in
// Angstrom, but XTP represents geometry internally in Bohr, so the
// finite-difference reference was dS/dR per-Angstrom while the analytic
// result here is dS/dR per-Bohr. Every nonzero entry of the two matrices
// differed by the exact same factor (~1.8897 = 1/0.529177, the
// Bohr-to-Angstrom conversion) -- a uniform multiplicative discrepancy
// across all entries is the signature of a unit mismatch, not a
// structural bug (wrong buffer/index/sign would not produce one uniform
// ratio). After converting the finite-difference reference to per-Bohr
// units in test_aoderivatives.cc, the two agree to ~1 part in 40000,
// consistent with expected O(h^2) truncation error at h=1e-4.
//
// Also independently consistent throughout: the translational-invariance
// symmetry check (deriv[0][2] == -deriv[1][2]) passed even before the
// units bug was found, since that check never involves finite differences
// -- it is a real, unit-independent confirmation that the analytic code
// satisfies the exact mathematical constraint it must. NOTE: this check
// is ALSO blind to the row-major/column-major bug above for an all-s-shell
// system -- d[0][2]==-d[1][2] holds regardless of storage order for 1x1
// blocks, so it could not have caught this either.
// ===========================================================================
//
// Build/run history on this file:
//   - Originally drafted with no local libint2 install available, so it
//     could not be compiled or run at all when first written.
//   - First real build attempt failed at link time (duplicate symbols for
//     libint2's internal static tables) -- fixed by removing a duplicate
//     inclusion of <libint2/statics_definition.h>, which must appear in
//     exactly one translation unit per library.
//   - First real run attempt crashed (null pointer, address 0x0) at test
//     entry, single-threaded. Multi-threaded, the same underlying issue
//     showed up as a repeating C++ exception across worker threads rather
//     than a clean crash, which is a separate, real defect in its own
//     right: the parallel loop below never wrapped its body in a
//     try/catch, so an exception escaping an OpenMP worker thread is
//     undefined behavior. That is not yet fixed either -- still worth
//     adding once the buffer-count question is settled, so a real error
//     in production use terminates cleanly instead of hanging.
//
// CURRENT HYPOTHESIS (revised from the original 6-buffer assumption):
//   For a two-center one-electron integral, libint2 most likely returns
//   only 3 derivative buffers (d/dx,dy,dz with respect to shell1's atom),
//   not 6, exploiting the translational-invariance sum rule
//   d(integral)/dR_atom1 + d(integral)/dR_atom2 = 0 (valid because this
//   integral depends only on the relative position of the two atoms) to
//   avoid computing shell2's derivative explicitly. The original
//   assert(buf.size() == 6) did not catch this, most likely because
//   target_ptr_vec's size() reports a fixed nominal capacity rather than
//   the number of meaningfully-populated buffers -- so an out-of-bounds
//   read at buf[3+xyz] returned a null pointer rather than tripping a
//   bounds check. The code below has been revised to derive shell2's
//   derivative as the negative of shell1's, per shell pair, rather than
//   reading buf[3+xyz] at all. A diagnostic print has been left in to
//   directly confirm the real buf.size() and which entries are non-null
//   on the next run -- please check that output before trusting this fix,
//   and remove the diagnostic once confirmed either way.
// ===========================================================================

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/qmmolecule.h"
#include <atomic>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

// include libint last otherwise it overrides eigen
#include "votca/xtp/make_libint_work.h"
#define LIBINT2_CONSTEXPR_STATICS 0
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wcpp"
#endif
#include <libint2.hpp>
// NOTE: deliberately NOT including <libint2/statics_definition.h> here.
// That header *defines* storage for libint2's internal static tables
// (e.g. FmEval_Chebyshev7<double>::cheb_table) and must appear in exactly
// one translation unit per library -- libint2_calls.cc already provides
// it for the whole votca_xtp library. Including it here too caused a
// duplicate-symbol link error (caught when actually building against a
// real libint2 install, which wasn't possible in the environment this
// file was originally drafted in -- see file status note above).
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace votca {
namespace xtp {

// Defined in libint2_calls.cc, not declared in any header; forward
// declared here for use by ComputeThreeCenterIntegrals below, to avoid
// re-deriving the same shell-iteration logic a third time.
std::vector<Eigen::MatrixXd> ComputeAO3cBlock(const libint2::Shell& auxshell,
                                               const AOBasis& dftbasis,
                                               libint2::Engine& engine);

// Result type: one entry per atom, each holding the three Cartesian
// derivative matrices (d/dx, d/dy, d/dz) of the full AO matrix with respect
// to that atom's nuclear coordinate. Matches the AO matrix convention used
// elsewhere in this file (dense Eigen::MatrixXd, size AOBasisSize^2).
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;

// BUG FOUND (real pyxtp run on CO, def2-tzvp -- see conversation history
// for the full diagnostic trail): libint2 returns each shell-pair
// derivative buffer in ROW-MAJOR order (n1 rows x n2 columns, matching
// libint2_calls.cc's own already-validated MatrixLibInt convention used
// throughout the energy-level, deriv_order=0 integral code). Every
// function below originally mapped these buffers with plain
// Eigen::Map<const Eigen::MatrixXd>, which defaults to Eigen's
// COLUMN-MAJOR storage order -- silently transposing each block within
// a shell pair whenever either shell has more than one function (any
// p/d/f angular momentum, not pure s-type). This was invisible in
// EVERY finite-difference test in this branch, all of which used H2
// with a minimal, all-s-shell basis (3-21G) -- row-major and
// column-major are identical for a 1x1 block. It showed up as a
// genuine, large (~2-4%) rotational-covariance violation for a real
// CO/def2-tzvp calculation (both O and C have p and d shells).
using MatrixLibInt =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Computes nuclear-coordinate derivatives of a two-center one-electron
// integral matrix (overlap, kinetic) via libint2's deriv_order=1 engine.
// obtype must be a two-center operator (libint2::Operator::overlap or
// libint2::Operator::kinetic); this function has not been checked against
// operators with a params argument (e.g. emultipole) or against operators
// with more than two centers (e.g. coulomb/three-center RI integrals) --
// those need separate handling for the additional center(s)' derivatives
// and are deliberately out of scope for this first pass.
//
// COMPILE GUARD: many pre-packaged libint2 builds (Homebrew, Ubuntu
// apt, etc.) are NOT built with derivative support at all -- this
// branch's own build hit exactly this issue early on (traced to
// Homebrew resolving formulae from its hosted API rather than a
// locally-edited tap file; see the file-level history comment above).
// libint2 exposes its own build-time maximum derivative order as the
// LIBINT2_MAX_DERIV_ORDER macro (confirmed directly: libint2's own
// Engine::initialize asserts deriv_order_ <= LIBINT2_MAX_DERIV_ORDER).
// Every function below this point (up to the matching #endif) uses
// deriv_order=1 engines -- if LIBINT2_MAX_DERIV_ORDER is 0, calling any
// of them would hit that assertion (a hard crash), or worse, undefined
// behavior if NDEBUG has compiled the assertion out. Guarding at
// compile time means a build against a derivative-less libint2 simply
// never generates code that could do this, rather than relying on
// every call site remembering to check first.
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1
template <libint2::Operator obtype>
std::vector<AOMatrixDerivative> computeOneBodyIntegralDerivatives(
    const AOBasis& aobasis) {

  static_assert(
      libint2::operator_traits<obtype>::nopers == 1,
      "computeOneBodyIntegralDerivatives only supports single-operator "
      "types (overlap, kinetic) in this first pass; multipole-type "
      "operators with multiple result matrices are not yet handled here.");

  // AOBasis has no direct getNumofAtoms(); getFuncPerAtom() is sized by
  // number of atoms and IS a confirmed real accessor (checked directly
  // against aobasis.h), so its size is used here instead.
  Index natoms = static_cast<Index>(aobasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  std::vector<std::vector<Index>> shellpair_list = aobasis.ComputeShellPairs();
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  // shell -> atom index, needed to know which atom's derivative each
  // returned buffer belongs to.
  std::vector<Index> shell2atom;
  shell2atom.reserve(shells.size());
  for (Index s = 0; s < aobasis.getNumofShells(); ++s) {
    // AOShell::getAtomIndex() and AOBasis::getShell() both confirmed
    // directly against the real headers (aoshell.h, aobasis.h).
    shell2atom.push_back(aobasis.getShell(s).getAtomIndex());
  }

  Index nbf = aobasis.AOBasisSize();
  std::vector<AOMatrixDerivative> result(natoms);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      result[a][xyz] = Eigen::MatrixXd::Zero(nbf, nbf);
    }
  }

  std::vector<libint2::Engine> engines(nthreads);
  // deriv_order = 1 instead of 0 is the one substantive change relative to
  // computeOneBodyIntegrals(); everything else follows the same pattern.
  engines[0] = libint2::Engine(obtype, aobasis.getMaxNprim(),
                               static_cast<int>(aobasis.getMaxL()), 1);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::exception_ptr eptr = nullptr;
  std::atomic<bool> any_nonnull_buffer{false};
#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
   try {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();
    Index atom1 = shell2atom[s1];

    for (Index s2 : shellpair_list[s1]) {

      // TARGETED DIAGNOSTIC: confirm whether the crash happens INSIDE
      // engine.compute() itself (libint2's own generated code), before
      // this function's own buffer-checking logic (added in the
      // previous fix) ever gets a chance to run at all -- if so, that
      // fix could never have addressed this, and the print for the
      // CURRENT shell pair would be the last output seen before the
      // crash. std::endl (not just "\n") to force a flush immediately,
      // so this survives even a crash on the very next line.
      std::cerr << "[PRE-COMPUTE DIAGNOSTIC] about to call engine.compute "
                   "for shell pair s1="
                << s1 << " s2=" << s2 << std::endl;
      engine.compute(shells[s1], shells[s2]);
      std::cerr << "[PRE-COMPUTE DIAGNOSTIC] engine.compute RETURNED for "
                   "shell pair s1="
                << s1 << " s2=" << s2 << std::endl;
      // CRITICAL FIX: buf[0]==nullptr alone is not a sufficient guard.
      // A real, hard crash (SIGSEGV, "no mapping at fault address 0x0")
      // was observed on a different CI architecture than the one that
      // showed the "silently all zero" failure mode -- consistent with
      // buf[0] coming back non-null (passing this check) while
      // buf[3+xyz] (dereferenced below, unconditionally, for shell2's/
      // atom2's derivative) is still null. This exact possibility was
      // flagged in this file's own history from early in this branch
      // ("buf[3..5] being unpopulated even though buf.size() reports
      // 6") but the check was never actually added until now. A
      // segfault is a hardware signal, not a C++ exception -- no
      // try/catch anywhere can catch it, so this must be prevented
      // before the dereference, not handled after.
      if (buf[0] == nullptr || buf[3] == nullptr) {
        continue;  // integrals screened out, or this operator's
                   // derivative support is genuinely absent (caught
                   // below by any_nonnull_buffer/the zero-result check
                   // once every shell pair has been skipped this way).
      }
      any_nonnull_buffer.store(true, std::memory_order_relaxed);

      // See HIGHEST-RISK ASSUMPTION note at the top of this file -- the
      // original assert(buf.size()==6) here did not catch a real problem:
      // a null-pointer crash (address 0x0) was observed at runtime,
      // consistent with buf[3..5] being unpopulated even though buf.size()
      // reports 6 (likely a fixed nominal capacity rather than the number
      // of meaningfully-populated buffers). Diagnostic left in place below
      // to confirm this directly on the next run rather than guess again.
      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();
      Index atom2 = shell2atom[s2];

      // CORRECTED (see file header): libint2::Engine provides all centers'
      // derivatives explicitly -- confirmed against Psi4's own
      // integral-programming documentation, which contrasts this directly
      // against an older, legacy one-electron implementation that DID rely
      // on translational invariance to omit one center. Reading buf[3+xyz]
      // directly (not deriving it via negation) is the corrected approach.
      for (Index xyz = 0; xyz < 3; ++xyz) {
        Eigen::Map<const MatrixLibInt> buf_mat1(buf[xyz], n1, n2);
        result[atom1][xyz].block(bf1, bf2, n1, n2) += buf_mat1;
        if (s1 != s2) {
          result[atom1][xyz].block(bf2, bf1, n2, n1) += buf_mat1.transpose();
        }

        Eigen::Map<const MatrixLibInt> buf_mat2(buf[3 + xyz], n1, n2);
        result[atom2][xyz].block(bf1, bf2, n1, n2) += buf_mat2;
        if (s1 != s2) {
          result[atom2][xyz].block(bf2, bf1, n2, n1) += buf_mat2.transpose();
        }
      }
    }
   } catch (...) {
     // An exception escaping an OpenMP worker thread uncaught is
     // undefined behavior -- observed in this exact function's own
     // early history (see file header) to manifest as a hang/repeating
     // exception across worker threads rather than a clean crash, not
     // just theoretically risky. Capture it here, inside the thread,
     // and rethrow once, safely, in the sequential context after the
     // parallel region ends.
#pragma omp critical
     {
       if (!eptr) {
         eptr = std::current_exception();
       }
     }
   }
  }
  if (eptr) {
    std::rethrow_exception(eptr);
  }
  if (!any_nonnull_buffer.load() && aobasis.getNumofShells() > 0) {
    // See the detailed explanation on ComputeNuclearAttractionDerivatives'
    // own version of this check -- distinct from the compile-time
    // LIBINT2_MAX_DERIV_ORDER guard, since individual operators can each
    // have their own build configuration; this specific operator
    // (obtype -- overlap or kinetic, depending on which public entry
    // point called this template) never produced a single non-null
    // buffer across the entire computation.
    throw std::runtime_error(
        "computeOneBodyIntegralDerivatives: engine.results() returned a "
        "null buffer for EVERY shell pair -- this libint2 build does not "
        "actually support this operator's derivative integrals at "
        "runtime, even though it may report LIBINT2_MAX_DERIV_ORDER >= 1 "
        "for other operators. Rebuild libint2 with this operator's "
        "derivative support enabled (--enable-1body=1) to use this "
        "feature.");
  }
  return result;
}

// Public entry points, mirroring AOOverlap::Fill / AOKinetic::Fill naming.
// Deliberately free functions rather than new methods on AOOverlap/AOKinetic
// for this first pass, to avoid touching those classes' existing (working,
// deriv_order=0) code paths at all.
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis) {
  return computeOneBodyIntegralDerivatives<libint2::Operator::overlap>(
      aobasis);
}

std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis) {
  return computeOneBodyIntegralDerivatives<libint2::Operator::kinetic>(
      aobasis);
}

// ===========================================================================
// Two-center Coulomb metric (P|Q) derivatives, needed for the RI-J/RI-K
// gradient assembly: d(V_PQ)/dR, contracted per the formula derived
// separately (dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR).
//
// STATUS: written following the same validated pattern as the one-body
// case above (computeOneBodyIntegralDerivatives), but NOT yet run or
// tested. Two things specifically are hypotheses, not confirmed facts,
// unlike the one-body case where both are now confirmed:
//
// 1. Buffer count/content: this integral is represented internally as a
//    degenerate 4-center integral via two dummy "unit" shells (s-type,
//    no physical center) -- see AOCoulomb::computeCoulombIntegrals for
//    the deriv_order=0 version this mirrors. The hypothesis here is that
//    at deriv_order=1, libint2 still returns derivatives for only the TWO
//    REAL centers (shells[s1]'s atom, shells[s2]'s atom), i.e. 6 buffers,
//    the same count as the one-body case, since the dummy unit shells
//    have no nuclear position to differentiate with respect to. This is
//    a reasonable extrapolation from the one-body case (confirmed: all
//    real centers' derivatives are returned explicitly, not omitted via
//    translational invariance) but has NOT been checked against this
//    specific operator/braket combination. If this is wrong, expect
//    either a wrong buffer count (should show up immediately as an
//    assertion or crash, same as the one-body case's build-configuration
//    issue did) or, more subtly, a right count with wrong center
//    ordering (would show up as a numerical discrepancy in a
//    finite-difference test, same as the units bug did for one-body).
//
// 2. libint2 build configuration: needs ENABLE_ERI2=2 (per Psi4's
//    documented pattern: "In order for gradient and Hessian tests to not
//    segfault, need ENABLE_ERI and ENABLE_ERI3 and ENABLE_ERI2 =2"),
//    which is a DIFFERENT flag from the ENABLE_ONEBODY=2 already
//    confirmed necessary and sufficient for the one-body case. This has
//    NOT been separately verified working -- if this crashes with the
//    same "null build function ptr" assertion the one-body case hit
//    before its build was fixed, check this flag/value first, not the
//    buffer-indexing code below.
// ===========================================================================
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis& aobasis) {
  Index natoms = static_cast<Index>(aobasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  std::vector<Index> shell2atom;
  shell2atom.reserve(aobasis.getNumofShells());
  for (Index s = 0; s < aobasis.getNumofShells(); ++s) {
    shell2atom.push_back(aobasis.getShell(s).getAtomIndex());
  }

  Index nbf = aobasis.AOBasisSize();
  std::vector<AOMatrixDerivative> result(natoms);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      result[a][xyz] = Eigen::MatrixXd::Zero(nbf, nbf);
    }
  }

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(libint2::Operator::coulomb,
                                aobasis.getMaxNprim(),
                                static_cast<int>(aobasis.getMaxL()), 1);
  engines[0].set(libint2::BraKet::xs_xs);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::exception_ptr eptr_coulmetric = nullptr;
  std::atomic<bool> any_nonnull_buffer_coulmetric{false};
#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
   try {
    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();
    Index atom1 = shell2atom[s1];

    // Same note as AOCoulomb::computeCoulombIntegrals: cannot use
    // shellpairs here, this is a two-center integral and overlap
    // screening would give the wrong result.
    for (Index s2 = 0; s2 <= s1; ++s2) {
      engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 1>(
          shells[s1], libint2::Shell::unit(), shells[s2],
          libint2::Shell::unit());

      // See the detailed explanation on the analogous check in
      // computeOneBodyIntegralDerivatives above -- buf[0]==nullptr
      // alone is not sufficient; buf[3+xyz] is dereferenced
      // unconditionally below too.
      if (buf[0] == nullptr || buf[3] == nullptr) {
        continue;
      }
      any_nonnull_buffer_coulmetric.store(true, std::memory_order_relaxed);

      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();
      Index atom2 = shell2atom[s2];

      for (Index xyz = 0; xyz < 3; ++xyz) {
        Eigen::Map<const MatrixLibInt> buf_mat1(buf[xyz], n1, n2);
        result[atom1][xyz].block(bf1, bf2, n1, n2) += buf_mat1;
        if (s1 != s2) {
          result[atom1][xyz].block(bf2, bf1, n2, n1) += buf_mat1.transpose();
        }

        Eigen::Map<const MatrixLibInt> buf_mat2(buf[3 + xyz], n1, n2);
        result[atom2][xyz].block(bf1, bf2, n1, n2) += buf_mat2;
        if (s1 != s2) {
          result[atom2][xyz].block(bf2, bf1, n2, n1) += buf_mat2.transpose();
        }
      }
    }
   } catch (...) {
#pragma omp critical
     {
       if (!eptr_coulmetric) {
         eptr_coulmetric = std::current_exception();
       }
     }
   }
  }
  if (eptr_coulmetric) {
    std::rethrow_exception(eptr_coulmetric);
  }
  if (!any_nonnull_buffer_coulmetric.load() && aobasis.getNumofShells() > 0) {
    // See the detailed explanation on ComputeNuclearAttractionDerivatives'
    // own version of this check.
    throw std::runtime_error(
        "ComputeCoulombMetricDerivatives: engine.results() returned a "
        "null buffer for EVERY shell pair -- this libint2 build does not "
        "actually support this operator's derivative integrals at "
        "runtime, even though it may report LIBINT2_MAX_DERIV_ORDER >= 1 "
        "for other operators. Rebuild libint2 with this operator's "
        "derivative support enabled (--enable-eri2=1) to use this "
        "feature.");
  }
  return result;
}

// ===========================================================================
// Three-center RI integral derivatives, d(mu,nu|P)/dR, the remaining piece
// needed (alongside the two-center metric above) for the RI-J/RI-K
// gradient assembly.
//
// STATUS: written but NOT yet run or tested. This case is genuinely new
// relative to everything validated so far (one-body: 2 real centers, 6
// buffers, confirmed; two-center metric: 2 real centers via dummy unit
// shells, 6 buffers, confirmed): THIS integral has THREE real centers
// (the aux shell's atom, and both dft shells' atoms -- see
// ComputeAO3cBlock above, which this mirrors at deriv_order=1 instead of
// 0, using the same auxshell/unit()/shell_col/shell_row shell ordering).
// The hypothesis, extrapolating from the confirmed "libint2 returns every
// real center's derivative explicitly, never omits one via translational
// invariance" pattern, is that this returns 9 buffers (3 real centers x
// 3 Cartesian), NOT yet checked. If the buffer count is actually
// different, expect either an immediate crash/assertion (same failure
// mode as the wrong libint2 build-configuration flag earlier) or a
// right-count-wrong-ordering bug that would only show up as a numerical
// discrepancy in a finite-difference test (same as the units bug did for
// the one-body case) -- do not trust this without running the
// corresponding test first.
//
// Also NOT yet verified: which libint2 build flag(s) this needs. Given
// Psi4's documented pattern ("need ENABLE_ERI and ENABLE_ERI3 and
// ENABLE_ERI2 =2" for gradients), ENABLE_ERI3=2 is the most likely
// requirement, separate from both ENABLE_ONEBODY=2 (confirmed necessary
// for the one-body case) and whatever ENABLE_ERI2 level the two-center
// metric case above turned out to need.
//
// Memory/design note, not a correctness concern for this validation-scale
// function but worth flagging before this gets used in a real gradient
// assembly on a real molecule: this materializes the FULL derivative
// tensor (all aux functions x all dft functions x all dft functions x
// all atoms x 3) at once, which scales as O(Naux * Ndft^2 * Natoms * 3) --
// fine for a small test system, but a real RI-J gradient assembly should
// almost certainly contract with the density matrix and fitting
// coefficients on the fly, per shell, rather than ever materializing this
// full tensor. Not addressed here; this function exists to validate the
// underlying integral derivative in isolation first, same as the
// one-body and two-center cases above.
// ===========================================================================

// One entry per Cartesian direction (x,y,z); each holds one matrix per
// aux basis function (indexed by absolute aux AO index), each matrix
// being (dftbasis size x dftbasis size).
using ThreeCenterDerivative = std::array<std::vector<Eigen::MatrixXd>, 3>;

std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis& auxbasis, const AOBasis& dftbasis) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  // NOTE: assumes auxbasis and dftbasis are defined on the SAME molecule
  // (same atoms, same count) -- true for every real use case (RI fitting
  // basis and DFT basis both live on the same molecular geometry), not
  // separately checked here.

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> dftshells = dftbasis.GenerateLibintBasis();
  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = dftbasis.getMapToBasisFunctions();
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

  std::vector<Index> dftshell2atom;
  dftshell2atom.reserve(dftbasis.getNumofShells());
  for (Index s = 0; s < dftbasis.getNumofShells(); ++s) {
    dftshell2atom.push_back(dftbasis.getShell(s).getAtomIndex());
  }
  std::vector<Index> auxshell2atom;
  auxshell2atom.reserve(auxbasis.getNumofShells());
  for (Index s = 0; s < auxbasis.getNumofShells(); ++s) {
    auxshell2atom.push_back(auxbasis.getShell(s).getAtomIndex());
  }

  Index n_dft_bf = dftbasis.AOBasisSize();
  Index n_aux_bf = auxbasis.AOBasisSize();

  std::vector<ThreeCenterDerivative> result(natoms);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      result[a][xyz] = std::vector<Eigen::MatrixXd>(
          n_aux_bf, Eigen::MatrixXd::Zero(n_dft_bf, n_dft_bf));
    }
  }

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 1);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::exception_ptr eptr_3c = nullptr;
  std::atomic<bool> any_nonnull_buffer_3c{false};
#pragma omp parallel for schedule(dynamic)
  for (Index aux = 0; aux < auxbasis.getNumofShells(); ++aux) {
   try {
    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    const libint2::Shell& auxshell = auxshells[aux];
    Index aux_start = auxshell2bf[aux];
    Index atom_aux = auxshell2atom[aux];

    for (Index row = 0; row < dftbasis.getNumofShells(); ++row) {
      const libint2::Shell& shell_row = dftshells[row];
      Index row_start = shell2bf[row];
      Index atom_row = dftshell2atom[row];

      // NOTE: unlike ComputeAO3cBlock (deriv_order=0), NOT restricting to
      // col <= row here even though (mu,nu|P) is symmetric in mu,nu --
      // keeping the full loop for this validation-scale function to
      // avoid also having to get a derivative-specific symmetrization
      // step right on the first pass. Revisit once this is confirmed
      // correct; the deriv_order=0 code's triangular-then-mirror
      // approach should still apply to the derivative case in principle
      // (differentiating a symmetric quantity preserves the symmetry),
      // but that itself is an assumption worth checking separately
      // rather than combining with the buffer-count question here.
      for (Index col = 0; col < dftbasis.getNumofShells(); ++col) {
        const libint2::Shell& shell_col = dftshells[col];
        Index col_start = shell2bf[col];
        Index atom_col = dftshell2atom[col];

        engine
            .compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(
                auxshell, libint2::Shell::unit(), shell_col, shell_row);

        // See the detailed explanation on the analogous check in
        // computeOneBodyIntegralDerivatives above -- this operator
        // dereferences buf[3+xyz] AND buf[6+xyz] unconditionally below
        // too (3 real centers, not 2), so both need checking here.
        if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
          continue;
        }
        any_nonnull_buffer_3c.store(true, std::memory_order_relaxed);

        // See HIGHEST-RISK note above: 9 buffers assumed (3 real centers
        // x 3 Cartesian), ordered [aux center xyz][col-shell's atom
        // xyz][row-shell's atom xyz] -- matching the shell ARGUMENT order
        // passed to compute2 above (auxshell, unit, shell_col,
        // shell_row), by direct analogy with the two-center case where
        // buffer order matched argument order. NOT independently
        // confirmed for this 3-real-center case.
        for (Index xyz = 0; xyz < 3; ++xyz) {
          Eigen::TensorMap<
              Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
              result_aux(buf[xyz], auxshell.size(), shell_col.size(),
                         shell_row.size());
          Eigen::TensorMap<
              Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
              result_col(buf[3 + xyz], auxshell.size(), shell_col.size(),
                         shell_row.size());
          Eigen::TensorMap<
              Eigen::Tensor<const double, 3, Eigen::RowMajor> const>
              result_row(buf[6 + xyz], auxshell.size(), shell_col.size(),
                         shell_row.size());

          for (size_t aux_c = 0; aux_c < auxshell.size(); ++aux_c) {
            Index global_aux = aux_start + static_cast<Index>(aux_c);
            for (size_t col_c = 0; col_c < shell_col.size(); ++col_c) {
              for (size_t row_c = 0; row_c < shell_row.size(); ++row_c) {
                double val_aux = result_aux(aux_c, col_c, row_c);
                double val_col = result_col(aux_c, col_c, row_c);
                double val_row = result_row(aux_c, col_c, row_c);
                Index r = row_start + static_cast<Index>(row_c);
                Index c = col_start + static_cast<Index>(col_c);
                result[atom_aux][xyz][global_aux](r, c) += val_aux;
                result[atom_col][xyz][global_aux](r, c) += val_col;
                result[atom_row][xyz][global_aux](r, c) += val_row;
              }
            }
          }
        }
      }
    }
   } catch (...) {
#pragma omp critical
     {
       if (!eptr_3c) {
         eptr_3c = std::current_exception();
       }
     }
   }
  }
  if (eptr_3c) {
    std::rethrow_exception(eptr_3c);
  }
  if (!any_nonnull_buffer_3c.load() && auxbasis.getNumofShells() > 0 &&
      dftbasis.getNumofShells() > 0) {
    // See the detailed explanation on ComputeNuclearAttractionDerivatives'
    // own version of this check.
    throw std::runtime_error(
        "ComputeThreeCenterDerivatives: engine.results() returned a null "
        "buffer for EVERY shell triple -- this libint2 build does not "
        "actually support this operator's derivative integrals at "
        "runtime, even though it may report LIBINT2_MAX_DERIV_ORDER >= 1 "
        "for other operators. Rebuild libint2 with this operator's "
        "derivative support enabled (--enable-eri3=1) to use this "
        "feature.");
  }
  return result;
}
#else   // !(LIBINT2_MAX_DERIV_ORDER >= 1)
// Stubs for a libint2 build without derivative support (see the compile
// guard's own explanatory comment above computeOneBodyIntegralDerivatives).
// These exist only so the rest of the codebase (dftgradient.cc,
// dftengine.cc) still links; DFTEngine::ComputeAndStoreForces(UKS) has
// its own, separate compile-time check and never actually calls these,
// so in normal use this throw is never reached at all -- it is a
// last-resort safety net, not the primary way this limitation is
// communicated to a user.
namespace {
[[noreturn]] void ThrowNoDerivativeSupport(const char* function_name) {
  throw std::runtime_error(
      std::string(function_name) +
      ": this libint2 build was compiled without derivative integral "
      "support (LIBINT2_MAX_DERIV_ORDER < 1) -- analytic nuclear forces "
      "are unavailable. Rebuild libint2 with derivative support enabled "
      "(e.g. --enable-1body=1 and, for RI-J/RI-K forces, --enable-eri2=1 "
      "at libint2 configure time) to use this feature. Many pre-packaged "
      "libint2 builds (Homebrew, Ubuntu apt, etc.) do not enable this by "
      "default.");
}
}  // namespace

std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeOverlapDerivatives");
}
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeKineticDerivatives");
}
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeCoulombMetricDerivatives");
}
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis&, const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeThreeCenterDerivatives");
}
#endif  // LIBINT2_MAX_DERIV_ORDER >= 1

// Energy-level (deriv_order=0) three-center integral (mu,nu|P), kept here
// (rather than in dftgradient.cc) so all direct libint2 API usage stays
// confined to this file and libint2_calls.cc -- dftgradient.cc, which
// needs this for the RI-J gradient assembly, consumes only the returned
// matrices, never touching libint2 types directly. Structurally the same
// shell-loop as ComputeThreeCenterDerivatives above, minus the
// atom/Cartesian bookkeeping that only applies at deriv_order=1, and
// reusing the already-existing ComputeAO3cBlock helper from
// libint2_calls.cc (declared there, not in any header -- forward
// declared here) rather than re-deriving the same shell-iteration logic
// a third time. UNGUARDED (unlike the functions above): deriv_order=0
// integrals are always supported, by any libint2 build.
std::vector<Eigen::MatrixXd> ComputeThreeCenterIntegrals(
    const AOBasis& auxbasis, const AOBasis& dftbasis) {
  std::vector<libint2::Shell> auxshells = auxbasis.GenerateLibintBasis();
  std::vector<Index> auxshell2bf = auxbasis.getMapToBasisFunctions();

  std::vector<Eigen::MatrixXd> tensor(
      auxbasis.AOBasisSize(),
      Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));

  libint2::Engine engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 0);
  engine.set(libint2::BraKet::xs_xx);

  for (Index aux = 0; aux < auxbasis.getNumofShells(); ++aux) {
    std::vector<Eigen::MatrixXd> block =
        ComputeAO3cBlock(auxshells[aux], dftbasis, engine);
    Index aux_start = auxshell2bf[aux];
    for (size_t i = 0; i < block.size(); ++i) {
      tensor[aux_start + static_cast<Index>(i)] = block[i];
    }
  }
  return tensor;
}

// ===========================================================================
// Nuclear attraction (electron-nucleus Coulomb) integral derivative --
// d(V_ne)/dR, the piece discovered MISSING from the total gradient
// assembly by the first genuine end-to-end SCF+forces test
// (test_dftengine_forces.cc). Every earlier gradient test in this branch
// validated individual terms against FIXED density matrices, never the
// complete picture against a real total SCF energy -- which is exactly
// why this gap went undetected until now: kinetic derivatives were
// validated at the very start of this session and then never actually
// used in any gradient assembly, and this nuclear attraction piece was
// never implemented at all.
//
// V_ne_munu = -sum_A Z_A * integral[ chi_mu(r) (1/|r-R_A|) chi_nu(r) dr ]
//
// Genuinely a 3-real-center integral per point charge: the two basis
// function centers, PLUS the point charge's own position -- structurally
// identical to the already-validated 3-center RI derivative case (2
// basis centers + 1 external center = 3 total), just with a different
// operator (libint2::Operator::nuclear instead of coulomb/xs_xx) and
// looping over EACH atom as a single point charge at a time (via
// engine.set_params(make_point_charges(...)) per atom, matching the
// pattern in libint2's own reference example --
// compute_1body_ints_deriv<Operator::nuclear> in
// libint2/include/libint2/lcao/1body.h, used directly to build real HF
// forces there) rather than all atoms simultaneously.
//
// SIGN: this codebase's own existing, already-validated energy-level
// code confirms V_ne (AOMultipole::FillPotential's output) is negative
// (attractive) -- AOMultipole computes a POSITIVE-charge Fill(aobasis)
// per atom then does `aopotential_ -= Fill(aobasis)`, and
// DFTEngine::SetupH0 does `H0 = kinetic + dftAOESP` (a plain addition,
// no extra negation by the caller). Initially reasoned (incorrectly)
// that libint2::Operator::nuclear must therefore need an explicit
// negation to match -- a finite-difference test (exact magnitude match,
// opposite sign in every single element) confirmed the OPPOSITE:
// libint2::Operator::nuclear already returns the attractive (negative)
// convention directly, no negation needed. Fixed by using += throughout
// rather than -=. Worth remembering as a concrete example of why this
// branch insists on actually running things rather than trusting
// plausible-sounding sign reasoning alone.
//
// STATUS: buffer count/ordering (9 buffers per shell-pair-per-point-
// charge, [shell1 atom][shell2 atom][point-charge atom]) CONFIRMED
// correct by the same finite-difference test -- only the sign was
// wrong, now fixed. Re-run pending to confirm the fix.
// ===========================================================================
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis& aobasis, const QMMolecule& mol) {
  Index natoms = mol.size();
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  std::vector<Index> shell2atom;
  shell2atom.reserve(aobasis.getNumofShells());
  for (Index s = 0; s < aobasis.getNumofShells(); ++s) {
    shell2atom.push_back(aobasis.getShell(s).getAtomIndex());
  }

  Index nbf = aobasis.AOBasisSize();
  std::vector<std::vector<AOMatrixDerivative>> result_thread(nthreads);
  for (Index t = 0; t < nthreads; ++t) {
    result_thread[t].resize(natoms);
    for (Index a = 0; a < natoms; ++a) {
      for (Index xyz = 0; xyz < 3; ++xyz) {
        result_thread[t][a][xyz] = Eigen::MatrixXd::Zero(nbf, nbf);
      }
    }
  }

  std::vector<libint2::Engine> engines(nthreads);
  for (Index i = 0; i < nthreads; ++i) {
    engines[i] = libint2::Engine(libint2::Operator::nuclear,
                                 aobasis.getMaxNprim(),
                                 static_cast<int>(aobasis.getMaxL()), 1);
  }

  // Parallelized over the OUTER atom (point-charge) loop, not the inner
  // shell-pair loop as in every other function in this file -- each
  // atom needs its own engine.set_params() call (mutating shared engine
  // state), so a per-thread engine set up fresh per atom is the natural
  // parallelization axis here, not shell pairs within one fixed engine
  // configuration.
  std::exception_ptr eptr_nucattr = nullptr;
  std::atomic<bool> any_nonnull_buffer_nucattr{false};
#pragma omp parallel for schedule(dynamic)
  for (Index A = 0; A < natoms; ++A) {
   try {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];

    std::vector<libint2::Atom> single_atom(1);
    single_atom[0].atomic_number =
        static_cast<int>(mol[A].getNuccharge());
    single_atom[0].x = mol[A].getPos().x();
    single_atom[0].y = mol[A].getPos().y();
    single_atom[0].z = mol[A].getPos().z();
    engine.set_params(libint2::make_point_charges(single_atom));

    const libint2::Engine::target_ptr_vec& buf = engine.results();

    for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
      Index bf1 = shell2bf[s1];
      Index n1 = shells[s1].size();
      Index atom1 = shell2atom[s1];

      for (Index s2 = 0; s2 <= s1; ++s2) {
        engine.compute(shells[s1], shells[s2]);

        // See the detailed explanation on the analogous check in
        // computeOneBodyIntegralDerivatives above -- this operator
        // dereferences buf[3+xyz] AND buf[6+xyz] unconditionally below
        // too (3 real centers: shell1's atom, shell2's atom, the point
        // charge), so both need checking here. This exact gap is what
        // produced BOTH observed failure modes on different CI
        // architectures for this specific function: a hard segfault
        // (buf[3] or buf[6] actually null) on one, and a silently
        // all-zero result (buf[3]/buf[6] non-null but zero-filled,
        // caught instead by the separate norm check further below) on
        // another -- same underlying root cause, different libint2
        // build behavior for an unsupported operator.
        if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
          continue;
        }
        any_nonnull_buffer_nucattr.store(true, std::memory_order_relaxed);

        Index bf2 = shell2bf[s2];
        Index n2 = shells[s2].size();
        Index atom2 = shell2atom[s2];

        for (Index xyz = 0; xyz < 3; ++xyz) {
          // SIGN FIX: confirmed via finite-difference test (exact
          // magnitude match, opposite sign everywhere) that
          // libint2::Operator::nuclear already returns the attractive
          // (negative) V_ne convention directly -- the explicit
          // negation originally here (reasoned from AOMultipole's own
          // "positive charge in, subtract" pattern) was backwards,
          // double-negating an already-correctly-signed quantity. Using
          // += throughout now, not -=.
          Eigen::Map<const MatrixLibInt> buf_mat1(buf[xyz], n1, n2);
          result_thread[thread_id][atom1][xyz].block(bf1, bf2, n1, n2) +=
              buf_mat1;
          if (s1 != s2) {
            result_thread[thread_id][atom1][xyz].block(bf2, bf1, n2, n1) +=
                buf_mat1.transpose();
          }

          Eigen::Map<const MatrixLibInt> buf_mat2(buf[3 + xyz], n1, n2);
          result_thread[thread_id][atom2][xyz].block(bf1, bf2, n1, n2) +=
              buf_mat2;
          if (s1 != s2) {
            result_thread[thread_id][atom2][xyz].block(bf2, bf1, n2, n1) +=
                buf_mat2.transpose();
          }

          Eigen::Map<const MatrixLibInt> buf_mat3(buf[6 + xyz], n1, n2);
          result_thread[thread_id][A][xyz].block(bf1, bf2, n1, n2) +=
              buf_mat3;
          if (s1 != s2) {
            result_thread[thread_id][A][xyz].block(bf2, bf1, n2, n1) +=
                buf_mat3.transpose();
          }
        }
      }
    }
   } catch (...) {
#pragma omp critical
     {
       if (!eptr_nucattr) {
         eptr_nucattr = std::current_exception();
       }
     }
   }
  }
  if (eptr_nucattr) {
    std::rethrow_exception(eptr_nucattr);
  }
  if (!any_nonnull_buffer_nucattr.load() && natoms > 0 &&
      aobasis.getNumofShells() > 0) {
    // Distinct from the LIBINT2_MAX_DERIV_ORDER compile-time guard
    // above: that macro is library-wide (the maximum derivative order
    // ANY operator supports), but individual operators can each have
    // their own, separate build configuration -- a libint2 build can
    // report LIBINT2_MAX_DERIV_ORDER>=1 (enough for e.g. overlap/
    // kinetic) while this specific operator (nuclear attraction with
    // point-charge derivatives) was never actually generated, in which
    // case engine.results() silently returns null buffers for EVERY
    // shell pair rather than throwing -- exactly the failure mode a CI
    // run against such a libint2 hit: this function ran to completion
    // and returned an all-zero matrix, not an exception, which a
    // finite-difference test then correctly reports as a large,
    // confusing discrepancy rather than a clear "not supported"
    // message. Catching it here, at the one place that already knows
    // whether a single buffer was ever populated across the entire
    // computation, turns that into an explicit, actionable error.
    throw std::runtime_error(
        "ComputeNuclearAttractionDerivatives: engine.results() returned "
        "a null buffer for EVERY shell pair and point charge -- this "
        "libint2 build does not actually support nuclear-attraction "
        "derivative integrals at runtime, even though it may report "
        "LIBINT2_MAX_DERIV_ORDER >= 1 for other operators. Rebuild "
        "libint2 with this operator's derivative support enabled "
        "(--enable-1body=1, which nuclear attraction is part of) to "
        "use this feature.");
  }

  // SECOND, COMPLEMENTARY CHECK: the check above only catches every
  // buffer coming back NULL. A real CI run on a different architecture
  // showed that check did NOT fire, yet the assembled result was still
  // exactly, precisely all-zero -- some libint2 builds apparently
  // return valid, non-null buffer POINTERS for this operator, but the
  // underlying computation for it was never actually generated, so the
  // pointed-to memory is zero-filled rather than containing a real
  // result. A genuine nuclear attraction derivative for any real
  // molecule with nonzero nuclear charges and a non-empty basis can
  // never be exactly zero everywhere.
  double total_norm_sq_nucattr = 0.0;
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      for (Index t = 0; t < nthreads; ++t) {
        total_norm_sq_nucattr += result_thread[t][a][xyz].squaredNorm();
      }
    }
  }
  if (total_norm_sq_nucattr < 1.e-20 && natoms > 0 &&
      aobasis.getNumofShells() > 0) {
    throw std::runtime_error(
        "ComputeNuclearAttractionDerivatives: the assembled result is "
        "exactly zero everywhere, which is physically impossible for a "
        "real molecule -- this libint2 build likely returns valid but "
        "zero-filled buffers for this operator (rather than either "
        "computing it correctly or returning null, which the separate, "
        "earlier check in this function already handles). Rebuild "
        "libint2 with this operator's derivative support enabled "
        "(--enable-1body=1) to use this feature.");
  }

  std::vector<AOMatrixDerivative> result(natoms);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      result[a][xyz] = Eigen::MatrixXd::Zero(nbf, nbf);
      for (Index t = 0; t < nthreads; ++t) {
        result[a][xyz] += result_thread[t][a][xyz];
      }
    }
  }
  return result;
}
#else   // !(LIBINT2_MAX_DERIV_ORDER >= 1)
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis&, const QMMolecule&) {
  throw std::runtime_error(
      "ComputeNuclearAttractionDerivatives: this libint2 build was "
      "compiled without derivative integral support "
      "(LIBINT2_MAX_DERIV_ORDER < 1) -- analytic nuclear forces are "
      "unavailable. Rebuild libint2 with derivative support enabled to "
      "use this feature. Many pre-packaged libint2 builds (Homebrew, "
      "Ubuntu apt, etc.) do not enable this by default.");
}
#endif  // LIBINT2_MAX_DERIV_ORDER >= 1

// Small, deliberately libint2-header-free helper: lets other files (e.g.
// dftengine.cc, which does not otherwise include any libint2 header at
// all) check whether this build's linked libint2 supports derivative
// integrals, WITHOUT needing to include <libint2.hpp> themselves just
// for this one macro. Used by DFTEngine::ComputeAndStoreForces(UKS) to
// skip cleanly (log and return) rather than ever calling the throwing
// stubs above at all.
bool HasLibint2DerivativeSupport() {
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1
  return true;
#else
  return false;
#endif
}

}  // namespace xtp
}  // namespace votca
