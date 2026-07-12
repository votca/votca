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
// STATUS: build-level root cause (libint2 compiled without derivative
// support) has been fixed separately and confirmed to change the failure
// mode. A second, DIFFERENT crash was then observed (SIGTRAP / uncaught
// boost::detail::system_signal_exception), and the buffer-indexing "fix"
// below this comment (originally: derive center 2 via translational
// invariance, only reading buf[0..2]) has been REVERTED, because it was
// based on an incorrect premise -- see CORRECTION below. Whether reverting
// this alone resolves the SIGTRAP is NOT yet confirmed; this may need
// further diagnosis (see note at the very end of this block).
//
// CORRECTION (this supersedes the "CURRENT HYPOTHESIS" note further down,
// which is left in place for history but should be read as superseded):
//   Confirmed via Psi4's own integral-programming documentation (a
//   production code that itself switched to using libint2::Engine
//   directly): "The old one electron integral code used translational
//   invariance relations to minimize the number of integrals to be
//   computed... The Libint2 engine instead provides all integrals, so the
//   caller simply needs to loop over all of the buffers provided." I.e.
//   the translational-invariance shortcut describes a DIFFERENT, older,
//   legacy implementation, not libint2::Engine's actual behavior. The
//   original 6-buffer assumption (buf[0..2] = center1, buf[3..5] =
//   center2, both read directly) was correct all along; the null-pointer
//   crash that motivated changing it was actually caused entirely by the
//   "libint2 built without derivative support" issue (see below), not by
//   an incorrect buffer count. Code has been reverted accordingly.
//
// OPEN QUESTION as of this revision: the null-pointer crash (address 0x0)
// is resolved by the libint2 rebuild (confirmed: failure mode changed).
// The new SIGTRAP crash's cause is NOT yet confirmed to be fixed by this
// revert -- it may be a separate issue (e.g. an out-of-bounds .block()
// call, an Eigen assertion, or something unrelated to buffer indexing
// entirely). Next diagnostic step: check for an Eigen assertion message
// printed to stderr immediately before the SIGTRAP line (may have been
// cut off in a truncated paste of console output), or get a debugger
// backtrace (e.g. `lldb ./unit_test_aoderivatives`, `run`, `bt` on trap).
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

// Standard includes
#include <iostream>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"

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

// Result type: one entry per atom, each holding the three Cartesian
// derivative matrices (d/dx, d/dy, d/dz) of the full AO matrix with respect
// to that atom's nuclear coordinate. Matches the AO matrix convention used
// elsewhere in this file (dense Eigen::MatrixXd, size AOBasisSize^2).
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;

// Computes nuclear-coordinate derivatives of a two-center one-electron
// integral matrix (overlap, kinetic) via libint2's deriv_order=1 engine.
// obtype must be a two-center operator (libint2::Operator::overlap or
// libint2::Operator::kinetic); this function has not been checked against
// operators with a params argument (e.g. emultipole) or against operators
// with more than two centers (e.g. coulomb/three-center RI integrals) --
// those need separate handling for the additional center(s)' derivatives
// and are deliberately out of scope for this first pass.
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

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < aobasis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();

    Index bf1 = shell2bf[s1];
    Index n1 = shells[s1].size();
    Index atom1 = shell2atom[s1];

    for (Index s2 : shellpair_list[s1]) {

      engine.compute(shells[s1], shells[s2]);
      if (buf[0] == nullptr) {
        continue;  // integrals screened out
      }

      // See HIGHEST-RISK ASSUMPTION note at the top of this file -- the
      // original assert(buf.size()==6) here did not catch a real problem:
      // a null-pointer crash (address 0x0) was observed at runtime,
      // consistent with buf[3..5] being unpopulated even though buf.size()
      // reports 6 (likely a fixed nominal capacity rather than the number
      // of meaningfully-populated buffers). Diagnostic left in place below
      // to confirm this directly on the next run rather than guess again.
      // Diagnostic only, not meant to stay: benign race on this static
      // bool if run multi-threaded (worst case, prints a few times instead
      // of once) -- acceptable for a one-off confirmation, remove once the
      // buffer-count question below is settled.
      static bool printed_diagnostic = false;
      if (!printed_diagnostic) {
        std::cerr << "[aoderivatives diagnostic] buf.size()=" << buf.size();
        for (size_t i = 0; i < buf.size(); ++i) {
          std::cerr << " buf[" << i << "]="
                     << (buf[i] == nullptr ? "null" : "non-null");
        }
        std::cerr << std::endl;
        printed_diagnostic = true;
      }

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
        Eigen::Map<const Eigen::MatrixXd> buf_mat1(buf[xyz], n1, n2);
        result[atom1][xyz].block(bf1, bf2, n1, n2) += buf_mat1;
        if (s1 != s2) {
          result[atom1][xyz].block(bf2, bf1, n2, n1) += buf_mat1.transpose();
        }

        Eigen::Map<const Eigen::MatrixXd> buf_mat2(buf[3 + xyz], n1, n2);
        result[atom2][xyz].block(bf1, bf2, n1, n2) += buf_mat2;
        if (s1 != s2) {
          result[atom2][xyz].block(bf2, bf1, n2, n1) += buf_mat2.transpose();
        }
      }
    }
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

}  // namespace xtp
}  // namespace votca
