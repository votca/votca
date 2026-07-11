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
// STATUS: compiles and links cleanly as of the fix below (confirmed against
// a real build attempt). NOT YET RUNTIME-TESTED -- the finite-difference
// check in test_aoderivatives.cc has not yet been run to confirm the
// numbers coming out of this code are actually correct, only that it
// builds. Run that test next, before trusting any number this produces.
//
// Build history on this file:
//   - Originally drafted with no local libint2 install available, so it
//     could not be compiled or run at all when first written.
//   - First real build attempt failed at link time: duplicate symbols for
//     libint2's internal static tables (FmEval_Chebyshev7<double>::
//     cheb_table, TennoGmEval<double>::cheb_table), because this file
//     originally included <libint2/statics_definition.h> a second time --
//     that header defines storage and must appear in exactly one
//     translation unit per library. libint2_calls.cc already provides it
//     for the whole votca_xtp library; the duplicate include here has
//     been removed (see include block below).
//
// HIGHEST-RISK REMAINING ASSUMPTION, to check via the finite-difference
// test before trusting numerical results:
//   For a two-center one-electron integral engine constructed with
//   deriv_order = 1, this code assumes engine.results() returns exactly 6
//   buffers per shell pair: indices 0,1,2 are d/dx,dy,dz with respect to
//   shell1's atom center, and indices 3,4,5 are d/dx,dy,dz with respect to
//   shell2's atom center. Compiling successfully confirms the API *shape*
//   used here (constructor arguments, .results(), indexing) matches this
//   libint2 version -- it says nothing about whether 6 is actually the
//   right buffer count, since that assert is a runtime check. This is
//   only genuinely confirmed once the test binary has actually been run.
// ===========================================================================

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

      // See HIGHEST-RISK ASSUMPTION note at the top of this file.
      assert(buf.size() == 6 &&
             "Unexpected number of derivative buffers for a two-center "
             "one-electron integral at deriv_order=1 -- check this "
             "libint2 version's buffer-ordering convention before "
             "proceeding (see file header comment).");

      Index bf2 = shell2bf[s2];
      Index n2 = shells[s2].size();
      Index atom2 = shell2atom[s2];

      for (Index xyz = 0; xyz < 3; ++xyz) {
        // derivative with respect to shell1's atom (buf indices 0,1,2)
        {
          Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[xyz], n1, n2);
          result[atom1][xyz].block(bf1, bf2, n1, n2) += buf_mat;
          if (s1 != s2) {
            result[atom1][xyz].block(bf2, bf1, n2, n1) += buf_mat.transpose();
          }
        }
        // derivative with respect to shell2's atom (buf indices 3,4,5)
        {
          Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[3 + xyz], n1, n2);
          result[atom2][xyz].block(bf1, bf2, n1, n2) += buf_mat;
          if (s1 != s2) {
            result[atom2][xyz].block(bf2, bf1, n2, n1) += buf_mat.transpose();
          }
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
