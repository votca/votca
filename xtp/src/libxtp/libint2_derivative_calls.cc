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

// Nuclear-coordinate derivatives of AO integrals (overlap, kinetic, nuclear
// attraction, two-centre Coulomb metric, three-centre RI) via libint2 engines
// with deriv_order = 1, for the analytic DFT gradient. Derivatives are with
// respect to atomic positions in Bohr.
//
// Conventions:
// - libint2::Engine returns the derivatives with respect to every centre
//   explicitly (buffers ordered by centre, then x, y, z); none is omitted via
//   translational invariance.
// - Shell-pair buffers are row-major and are mapped with MatrixLibInt, as in
//   libint2_calls.cc. A column-major map is indistinguishable for s shells
//   only.
// - Each derivative category needs its own libint2 configure option (see the
//   compile guards below); missing support makes the functions throw.
// - Exceptions inside OpenMP regions are captured and rethrown afterwards.

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/votca_xtp_config.h"
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
// <libint2/statics_definition.h> defines libint2's static tables and must be
// included in exactly one translation unit; libint2_calls.cc does that.
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

// libint2 returns each shell-pair (derivative) buffer in row-major order,
// n1 rows x n2 columns, as in libint2_calls.cc. Mapping it as column-major
// transposes every block with p or higher shells; for s shells (1x1 blocks)
// both agree, so tests on all-s bases cannot detect it.
using MatrixLibInt =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

namespace {
// Thrown by the stubs compiled when a derivative category (ONEBODY, ERI2,
// ERI3; independent libint2 configure options) is unavailable.
// [[maybe_unused]] because no stub is compiled if all three are supported
// (-Werror -Wunused-function).
[[noreturn, maybe_unused]] void ThrowNoDerivativeSupport(
    const char* function_name, const char* configure_flag) {
  throw std::runtime_error(
      std::string(function_name) +
      ": this libint2 build was compiled without derivative integral "
      "support for this operator category, so analytic nuclear forces "
      "are unavailable. Rebuild libint2 with " +
      configure_flag +
      " to use this feature. Many packaged libint2 builds (Homebrew, "
      "Ubuntu apt, etc.) do not enable it; one-body, two-center Coulomb "
      "metric and three-center RI integrals each have their own "
      "configure option.");
}
}  // namespace

// Nuclear-coordinate derivatives of a two-centre one-electron integral
// matrix for obtype = libint2::Operator::overlap or kinetic (operators with
// parameters or more centres are not handled).
//
// Compile guard: LIBINT2_MAX_DERIV_ORDER is the maximum over all operator
// categories, so it is not sufficient; a category built without derivatives
// fails inside engine.compute() (null buffers or an uncatchable crash). The
// per-category macros LIBINT_INCLUDE_ONEBODY/_ERI2/_ERI3 hold the configured
// derivative order (LIBINT2_ENABLE_ONEBODY/_ERI2/_ERI3, -1 = off), so their
// value must be checked, not just whether they are defined: ONEBODY defaults
// to 0 (energies only), e.g. in the Homebrew formula, which only enables
// ERI/ERI2/ERI3, while Psi4's default enables ONEBODY but not ERI2/ERI3.
//
// ONEBODY: computeOneBodyIntegralDerivatives,
// ComputeNuclearAttractionDerivatives; ERI2: ComputeCoulombMetricDerivatives;
// ERI3: the three-centre functions.
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1 && \
    (LIBINT_INCLUDE_ONEBODY >= 1)
template <libint2::Operator obtype>
std::vector<AOMatrixDerivative> computeOneBodyIntegralDerivatives(
    const AOBasis& aobasis) {

  static_assert(
      libint2::operator_traits<obtype>::nopers == 1,
      "computeOneBodyIntegralDerivatives only supports single-operator "
      "types (overlap, kinetic); multipole-type operators with multiple "
      "result matrices are not handled.");

  // AOBasis has no getNumofAtoms(); getFuncPerAtom() has one entry per atom.
  Index natoms = static_cast<Index>(aobasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Shell> shells = aobasis.GenerateLibintBasis();
  std::vector<std::vector<Index>> shellpair_list = aobasis.ComputeShellPairs();
  std::vector<Index> shell2bf = aobasis.getMapToBasisFunctions();

  // shell -> atom index, to assign each derivative buffer to its atom
  std::vector<Index> shell2atom;
  shell2atom.reserve(shells.size());
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
  // Same as computeOneBodyIntegrals(), but with deriv_order = 1.
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

        engine.compute(shells[s1], shells[s2]);
        // buf[0..2]: derivatives w.r.t. the centre of s1, buf[3..5]: w.r.t.
        // the centre of s2. Both must be checked before dereferencing, since
        // buf.size() reports 6 even if buf[3..5] are unpopulated.
        if (buf[0] == nullptr || buf[3] == nullptr) {
          continue;  // screened out, or no derivative support (checked
                     // via any_nonnull_buffer after the loop)
        }
        any_nonnull_buffer.store(true, std::memory_order_relaxed);

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
      // An exception escaping an OpenMP thread is undefined behaviour;
      // keep the first one and rethrow after the parallel region.
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
    // Runtime counterpart of the compile guard: obtype never produced a
    // non-null buffer.
    throw std::runtime_error(
        "computeOneBodyIntegralDerivatives: engine.results() returned a "
        "null buffer for every shell pair; this libint2 build does not "
        "support this operator's derivative integrals at runtime, even "
        "though it may report LIBINT2_MAX_DERIV_ORDER >= 1 for other "
        "operators. Rebuild libint2 with --enable-1body=1 to use this "
        "feature.");
  }
  return result;
}

// Free-function counterparts of AOOverlap::Fill / AOKinetic::Fill.
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis) {
  return computeOneBodyIntegralDerivatives<libint2::Operator::overlap>(aobasis);
}

std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis) {
  return computeOneBodyIntegralDerivatives<libint2::Operator::kinetic>(aobasis);
}
#else   // !(LIBINT_INCLUDE_ONEBODY)
// Stubs for overlap and kinetic without ONEBODY derivative support; ERI2 and
// ERI3 are guarded separately below.
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeOverlapDerivatives", "--enable-1body=1");
}
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeKineticDerivatives", "--enable-1body=1");
}
#endif  // LIBINT_INCLUDE_ONEBODY

// ===========================================================================
// Two-centre Coulomb metric derivatives d(P|Q)/dR for the RI-J/RI-K
// gradient, dE_J/dR = sum_P c_P d(d_P)/dR - 1/2 sum_PQ c_P c_Q d(V_PQ)/dR.
//
// As in AOCoulomb::computeCoulombIntegrals, (P|Q) is computed as (P 1|Q 1)
// with unit shells (BraKet::xs_xs), which have no centre; libint2 returns 6
// buffers, the derivatives w.r.t. the centres of P and Q.
//
// Requires LIBINT_INCLUDE_ERI2 >= 1 (libint2 --enable-eri2), independent of
// ONEBODY and ERI3.
#if (LIBINT_INCLUDE_ERI2 >= 1)
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
  engines[0] =
      libint2::Engine(libint2::Operator::coulomb, aobasis.getMaxNprim(),
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

        // see computeOneBodyIntegralDerivatives
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
    throw std::runtime_error(
        "ComputeCoulombMetricDerivatives: engine.results() returned a "
        "null buffer for every shell pair; this libint2 build does not "
        "support this operator's derivative integrals at runtime, even "
        "though it may report LIBINT2_MAX_DERIV_ORDER >= 1 for other "
        "operators. Rebuild libint2 with --enable-eri2=1 to use this "
        "feature.");
  }
  return result;
}
#else   // !(LIBINT_INCLUDE_ERI2)
std::vector<AOMatrixDerivative> ComputeCoulombMetricDerivatives(
    const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeCoulombMetricDerivatives",
                           "--enable-eri2=1");
}
#endif  // LIBINT_INCLUDE_ERI2

// ===========================================================================
// Three-centre RI integral derivatives d(mu,nu|P)/dR. Same shell order as
// ComputeAO3cBlock (auxshell, unit, shell_col, shell_row, BraKet::xs_xx);
// libint2 returns 9 buffers, the derivatives w.r.t. the centres of the aux,
// col and row shells, each x, y, z.
//
// ComputeThreeCenterDerivatives stores the full tensor,
// O(natoms * 3 * Naux * Nao^2) memory, and is only meant for tests; the
// gradients use the contracted variants below.
//
// Requires LIBINT_INCLUDE_ERI3 >= 1 (libint2 --enable-eri3), independent of
// ONEBODY and ERI2.
// ===========================================================================

// One entry per Cartesian direction (x,y,z); each holds one matrix per
// aux basis function (indexed by absolute aux AO index), each matrix
// being (dftbasis size x dftbasis size).
using ThreeCenterDerivative = std::array<std::vector<Eigen::MatrixXd>, 3>;

#if (LIBINT_INCLUDE_ERI3 >= 1)
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis& auxbasis, const AOBasis& dftbasis) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  // assumes auxbasis and dftbasis are defined on the same molecule

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

        // Full col loop: the mu,nu symmetry of (mu,nu|P) is not exploited
        // here, unlike in ComputeAO3cBlock.
        for (Index col = 0; col < dftbasis.getNumofShells(); ++col) {
          const libint2::Shell& shell_col = dftshells[col];
          Index col_start = shell2bf[col];
          Index atom_col = dftshell2atom[col];

          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(
                  auxshell, libint2::Shell::unit(), shell_col, shell_row);

          // see computeOneBodyIntegralDerivatives; three centres here
          if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
            continue;
          }
          any_nonnull_buffer_3c.store(true, std::memory_order_relaxed);

          // buffers: [aux centre xyz][col centre xyz][row centre xyz], in
          // the shell argument order of compute2
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
    throw std::runtime_error(
        "ComputeThreeCenterDerivatives: engine.results() returned a null "
        "buffer for every shell triple; this libint2 build does not "
        "support this operator's derivative integrals at runtime, even "
        "though it may report LIBINT2_MAX_DERIV_ORDER >= 1 for other "
        "operators. Rebuild libint2 with --enable-eri3=1 to use this "
        "feature.");
  }
  return result;
}
#else   // !(LIBINT_INCLUDE_ERI3)
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivatives(
    const AOBasis&, const AOBasis&) {
  ThrowNoDerivativeSupport("ComputeThreeCenterDerivatives", "--enable-eri3=1");
}
#endif  // LIBINT_INCLUDE_ERI3

#if (LIBINT_INCLUDE_ERI3 >= 1)
// Density-contracted three-centre derivatives for DFTGradient::RIJGradient:
//   result[a](xyz, P) = sum_{mu,nu} D_{mu nu} d(mu,nu|P)/dR_a,xyz.
// Each shell triple is contracted as soon as it is computed, so memory is
// O(natoms * 3 * Naux) instead of the O(natoms * 3 * Naux * Nao^2) of
// ComputeThreeCenterDerivatives (e.g. 53 atoms, 943 AOs, 2320 aux functions
// would need ~2.6 PB). Loop structure and buffer layout as in
// ComputeThreeCenterDerivatives.
std::vector<Eigen::MatrixXd> ComputeThreeCenterDerivativeContraction(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& density) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());

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

  Index n_aux_bf = auxbasis.AOBasisSize();

  // result[a]: (3, n_aux_bf), rows xyz, columns aux function p
  std::vector<Eigen::MatrixXd> result(natoms,
                                      Eigen::MatrixXd::Zero(3, n_aux_bf));

  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(
      libint2::Operator::coulomb,
      std::max(dftbasis.getMaxNprim(), auxbasis.getMaxNprim()),
      static_cast<int>(std::max(dftbasis.getMaxL(), auxbasis.getMaxL())), 1);
  engines[0].set(libint2::BraKet::xs_xx);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  // No synchronization needed: the loop is parallel over aux shells, so
  // different threads write to different columns (global_aux) of
  // result[atom].
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

        for (Index col = 0; col < dftbasis.getNumofShells(); ++col) {
          const libint2::Shell& shell_col = dftshells[col];
          Index col_start = shell2bf[col];
          Index atom_col = dftshell2atom[col];

          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(
                  auxshell, libint2::Shell::unit(), shell_col, shell_row);

          if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
            continue;
          }
          any_nonnull_buffer_3c.store(true, std::memory_order_relaxed);

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
                  // contract with the density immediately
                  double d_rc = density(r, c);
                  result[atom_aux](xyz, global_aux) += d_rc * val_aux;
                  result[atom_col](xyz, global_aux) += d_rc * val_col;
                  result[atom_row](xyz, global_aux) += d_rc * val_row;
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
    throw std::runtime_error(
        "ComputeThreeCenterDerivativeContraction: engine.results() "
        "returned a null buffer for every shell triple; this libint2 "
        "build does not support this operator's derivative integrals at "
        "runtime, even though it may report LIBINT2_MAX_DERIV_ORDER >= 1 "
        "for other operators. Rebuild libint2 with --enable-eri3=1 to "
        "use this feature.");
  }
  return result;
}
#else   // !(LIBINT_INCLUDE_ERI3)
std::vector<Eigen::MatrixXd> ComputeThreeCenterDerivativeContraction(
    const AOBasis&, const AOBasis&, const Eigen::MatrixXd&) {
  ThrowNoDerivativeSupport("ComputeThreeCenterDerivativeContraction",
                           "--enable-eri3=1");
}
#endif  // LIBINT_INCLUDE_ERI3

#if (LIBINT_INCLUDE_ERI3 >= 1)
// Full (Nao x Nao) three-centre derivative matrices for one atom only,
// (3, n_aux_bf) of them. Shell triples with no shell on target_atom are
// skipped, but a loop over all atoms still evaluates each triple up to
// three times. Superseded for RIKGradient by
// ComputeThreeCenterDerivativesMOTransformed, which needs a single pass.
ThreeCenterDerivative ComputeThreeCenterDerivativesForAtom(
    const AOBasis& auxbasis, const AOBasis& dftbasis, Index target_atom) {
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
  Index nthreads = OPENMP::getMaxThreads();

  // (3, n_aux_bf) matrices of size nao x nao for target_atom only
  ThreeCenterDerivative result;
  for (Index xyz = 0; xyz < 3; ++xyz) {
    result[xyz] = std::vector<Eigen::MatrixXd>(
        n_aux_bf, Eigen::MatrixXd::Zero(n_dft_bf, n_dft_bf));
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

        for (Index col = 0; col < dftbasis.getNumofShells(); ++col) {
          const libint2::Shell& shell_col = dftshells[col];
          Index col_start = shell2bf[col];
          Index atom_col = dftshell2atom[col];

          // no shell on target_atom: no contribution
          if (atom_aux != target_atom && atom_col != target_atom &&
              atom_row != target_atom) {
            continue;
          }

          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(
                  auxshell, libint2::Shell::unit(), shell_col, shell_row);

          if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
            continue;
          }
          any_nonnull_buffer_3c.store(true, std::memory_order_relaxed);

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
                  Index r = row_start + static_cast<Index>(row_c);
                  Index c = col_start + static_cast<Index>(col_c);
                  // only the centres located on target_atom contribute
                  if (atom_aux == target_atom) {
                    result[xyz][global_aux](r, c) +=
                        result_aux(aux_c, col_c, row_c);
                  }
                  if (atom_col == target_atom) {
                    result[xyz][global_aux](r, c) +=
                        result_col(aux_c, col_c, row_c);
                  }
                  if (atom_row == target_atom) {
                    result[xyz][global_aux](r, c) +=
                        result_row(aux_c, col_c, row_c);
                  }
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
    throw std::runtime_error(
        "ComputeThreeCenterDerivativesForAtom: engine.results() returned "
        "a null buffer for every shell triple touching this atom; this "
        "libint2 build does not support this operator's derivative "
        "integrals at runtime, even though it may report "
        "LIBINT2_MAX_DERIV_ORDER >= 1 for other operators. Rebuild "
        "libint2 with --enable-eri3=1 to use this feature.");
  }
  return result;
}
#else   // !(LIBINT_INCLUDE_ERI3)
ThreeCenterDerivative ComputeThreeCenterDerivativesForAtom(const AOBasis&,
                                                           const AOBasis&,
                                                           Index) {
  ThrowNoDerivativeSupport("ComputeThreeCenterDerivativesForAtom",
                           "--enable-eri3=1");
}
#endif  // LIBINT_INCLUDE_ERI3

#if (LIBINT_INCLUDE_ERI3 >= 1)
// Three-centre derivatives transformed to occupied MOs, for
// DFTGradient::RIKGradient, in a single pass over the shell triples:
//   result[a][xyz][P](i, j) = sum_{mu,nu} C_mu,i C_nu,j d(mu,nu|P)/dR_a,xyz
// with C = occ_mo_coeffs, i.e. (nocc x nocc) matrices (ThreeCenterDerivative
// is shape-agnostic). Memory natoms * 3 * Naux * nocc^2 doubles (~24 GB for
// 53 atoms, 2320 aux functions, 93 occupied orbitals). Each shell pair and
// aux function contributes C_row^T B C_col with the (nrow x ncol) integral
// block B.
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivativesMOTransformed(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& occ_mo_coeffs) {
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index nocc = occ_mo_coeffs.cols();

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

  Index n_aux_bf = auxbasis.AOBasisSize();

  // result[a][xyz][p]: (nocc, nocc)
  std::vector<ThreeCenterDerivative> result(natoms);
  for (Index a = 0; a < natoms; ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      result[a][xyz] = std::vector<Eigen::MatrixXd>(
          n_aux_bf, Eigen::MatrixXd::Zero(nocc, nocc));
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

  // No synchronization needed: threads write to different global_aux
  // entries of result[atom][xyz].
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

        for (Index col = 0; col < dftbasis.getNumofShells(); ++col) {
          const libint2::Shell& shell_col = dftshells[col];
          Index col_start = shell2bf[col];
          Index atom_col = dftshell2atom[col];

          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(
                  auxshell, libint2::Shell::unit(), shell_col, shell_row);

          if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
            continue;
          }
          any_nonnull_buffer_3c.store(true, std::memory_order_relaxed);

          // MO coefficient rows of this shell pair
          Index nrow = static_cast<Index>(shell_row.size());
          Index ncol = static_cast<Index>(shell_col.size());
          Eigen::MatrixXd mo_row_block =
              occ_mo_coeffs.middleRows(row_start, nrow);
          Eigen::MatrixXd mo_col_block =
              occ_mo_coeffs.middleRows(col_start, ncol);

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

              // Copy the (nrow, ncol) block of this aux function and
              // transform it with two small matrix products; one rank-1
              // update per AO element would cost nao^2 * Naux * nocc^2.
              Eigen::MatrixXd block_aux(nrow, ncol);
              Eigen::MatrixXd block_col(nrow, ncol);
              Eigen::MatrixXd block_row(nrow, ncol);
              for (Index col_c = 0; col_c < ncol; ++col_c) {
                for (Index row_c = 0; row_c < nrow; ++row_c) {
                  block_aux(row_c, col_c) =
                      result_aux(aux_c, static_cast<size_t>(col_c),
                                 static_cast<size_t>(row_c));
                  block_col(row_c, col_c) =
                      result_col(aux_c, static_cast<size_t>(col_c),
                                 static_cast<size_t>(row_c));
                  block_row(row_c, col_c) =
                      result_row(aux_c, static_cast<size_t>(col_c),
                                 static_cast<size_t>(row_c));
                }
              }
              result[atom_aux][xyz][global_aux].noalias() +=
                  mo_row_block.transpose() * block_aux * mo_col_block;
              result[atom_col][xyz][global_aux].noalias() +=
                  mo_row_block.transpose() * block_col * mo_col_block;
              result[atom_row][xyz][global_aux].noalias() +=
                  mo_row_block.transpose() * block_row * mo_col_block;
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
    throw std::runtime_error(
        "ComputeThreeCenterDerivativesMOTransformed: engine.results() "
        "returned a null buffer for every shell triple; this libint2 "
        "build does not support this operator's derivative integrals at "
        "runtime, even though it may report LIBINT2_MAX_DERIV_ORDER >= 1 "
        "for other operators. Rebuild libint2 with --enable-eri3=1 to "
        "use this feature.");
  }
  return result;
}
#else   // !(LIBINT_INCLUDE_ERI3)
std::vector<ThreeCenterDerivative> ComputeThreeCenterDerivativesMOTransformed(
    const AOBasis&, const AOBasis&, const Eigen::MatrixXd&) {
  ThrowNoDerivativeSupport("ComputeThreeCenterDerivativesMOTransformed",
                           "--enable-eri3=1");
}
#endif  // LIBINT_INCLUDE_ERI3

// Three-centre integrals (mu,nu|P) (deriv_order = 0) for the RI gradients,
// one Nao x Nao matrix per aux function, built with ComputeAO3cBlock. Kept
// here so that libint2 is only used in this file and libint2_calls.cc. No
// compile guard needed.
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
// Nuclear attraction integral derivatives d(V_ne)/dR,
//   V_ne,munu = -sum_A Z_A int chi_mu(r) 1/|r - R_A| chi_nu(r) dr.
// For each point charge A this is a three-centre derivative (both basis
// function centres and R_A); as in libint2's
// compute_1body_ints_deriv<Operator::nuclear> (lcao/1body.h), each atom is
// set as a single point charge, giving 9 buffers ordered
// [shell1 centre][shell2 centre][R_A], each x, y, z.
//
// libint2::Operator::nuclear already includes the minus sign (attractive),
// consistent with AOMultipole::FillPotential/SetupH0, so the buffers are
// added without negation.
// ===========================================================================
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1 && \
    (LIBINT_INCLUDE_ONEBODY >= 1)
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
    engines[i] =
        libint2::Engine(libint2::Operator::nuclear, aobasis.getMaxNprim(),
                        static_cast<int>(aobasis.getMaxL()), 1);
  }

  // Parallel over point charges (each needs engine.set_params()), with
  // per-thread result matrices summed afterwards.
  std::exception_ptr eptr_nucattr = nullptr;
  std::atomic<bool> any_nonnull_buffer_nucattr{false};
#pragma omp parallel for schedule(dynamic)
  for (Index A = 0; A < natoms; ++A) {
    try {
      Index thread_id = OPENMP::getThreadId();
      libint2::Engine& engine = engines[thread_id];

      std::vector<libint2::Atom> single_atom(1);
      single_atom[0].atomic_number = static_cast<int>(mol[A].getNuccharge());
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

          // see computeOneBodyIntegralDerivatives; three centres here.
          // Zero-filled (non-null) buffers are caught by the norm check
          // after the loop.
          if (buf[0] == nullptr || buf[3] == nullptr || buf[6] == nullptr) {
            continue;
          }
          any_nonnull_buffer_nucattr.store(true, std::memory_order_relaxed);

          Index bf2 = shell2bf[s2];
          Index n2 = shells[s2].size();
          Index atom2 = shell2atom[s2];

          for (Index xyz = 0; xyz < 3; ++xyz) {
            // no negation, see above
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
    // A build without derivative support for this operator can return
    // null buffers instead of throwing, which would silently give zero.
    throw std::runtime_error(
        "ComputeNuclearAttractionDerivatives: engine.results() returned "
        "a null buffer for every shell pair and point charge; this "
        "libint2 build does not support nuclear-attraction derivative "
        "integrals at runtime, even though it may report "
        "LIBINT2_MAX_DERIV_ORDER >= 1 for other operators. Rebuild "
        "libint2 with --enable-1body=1 to use this feature.");
  }

  // Some builds instead return non-null but zero-filled buffers. The
  // derivative cannot vanish identically for a molecule with nonzero
  // nuclear charges.
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
        "exactly zero everywhere, which is impossible for a real "
        "molecule; this libint2 build likely returns zero-filled buffers "
        "for this operator. Rebuild libint2 with --enable-1body=1 to use "
        "this feature.");
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
#else   // !(LIBINT_INCLUDE_ONEBODY)
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis&, const QMMolecule&) {
  ThrowNoDerivativeSupport("ComputeNuclearAttractionDerivatives",
                           "--enable-1body=1");
}
#endif  // LIBINT_INCLUDE_ONEBODY

// Lets code that does not include libint2 headers (e.g. dftengine.cc) check
// for derivative support. The forces need all three categories: ONEBODY
// (overlap, kinetic, nuclear attraction), ERI2 (Coulomb metric) and ERI3
// (three-centre RI).
bool HasLibint2DerivativeSupport() {
#if defined(LIBINT2_MAX_DERIV_ORDER) && LIBINT2_MAX_DERIV_ORDER >= 1 && \
    (LIBINT_INCLUDE_ONEBODY >= 1) && (LIBINT_INCLUDE_ERI2 >= 1) &&      \
    (LIBINT_INCLUDE_ERI3 >= 1)
  return true;
#else
  return false;
#endif
}

}  // namespace xtp
}  // namespace votca
