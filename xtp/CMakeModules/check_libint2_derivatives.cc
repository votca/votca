// Standalone CMake-configure-time check: does this libint2 installation
// actually WORK for deriv_order=1 (first-derivative) integrals, not just
// claim to via LIBINT2_MAX_DERIV_ORDER?
//
// Background: LIBINT2_MAX_DERIV_ORDER is a library-wide macro (the
// maximum derivative order ANY operator supports), but individual
// libint2 builds can enable derivative support for some operators while
// silently failing for others -- observed directly, on real CI runs, as
// TWO different failure modes for the exact same underlying gap:
//   1. engine.results() returning valid-but-zero-filled buffers (a
//      wrong, but at least non-crashing, silent answer), and
//   2. a hard segfault (null pointer dereference) inside engine.compute()
//      itself, which cannot be caught, detected, or guarded against from
//      ANY application-level code -- no try/catch, no null-checking the
//      returned buffers afterward, nothing -- since the crash happens
//      before control ever returns to the caller.
// Given (2) specifically, the ONLY reliable way to detect this is to
// actually run a real deriv_order=1 computation before the main build
// exists at all -- if THIS standalone program itself crashes, CMake's
// try_run will detect that as a build-time failure and we can disable
// the corresponding functionality entirely, rather than discovering it
// via a segfault deep inside a real user's SCF run.
//
// Deliberately minimal and self-contained: a two-atom, single-primitive-
// per-shell system, using libint2's own native Shell/Engine API
// directly -- NOT any VOTCA wrapper (AOBasis, QMMolecule, etc.), none of
// which exist yet at this point in the build.
//
// Exit code 0: every operator this branch actually uses for derivative
// integrals (overlap, kinetic, nuclear attraction with a point charge,
// and the two-center Coulomb metric used for RI) computed a genuinely
// non-zero, non-null result. Any other outcome (non-zero exit, or the
// process not surviving at all, which is exactly what a crash looks
// like to try_run) means this libint2 cannot be trusted for derivative
// integrals.

#include <libint2.hpp>

#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

libint2::Shell MakeSShell(double exponent, std::array<double, 3> center) {
  libint2::Shell::do_enforce_unit_normalization(false);
  return libint2::Shell{{exponent},
                        {{0, false, {1.0}}},
                        {{center[0], center[1], center[2]}}};
}

// Returns true if buf has at least one non-null pointer among the
// first n_centers*3 slots AND the values behind those pointers are not
// all exactly (or suspiciously near) zero -- the direct check for the
// "silently wrong" failure mode. A crash inside engine.compute() itself
// (the OTHER failure mode) will simply end this whole program before
// this function is ever reached, which is exactly what try_run is
// meant to detect at the CMake level.
bool BuffersLookReal(const libint2::Engine::target_ptr_vec& buf,
                     int n_centers, std::size_t n1, std::size_t n2) {
  bool saw_nonnull = false;
  double norm_sq = 0.0;
  for (int c = 0; c < n_centers; ++c) {
    for (int xyz = 0; xyz < 3; ++xyz) {
      const double* p = buf[static_cast<std::size_t>(c * 3 + xyz)];
      if (p == nullptr) {
        continue;
      }
      saw_nonnull = true;
      for (std::size_t i = 0; i < n1 * n2; ++i) {
        norm_sq += p[i] * p[i];
      }
    }
  }
  return saw_nonnull && norm_sq > 1.e-20;
}

}  // namespace

int main() {
  libint2::initialize();
  int failures = 0;

  const std::array<double, 3> center0{0.0, 0.0, 0.0};
  const std::array<double, 3> center1{0.0, 0.0, 1.4};
  libint2::Shell s0 = MakeSShell(1.0, center0);
  libint2::Shell s1 = MakeSShell(1.0, center1);
  const std::size_t n = 1;  // single s-function per shell

  // --- overlap and kinetic (two-center, deriv_order=1) ---
  for (libint2::Operator op :
      {libint2::Operator::overlap, libint2::Operator::kinetic}) {
    libint2::Engine engine(op, 1, 0, 1);
    engine.compute(s0, s1);
    const auto& buf = engine.results();
    if (!BuffersLookReal(buf, 2, n, n)) {
      std::cerr << "FAIL: operator " << static_cast<int>(op)
                << " (overlap=0, kinetic=1) produced no real "
                   "derivative buffers.\n";
      ++failures;
    }
  }

  // --- nuclear attraction with an explicit point charge (three-center-
  //     like: shell1's atom, shell2's atom, the point charge itself;
  //     deriv_order=1) ---
  {
    libint2::Engine engine(libint2::Operator::nuclear, 1, 0, 1);
    std::vector<libint2::Atom> point_charge(1);
    point_charge[0].atomic_number = 1;
    point_charge[0].x = 0.0;
    point_charge[0].y = 0.0;
    point_charge[0].z = 0.7;
    engine.set_params(libint2::make_point_charges(point_charge));
    engine.compute(s0, s1);
    const auto& buf = engine.results();
    if (!BuffersLookReal(buf, 3, n, n)) {
      std::cerr << "FAIL: nuclear attraction (point charge) produced no "
                   "real derivative buffers.\n";
      ++failures;
    }
  }

  // --- two-center Coulomb metric, via unit shells (deriv_order=1) ---
  {
    libint2::Engine engine(libint2::Operator::coulomb, 1, 0, 1);
    engine.set(libint2::BraKet::xs_xs);
    engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 1>(
        s0, libint2::Shell::unit(), s1, libint2::Shell::unit());
    const auto& buf = engine.results();
    if (!BuffersLookReal(buf, 2, n, n)) {
      std::cerr << "FAIL: two-center Coulomb metric produced no real "
                   "derivative buffers.\n";
      ++failures;
    }
  }

  libint2::finalize();

  if (failures > 0) {
    std::cerr << failures << " operator(s) failed the derivative-integral "
                            "check.\n";
    return 1;
  }
  std::cout << "All derivative-integral checks passed.\n";
  return 0;
}
