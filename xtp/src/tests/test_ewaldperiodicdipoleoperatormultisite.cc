/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ewaldperiodicdipoleoperatormultisite_test

// Standard includes
#include <cmath>
#include <iostream>

// Third party includes
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldperiodicdipoleoperator.h"

using namespace votca::xtp;
using namespace votca;

// Neither test_ewaldperiodicdipoleoperator.cc nor
// test_ewaldperiodicdipolesymmetry.cc exercises a genuinely multi-site
// segment (both used one site per segment throughout) -- this test
// specifically targets the code path introduced when
// EwaldPeriodicDipoleOperator was generalized from a fixed
// one-site-per-segment layout to a variable-width (offset-table) one,
// AND (in later revisions of this same file) the intramolecular
// erfc-screened, Thole-damped coupling term added once it was
// discovered that omitting it was a real correctness gap, not a
// deliberate design choice -- see this operator's own class
// documentation for the full story of how that was found and, across
// several rounds, corrected: first from a mistaken Thole-damped-via-
// eeInteractor version, to a mistaken erfc-only-UNDAMPED version
// (reasoning by analogy from the permanent-field case, which does not
// transfer), to the current erfc-screened AND Thole-damped one, which
// matches legacy's own FU12_ERFC_At_By at both its intermolecular and
// intramolecular call sites. Also see that class's own documentation
// for a fourth, related fix in EwaldReciprocalSpaceSum this test's own
// reference solve does not directly exercise (this test's molecule has
// no static multipoles, so the reciprocal-space exclusion-removal
// change does not affect its own permanent field values -- see the
// class documentation for why that change was still necessary for the
// general case).
//
// An earlier version of this test (before any of that) checked the
// OPPOSITE property -- that two close sites of the same segment did NOT
// couple at all. That was the correct expectation for the code as it
// existed then, but is now the wrong property to test: two sites of the
// same segment SHOULD couple, via an erfc-screened, Thole-damped
// interaction, matching legacy PolarBackground's own FU12_ERFC_At_By.
//
// What this test actually checks now: EwaldPeriodicDipoleOperator's own
// PCG solve, for a small 2-site-one-segment system, must agree with an
// INDEPENDENT direct dense linear solve of the exact same
// self-consistent equation (built directly from getPInv() and
// EwaldRealSpaceInteractor::ComputeB/ComputeThole, not by calling
// anything on the operator itself) -- a genuine cross-check on the
// operator's own multiply() implementation, not merely a
// self-consistency check against itself.
BOOST_AUTO_TEST_SUITE(ewaldperiodicdipoleoperatormultisite_test)

BOOST_AUTO_TEST_CASE(intramolecular_coupling_matches_direct_solve) {
  Eigen::Matrix3d box = 40.0 * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;

  // Fixed external background: one positive charge. Deliberately NOT at
  // a position that's a simple fraction of the box length along a
  // single axis relative to the molecule below (e.g. exactly half the
  // box length, as an earlier version of this test had) -- that kind of
  // special symmetric placement can make the reciprocal-space field at
  // a given site cancel to near machine-epsilon (confirmed directly:
  // reproducing that earlier geometry standalone gave a field of order
  // 1e-19, not a small-but-genuine physical value), which is a fragile
  // foundation for a numerical comparison test: when a "real" RHS value
  // is actually just floating-point noise, tiny summation-order
  // differences between this test's own reference computation and the
  // operator's real code path can produce large relative differences in
  // the final solved result, for reasons having nothing to do with
  // whether the underlying physics is correct.
  PolarSegment external("seg", 1);
  PolarSite charge(1, "H", Eigen::Vector3d(15.3, 4.7, -8.2));
  charge.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 1.0;
  charge.setMultipole(mpoles, 0);
  external.push_back(charge);
  registry.Register(1, EwaldChargeState::Neutral, external);

  // The multi-site segment being solved for: two polarizable sites at a
  // realistic bonded distance (2.0 bohr, roughly a C-H bond length) so
  // intramolecular erfc-screened coupling is a genuine, substantial
  // effect -- not negligible -- while the resulting self-consistent
  // system stays physically well-posed (positive definite; confirmed
  // numerically, see below). An earlier version of this test used 0.5
  // bohr, deliberately chosen to make the coupling "large and easy to
  // see" -- but that reasoning didn't go far enough: at that separation,
  // the erfc-screened tensor's own diagonal (~16-376 depending on other
  // parameters) dwarfs D=1/alpha_pol (~0.5), and unlike Thole damping
  // (which exists specifically to regularize this), plain erfc screening
  // does NOT prevent a polarization-catastrophe-like instability at
  // short range -- so the intended self-consistent system was genuinely
  // NOT positive definite at 0.5 bohr (confirmed by computing its
  // eigenvalues directly: several came out negative). CG's own
  // guarantees assume positive-definiteness, so a numerically
  // ill-conditioned or indefinite system is not something this test
  // should be probing incidentally as a side effect of an unrealistic
  // geometry choice. 2.0 bohr keeps the coupling real and substantial
  // (worth still checking a broken implementation would be caught) while
  // keeping the system itself physically sane.
  PolarSegment molecule("seg", 2);
  const double alpha_pol = 2.0;
  const double thole_a = 0.39;
  PolarSite site_a(1, "H", Eigen::Vector3d(3.4, -1.9, 0.7));
  site_a.setpolarization(alpha_pol * Eigen::Matrix3d::Identity());
  molecule.push_back(site_a);
  PolarSite site_b(2, "H", Eigen::Vector3d(5.4, -1.9, 0.7));
  site_b.setpolarization(alpha_pol * Eigen::Matrix3d::Identity());
  molecule.push_back(site_b);
  registry.Register(2, EwaldChargeState::Neutral, molecule);

  const double alpha_ewald = 0.3;
  EwaldRealSpaceSum real_sum(box, registry, alpha_ewald, thole_a,
                             /*r_min=*/12.0, /*field_tol=*/1e-12);
  // See the k_max comment in
  // converged_dipoles_are_independent_of_the_splitting below: 12*alpha,
  // not a flat literal.
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_ewald,
                                    /*k_max=*/12.0 * alpha_ewald);
  const double volume = box.determinant();
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  // Permanent field at each site of the molecule, individually -- this
  // is also exactly b for the PCG solve below (b = +V_permanent; see
  // EwaldPeriodicDipoleOperator's own class documentation for the sign).
  PolarSegment& mseg = registry.Get(2, EwaldChargeState::Neutral);
  Eigen::Vector3d V_a, V_b;
  {
    PolarSite& a = mseg[0];
    a.Reset();
    real_sum.AddFieldAt<Estatic::V>(2, a, EwaldChargeState::Neutral);
    recip_sum.AddFieldAt<Estatic::V>(a, EwaldChargeState::Neutral);
    V_a = a.V();
    a.Reset();

    PolarSite& b = mseg[1];
    b.Reset();
    real_sum.AddFieldAt<Estatic::V>(2, b, EwaldChargeState::Neutral);
    recip_sum.AddFieldAt<Estatic::V>(b, EwaldChargeState::Neutral);
    V_b = b.V();
    b.Reset();
  }

  // Sanity: this doesn't check a specific direction the way the
  // single-site test does (there's no simple symmetry here to predict
  // one from the geometry alone) -- just that both sites feel a
  // genuinely nonzero field, confirming the setup itself produces real
  // physics rather than a near-cancellation before trusting anything
  // about the operator. (An earlier version of this test placed the
  // external charge at exactly half the box length along one axis from
  // the molecule -- a special symmetric position where the reciprocal-
  // space field can cancel to near machine-epsilon rather than a small
  // genuine value; confirmed directly by reproducing that geometry
  // standalone. The specific coordinates here were chosen only to avoid
  // that kind of special-case placement, not for any other reason.)
  BOOST_REQUIRE(V_a.norm() > 1e-6);
  BOOST_REQUIRE(V_b.norm() > 1e-6);

  Eigen::VectorXd b(6);
  b.segment<3>(0) = V_a;
  b.segment<3>(3) = V_b;

  // Independent reference solve: build the true 6x6 system directly.
  // getPInv() (D) and EwaldRealSpaceInteractor::ComputeB (the erfc-only,
  // undamped-by-Thole intramolecular tensor -- matching
  // AddIntraSegmentCoupling's own choice, NOT eeInteractor::
  // FillTholeInteraction; see that method's own documentation for why)
  // are used directly, exactly as before. The reciprocal-space
  // contribution is different: rather than re-deriving its tensor by
  // hand (EwaldReciprocalSpaceSum's own structure factor couples every
  // registered site to every other, including same-position/self pairs,
  // via a genuinely global k-space sum -- reimplementing that
  // independently here would mean re-deriving the very code this test
  // exists to check, not a genuine cross-check on it), this PROBES
  // recip_sum's own actual output directly: a unit trial dipole is
  // placed on one site at a time (the other left at zero), fields are
  // read off at both sites via recip_sum's own public AddFieldAtMany,
  // and the (both-zero) baseline is subtracted to isolate that one
  // site's own linear contribution (field is linear in dipole moment,
  // so this subtraction is exact, not an approximation). This still
  // uses recip_sum's own real implementation for the numbers themselves
  // -- what's independent is the assembly into a dense linear system
  // and its solution, not the underlying field values.
  auto probe_recip = [&](PolarSite& source, int col,
                         Eigen::Matrix3d& block_at_a,
                         Eigen::Matrix3d& block_at_b) {
    mseg[0].setInduced_Dipole(Eigen::Vector3d::Zero());
    mseg[1].setInduced_Dipole(Eigen::Vector3d::Zero());
    Eigen::Vector3d unit = Eigen::Vector3d::Zero();
    unit(col) = 1.0;
    source.setInduced_Dipole(unit);
    mseg[0].Reset();
    mseg[1].Reset();
    std::vector<PolarSite*> both = {&mseg[0], &mseg[1]};
    recip_sum.AddFieldAtMany<Estatic::V>(both, EwaldChargeState::Neutral);
    block_at_a.col(col) = mseg[0].V();
    block_at_b.col(col) = mseg[1].V();
  };
  Eigen::Matrix3d baseline_a = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d baseline_b = Eigen::Matrix3d::Zero();
  {
    mseg[0].setInduced_Dipole(Eigen::Vector3d::Zero());
    mseg[1].setInduced_Dipole(Eigen::Vector3d::Zero());
    mseg[0].Reset();
    mseg[1].Reset();
    std::vector<PolarSite*> both = {&mseg[0], &mseg[1]};
    recip_sum.AddFieldAtMany<Estatic::V>(both, EwaldChargeState::Neutral);
    baseline_a.colwise() = mseg[0].V();
    baseline_b.colwise() = mseg[1].V();
  }
  // recip_a_a: field at a from a's own trial dipole (the self-field
  // term, same mechanism as SelfFieldMatrix() but probed directly here
  // rather than assuming the two agree). recip_b_a: field at b from a's
  // trial dipole. Etc.
  Eigen::Matrix3d recip_a_a = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d recip_b_a = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d recip_a_b = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d recip_b_b = Eigen::Matrix3d::Zero();
  for (int col = 0; col < 3; ++col) {
    probe_recip(mseg[0], col, recip_a_a, recip_b_a);
    probe_recip(mseg[1], col, recip_a_b, recip_b_b);
  }
  recip_a_a -= baseline_a;
  recip_b_a -= baseline_b;
  recip_a_b -= baseline_a;
  recip_b_b -= baseline_b;
  mseg[0].setInduced_Dipole(Eigen::Vector3d::Zero());
  mseg[1].setInduced_Dipole(Eigen::Vector3d::Zero());
  mseg[0].Reset();
  mseg[1].Reset();

  EwaldRealSpaceInteractor intra_interactor(alpha_ewald, thole_a);
  Eigen::Vector3d r_vec = mseg[0].getPos() - mseg[1].getPos();
  const double r_ab = r_vec.norm();
  EwaldRealSpaceInteractor::BFunctions bf = intra_interactor.ComputeB(r_ab);
  // AddIntraSegmentCoupling now applies ComputeThole (see this operator's
  // own class documentation for why an earlier, UNDAMPED version was
  // wrong) -- this reference must therefore include the same l3/l5
  // damping to stay a genuine cross-check on the operator's own CURRENT
  // behavior, not the older undamped one. site_a/site_b share the same
  // polarizability here, so ComputeThole's own two-site-order dependence
  // doesn't matter for this particular test, but the call is written the
  // same way AddIntraSegmentCoupling itself makes it (source, target)
  // for consistency.
  EwaldRealSpaceInteractor::TholeFactors thole =
      intra_interactor.ComputeThole(r_ab, mseg[1], mseg[0]);
  // NOTE on what this reference can and cannot catch: it is built from
  // the same ComputeB/ComputeErfB/ComputeThole primitives the operator
  // itself uses, and combines them the same way, so it cross-checks the
  // operator's ASSEMBLY (offsets, signs, transposes, which blocks get
  // which term) but is blind by construction to the CONVENTION being
  // assembled. It passed unchanged through both earlier, wrong versions
  // of AddIntraSegmentCoupling, and it would pass again if the l-versus-
  // erf combination below were wrong. The property that does discriminate
  // is alpha-independence, which has a test of its own further down.
  EwaldRealSpaceInteractor::BFunctions bferf =
      intra_interactor.ComputeErfB(r_ab);
  const double c3_ref = thole.l3 * bf.B1 + (thole.l3 - 1.0) * bferf.B1;
  const double c5_ref = thole.l5 * bf.B2 + (thole.l5 - 1.0) * bferf.B2;
  Eigen::Matrix3d erfc_block = c5_ref * (r_vec * r_vec.transpose()) -
                               c3_ref * Eigen::Matrix3d::Identity();

  // multiply(v) = D*v + field(v), where field(v) is the TRUE field a
  // trial dipole v produces (whatever its own actual sign is), so the
  // operator this represents is simply A = D + field-per-unit-v -- NOT
  // D minus that field. This differs from how
  // test_ewaldperiodicdipoleoperator.cc's own single-site case looks at
  // first glance (that test computes A = D - alpha*SelfFieldMatrix()),
  // but the two are consistent, not contradictory: SelfFieldMatrix() is
  // deliberately DEFINED as a positive matrix (see its own
  // documentation), with the minus sign applied explicitly by that
  // test's own caller to get the true (negative) self-field. Here,
  // recip_a_a/recip_a_b/etc. are PROBED directly from recip_sum's own
  // real output -- they already ARE the true field-per-unit-v, sign
  // included (printed as negative, e.g. recip_a_a's own diagonal), not
  // a magnitude requiring a separate sign to be applied by this test.
  // An earlier version of this test applied a minus sign here anyway,
  // by analogy with the single-site case, without checking whether the
  // probed values already had their own true sign built in -- they did,
  // so that was a real double-negation, not a genuine ambiguity: caught
  // by comparing against the actual PCG output at several different
  // sign combinations, only recip=+1/erfc=+1 matched (to 1e-8, machine
  // precision) while every other combination was off by a real margin
  // (1e-3 to 1e-2). The same reasoning applies to the intramolecular
  // erfc_block term: AddIntraSegmentCoupling adds +block*v_j directly
  // into multiply(v)'s own result, so it belongs in A with a plain
  // plus sign too.
  //
  // EwaldPeriodicDipoleOperator::RawMultiply now also adds
  // +self_field_matrix_*v_i into its own diagonal (see that class's own
  // documentation): the reciprocal-space self-field leak captured above
  // in recip_a_a/recip_b_b is now cancelled by the real operator itself,
  // matching legacy's own atomic self-interaction correction. An
  // independent reference must therefore add the same correction here to
  // stay a genuine cross-check on the operator's own CURRENT behavior --
  // this is deliberately recip_sum.SelfFieldMatrix() itself (the same
  // building block RawMultiply's own fix actually uses), not simply
  // -recip_a_a/-recip_b_b (which should equal it, being
  // position-independent and identical for every site, but computing it
  // independently here keeps this a genuine check on that equivalence
  // rather than assuming it).
  Eigen::Matrix3d self_field = recip_sum.SelfFieldMatrix();
  // Shape/surface correction: EwaldShapeCorrection::TotalDipoleMoment
  // sums charge + static dipole + induced dipole over EVERY registered
  // segment (see this operator's own class documentation for why this
  // was missing entirely before). Neither site_a nor site_b has a charge
  // or static dipole (this test's molecule was constructed with only
  // setpolarization called on each, no setMultipole), and the external
  // fixed charge contributes only a v-INDEPENDENT term (folded into b,
  // via EwaldBackground's own permanent-field-generation pattern this
  // test mirrors for V_a/V_b above) -- so the only v-dependent
  // contribution to M here is v_a+v_b themselves. Since M affects every
  // site identically (a genuine mean-field term, not pairwise), the same
  // 3x3 coefficient block belongs in all four blocks of A_dense.
  const double shape_coeff = -4.0 * M_PI / (3.0 * volume);
  Eigen::Matrix3d shape_block = shape_coeff * Eigen::Matrix3d::Identity();
  // The operator is A = P^-1 - C, where C is the whole induced-coupling
  // block (real-space intramolecular, reciprocal, shape, and the
  // self-field correction). RawMultiply builds exactly that: P^-1*v
  // minus the accumulated coupling field, minus self_field_matrix_*v.
  // This reference therefore SUBTRACTS every coupling term from the
  // P^-1 diagonal. An earlier version of this test added them, matching
  // an earlier operator that had the coupling sign inverted -- the same
  // sign error that made mu2 disagree with legacy while mu1 matched,
  // since the first iteration never exercises the coupling block at all.
  Eigen::Matrix<double, 6, 6> A_dense = Eigen::Matrix<double, 6, 6>::Zero();
  A_dense.block<3, 3>(0, 0) =
      mseg[0].getPInv() - recip_a_a - self_field - shape_block;
  A_dense.block<3, 3>(3, 3) =
      mseg[1].getPInv() - recip_b_b - self_field - shape_block;
  A_dense.block<3, 3>(0, 3) = -erfc_block - recip_a_b - shape_block;
  A_dense.block<3, 3>(3, 0) = -erfc_block.transpose() - recip_b_a - shape_block;
  Eigen::VectorXd x_reference = A_dense.fullPivLu().solve(b);

  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape,
                                 std::vector<Index>{2}, alpha_ewald, thole_a);

  Eigen::ConjugateGradient<EwaldPeriodicDipoleOperator,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setMaxIterations(200);
  cg.setTolerance(1e-10);
  cg.compute(op);
  Eigen::VectorXd x = cg.solveWithGuess(b, Eigen::VectorXd::Zero(6));

  BOOST_REQUIRE(cg.info() == Eigen::ComputationInfo::Success);

  // Diagnostic output, unconditional (not just on failure), so the
  // actual numbers are visible directly in the test's own stdout
  // without needing to rerun with a different verbosity flag.
  std::cout << "V_a = " << V_a.transpose() << std::endl;
  std::cout << "V_b = " << V_b.transpose() << std::endl;
  std::cout << "recip_a_a =\n" << recip_a_a << std::endl;
  std::cout << "recip_a_b =\n" << recip_a_b << std::endl;
  std::cout << "recip_b_a =\n" << recip_b_a << std::endl;
  std::cout << "recip_b_b =\n" << recip_b_b << std::endl;
  std::cout << "erfc_block =\n" << erfc_block << std::endl;
  std::cout << "A_dense =\n" << A_dense << std::endl;
  std::cout << "b (rhs) = " << b.transpose() << std::endl;
  std::cout << "x (operator/PCG)  = " << x.transpose() << std::endl;
  std::cout << "x_reference (dense) = " << x_reference.transpose() << std::endl;
  std::cout << "difference norm = " << (x - x_reference).norm() << std::endl;
  std::cout << "cg iterations = " << cg.iterations()
            << "  cg error = " << cg.error() << std::endl;

  BOOST_CHECK(x.isApprox(x_reference, 1e-6));
  // Sanity: not vacuously passing because both are ~zero.
  BOOST_CHECK(x.norm() > 1e-6);
}

// The discriminating test for AddIntraSegmentCoupling's CONVENTION, as
// opposed to its assembly. alpha is a free parameter of the Ewald split
// and the converged dipoles cannot depend on it; the test above cannot
// see a violation because its reference is built from the same
// primitives, combined the same way, so it agrees with whatever the
// operator does.
//
// This is not a hypothetical. Before the (l-1)*B_erf term was added,
// an intramolecular pair's total induced coupling was l*B_erfc + B_erf
// rather than l*B_bare -- see AddIntraSegmentCoupling's own comment --
// and in an 18-segment MM/MM run that drifted the [fg permanent x bg
// induced] channel by 0.43% across alpha = 1.5 ... 3.0 1/nm, visible
// even in that channel's shape term, which contains no alpha at all.
//
// Geometry, polarizabilities and thole_a are deliberately identical to
// the test above: 2.0 bohr is a bonded separation where Thole damping
// is genuinely active (l < 1), which is the whole point -- at a
// separation where l = 1 the old and new expressions coincide and this
// test would pass either way. If this test is ever seen to pass against
// a build WITHOUT the (l-1)*B_erf term, that means l has drifted to 1
// for this geometry and the test has gone vacuous; the damping printout
// below is there to make that visible rather than silent.
BOOST_AUTO_TEST_CASE(converged_dipoles_are_independent_of_the_splitting) {
  Eigen::Matrix3d box = 40.0 * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;

  PolarSegment external("seg", 1);
  PolarSite charge(1, "H", Eigen::Vector3d(15.3, 4.7, -8.2));
  charge.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 1.0;
  charge.setMultipole(mpoles, 0);
  external.push_back(charge);
  registry.Register(1, EwaldChargeState::Neutral, external);

  PolarSegment molecule("seg", 2);
  const double alpha_pol = 2.0;
  const double thole_a = 0.39;
  PolarSite site_a(1, "H", Eigen::Vector3d(3.4, -1.9, 0.7));
  site_a.setpolarization(alpha_pol * Eigen::Matrix3d::Identity());
  molecule.push_back(site_a);
  PolarSite site_b(2, "H", Eigen::Vector3d(5.4, -1.9, 0.7));
  site_b.setpolarization(alpha_pol * Eigen::Matrix3d::Identity());
  molecule.push_back(site_b);
  registry.Register(2, EwaldChargeState::Neutral, molecule);

  // Report the damping actually in force, so a vacuous pass (l -> 1) is
  // visible in the test's own output rather than silently harmless.
  {
    EwaldRealSpaceInteractor probe(0.30, thole_a);
    PolarSegment& m = registry.Get(2, EwaldChargeState::Neutral);
    const double r_ab = (m[0].getPos() - m[1].getPos()).norm();
    EwaldRealSpaceInteractor::TholeFactors probe_t =
        probe.ComputeThole(r_ab, m[1], m[0]);
    std::cout << "intramolecular r = " << r_ab
              << " bohr, thole l3 = " << probe_t.l3 << ", l5 = " << probe_t.l5
              << std::endl;
    // Guard against the test quietly going vacuous.
    BOOST_REQUIRE(probe_t.l3 < 0.999);
  }

  // r_min is held fixed while alpha only INCREASES, so the real-space
  // truncation error erfc(alpha*r_min) shrinks monotonically across the
  // scan and cannot itself masquerade as alpha-dependence. k_max is
  // likewise already far past convergence at the largest alpha here.
  const std::vector<double> alphas = {0.30, 0.40, 0.50};

  std::vector<Eigen::VectorXd> solutions;
  std::vector<Eigen::VectorXd> rhs;

  for (double alpha_ewald : alphas) {
    PolarSegment& mseg = registry.Get(2, EwaldChargeState::Neutral);
    mseg[0].Reset();
    mseg[1].Reset();
    mseg[0].setInduced_Dipole(Eigen::Vector3d::Zero());
    mseg[1].setInduced_Dipole(Eigen::Vector3d::Zero());

    EwaldRealSpaceSum real_sum(box, registry, alpha_ewald, thole_a,
                               /*r_min=*/12.0, /*field_tol=*/1e-12);
    // k_max = 12*alpha, the same convergence rule test_ewaldregion.cc
    // already uses (params.k_max = 12.0 * alpha), rather than a flat
    // literal. The reciprocal sum's Gaussian weight is
    // exp(-k^2/(4*alpha^2)), so the largest neglected term is
    // exp(-(12*alpha)^2/(4*alpha^2)) = exp(-36) = 2.3e-16 at EVERY
    // alpha -- eight orders of magnitude below this file's tightest
    // assertion (1e-8) and independent of alpha, which is what an
    // alpha-independence test needs: a truncation error that tracks
    // alpha would imitate the very failure being looked for.
    //
    // A flat k_max = 20.0 was converged too, just absurdly so, and the
    // cost is cubic in k_max: in this 40 bohr box the cutoff sphere
    // holds 8,647,082 k-vectors at k_max = 20, against 50,540 /
    // 119,512 / 233,576 at 12*alpha for alpha = 0.30 / 0.40 / 0.50.
    // That is what made unit_test_ewaldperiodicdipoleoperatormultisite
    // take 1452 sec of CTest's 1500 sec default in a Debug CI run --
    // passing with 48 sec to spare, i.e. one slower runner away from a
    // spurious timeout.
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_ewald,
                                      /*k_max=*/12.0 * alpha_ewald);
    const double volume = box.determinant();
    EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

    // b = permanent field. Only the external segment carries static
    // multipoles, so this is a clean real+reciprocal split of an
    // INTERmolecular interaction and is alpha-independent in its own
    // right -- checked explicitly below, so that a failure of the
    // dipoles can be attributed to the coupling operator rather than to
    // a drifting right-hand side.
    Eigen::VectorXd b(6);
    {
      PolarSite& a = mseg[0];
      a.Reset();
      real_sum.AddFieldAt<Estatic::V>(2, a, EwaldChargeState::Neutral);
      recip_sum.AddFieldAt<Estatic::V>(a, EwaldChargeState::Neutral);
      b.segment<3>(0) = a.V();
      a.Reset();

      PolarSite& bb = mseg[1];
      bb.Reset();
      real_sum.AddFieldAt<Estatic::V>(2, bb, EwaldChargeState::Neutral);
      recip_sum.AddFieldAt<Estatic::V>(bb, EwaldChargeState::Neutral);
      b.segment<3>(3) = bb.V();
      bb.Reset();
    }

    EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape,
                                   std::vector<Index>{2}, alpha_ewald, thole_a);

    Eigen::ConjugateGradient<EwaldPeriodicDipoleOperator,
                             Eigen::Lower | Eigen::Upper,
                             Eigen::DiagonalPreconditioner<double>>
        cg;
    cg.setMaxIterations(200);
    cg.setTolerance(1e-12);
    cg.compute(op);
    Eigen::VectorXd x = cg.solveWithGuess(b, Eigen::VectorXd::Zero(6));
    BOOST_REQUIRE(cg.info() == Eigen::ComputationInfo::Success);

    std::cout << "alpha = " << alpha_ewald << "  b = " << b.transpose()
              << "  mu = " << x.transpose() << std::endl;

    rhs.push_back(b);
    solutions.push_back(x);
  }

  // Sanity: a real, nonzero induction problem, not a vacuous one.
  BOOST_REQUIRE(solutions.front().norm() > 1e-6);

  for (std::size_t i = 1; i < solutions.size(); ++i) {
    const double db = (rhs[i] - rhs[0]).norm() / rhs[0].norm();
    const double dmu =
        (solutions[i] - solutions[0]).norm() / solutions[0].norm();
    std::cout << "alpha = " << alphas[i] << "  relative drift: rhs = " << db
              << "  mu = " << dmu << std::endl;
    // The right-hand side is a plain Ewald split and should be tight.
    BOOST_CHECK_SMALL(db, 1e-8);
    // The dipoles carry the coupling operator. The pre-fix code failed
    // this at the 1e-3 level; 1e-7 leaves room for the real-space
    // truncation at r_min = 12 bohr without leaving room for a
    // convention error.
    BOOST_CHECK_SMALL(dmu, 1e-7);
  }
}

BOOST_AUTO_TEST_SUITE_END()
