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

#define BOOST_TEST_MODULE ewaldperiodicdipoleoperator_test

// Standard includes
#include <cmath>

// Third party includes
#include <Eigen/IterativeLinearSolvers>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/ewaldperiodicdipoleoperator.h"

using namespace votca::xtp;
using namespace votca;

// This test deliberately avoids hand-deriving the sign convention that
// connects EwaldPeriodicDipoleOperator's own A*x=b equation back to the
// physical induced dipole -- doing that by hand is exactly the kind of
// manual derivation that has introduced real bugs elsewhere in this
// project (see the factor-of-4 bug caught by
// test_ewaldrealspaceinteractor.cc's own dipole-dipole test). Instead,
// this checks a case simple enough that physical intuition alone gives
// an unambiguous, hard-to-get-wrong QUALITATIVE expectation: an
// isotropically polarizable site sitting in the field of one fixed
// external positive charge must develop an induced dipole moment
// *aligned with* that external field, not opposed to it. With a single
// site in `ids`, there is no periodic self-interaction between `ids`
// members to also get right -- this isolates the self-term/RHS sign
// convention specifically.
//
// The QUANTITATIVE check below is a separate matter: an earlier version
// of this test expected the naive mu=alpha*E exactly, which held when
// this was written but is no longer quite right (see the check itself
// for why -- EwaldReciprocalSpaceSum's own reciprocal-space self-field
// term, once its own class documentation is understood, means even this
// single-site case has a genuine, computable correction on top of the
// naive expectation). Getting the qualitative check right first, before
// trusting a hand-derived quantitative one, is exactly the discipline
// this test was originally written to encourage -- worth keeping in mind
// if this quantitative check ever needs revisiting again.
BOOST_AUTO_TEST_SUITE(ewaldperiodicdipoleoperator_test)

BOOST_AUTO_TEST_CASE(induced_dipole_aligns_with_external_field) {
  Eigen::Matrix3d box = 40.0 * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;

  // Fixed external background: a single positive charge, not part of
  // `ids`, so it is never itself solved for -- exactly the "fixed,
  // non-reacting background" role described in the class documentation.
  PolarSegment external("seg", 1);
  PolarSite charge(1, "H", Eigen::Vector3d::Zero());
  charge.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 1.0;
  charge.setMultipole(mpoles, 0);
  external.push_back(charge);
  registry.Register(1, EwaldChargeState::Neutral, external);

  // The one polarizable site being solved for.
  PolarSegment polarizable("seg", 2);
  const double alpha_pol = 2.0;
  PolarSite site(2, "H", Eigen::Vector3d(3.0, 0.0, 0.0));
  site.setpolarization(alpha_pol * Eigen::Matrix3d::Identity());
  polarizable.push_back(site);
  registry.Register(2, EwaldChargeState::Neutral, polarizable);

  const double alpha_ewald = 0.3;
  EwaldRealSpaceSum real_sum(box, registry, alpha_ewald, /*thole_a=*/0.39,
                             /*r_min=*/12.0, /*field_tol=*/1e-12);
  // 12*alpha, the convergence rule used across the Ewald tests: the
  // largest neglected Gaussian weight is exp(-36) = 2.3e-16 at any
  // alpha. In this 40 bohr box that is 50,540 k-vectors instead of the
  // 8,647,082 a flat k_max = 20.0 generates, for the same answer.
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_ewald,
                                    /*k_max=*/12.0 * alpha_ewald);
  const double volume = box.determinant();
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  // Right-hand side: the permanent field from the fixed external charge
  // at the polarizable site. IMPORTANT: this is b = +V_permanent, not
  // b = -V as CalcInducedDipolesViaPCG itself uses -- see the class
  // documentation on why the sign differs from that legacy convention
  // (V()'s own sign convention differs between eeInteractor and
  // EwaldRealSpaceSum/EwaldReciprocalSpaceSum).
  PolarSegment& poseg = registry.Get(2, EwaldChargeState::Neutral);
  PolarSite& target = poseg[0];
  target.Reset();
  real_sum.AddFieldAt<Estatic::V>(2, target, EwaldChargeState::Neutral);
  recip_sum.AddFieldAt<Estatic::V>(target, EwaldChargeState::Neutral);
  Eigen::Vector3d V_permanent = target.V();
  target.Reset();

  // Sanity: the field from a positive charge at the origin, felt at
  // (3,0,0), must point in +x (away from the source) -- confirmed by
  // test_ewaldrealspaceinteractor.cc's own already-validated sign
  // convention. If this fails, the setup itself is wrong, independent of
  // anything about EwaldPeriodicDipoleOperator.
  BOOST_REQUIRE(V_permanent.x() > 0.0);

  Eigen::VectorXd b(3);
  b = V_permanent;

  EwaldPeriodicDipoleOperator op(registry, real_sum, recip_sum, shape,
                                 std::vector<Index>{2}, alpha_ewald,
                                 /*thole_a=*/0.39);

  Eigen::ConjugateGradient<EwaldPeriodicDipoleOperator,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg;
  cg.setMaxIterations(200);
  cg.setTolerance(1e-10);
  cg.compute(op);
  Eigen::VectorXd x = cg.solveWithGuess(b, Eigen::VectorXd::Zero(3));

  BOOST_REQUIRE(cg.info() == Eigen::ComputationInfo::Success);

  // The physical expectation, independent of any sign convention inside
  // the operator's own A*x=b equation: mu = alpha*E, aligned with the
  // external field, i.e. also pointing in +x here.
  BOOST_CHECK(x(0) > 0.0);
  BOOST_CHECK_SMALL(x(1), 1e-8);
  BOOST_CHECK_SMALL(x(2), 1e-8);

  // With no other `ids` member to interact with intermolecularly, this
  // reduces to a single-site self-consistent equation. Two artifacts
  // that would otherwise complicate this are already accounted for:
  //
  //   - EwaldReciprocalSpaceSum's own structure factor never excludes
  //     any site, not even the one a field is being computed at (see
  //     that class's own class documentation), so this site's own trial
  //     dipole would otherwise feed back into its own field via a
  //     spurious reciprocal-space "self-field" term, -M*mu (M =
  //     EwaldReciprocalSpaceSum::SelfFieldMatrix()). That term is
  //     cancelled directly inside EwaldPeriodicDipoleOperator::
  //     RawMultiply itself (see its own self_field_matrix_ member).
  //
  //   - EwaldShapeCorrection's own mean-field surface term is now also
  //     applied inside RawMultiply (see this operator's own class
  //     documentation for why this was missing entirely before, and how
  //     it was found). Unlike the self-field term above, this one does
  //     NOT cancel to zero here: EwaldShapeCorrection::TotalDipoleMoment
  //     sums every registered segment's charge + static dipole + induced
  //     dipole, so even in this single-site case it picks up the fixed
  //     external charge's own q*r (folded into V_permanent's own value
  //     the same way as any other permanent-field contribution, so it
  //     does not appear separately here) AND this site's own induced
  //     dipole mu -- a genuine coefficient
  //     shape_coeff = -4*pi/(3*volume) multiplying mu, added to the
  //     diagonal alongside D = getPInv() = (1/alpha_pol)*I.
  //
  // So for this single-site case the equation is simply
  // D*mu = V_permanent, with D = getPInv() = (1/alpha_pol)*I -- the
  // plain mu = alpha*E result after all, but for a non-obvious reason
  // worth recording, since this expectation has now been revised twice
  // in opposite directions.
  //
  // The three v-dependent terms the operator contributes here cancel
  // exactly:
  //
  //     reciprocal self-term (the site's own dipole, via the
  //       unconditional k-sum)
  //   + shape/surface term (shape_coeff * mu)
  //   + analytic self-field correction (subtracted by RawMultiply)
  //   = 0
  //
  // Measured on this geometry: -0.0202453751608 + (-6.54498469e-05)
  // = -0.0203108250077, which is exactly minus SelfFieldMatrix()'s own
  // +0.0203108250077. So A comes out as (1/alpha_pol)*I to machine
  // precision and no shape coefficient survives.
  //
  // That cancellation is a property of the ANALYTIC self-field term.
  // An earlier version of this test added shape_coeff to the diagonal,
  // which was right back when SelfFieldMatrix() returned the discrete
  // k-lattice sum: that sum is ~1.6% off the analytic value, so the
  // terms did not cancel and the remainder looked like a shape
  // coefficient. This assertion therefore doubles as a check on that
  // term -- if SelfFieldMatrix() ever reverts to a lattice sum, the
  // cancellation breaks and this fails.
  Eigen::Matrix3d A = Eigen::Matrix3d::Identity() / alpha_pol;
  Eigen::Vector3d expected = A.fullPivLu().solve(V_permanent);
  BOOST_CHECK(x.isApprox(expected, 1e-6));
  // Sanity: not vacuously passing because both are ~zero.
  BOOST_CHECK(x.norm() > 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
