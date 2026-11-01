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

// Tests EwaldPeriodicDipoleOperator for segments with more than one site:
// the variable-width vector layout and the intramolecular Thole-damped
// induced-dipole coupling (AddIntraSegmentCoupling), which must be
// present, as in legacy PolarBackground (FU12_ERFC_At_By).
//
// The PCG solve for one 2-site segment is compared with a dense direct
// solve of the same equation, assembled from getPInv(),
// EwaldRealSpaceInteractor::ComputeB/ComputeErfB/ComputeThole, and
// probed reciprocal-space fields, without calling the operator itself.
BOOST_AUTO_TEST_SUITE(ewaldperiodicdipoleoperatormultisite_test)

BOOST_AUTO_TEST_CASE(intramolecular_coupling_matches_direct_solve) {
  Eigen::Matrix3d box = 40.0 * Eigen::Matrix3d::Identity();

  EwaldRegistry registry;

  // Fixed external background: one positive charge, placed off any
  // symmetric position relative to the molecule. At e.g. half a box
  // length the reciprocal-space field cancels to ~1e-19, and the
  // comparison would be dominated by round-off.
  PolarSegment external("seg", 1);
  PolarSite charge(1, "H", Eigen::Vector3d(15.3, 4.7, -8.2));
  charge.setpolarization(Eigen::Matrix3d::Identity());
  Vector9d mpoles = Vector9d::Zero();
  mpoles(0) = 1.0;
  charge.setMultipole(mpoles, 0);
  external.push_back(charge);
  registry.Register(1, EwaldChargeState::Neutral, external);

  // The segment solved for: two polarizable sites at a bonded distance
  // (2.0 bohr, ~C-H), where the intramolecular coupling is substantial
  // but the system is still positive definite. At 0.5 bohr the
  // erfc-screened coupling dominates 1/alpha_pol and the system becomes
  // indefinite (polarization catastrophe), which this test is not meant
  // to probe.
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
  // k_max = 12*alpha, see converged_dipoles_are_independent_of_the_splitting.
  EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_ewald,
                                    /*k_max=*/12.0 * alpha_ewald);
  const double volume = box.determinant();
  EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

  // Permanent field at each site; this is b = +V_permanent for the solve
  // (see EwaldPeriodicDipoleOperator for the sign).
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

  // Sanity: both sites feel a nonzero field (no symmetric cancellation).
  BOOST_REQUIRE(V_a.norm() > 1e-6);
  BOOST_REQUIRE(V_b.norm() > 1e-6);

  Eigen::VectorXd b(6);
  b.segment<3>(0) = V_a;
  b.segment<3>(3) = V_b;

  // Reference solve: assemble the 6x6 system directly. The
  // reciprocal-space blocks are not re-derived (that would duplicate the
  // code under test) but probed from recip_sum: a unit trial dipole is
  // placed on one site at a time, the fields at both sites are read via
  // AddFieldAtMany, and the zero-dipole baseline is subtracted (exact,
  // since the field is linear in the dipoles). The independent part is
  // the assembly and solution of the dense system.
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
  // recip_a_a: field at a from a's trial dipole (includes the
  // reciprocal-space self-field); recip_b_a: field at b from a's trial
  // dipole; etc.
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
  // Thole damping as in AddIntraSegmentCoupling, with the same (source,
  // target) argument order.
  EwaldRealSpaceInteractor::TholeFactors thole =
      intra_interactor.ComputeThole(r_ab, mseg[1], mseg[0]);
  // This reference uses the same primitives, combined the same way, as
  // the operator. It therefore checks the ASSEMBLY (offsets, signs,
  // transposes) but not the CONVENTION (the l versus erf combination);
  // that is checked by the alpha-independence test below.
  EwaldRealSpaceInteractor::BFunctions bferf =
      intra_interactor.ComputeErfB(r_ab);
  const double c3_ref = thole.l3 * bf.B1 + (thole.l3 - 1.0) * bferf.B1;
  const double c5_ref = thole.l5 * bf.B2 + (thole.l5 - 1.0) * bferf.B2;
  Eigen::Matrix3d erfc_block = c5_ref * (r_vec * r_vec.transpose()) -
                               c3_ref * Eigen::Matrix3d::Identity();

  // The probed recip_* blocks are the field per unit dipole with its
  // sign included. RawMultiply also cancels the reciprocal-space
  // self-field with SelfFieldMatrix(); the reference uses that same
  // matrix rather than -recip_a_a.
  Eigen::Matrix3d self_field = recip_sum.SelfFieldMatrix();
  // Shape correction: the only v-dependent part of the total dipole M is
  // v_a + v_b (the molecule has no permanent moments; the external
  // charge's part is v-independent and cancels in the operator). The
  // field -(4 pi/3V) M is the same at every site, so the same block
  // enters all four blocks of A_dense.
  const double shape_coeff = -4.0 * M_PI / (3.0 * volume);
  Eigen::Matrix3d shape_block = shape_coeff * Eigen::Matrix3d::Identity();
  // A = P^-1 - C, with C the full induced-coupling block (intramolecular
  // real space, reciprocal space, shape and self-field correction), so
  // every coupling term is subtracted.
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

  // Diagnostic output, printed unconditionally.
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

// Tests the CONVENTION of AddIntraSegmentCoupling (the test above only
// checks its assembly): the converged dipoles must not depend on the
// Ewald splitting parameter alpha. With l*B_erfc instead of
// l*B_erfc + (l-1)*B_erf, an 18-segment MM/MM run showed a 0.43% drift
// of the [fg permanent x bg induced] channel over alpha = 1.5 ... 3.0 1/nm.
//
// Same geometry as above: at 2.0 bohr Thole damping is active (l < 1).
// Where l = 1 both expressions coincide and the test would be vacuous;
// the printed damping factors make that visible.
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

  // Report the damping factors, so that a vacuous pass (l -> 1) is
  // visible.
  {
    EwaldRealSpaceInteractor probe(0.30, thole_a);
    PolarSegment& m = registry.Get(2, EwaldChargeState::Neutral);
    const double r_ab = (m[0].getPos() - m[1].getPos()).norm();
    EwaldRealSpaceInteractor::TholeFactors probe_t =
        probe.ComputeThole(r_ab, m[1], m[0]);
    std::cout << "intramolecular r = " << r_ab
              << " bohr, thole l3 = " << probe_t.l3 << ", l5 = " << probe_t.l5
              << std::endl;
    // Guard against a vacuous test.
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
    // k_max = 12*alpha (as in test_ewaldregion.cc): the largest neglected
    // reciprocal term is exp(-(12 alpha)^2/(4 alpha^2)) = exp(-36) =
    // 2.3e-16 for every alpha, so the truncation error does not itself
    // depend on alpha. A fixed k_max = 20 needs 8,647,082 k-vectors in
    // this 40 bohr box (vs 50,540 to 233,576 here) and brought the Debug
    // CI runtime to 1452 s of the 1500 s CTest limit.
    EwaldReciprocalSpaceSum recip_sum(box, registry, alpha_ewald,
                                      /*k_max=*/12.0 * alpha_ewald);
    const double volume = box.determinant();
    EwaldShapeCorrection shape(volume, registry, EwaldShape::Cube);

    // b = permanent field. Only the external segment has permanent
    // moments, so b is a purely intermolecular Ewald sum and must itself
    // be alpha-independent; checked below, so that a failure of the
    // dipoles can be attributed to the coupling operator.
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

  // Sanity: a nonzero induction problem.
  BOOST_REQUIRE(solutions.front().norm() > 1e-6);

  for (std::size_t i = 1; i < solutions.size(); ++i) {
    const double db = (rhs[i] - rhs[0]).norm() / rhs[0].norm();
    const double dmu =
        (solutions[i] - solutions[0]).norm() / solutions[0].norm();
    std::cout << "alpha = " << alphas[i] << "  relative drift: rhs = " << db
              << "  mu = " << dmu << std::endl;
    // The right-hand side is a plain Ewald split and should be tight.
    BOOST_CHECK_SMALL(db, 1e-8);
    // The dipoles test the coupling operator. The wrong convention fails
    // at the 1e-3 level; 1e-7 allows for the real-space truncation at
    // r_min = 12 bohr.
    BOOST_CHECK_SMALL(dmu, 1e-7);
  }
}

BOOST_AUTO_TEST_SUITE_END()
