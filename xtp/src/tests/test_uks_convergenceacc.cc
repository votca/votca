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
//
// No dedicated unit test file for UKSConvergenceAcc existed before this one
// -- everything it does was previously exercised only indirectly, via
// integration-style tests elsewhere (test_dftengine.cc,
// test_dftengine_cdft.cc, etc.) that run full SCF/CDFT calculations.
// Starting here with the cheapest, most directly testable pieces: the
// pure-math rotation-vector <-> antisymmetric-matrix conversion functions,
// which have no dependency on any Fock-builder callback, integral engine,
// or converged SCF state at all.
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE uks_convergenceacc_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/uks_convergenceacc.h"
#include "xtp_libint2.h"

#include <cmath>
#include <fstream>
#include <tuple>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(uks_convergenceacc_test)

BOOST_AUTO_TEST_CASE(unflatten_rotation_test) {
  UKSConvergenceAcc conv;

  Index nao = 5;
  Index nocclevels = 2;
  Index nvirt = nao - nocclevels;  // 3
  Eigen::VectorXd v_ov(nocclevels * nvirt);
  // Deliberately distinct values at every (i, a) position, so a mistake
  // in the index mapping (e.g. row-major vs column-major, or an
  // off-by-one in the occ/virt split) would show up as a value landing
  // in the wrong matrix element, not just a wrong-but-plausible one.
  for (Index k = 0; k < v_ov.size(); ++k) {
    v_ov(k) = 10.0 + double(k);
  }

  Eigen::MatrixXd kappa = conv.UnflattenRotation(v_ov, nao, nocclevels);

  BOOST_REQUIRE_EQUAL(kappa.rows(), nao);
  BOOST_REQUIRE_EQUAL(kappa.cols(), nao);

  // Occ-virt block matches v_ov exactly, at the expected (i, nocclevels+a)
  // positions, in row-major (i outer, a inner) order.
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      double expected = v_ov(i * nvirt + a);
      BOOST_CHECK_CLOSE(kappa(i, nocclevels + a), expected, 1e-10);
      // Antisymmetric: the virt-occ entry is the exact negative.
      BOOST_CHECK_CLOSE(kappa(nocclevels + a, i), -expected, 1e-10);
    }
  }

  // Occ-occ and virt-virt blocks are untouched (exactly zero) -- a
  // rotation between two occupied (or two virtual) orbitals is not part
  // of what this function is meant to construct at all.
  for (Index i = 0; i < nocclevels; ++i) {
    for (Index j = 0; j < nocclevels; ++j) {
      BOOST_CHECK_SMALL(kappa(i, j), 1e-12);
    }
  }
  for (Index a = 0; a < nvirt; ++a) {
    for (Index b = 0; b < nvirt; ++b) {
      BOOST_CHECK_SMALL(kappa(nocclevels + a, nocclevels + b), 1e-12);
    }
  }
}

BOOST_AUTO_TEST_CASE(unflatten_rotation_zero_test) {
  // Edge case: an all-zero input should give an all-zero, but still
  // correctly-SIZED, output -- not, say, an empty or mis-shaped matrix.
  UKSConvergenceAcc conv;
  Index nao = 4;
  Index nocclevels = 1;
  Eigen::VectorXd v_ov = Eigen::VectorXd::Zero(nocclevels * (nao - nocclevels));
  Eigen::MatrixXd kappa = conv.UnflattenRotation(v_ov, nao, nocclevels);
  BOOST_REQUIRE_EQUAL(kappa.rows(), nao);
  BOOST_REQUIRE_EQUAL(kappa.cols(), nao);
  BOOST_CHECK_SMALL(kappa.norm(), 1e-12);
}

BOOST_AUTO_TEST_CASE(unflatten_coupled_rotation_test) {
  // The coupled version splits ONE combined vector (alpha's own block
  // first, beta's immediately after) into two SEPARATE antisymmetric
  // matrices, one per channel, each with its own, potentially different
  // nao/nocclevels -- confirmed here by using different sizes for alpha
  // and beta, rather than two identical channels that could accidentally
  // pass even with a channel-swapping bug.
  UKSConvergenceAcc conv;

  Index nao_alpha = 5;
  Index nocclevels_alpha = 2;
  Index nvirt_alpha = nao_alpha - nocclevels_alpha;   // 3
  Index n_ov_alpha = nocclevels_alpha * nvirt_alpha;  // 6

  Index nao_beta = 4;
  Index nocclevels_beta = 1;
  Index nvirt_beta = nao_beta - nocclevels_beta;   // 3
  Index n_ov_beta = nocclevels_beta * nvirt_beta;  // 3

  Eigen::VectorXd v_ov_alpha(n_ov_alpha);
  for (Index k = 0; k < n_ov_alpha; ++k) {
    v_ov_alpha(k) = 1.0 + double(k);
  }
  Eigen::VectorXd v_ov_beta(n_ov_beta);
  for (Index k = 0; k < n_ov_beta; ++k) {
    // Deliberately a disjoint value range from alpha's own (100+ vs
    // 1-6), so a bug that mixed up which half of the combined vector
    // belongs to which channel would produce values clearly out of
    // the expected range for that channel, not just a subtly wrong
    // number within a plausible range.
    v_ov_beta(k) = 100.0 + double(k);
  }

  Eigen::VectorXd v_combined(n_ov_alpha + n_ov_beta);
  v_combined.head(n_ov_alpha) = v_ov_alpha;
  v_combined.tail(n_ov_beta) = v_ov_beta;

  auto [kappa_alpha, kappa_beta] = conv.UnflattenCoupledRotation(
      v_combined, nao_alpha, nocclevels_alpha, nao_beta, nocclevels_beta);

  BOOST_REQUIRE_EQUAL(kappa_alpha.rows(), nao_alpha);
  BOOST_REQUIRE_EQUAL(kappa_alpha.cols(), nao_alpha);
  BOOST_REQUIRE_EQUAL(kappa_beta.rows(), nao_beta);
  BOOST_REQUIRE_EQUAL(kappa_beta.cols(), nao_beta);

  // Each channel's own result must match calling UnflattenRotation
  // directly on that channel's own sub-vector -- the coupled version is
  // meant to be exactly equivalent to two separate calls, just sharing
  // one, concatenated input vector.
  Eigen::MatrixXd kappa_alpha_expected =
      conv.UnflattenRotation(v_ov_alpha, nao_alpha, nocclevels_alpha);
  Eigen::MatrixXd kappa_beta_expected =
      conv.UnflattenRotation(v_ov_beta, nao_beta, nocclevels_beta);

  BOOST_CHECK_SMALL((kappa_alpha - kappa_alpha_expected).norm(), 1e-10);
  BOOST_CHECK_SMALL((kappa_beta - kappa_beta_expected).norm(), 1e-10);
}

namespace {
// Toy unrestricted Hartree-Fock model on a real AO basis (methane,
// 3-21G, 17 functions): random one-electron h and a positive
// semidefinite, RI-like two-electron interaction
// (mn|ls) = sum_P B_P(m,n) B_P(l,s). It provides the energy and the
// coupled Fock builder, so that the second-order step can be checked
// against energies computed independently of UKSConvergenceAcc.
struct ToyUHF {
  AOOverlap overlap;
  Eigen::MatrixXd S;
  Eigen::MatrixXd h;
  std::vector<Eigen::MatrixXd> B;
  Index nocc_alpha = 5;
  Index nocc_beta = 4;

  Eigen::MatrixXd J(const Eigen::MatrixXd& D) const {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(D.rows(), D.cols());
    for (const Eigen::MatrixXd& b : B) {
      result += b.cwiseProduct(D).sum() * b;
    }
    return result;
  }
  Eigen::MatrixXd K(const Eigen::MatrixXd& D) const {
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(D.rows(), D.cols());
    for (const Eigen::MatrixXd& b : B) {
      result += b * D * b;
    }
    return result;
  }
  UKSConvergenceAcc::SpinFock Fock(const Eigen::MatrixXd& Da,
                                   const Eigen::MatrixXd& Db) const {
    Eigen::MatrixXd Jt = J(Da + Db);
    return {h + Jt - K(Da), h + Jt - K(Db)};
  }
  double Energy(const Eigen::MatrixXd& Da, const Eigen::MatrixXd& Db) const {
    UKSConvergenceAcc::SpinFock F = Fock(Da, Db);
    return 0.5 * (Da.cwiseProduct(h + F.alpha).sum() +
                  Db.cwiseProduct(h + F.beta).sum());
  }
  static Eigen::MatrixXd Density(const Eigen::MatrixXd& C, Index nocc) {
    return C.leftCols(nocc) * C.leftCols(nocc).transpose();
  }
  double EnergyOfOrbitals(const Eigen::MatrixXd& Ca,
                          const Eigen::MatrixXd& Cb) const {
    return Energy(Density(Ca, nocc_alpha), Density(Cb, nocc_beta));
  }
  // C (1 + kappa), Lowdin re-orthonormalised: the convention of
  // UKSConvergenceAcc.
  Eigen::MatrixXd Rotate(const Eigen::MatrixXd& C,
                         const Eigen::MatrixXd& kappa) const {
    Eigen::MatrixXd Cn =
        C * (Eigen::MatrixXd::Identity(C.cols(), C.cols()) + kappa);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cn.transpose() * S * Cn);
    return Cn * es.operatorInverseSqrt();
  }
  tools::EigenSystem Diagonalize(const Eigen::MatrixXd& F) const {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(F, S);
    tools::EigenSystem result;
    result.eigenvalues() = es.eigenvalues();
    result.eigenvectors() = es.eigenvectors();
    return result;
  }
};

ToyUHF MakeToyUHF() {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/aomatrix/molecule.xyz");
  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basisset, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  libint2::finalize();

  ToyUHF toy;
  toy.overlap = overlap;
  toy.S = overlap.Matrix();
  const Index n = toy.S.rows();
  std::srand(7);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  toy.h = 0.5 * (A + A.transpose()) - 2.0 * toy.S;
  for (Index P = 0; P < 10; ++P) {
    Eigen::MatrixXd b = 0.3 * Eigen::MatrixXd::Random(n, n);
    toy.B.push_back(0.5 * (b + b.transpose()));
  }
  return toy;
}

// Orbitals away from self-consistency: three Roothaan steps from the core
// guess, then a fixed random occupied-virtual rotation in each channel.
std::pair<tools::EigenSystem, tools::EigenSystem> StartingOrbitals(
    const ToyUHF& toy, const UKSConvergenceAcc& conv) {
  tools::EigenSystem Ma = toy.Diagonalize(toy.h);
  tools::EigenSystem Mb = Ma;
  for (int it = 0; it < 3; ++it) {
    UKSConvergenceAcc::SpinFock F =
        toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
                 ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
    Ma = toy.Diagonalize(F.alpha);
    Mb = toy.Diagonalize(F.beta);
  }
  const Index n = toy.S.rows();
  std::srand(11);
  Eigen::VectorXd va =
      0.05 * Eigen::VectorXd::Random(toy.nocc_alpha * (n - toy.nocc_alpha));
  Eigen::VectorXd vb =
      0.05 * Eigen::VectorXd::Random(toy.nocc_beta * (n - toy.nocc_beta));
  Ma.eigenvectors() = toy.Rotate(Ma.eigenvectors(),
                                 conv.UnflattenRotation(va, n, toy.nocc_alpha));
  Mb.eigenvectors() = toy.Rotate(Mb.eigenvectors(),
                                 conv.UnflattenRotation(vb, n, toy.nocc_beta));
  // orbital energies as UKSConvergenceAcc sees them after a rotation
  UKSConvergenceAcc::SpinFock F =
      toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
               ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
  Ma.eigenvalues() =
      (Ma.eigenvectors().transpose() * F.alpha * Ma.eigenvectors()).diagonal();
  Mb.eigenvalues() =
      (Mb.eigenvectors().transpose() * F.beta * Mb.eigenvectors()).diagonal();
  return {Ma, Mb};
}

// Energy along x (layout of UnflattenCoupledRotation) from (Ma, Mb).
double EnergyAlong(const ToyUHF& toy, const UKSConvergenceAcc& conv,
                   const tools::EigenSystem& Ma, const tools::EigenSystem& Mb,
                   const Eigen::VectorXd& x) {
  const Index n = toy.S.rows();
  auto [ka, kb] =
      conv.UnflattenCoupledRotation(x, n, toy.nocc_alpha, n, toy.nocc_beta);
  return toy.EnergyOfOrbitals(toy.Rotate(Ma.eigenvectors(), ka),
                              toy.Rotate(Mb.eigenvectors(), kb));
}

// Rotation actually applied: occupied-virtual block of C_old^T S C_new
// (equal to -x to first order, see UnflattenRotation).
double AppliedStepNorm(const ToyUHF& toy, const Eigen::MatrixXd& C_old,
                       const Eigen::MatrixXd& C_new, Index nocc) {
  const Index n = C_old.cols();
  Eigen::MatrixXd X = C_old.transpose() * toy.S * C_new;
  return X.block(nocc, 0, n - nocc, nocc).norm();
}
}  // namespace

BOOST_AUTO_TEST_CASE(energy_gradient_and_hessian_scale_test) {
  // dE/dx_ia = EnergyGradientScale * F_ia, and BuildCoupledSigmaVector is
  // the energy Hessian applied to v: both against finite differences of
  // the energy itself.
  ToyUHF toy = MakeToyUHF();
  Logger log;
  UKSConvergenceAcc conv;
  conv.setLogger(&log);
  conv.setOverlap(toy.overlap, 1e-8);
  auto [Ma, Mb] = StartingOrbitals(toy, conv);
  const Index n = toy.S.rows();
  const Index nva = n - toy.nocc_alpha;
  const Index nvb = n - toy.nocc_beta;
  const Index n_ov = toy.nocc_alpha * nva + toy.nocc_beta * nvb;

  UKSConvergenceAcc::SpinFock F =
      toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
               ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
  Eigen::MatrixXd Fa =
      Ma.eigenvectors().transpose() * F.alpha * Ma.eigenvectors();
  Eigen::MatrixXd Fb =
      Mb.eigenvectors().transpose() * F.beta * Mb.eigenvectors();
  Eigen::VectorXd g(n_ov);
  for (Index i = 0; i < toy.nocc_alpha; ++i) {
    for (Index a = 0; a < nva; ++a) {
      g(i * nva + a) =
          UKSConvergenceAcc::EnergyGradientScale * Fa(i, toy.nocc_alpha + a);
    }
  }
  for (Index i = 0; i < toy.nocc_beta; ++i) {
    for (Index a = 0; a < nvb; ++a) {
      g(toy.nocc_alpha * nva + i * nvb + a) =
          UKSConvergenceAcc::EnergyGradientScale * Fb(i, toy.nocc_beta + a);
    }
  }

  std::srand(3);
  Eigen::VectorXd v = Eigen::VectorXd::Random(n_ov).normalized();
  const double t = 1e-3;
  const double e0 = EnergyAlong(toy, conv, Ma, Mb, Eigen::VectorXd::Zero(n_ov));
  const double ep = EnergyAlong(toy, conv, Ma, Mb, t * v);
  const double em = EnergyAlong(toy, conv, Ma, Mb, -t * v);
  const double dE = (ep - em) / (2 * t);
  const double d2E = (ep - 2 * e0 + em) / (t * t);

  auto builder = [&toy](const Eigen::MatrixXd& Da, const Eigen::MatrixXd& Db) {
    return toy.Fock(Da, Db);
  };
  Eigen::VectorXd sigma =
      conv.BuildCoupledSigmaVector(v, Ma.eigenvectors(), toy.nocc_alpha,
                                   Mb.eigenvectors(), toy.nocc_beta, builder);

  BOOST_CHECK_CLOSE(g.dot(v), dE, 1e-3);
  BOOST_CHECK_CLOSE(v.dot(sigma), d2E, 0.1);
}

BOOST_AUTO_TEST_CASE(coupled_augmented_hessian_step_test) {
  // The step must stay inside the trust radius, lower the energy and
  // predict the energy change of the toy model.
  ToyUHF toy = MakeToyUHF();
  Logger log;
  UKSConvergenceAcc conv;
  conv.setLogger(&log);
  conv.setOverlap(toy.overlap, 1e-8);

  auto [Ma, Mb] = StartingOrbitals(toy, conv);
  UKSConvergenceAcc::SpinFock F =
      toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
               ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
  const double e0 = toy.EnergyOfOrbitals(Ma.eigenvectors(), Mb.eigenvectors());
  auto builder = [&toy](const Eigen::MatrixXd& Da, const Eigen::MatrixXd& Db) {
    return toy.Fock(Da, Db);
  };

  for (double trust : {0.02, 0.1}) {
    double predicted = 0.0;
    auto [Ca, Cb] = conv.CoupledAugmentedHessianStep(
        F.alpha, Ma, toy.nocc_alpha, F.beta, Mb, toy.nocc_beta, builder, trust,
        predicted);
    const double actual = toy.EnergyOfOrbitals(Ca, Cb) - e0;
    const double step = std::sqrt(
        std::pow(AppliedStepNorm(toy, Ma.eigenvectors(), Ca, toy.nocc_alpha),
                 2) +
        std::pow(AppliedStepNorm(toy, Mb.eigenvectors(), Cb, toy.nocc_beta),
                 2));
    BOOST_TEST_MESSAGE("trust " << trust << ": step " << step << " predicted "
                                << predicted << " actual " << actual);
    BOOST_CHECK_LE(step, 1.01 * trust);
    BOOST_CHECK_GE(step, 0.85 * trust);
    BOOST_CHECK_LT(predicted, 0.0);
    BOOST_CHECK_LT(actual, 0.0);
    BOOST_CHECK_CLOSE(actual, predicted, 5.0);
  }
}

BOOST_AUTO_TEST_CASE(direct_minimization_rotation_test) {
  // The diagonal fallback must lower the energy of the channel it rotates.
  ToyUHF toy = MakeToyUHF();
  Logger log;
  UKSConvergenceAcc conv;
  conv.setLogger(&log);
  conv.setOverlap(toy.overlap, 1e-8);

  auto [Ma, Mb] = StartingOrbitals(toy, conv);
  UKSConvergenceAcc::SpinFock F =
      toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
               ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
  const double e0 = toy.EnergyOfOrbitals(Ma.eigenvectors(), Mb.eigenvectors());
  double predicted = 0.0;
  Eigen::MatrixXd Ca =
      conv.DirectMinimizationRotation(F.alpha, Ma, toy.nocc_alpha, predicted);
  const double actual = toy.EnergyOfOrbitals(Ca, Mb.eigenvectors()) - e0;
  BOOST_TEST_MESSAGE("diagonal step: predicted " << predicted << " actual "
                                                 << actual);
  BOOST_CHECK_LT(predicted, 0.0);
  BOOST_CHECK_LT(actual, 0.0);
}

BOOST_AUTO_TEST_CASE(augmented_hessian_minimizes_toy_uhf_test) {
  // Repeated steps with a Fletcher trust-radius update, as Iterate does,
  // must lower the energy monotonically and converge quadratically once
  // close: from a gradient norm of ~200 to below 1e-3 within 15 steps.
  // (Further down, the loose Davidson tolerance limits the steps; there
  // Iterate hands over to DIIS anyway.)
  ToyUHF toy = MakeToyUHF();
  Logger log;
  UKSConvergenceAcc conv;
  conv.setLogger(&log);
  conv.setOverlap(toy.overlap, 1e-8);
  auto [Ma, Mb] = StartingOrbitals(toy, conv);
  auto builder = [&toy](const Eigen::MatrixXd& Da, const Eigen::MatrixXd& Db) {
    return toy.Fock(Da, Db);
  };
  auto gradient_norm = [&](const tools::EigenSystem& A,
                           const tools::EigenSystem& B) {
    UKSConvergenceAcc::SpinFock F =
        toy.Fock(ToyUHF::Density(A.eigenvectors(), toy.nocc_alpha),
                 ToyUHF::Density(B.eigenvectors(), toy.nocc_beta));
    Eigen::MatrixXd Fa =
        A.eigenvectors().transpose() * F.alpha * A.eigenvectors();
    Eigen::MatrixXd Fb =
        B.eigenvectors().transpose() * F.beta * B.eigenvectors();
    const Index n = toy.S.rows();
    return std::sqrt(
        Fa.block(0, toy.nocc_alpha, toy.nocc_alpha, n - toy.nocc_alpha)
            .squaredNorm() +
        Fb.block(0, toy.nocc_beta, toy.nocc_beta, n - toy.nocc_beta)
            .squaredNorm());
  };

  double energy = toy.EnergyOfOrbitals(Ma.eigenvectors(), Mb.eigenvectors());
  double trust = 0.2;
  bool monotonic = true;
  Index steps = 0;
  for (; steps < 15 && gradient_norm(Ma, Mb) > 1e-3; ++steps) {
    UKSConvergenceAcc::SpinFock F =
        toy.Fock(ToyUHF::Density(Ma.eigenvectors(), toy.nocc_alpha),
                 ToyUHF::Density(Mb.eigenvectors(), toy.nocc_beta));
    double predicted = 0.0;
    auto [Ca, Cb] = conv.CoupledAugmentedHessianStep(
        F.alpha, Ma, toy.nocc_alpha, F.beta, Mb, toy.nocc_beta, builder, trust,
        predicted);
    const double new_energy = toy.EnergyOfOrbitals(Ca, Cb);
    const double r = (new_energy - energy) / predicted;
    if (r < 0.0) {
      trust *= 0.7;
      continue;
    }
    if (new_energy > energy) {
      monotonic = false;
    }
    trust = (r > 0.75) ? 1.2 * trust : (r <= 0.25 ? 0.7 * trust : trust);
    energy = new_energy;
    UKSConvergenceAcc::SpinFock Fn =
        toy.Fock(ToyUHF::Density(Ca, toy.nocc_alpha),
                 ToyUHF::Density(Cb, toy.nocc_beta));
    Ma.eigenvectors() = Ca;
    Mb.eigenvectors() = Cb;
    Ma.eigenvalues() = (Ca.transpose() * Fn.alpha * Ca).diagonal();
    Mb.eigenvalues() = (Cb.transpose() * Fn.beta * Cb).diagonal();
  }
  BOOST_TEST_MESSAGE("toy UHF: " << steps << " steps, gradient "
                                 << gradient_norm(Ma, Mb) << ", energy "
                                 << energy);
  BOOST_CHECK(monotonic);
  BOOST_CHECK_LT(gradient_norm(Ma, Mb), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
