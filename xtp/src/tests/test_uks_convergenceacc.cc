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
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/uks_convergenceacc.h"
#include "xtp_libint2.h"

#include <fstream>
#include <tuple>
#include <cmath>

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
  Index nvirt_alpha = nao_alpha - nocclevels_alpha;  // 3
  Index n_ov_alpha = nocclevels_alpha * nvirt_alpha;  // 6

  Index nao_beta = 4;
  Index nocclevels_beta = 1;
  Index nvirt_beta = nao_beta - nocclevels_beta;  // 3
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


BOOST_AUTO_TEST_SUITE_END()
