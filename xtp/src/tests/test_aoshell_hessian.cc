/*
 * Copyright 2009-2024 The VOTCA Development Team (http://www.votca.org)
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

// ===========================================================================
// STATUS: written but NOT yet run/tested. Validates AOShell::EvalAOspaceHessian
// (added for GGA XC-gradient support -- see aoshell.h/aoshell.cc for the
// full derivation) against a finite difference of the existing, already-
// validated first-derivative EvalAOspace output.
//
// Every individual d2P polynomial coefficient (the shell-specific new
// content) was already cross-checked symbolically against exact sympy
// differentiation before this test was written -- see conversation
// history. This test validates the FULL assembled Hessian (including the
// Gaussian-times-polynomial product-rule terms, the per-primitive
// accumulation for a genuine multi-primitive contraction, and the one
// scaling-convention bug already found and fixed by that symbolic check)
// against real numerical evidence, the same way every other derivative
// in this whole branch has been validated.
//
// Reuses aoshell/Al.xyz + aoshell/largeshell.xml from the existing,
// already-validated test_aoshell.cc -- a single SPDFG shell (S through G,
// every angular momentum EvalAOspaceHessian supports) with TWO Gaussian
// primitives (a genuine contraction, not a single-primitive edge case).
// ===========================================================================

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aoshell_hessian_test

// Standard includes
#include <cmath>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/orbitals.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(aoshell_hessian_test)

BOOST_AUTO_TEST_CASE(EvalAOspaceHessian_finite_difference) {
  libint2::initialize();

  QMMolecule mol = QMMolecule("", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aoshell/Al.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aoshell/largeshell.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  BOOST_REQUIRE_EQUAL(aobasis.getNumofShells(), 1);
  const AOShell& shell = aobasis.getShell(0);
  BOOST_REQUIRE_EQUAL(shell.getNumFunc(), 25);  // 1(S)+3(P)+5(D)+7(F)+9(G)

  // Deliberately NOT at the shell center (Al is at the origin) or on any
  // coordinate axis/plane -- an off-axis, generic point exercises every
  // term in every polynomial (many of the d2P entries vanish identically
  // on-axis, e.g. d2P/dxy for several functions is proportional to
  // x*y or similar, which would trivially vanish at x=0 or y=0 and not
  // actually test that term).
  Eigen::Vector3d point(0.83, -0.51, 1.27);
  double h = 1e-5;

  AOShell::AOValuesHessian AO = shell.EvalAOspaceHessian(point);

  double max_rel_error = 0.0;
  Index worst_func = -1;
  Index worst_i = -1, worst_j = -1;

  for (Index j = 0; j < 3; ++j) {
    Eigen::Vector3d point_plus = point;
    Eigen::Vector3d point_minus = point;
    point_plus(j) += h;
    point_minus(j) -= h;

    AOShell::AOValues AO_plus = shell.EvalAOspace(point_plus);
    AOShell::AOValues AO_minus = shell.EvalAOspace(point_minus);

    // Finite difference of the first derivative (a Matrix, size
    // nfunc x 3) w.r.t. displacement in direction j gives column j of
    // every function's Hessian: d(derivatives(:,i))/dx_j.
    Eigen::MatrixX3d fd_hessian_col_j =
        (AO_plus.derivatives - AO_minus.derivatives) / (2.0 * h);

    for (Index k = 0; k < shell.getNumFunc(); ++k) {
      for (Index i = 0; i < 3; ++i) {
        double analytic = AO.hessians[k](i, j);
        double numeric = fd_hessian_col_j(k, i);
        double scale = std::max({std::abs(analytic), std::abs(numeric), 1e-8});
        double rel_error = std::abs(analytic - numeric) / scale;
        if (rel_error > max_rel_error) {
          max_rel_error = rel_error;
          worst_func = k;
          worst_i = i;
          worst_j = j;
        }
      }
    }
  }

  if (max_rel_error > 1e-4) {
    std::cout << "Largest relative Hessian error: " << max_rel_error
              << " at function " << worst_func << ", component (" << worst_i
              << "," << worst_j << ")" << std::endl;
    std::cout << "Analytic Hessian for that function:\n"
              << AO.hessians[worst_func] << std::endl;
  }
  BOOST_CHECK_SMALL(max_rel_error, 1e-4);

  // Also check the symmetry of every Hessian directly (a real,
  // independent mathematical requirement -- H_ij must equal H_ji for
  // any twice-differentiable function -- and a cheap check that doesn't
  // depend on finite-difference precision at all).
  double max_asymmetry = 0.0;
  for (Index k = 0; k < shell.getNumFunc(); ++k) {
    double asymmetry = (AO.hessians[k] - AO.hessians[k].transpose())
                            .cwiseAbs()
                            .maxCoeff();
    max_asymmetry = std::max(max_asymmetry, asymmetry);
  }
  BOOST_CHECK_SMALL(max_asymmetry, 1e-10);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
