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

// Checks AOShell::EvalAOspaceHessian (second derivatives of the AO values
// w.r.t. the electron position, used by the GGA XC gradient) against
// central finite differences of the first derivatives from EvalAOspace, for
// an SPDFG shell with two primitives (aoshell/Al.xyz,
// aoshell/largeshell.xml) at a generic off-axis point; also checks that the
// Hessians are symmetric.

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

  BOOST_REQUIRE(aobasis.getNumofShells() >= 1);
  // AOBasis::Fill splits the SPDFG shell into one AOShell per l (5 shells
  // sharing the primitives), so count the functions over all shells.
  Index total_funcs = 0;
  for (Index s = 0; s < aobasis.getNumofShells(); ++s) {
    total_funcs += aobasis.getShell(s).getNumFunc();
  }
  BOOST_REQUIRE_EQUAL(total_funcs, 25);  // 1(S)+3(P)+5(D)+7(F)+9(G)

  // Generic point, not on an axis or plane through the centre (Al at the
  // origin), where many second-derivative terms vanish.
  Eigen::Vector3d point(0.83, -0.51, 1.27);
  double h = 1e-5;

  double max_rel_error = 0.0;
  Index worst_shell = -1, worst_func = -1;
  Index worst_i = -1, worst_j = -1;
  double max_asymmetry = 0.0;

  for (Index s = 0; s < aobasis.getNumofShells(); ++s) {
    const AOShell& shell = aobasis.getShell(s);
    AOShell::AOValuesHessian AO = shell.EvalAOspaceHessian(point);

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
          double scale =
              std::max({std::abs(analytic), std::abs(numeric), 1e-8});
          double rel_error = std::abs(analytic - numeric) / scale;
          if (rel_error > max_rel_error) {
            max_rel_error = rel_error;
            worst_shell = s;
            worst_func = k;
            worst_i = i;
            worst_j = j;
          }
        }
      }
    }

    // Hessian symmetry, independent of the finite differences
    for (Index k = 0; k < shell.getNumFunc(); ++k) {
      double asymmetry =
          (AO.hessians[k] - AO.hessians[k].transpose()).cwiseAbs().maxCoeff();
      max_asymmetry = std::max(max_asymmetry, asymmetry);
    }
  }

  if (max_rel_error > 1e-4) {
    std::cout << "Largest relative Hessian error: " << max_rel_error
              << " at shell " << worst_shell
              << " (L=" << EnumToString(aobasis.getShell(worst_shell).getL())
              << "), function " << worst_func << ", component (" << worst_i
              << "," << worst_j << ")" << std::endl;
  }
  BOOST_CHECK_SMALL(max_rel_error, 1e-4);
  BOOST_CHECK_SMALL(max_asymmetry, 1e-10);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
