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
// STATUS: DRAFT, UNCOMPILED, UNTESTED (see libint2_derivative_calls.cc).
//
// This is the intended first correctness check for ComputeOverlapDerivatives
// / ComputeKineticDerivatives: compare the analytic derivative against a
// central finite difference of the corresponding AOOverlap/AOKinetic
// matrices under small nuclear displacements. If the buffer-ordering
// assumption flagged in libint2_derivative_calls.cc is wrong, this test
// should fail clearly (large discrepancy, not a subtle one) rather than
// silently passing.
// ===========================================================================

#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aoderivatives_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

// Forward declarations of the new functions under test (declared in
// libint2_derivative_calls.cc; not yet added to a public header, since this
// is still a draft -- move to aomatrix.h or a new header once the API
// stabilizes).
namespace votca {
namespace xtp {
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);
}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_SUITE(aoderivatives_test)

// Builds a small test molecule (H2, two atoms, minimal basis) directly
// rather than loading from a data file, since no dedicated test geometry
// exists yet for this new test -- two atoms is also the minimal case that
// exercises both the "derivative w.r.t. shell1's atom" and "derivative
// w.r.t. shell2's atom" branches in the new code, including the s1 != s2
// transpose-block path.
QMMolecule BuildH2(double bond_length_angstrom) {
  QMMolecule mol(" ", 0);
  std::string xyz_content =
      "2\n\n"
      "H 0.0 0.0 0.0\n"
      "H 0.0 0.0 " +
      std::to_string(bond_length_angstrom) + "\n";
  // NOTE: QMMolecule::LoadFromFile reads from a file path in the existing
  // test pattern (test_threecenter_dft.cc); writing a temp file here rather
  // than assuming an in-memory constructor exists, since that API wasn't
  // confirmed against qmmolecule.h in this pass.
  std::string tmp_path = "/tmp/xtp_test_h2_tmp.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();
  mol.LoadFromFile(tmp_path);
  return mol;
}

BOOST_AUTO_TEST_CASE(overlap_derivative_finite_difference) {
  libint2::initialize();

  // NOTE: basis set XML path assumed available from existing test data
  // conventions (e.g. a minimal STO-3G-like set already used by other
  // tests) -- point this at whatever small basis XML already exists under
  // XTP_TEST_DATA_FOLDER once compiling; using a placeholder path here.
  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-4;            // finite-difference step, Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildH2(bond_length);
  AOBasis aobasis0;
  aobasis0.Fill(basisset, mol0);

  auto deriv = ComputeOverlapDerivatives(aobasis0);
  // deriv[0] = derivative w.r.t. atom 0 (first H), deriv[1] = w.r.t. atom 1
  // (second H); by symmetry for this diatomic along z, only the z-component
  // (index 2) should be non-negligible, and deriv[0][2] should equal
  // -deriv[1][2] since displacing atom 0 by +dz has the same effect on the
  // bond-length-dependent overlap as displacing atom 1 by -dz. This
  // symmetry is itself a cheap first sanity check, independent of the
  // finite-difference comparison below.
  BOOST_CHECK_SMALL(
      (deriv[0][2] + deriv[1][2]).cwiseAbs().maxCoeff(), 1e-8);

  // Central finite difference: displace atom 1 (index 1, the second H)
  // along z by +-h, recompute overlap, compare to analytic deriv[1][2].
  QMMolecule mol_plus = BuildH2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basisset, mol_plus);
  AOOverlap overlap_plus;
  overlap_plus.Fill(aobasis_plus);

  QMMolecule mol_minus = BuildH2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basisset, mol_minus);
  AOOverlap overlap_minus;
  overlap_minus.Fill(aobasis_minus);

  Eigen::MatrixXd finite_diff_deriv =
      (overlap_plus.Matrix() - overlap_minus.Matrix()) / (2.0 * h);

  bool matches = finite_diff_deriv.isApprox(deriv[1][2], 1e-5);
  if (!matches) {
    std::cout << "Analytic dS/dz(atom1):\n" << deriv[1][2] << std::endl;
    std::cout << "Finite-difference dS/dz(atom1):\n"
              << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails with a large (not small numerical) "
                 "discrepancy, check the buffer-ordering assumption flagged "
                 "in libint2_derivative_calls.cc first."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
