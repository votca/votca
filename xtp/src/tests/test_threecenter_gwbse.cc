/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE threecenter_gwbse_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA inlcudes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/threecenter.h"
#include "xtp_libint2.h"
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(threecenter_gwbse_test)
BOOST_AUTO_TEST_CASE(threecenter_gwbse) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
             "/threecenter_gwbse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  TCMatrix_gwbse tc;
  tc.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aobasis, aobasis, MOs);

  Eigen::MatrixXd ref0b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref0b.mm");

  bool check0_before = ref0b.isApprox(tc[0], 1e-5);
  if (!check0_before) {
    cout << "tc0" << endl;
    cout << tc[0] << endl;
    cout << "tc0_ref" << endl;
    cout << ref0b << endl;
  }
  BOOST_CHECK_EQUAL(check0_before, true);

  Eigen::MatrixXd ref2b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref2b.mm");

  bool check2_before = ref2b.isApprox(tc[2], 1e-5);
  if (!check2_before) {
    cout << "tc2" << endl;
    cout << tc[2] << endl;
    cout << "tc2_ref" << endl;
    cout << ref2b << endl;
  }

  BOOST_CHECK_EQUAL(check2_before, true);

  Eigen::MatrixXd ref4b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref4b.mm");

  bool check4_before = ref4b.isApprox(tc[4], 1e-5);
  if (!check4_before) {
    cout << "tc4" << endl;
    cout << tc[4] << endl;
    cout << "tc4_ref" << endl;
    cout << ref4b << endl;
  }

  BOOST_CHECK_EQUAL(check4_before, true);

  Eigen::MatrixXd auxmatrix =
      Eigen::MatrixXd::Identity(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  tc.MultiplyRightWithAuxMatrix(auxmatrix);

  bool check0_after = ref0b.isApprox(tc[0], 1e-5);
  if (!check0_after) {
    cout << "tc0" << endl;
    cout << tc[0] << endl;
    cout << "tc0_ref" << endl;
    cout << ref0b << endl;
  }
  BOOST_CHECK_EQUAL(check0_after, true);

  bool check2_after = ref2b.isApprox(tc[2], 1e-5);
  if (!check2_after) {
    cout << "tc2" << endl;
    cout << tc[2] << endl;
    cout << "tc2_ref" << endl;
    cout << ref2b << endl;
  }

  BOOST_CHECK_EQUAL(check2_after, true);

  bool check4_after = ref4b.isApprox(tc[4], 1e-5);
  if (!check4_after) {
    cout << "tc4" << endl;
    cout << tc[4] << endl;
    cout << "tc4_ref" << endl;
    cout << ref4b << endl;
  }

  BOOST_CHECK_EQUAL(check4_after, true);

  libint2::finalize();
}

// InvSqrt() must be the metric that was actually folded into the stored
// integrals, not merely some inverse square root of V. What everything
// downstream relies on is that M M^T reproduces the RI Coulomb interaction,
// which holds exactly when T^T V T is the projector onto the retained
// auxiliary functions: the identity if nothing was removed, and in general
// symmetric, idempotent, with trace = size - Removedfunctions().
//
// Checked against a V filled independently here, so the test cannot pass
// by comparing the stored matrix with itself.
BOOST_AUTO_TEST_CASE(stored_metric_is_the_one_folded_into_the_integrals) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
             "/threecenter_gwbse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  TCMatrix_gwbse tc;
  tc.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aobasis, aobasis, MOs);

  const Eigen::MatrixXd& T = tc.InvSqrt();
  const votca::Index n = aobasis.AOBasisSize();
  BOOST_REQUIRE_EQUAL(T.rows(), n);
  BOOST_REQUIRE_EQUAL(T.cols(), n);

  AOCoulomb V;
  V.Fill(aobasis);
  const Eigen::MatrixXd P = T.transpose() * V.Matrix() * T;

  BOOST_CHECK_SMALL((P - P.transpose()).cwiseAbs().maxCoeff(), 1e-10);
  BOOST_CHECK_SMALL((P * P - P).cwiseAbs().maxCoeff(), 1e-8);
  BOOST_CHECK_CLOSE(P.trace(), double(n - tc.Removedfunctions()), 1e-8);
  if (tc.Removedfunctions() == 0) {
    BOOST_CHECK_SMALL(
        (P - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff(), 1e-8);
  }

  libint2::finalize();
}
BOOST_AUTO_TEST_CASE(leading_rotation_touches_only_the_leading_block) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
             "/threecenter_gwbse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  TCMatrix_gwbse tc;
  tc.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aobasis, aobasis, MOs);
  std::vector<Eigen::MatrixXd> before;
  for (votca::Index m = 0; m < tc.msize(); ++m) {
    before.push_back(tc[m]);
  }
  const votca::Index n = tc.auxsize();
  const Eigen::MatrixXd U =
      Eigen::MatrixXd::Random(n, n).householderQr().householderQ();

  const votca::Index full = 2;
  const votca::Index lead = 4;
  tc.MultiplyRightWithAuxMatrixLeading(U, full, lead);
  BOOST_CHECK(tc.PartiallyRotated());
  for (votca::Index m = 0; m < tc.msize(); ++m) {
    Eigen::MatrixXd expected = before[std::size_t(m)];
    if (m < full) {
      expected = before[std::size_t(m)] * U;
    } else if (m < lead) {
      expected.topRows(lead) = before[std::size_t(m)].topRows(lead) * U;
    }
    BOOST_CHECK_LE((tc[m] - expected).cwiseAbs().maxCoeff(), 1e-12);
  }
  // the frame of the rotated part is recorded; dressing needs the whole
  // tensor in one frame
  BOOST_CHECK_LE((tc.AuxFrame() - U).cwiseAbs().maxCoeff(), 0.0);
  BOOST_CHECK_THROW(tc.DressAuxIndex(Eigen::MatrixXd::Identity(n, n)),
                    std::runtime_error);

  // a rebuild starts over
  tc.Rebuild();
  BOOST_CHECK(!tc.PartiallyRotated());

  // leading block covering everything is a full rotation
  tc.MultiplyRightWithAuxMatrixLeading(U, tc.msize(), tc.nsize());
  BOOST_CHECK(!tc.PartiallyRotated());
  for (votca::Index m = 0; m < tc.msize(); ++m) {
    BOOST_CHECK_LE((tc[m] - before[std::size_t(m)] * U).cwiseAbs().maxCoeff(),
                   1e-12);
  }
  libint2::finalize();
}

// Rotating the orbitals of a window is the same as filling the integrals
// with the rotated MO coefficients, for slices and rows in and outside the
// window.
BOOST_AUTO_TEST_CASE(rotate_orbitals_equals_fill_with_rotated_mos) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
             "/threecenter_gwbse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  const votca::Index first = 1;
  const votca::Index last = 4;
  const votca::Index window = last - first + 1;
  std::srand(7);
  const Eigen::MatrixXd U =
      Eigen::MatrixXd::Random(window, window).householderQr().householderQ();
  Eigen::MatrixXd MOs_rotated = MOs;
  MOs_rotated.middleCols(first, window) = MOs.middleCols(first, window) * U;

  TCMatrix_gwbse tc;
  tc.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aobasis, aobasis, MOs);
  tc.RotateOrbitals(U, first, last);
  TCMatrix_gwbse ref;
  ref.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  ref.Fill(aobasis, aobasis, MOs_rotated);
  for (votca::Index m = 0; m < tc.msize(); ++m) {
    BOOST_CHECK_LE((tc[m] - ref[m]).cwiseAbs().maxCoeff(), 1e-12);
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
