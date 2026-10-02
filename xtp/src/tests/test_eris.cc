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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE eris_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/ERIs.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace std;
using votca::Index;

BOOST_AUTO_TEST_SUITE(eris_test)

BOOST_AUTO_TEST_CASE(fourcenter) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/dmat.mm");

  ERIs eris;
  eris.Initialize_4c(aobasis);
  Eigen::MatrixXd erissmall = eris.CalculateERIs_4c(dmat, 1e-20);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref.mm");

  bool eris_check = erissmall.isApprox(eris_ref, 0.00001);
  if (!eris_check) {
    std::cout << "result eri" << std::endl;
    std::cout << erissmall << std::endl;
    std::cout << "ref eri" << std::endl;
    std::cout << eris_ref << std::endl;
    std::cout << " quotient" << std::endl;
    std::cout << erissmall.cwiseQuotient(eris_ref);
  }
  BOOST_CHECK_EQUAL(eris_check, 1);

  std::array<Eigen::MatrixXd, 2> both = eris.CalculateERIs_EXX_4c(dmat, 1e-20);
  const Eigen::MatrixXd& exx_small = both[1];
  Eigen::MatrixXd exx_ref = -votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref.mm");

  Eigen::MatrixXd sum = both[0] + both[1];
  Eigen::MatrixXd sum_ref = exx_ref + eris_ref;
  bool sum_check = sum.isApprox(sum_ref, 1e-5);
  BOOST_CHECK_EQUAL(sum_check, 1);
  if (!sum_check) {
    std::cout << "result sum" << std::endl;
    std::cout << sum << std::endl;
    std::cout << "ref sum" << std::endl;
    std::cout << sum_ref << std::endl;
  }

  bool exxs_check = exx_small.isApprox(exx_ref, 0.00001);
  if (!eris_check) {
    std::cout << "result exx" << std::endl;
    std::cout << exx_small << std::endl;
    std::cout << "ref exx" << std::endl;
    std::cout << exx_ref << std::endl;
    std::cout << "quotient" << std::endl;
    std::cout << exx_small.cwiseQuotient(exx_ref);
  }
  BOOST_CHECK_EQUAL(exxs_check, 1);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(threecenter) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/dmat2.mm");

  Eigen::MatrixXd mos = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/mos.mm");

  ERIs eris;
  eris.Initialize(aobasis, aobasis);
  Eigen::MatrixXd exx_dmat =
      eris.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), dmat)[1];
  Eigen::MatrixXd exx_mo =
      eris.CalculateERIs_EXX_3c(mos.block(0, 0, 17, 4), dmat)[1];

  bool compare_exx = exx_mo.isApprox(exx_dmat, 1e-4);
  BOOST_CHECK_EQUAL(compare_exx, true);

  Eigen::MatrixXd exx_ref = -votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref2.mm");

  bool compare_exx_ref = exx_ref.isApprox(exx_mo, 1e-5);
  if (!compare_exx_ref) {
    std::cout << "result exx" << std::endl;
    std::cout << exx_mo << std::endl;
    std::cout << "ref exx" << std::endl;
    std::cout << exx_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_exx_ref, true);

  Eigen::MatrixXd eri = eris.CalculateERIs_3c(dmat);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref2.mm");

  bool compare_eris = eris_ref.isApprox(eri, 1e-5);
  if (!compare_eris) {
    std::cout << "result eris" << std::endl;
    std::cout << eri << std::endl;
    std::cout << "ref eris" << std::endl;
    std::cout << eris_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_eris, true);
}

namespace {
// The density-matrix exchange as it used to be computed: two N^3 products
// per aux function, straight from the 3c tensor.
Eigen::MatrixXd BruteForceExchange(const TCMatrix_dft& tc,
                                   const Eigen::MatrixXd& dmat) {
  Eigen::MatrixXd exx = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  for (Index i = 0; i < tc.size(); i++) {
    const Eigen::MatrixXd b = tc[i].FullMatrix();
    exx -= b * dmat * b;
  }
  return exx;
}
}  // namespace

// The factorised density-matrix route must reproduce K[D] = -sum_P B_P D B_P
// for any symmetric D, not just for densities of rank N_occ: an indefinite
// full-rank matrix (like an incremental density difference), a mixed density
// of two different idempotent ones, and 2 C C^T, where it must also agree
// with the occupied-MO route to rounding.
BOOST_AUTO_TEST_CASE(exchange_from_factorised_density) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  const Index n = aobasis.AOBasisSize();

  ERIs eris;
  eris.Initialize(aobasis, aobasis);
  TCMatrix_dft tc;
  tc.Fill(aobasis, aobasis);

  auto exchange = [&](const Eigen::MatrixXd& d) {
    return eris.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), d)[1];
  };
  auto check = [](const Eigen::MatrixXd& result, const Eigen::MatrixXd& ref,
                  const std::string& what) {
    const double err = (result - ref).cwiseAbs().maxCoeff();
    BOOST_TEST_MESSAGE(what << ": max deviation " << err);
    BOOST_CHECK_MESSAGE(err < 1e-10 * ref.cwiseAbs().maxCoeff(),
                        what << ": max deviation " << err);
  };

  std::srand(7);
  Eigen::MatrixXd r = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd indefinite = r + r.transpose();
  check(exchange(indefinite), BruteForceExchange(tc, indefinite),
        "indefinite full-rank matrix");

  const Eigen::MatrixXd mos = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/mos.mm");
  const Eigen::MatrixXd occ = mos.leftCols(4);
  const Eigen::MatrixXd d1 = 2.0 * occ * occ.transpose();
  // a different occupied set: swap HOMO for LUMO
  Eigen::MatrixXd occ2 = mos.leftCols(5);
  occ2.col(3) = occ2.col(4);
  occ2.conservativeResize(n, 4);
  const Eigen::MatrixXd d2 = 2.0 * occ2 * occ2.transpose();
  const Eigen::MatrixXd mixed = 0.7 * d1 + 0.3 * d2;
  check(exchange(mixed), BruteForceExchange(tc, mixed), "mixed density");

  check(exchange(d1), BruteForceExchange(tc, d1), "2 C C^T");
  check(exchange(d1), eris.CalculateERIs_EXX_3c(occ, d1)[1],
        "2 C C^T vs occupied-MO route");

  BOOST_CHECK_EQUAL(exchange(Eigen::MatrixXd::Zero(n, n)).cwiseAbs().maxCoeff(),
                    0.0);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
