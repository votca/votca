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
    const Eigen::MatrixXd b = tc.FullMatrix(i);
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

namespace {
// Two copies of the test molecule, the second one shifted by distance (bohr),
// so that part of the shell pairs between them fall below the threshold.
QMMolecule TwoCopies(double distance) {
  QMMolecule single("single", 0);
  single.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/eris/molecule.xyz");
  QMMolecule pair("pair", 0);
  Index id = 0;
  for (const Eigen::Vector3d& shift :
       {Eigen::Vector3d::Zero().eval(), Eigen::Vector3d(0.0, 0.0, distance)}) {
    for (const QMAtom& atom : single) {
      pair.push_back(QMAtom(id++, atom.getElement(), atom.getPos() + shift));
    }
  }
  return pair;
}
}  // namespace

// The pair-screened tensor equals the full one on the stored pairs, its
// dropped entries are below the threshold before the metric transform, and J
// and K change only at that level. Threshold 0 and any build chunking give
// the full tensor.
BOOST_AUTO_TEST_CASE(pair_screened_tensor) {
  libint2::initialize();
  const QMMolecule mol = TwoCopies(12.0);
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  const Index n = aobasis.AOBasisSize();

  TCMatrix_dft full;
  full.Fill(aobasis, aobasis, 0.0);
  BOOST_CHECK_EQUAL(full.StoredPairs(), full.AllPairs());

  // small build chunks: same tensor
  TCMatrix_dft chunked;
  chunked.Fill(aobasis, aobasis, 0.0, 8);
  TCMatrix_dft screened;
  screened.Fill(aobasis, aobasis, 1e-10);
  TCMatrix_dft screened_chunked;
  screened_chunked.Fill(aobasis, aobasis, 1e-10, 8);
  BOOST_TEST_MESSAGE("stored pairs " << screened.StoredPairs() << " of "
                                     << screened.AllPairs());
  BOOST_CHECK_LT(screened.StoredPairs(), screened.AllPairs());
  BOOST_CHECK_GT(screened.StoredPairs(), screened.AllPairs() / 2);
  BOOST_CHECK_EQUAL(screened_chunked.StoredPairs(), screened.StoredPairs());

  double scale = 0.0;
  double chunk_dev = 0.0;
  double kept_dev = 0.0;
  double dropped_max = 0.0;
  for (Index P = 0; P < full.size(); ++P) {
    const Eigen::MatrixXd b = full.FullMatrix(P);
    const Eigen::MatrixXd s = screened.FullMatrix(P);
    scale = std::max(scale, b.cwiseAbs().maxCoeff());
    chunk_dev =
        std::max(chunk_dev, (chunked.FullMatrix(P) - b).cwiseAbs().maxCoeff());
    chunk_dev = std::max(
        chunk_dev, (screened_chunked.FullMatrix(P) - s).cwiseAbs().maxCoeff());
    for (Index j = 0; j < n; ++j) {
      for (Index i = 0; i < n; ++i) {
        if (s(i, j) != 0.0) {
          kept_dev = std::max(kept_dev, std::abs(s(i, j) - b(i, j)));
        } else {
          dropped_max = std::max(dropped_max, std::abs(b(i, j)));
        }
      }
    }
  }
  BOOST_TEST_MESSAGE("largest dropped B " << dropped_max << ", kept deviation "
                                          << kept_dev);
  BOOST_CHECK_SMALL(chunk_dev, 1e-13 * scale);
  BOOST_CHECK_SMALL(kept_dev, 1e-13 * scale);
  BOOST_CHECK_SMALL(dropped_max, 1e-7);

  // packing and unpacking a symmetric matrix
  std::srand(3);
  Eigen::MatrixXd r = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd sym = r + r.transpose();
  const Eigen::MatrixXd b0 = screened.FullMatrix(0);
  BOOST_CHECK_CLOSE(screened.PackWeighted(sym).dot(screened.Data().col(0)),
                    sym.cwiseProduct(b0).sum(), 1e-10);
  const Eigen::MatrixXd roundtrip = screened.Unpack(screened.Data().col(0));
  BOOST_CHECK_SMALL((roundtrip - b0).cwiseAbs().maxCoeff(), 1e-15);

  // J and K from a density of the two molecules
  const Eigen::MatrixXd mos = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/mos.mm");
  const Index nsingle = mos.rows();
  Eigen::MatrixXd occ = Eigen::MatrixXd::Zero(n, 8);
  occ.block(0, 0, nsingle, 4) = mos.leftCols(4);
  occ.block(nsingle, 4, nsingle, 4) = mos.leftCols(4);
  const Eigen::MatrixXd dmat = 2.0 * occ * occ.transpose();

  ERIs eris_full;
  eris_full.Initialize(aobasis, aobasis, 0.0);
  ERIs eris_screened;
  eris_screened.Initialize(aobasis, aobasis, 1e-10);
  const auto jk_full =
      eris_full.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), dmat);
  const auto jk_screened =
      eris_screened.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), dmat);
  const double ej_full = 0.5 * ERIs::CalculateEnergy(dmat, jk_full[0]);
  const double ej_screened = 0.5 * ERIs::CalculateEnergy(dmat, jk_screened[0]);
  const double ek_full = 0.25 * ERIs::CalculateEnergy(dmat, jk_full[1]);
  const double ek_screened = 0.25 * ERIs::CalculateEnergy(dmat, jk_screened[1]);
  BOOST_TEST_MESSAGE("E_J " << ej_full << " deviation " << ej_screened - ej_full
                            << ", E_K " << ek_full << " deviation "
                            << ek_screened - ek_full);
  BOOST_CHECK_SMALL(ej_screened - ej_full, 1e-8);
  BOOST_CHECK_SMALL(ek_screened - ek_full, 1e-8);
  BOOST_CHECK_SMALL((jk_screened[0] - jk_full[0]).cwiseAbs().maxCoeff(), 1e-8);
  BOOST_CHECK_SMALL((jk_screened[1] - jk_full[1]).cwiseAbs().maxCoeff(), 1e-8);

  // the unscreened J agrees with the direct sum over aux functions
  Eigen::MatrixXd j_direct = Eigen::MatrixXd::Zero(n, n);
  for (Index P = 0; P < full.size(); ++P) {
    const Eigen::MatrixXd b = full.FullMatrix(P);
    j_direct += b.cwiseProduct(dmat).sum() * b;
  }
  BOOST_CHECK_SMALL((jk_full[0] - j_direct).cwiseAbs().maxCoeff(),
                    1e-12 * j_direct.cwiseAbs().maxCoeff());
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
