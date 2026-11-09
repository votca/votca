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
#include "xtp_libint2.h"

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE diabatwriter_test

// Standard includes
#include <cmath>
#include <regex>
#include <sstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/diabatwriter.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(diabatwriter_test)

namespace {

// methane/def2-TZVP (no MO data needed), 5 occupied, BSE window MO 1..8
Orbitals SetupOrbitals(bool tda) {
  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/gaussianwriter/methane.xyz");
  orb.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                    "/gaussianwriter/def2-tzvp.xml");
  orb.setNumberOfAlphaElectrons(5);
  orb.setNumberOfOccupiedLevels(5);
  orb.setChargeAndSpin(0, 1);
  orb.setQMEnergy(-40.5);
  orb.setBSEindices(1, 8);  // v = 1..4, c = 5..8 (0-based)
  orb.setTDAApprox(tda);

  const Index size = 16;
  Eigen::MatrixXd X(size, 2);
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(size, 2);
  for (Index i = 0; i < size; ++i) {
    X(i, 0) = 0.1 * double(i + 1);
    X(i, 1) = std::cos(double(i));
    if (!tda) {
      Y(i, 0) = 0.01 * double(i);
      Y(i, 1) = -0.02 * double(i % 3);
    }
  }
  for (Index s = 0; s < 2; ++s) {  // X^T X - Y^T Y = 1
    const double n = std::sqrt(X.col(s).squaredNorm() - Y.col(s).squaredNorm());
    X.col(s) /= n;
    Y.col(s) /= n;
  }
  orb.BSESinglets().eigenvalues() = Eigen::Vector2d(0.2, 0.3);
  orb.BSESinglets().eigenvectors() = X;
  if (!tda) {
    orb.BSESinglets().eigenvectors2() = Y;
  }
  return orb;
}

struct Amplitude {
  Index occ;
  Index virt;
  bool deexcitation;
  double value;
};

// amplitudes of the first state in the file
std::vector<Amplitude> FirstStateAmplitudes(const std::string& text) {
  std::vector<Amplitude> result;
  std::istringstream in(text);
  std::string line;
  const std::regex amp(R"(^\s+(\d+)\s+(-->|<--)\s+(\d+)\s+(\S+)\s+alpha\s*$)");
  Index state = 0;
  while (std::getline(in, line)) {
    if (line.find("State index:") != std::string::npos) {
      ++state;
    }
    std::smatch m;
    if (state == 1 && std::regex_match(line, m, amp)) {
      result.push_back(
          {std::stoi(m[1]), std::stoi(m[3]), m[2] == "<--", std::stod(m[4])});
    }
  }
  return result;
}

}  // namespace

BOOST_AUTO_TEST_CASE(full_bse_singlets) {
  libint2::initialize();
  Orbitals orb = SetupOrbitals(false);
  std::stringstream ss;
  DiabatWriter().Write(ss, orb, QMStateType(QMStateType::Singlet));
  const std::string text = ss.str();

  BOOST_CHECK_EQUAL(text.substr(0, text.find('\n')),
                    DiabatWriter::FormatIdentifier);
  BOOST_CHECK(text.find("\nMethod: TDDFT\n") != std::string::npos);
  BOOST_CHECK(text.find("Number of virtual orbitals (alpha): " +
                        std::to_string(orb.getBasisSetSize() - 5)) !=
              std::string::npos);
  BOOST_CHECK(text.find("  State index: 2\n") != std::string::npos);

  std::vector<Amplitude> amps = FirstStateAmplitudes(text);
  BOOST_REQUIRE_EQUAL(amps.size(), 32);
  double nx = 0.0;
  double ny = 0.0;
  const double s = 1.0 / std::sqrt(2.0);
  for (const Amplitude& a : amps) {
    // absolute 1-based indices, VOTCA compound index v * ctotal + c
    const Index v = a.occ - 2;
    const Index c = a.virt - 6;
    BOOST_REQUIRE(v >= 0 && v < 4 && c >= 0 && c < 4);
    const double ref = a.deexcitation
                           ? orb.BSESinglets().eigenvectors2()(v * 4 + c, 0)
                           : orb.BSESinglets().eigenvectors()(v * 4 + c, 0);
    BOOST_CHECK_CLOSE_FRACTION(a.value, s * ref, 1e-12);
    (a.deexcitation ? ny : nx) += a.value * a.value;
  }
  BOOST_CHECK_CLOSE_FRACTION(nx - ny, 0.5, 1e-12);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(tda_singlets_and_errors) {
  libint2::initialize();
  Orbitals orb = SetupOrbitals(true);
  std::stringstream ss;
  DiabatWriter().Write(ss, orb, QMStateType(QMStateType::Singlet), 1);
  const std::string text = ss.str();
  BOOST_CHECK(text.find("\nMethod: TDA\n") != std::string::npos);
  BOOST_CHECK(text.find("<--") == std::string::npos);
  BOOST_CHECK(text.find("State index: 2") == std::string::npos);
  BOOST_CHECK_EQUAL(FirstStateAmplitudes(text).size(), 16);

  // no triplets stored
  std::stringstream dummy;
  BOOST_CHECK_THROW(
      DiabatWriter().Write(dummy, orb, QMStateType(QMStateType::Triplet)),
      std::runtime_error);
  // open shell not supported
  orb.setChargeAndSpin(1, 2);
  BOOST_CHECK_THROW(
      DiabatWriter().Write(dummy, orb, QMStateType(QMStateType::Singlet)),
      std::runtime_error);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
