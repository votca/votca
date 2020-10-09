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

#define BOOST_TEST_MODULE threecenter_dft_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include <votca/xtp/aobasis.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/threecenter.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(threecenter_dft_test)

Eigen::MatrixXd ReadMatrixFromString(const std::string& matrix) {
  votca::tools::Tokenizer lines(matrix, "\n");

  std::vector<double> entries;
  Index cols = 0;
  Index rows = 0;
  for (auto line : lines) {
    if (line[0] == '#') {
      continue;
    }
    votca::tools::Tokenizer entries_tok(line, " ");
    std::vector<std::string> temp = entries_tok.ToVector();
    cols = Index(temp.size());
    rows++;
    for (const auto& s : temp) {
      entries.push_back(std::stod(s));
    }
  }

  return Eigen::Map<Eigen::MatrixXd>(entries.data(), rows, cols);
}

BOOST_AUTO_TEST_CASE(small_basis) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_dft/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  TCMatrix_dft threec;
  threec.Fill(aobasis, aobasis);

  TCMatrix_dft threec2;
  threec2.Fill2(aobasis, aobasis);

  Eigen::MatrixXd Ref0 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/Ref0.mm");

  Eigen::MatrixXd Ref4 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/Ref4.mm");

  bool check_three1 = Ref0.isApprox(threec[0].FullMatrix(), 0.00001);
  if (!check_three1) {
    std::cout << "Res0" << std::endl;
    std::cout << threec[0].FullMatrix() << std::endl;
    std::cout << "0_ref" << std::endl;
    std::cout << Ref0 << std::endl;
  }
  BOOST_CHECK_EQUAL(check_three1, true);
  bool check_three2 = Ref4.isApprox(threec[4].FullMatrix(), 0.00001);
  if (!check_three2) {
    std::cout << "Res4" << std::endl;
    std::cout << threec[4].FullMatrix() << std::endl;
    std::cout << "4_ref" << std::endl;
    std::cout << Ref4 << std::endl;
  }
  BOOST_CHECK_EQUAL(check_three2, true);
  libint2::finalize();
}

/*BOOST_AUTO_TEST_CASE(large_l_test) {

  QMMolecule mol("C", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_dft/C2.xyz");

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/G.xml");
  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);
  BasisSet auxbasisset;
  auxbasisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_dft/I.xml");
  AOBasis auxbasis;
  auxbasis.Fill(auxbasisset, mol);

  TCMatrix_dft threec;
  threec.Fill(auxbasis, dftbasis);
  // we only test half of it because it gets a bit big
  std::array<Eigen::MatrixXd, 4> ref;
  std::array<std::string, 4> ref_string;

  std::array<int, 4> indeces = {0, 1, 3, 12};
  for (Index i = 0; i < 4; i++) {
    ref[i] = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/RefList" +
        std::to_string(i) + ".mm");
  }

  for (Index i = 0; i < 4; i++) {
    bool check = ref[i].isApprox(threec[indeces[i]].FullMatrix(), 1e-5);
    BOOST_CHECK_EQUAL(check, true);
    if (!check) {
      std::cout << "ref " << indeces[i] << std::endl;
      std::cout << ref[i] << std::endl;
      std::cout << "result " << indeces[i] << std::endl;
      std::cout << threec[indeces[i]].FullMatrix() << std::endl;
    }
  }
} */
BOOST_AUTO_TEST_SUITE_END()
