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

#define BOOST_TEST_MODULE cubefilewriter_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/cubefile_writer.h"
#include "votca/tools/eigenio_matrixmarket.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(cubefilewriter_test)

Eigen::VectorXd Readcubefile(const std::string& filename) {

  std::ifstream in1;
  in1.open(filename, std::ios::in);

  std::string result = "";
  std::string s;
  getline(in1, s);
  getline(in1, s);
  getline(in1, s);
  std::vector<double> cube_values;
  do {
    votca::tools::Tokenizer tok(s, " ");
    std::vector<double> values;
    tok.ConvertToVector<double>(values);
    cube_values.insert(cube_values.end(), values.begin(), values.end());
  } while (getline(in1, s));
  return Eigen::Map<Eigen::VectorXd>(cube_values.data(), cube_values.size());
}

BOOST_AUTO_TEST_CASE(constructors_test) {

  Orbitals A;
  A.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/3-21G.xml");
  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/cubefile_writer/molecule.xyz");
  A.setBasisSetSize(17);
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/A_MOs.mm");

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/spsi_ref.mm");

  A.BSESinglets().eigenvectors() = spsi_ref;

  Eigen::Array<votca::Index, 3, 1> steps(3, 4, 7);
  Logger log;
  double padding = 0.5;
  CubeFile_Writer writer(steps, padding, log);
  QMState state("s1");
  writer.WriteFile("test_writer.cube", A, state, false);

  auto result1 = Readcubefile("test_writer.cube");

  Eigen::VectorXd values_ref1 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/values_ref1.mm");

  BOOST_CHECK_EQUAL(values_ref1.size(), result1.size());

  bool check_ref1 = values_ref1.isApprox(result1, 1e-4);
  BOOST_CHECK_EQUAL(check_ref1, true);
  if (!check_ref1) {
    std::cout << "ref" << std::endl;
    std::cout << values_ref1.transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << result1.transpose() << std::endl;
  }

  writer.WriteFile("test_writer2.cube", A, state, true);

  auto result2 = Readcubefile("test_writer2.cube");

  Eigen::VectorXd values_ref2 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/values_ref2.mm");

  bool check_ref2 = values_ref2.isApprox(result2, 1e-4);
  BOOST_CHECK_EQUAL(values_ref2.size(), result2.size());
  BOOST_CHECK_EQUAL(check_ref2, true);
  if (!check_ref2) {
    std::cout << "ref" << std::endl;
    std::cout << values_ref2.transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << result2.transpose() << std::endl;
  }

  Eigen::Array<votca::Index, 3, 1> steps2(3, 4, 8);
  Logger log2;
  double padding2 = 0.75;
  CubeFile_Writer writer2(steps2, padding2, log2);

  QMState state2("ks5");
  writer2.WriteFile("test_writer3.cube", A, state2, false);

  auto result3 = Readcubefile("test_writer3.cube");

  Eigen::VectorXd values_ref3 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/cubefile_writer/values_ref3.mm");

  

  BOOST_CHECK_EQUAL(values_ref3.size(), result3.size());
  bool check_ref3 = values_ref3.isApprox(result3, 1e-4);
  BOOST_CHECK_EQUAL(check_ref3, true);
  if (!check_ref3) {
    std::cout << "ref" << std::endl;
    std::cout << values_ref3.transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << result3.transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
