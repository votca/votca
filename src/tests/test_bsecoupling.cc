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

#define BOOST_TEST_MODULE bsecoupling_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bsecoupling.h"
#include "votca/tools/eigenio_matrixmarket.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(bsecoupling_test)

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

BOOST_AUTO_TEST_CASE(coupling_test) {

  Orbitals A;
  A.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bsecoupling/3-21G.xml");
  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bsecoupling/molecule.xyz");
  A.setBasisSetSize(17);
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/A_MOs.mm");

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/spsi_ref.mm");

  A.BSESinglets().eigenvectors() = spsi_ref;

  Orbitals B = A;
  B.QMAtoms().Translate(4 * Eigen::Vector3d::UnitX());

  Orbitals AB;
  AB.QMAtoms() = A.QMAtoms();
  AB.QMAtoms().AddContainer(B.QMAtoms());
  AB.MOs().eigenvalues().resize(34);
  AB.MOs().eigenvalues() << -10.1341, -10.1337, -0.808607, -0.665103, -0.474928,
      -0.455857, -0.455857, -0.365971, -0.365971, -0.263259, 0.140444, 0.154745,
      0.168775, 0.168775, 0.223948, 0.231217, 0.26323, 0.26323, 0.713478,
      0.713478, 0.793559, 0.885998, 0.944915, 0.944915, 1.01169, 1.04977,
      1.04977, 1.08863, 1.10318, 1.17822, 1.18094, 1.18094, 1.69037, 1.91046;
 
  AB.setBasisSetSize(34);
  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);
  AB.setDFTbasisName(A.getDFTbasisName());
  AB.setAuxbasisName(A.getDFTbasisName());
  AB.setRPAindices(0, 33);
  AB.setBSEindices(0, 33);
  AB.setGWindices(0, 33);
  AB.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/AB_MOs.mm");

  AB.setNumberOfAlphaElectrons(10);
  AB.setNumberOfOccupiedLevels(10);

  std::ofstream opt("bsecoupling.xml");
  opt << "<bsecoupling>" << std::endl;
  opt << "        <spin>singlet</spin>" << std::endl;
  opt << "       <moleculeA>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "        </moleculeA>" << std::endl;
  opt << "        <moleculeB>" << std::endl;
  opt << "                <states>1</states>" << std::endl;
  opt << "                <occLevels>3</occLevels>" << std::endl;
  opt << "                <unoccLevels>3</unoccLevels>" << std::endl;
  opt << "         </moleculeB>" << std::endl;
  opt << "</bsecoupling>" << std::endl;
  opt.close();
  votca::tools::Property prop;
  prop.LoadFromXML("bsecoupling.xml");
  BSECoupling coup;
  Logger log;
  log.setCommonPreface("\n... ...");
  coup.setLogger(&log);

  AB.QPdiag().eigenvalues().resize(34);
  AB.QPdiag().eigenvalues() << -10.504, -10.5038, -0.923616, -0.775673,
      -0.549084, -0.530193, -0.530193, -0.430293, -0.430293, -0.322766,
      0.267681, 0.307809, 0.326961, 0.326961, 0.36078, 0.381947, 0.414845,
      0.414845, 0.906609, 0.906609, 0.993798, 1.09114, 1.14639, 1.14639, 1.1966,
      1.25629, 1.25629, 1.27991, 1.29122, 1.35945, 1.36705, 1.36705, 1.93286,
      2.11739;

  
  AB.QPdiag().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bsecoupling/Hqp.mm");

  const Eigen::MatrixXd& qpcoeff = AB.QPdiag().eigenvectors();
  Eigen::MatrixXd Hqp =
      qpcoeff * AB.QPdiag().eigenvalues().asDiagonal() * qpcoeff.transpose();
  AB.RPAInputEnergies() = Hqp.diagonal();
  coup.Initialize(prop);
  log.setReportLevel(Log::error);
  coup.CalculateCouplings(A, B, AB);
  votca::tools::Property output;
  coup.Addoutput(output, A, B);
  double diag_J_ref = 32.67651;
  double pert_J_ref = 4.434018;

  double diag_j =
      output.get("bsecoupling.singlet.coupling").getAttribute<double>("j_diag");
  double pert_j =
      output.get("bsecoupling.singlet.coupling").getAttribute<double>("j_pert");

  BOOST_CHECK_CLOSE(diag_J_ref, diag_j, 1e-4);
  BOOST_CHECK_CLOSE(pert_J_ref, pert_j, 1e-4);
}
BOOST_AUTO_TEST_SUITE_END()
