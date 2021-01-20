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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE gaussianwriter_test

#include <fstream>
#include <sstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(gaussianwriter_test)

std::pair<std::string, std::vector<double>> readFchkFile(std::string filename,
                                                         Index skipLines = 0) {
  std::ifstream inFile(filename);
  if (inFile) {
    std::ostringstream ss;
    std::string line;
    for (int i = 0; i < skipLines; ++i) {
      std::getline(inFile, line);
    }
    while (!inFile.eof()) {
      std::getline(inFile, line);
      ss << line << "\n";
      if (line.find("Total SCF Density") != std::string::npos) {
        break;
      }
    }
    std::vector<double> density;
    double dummy;
    inFile >> dummy;
    while (!inFile.eof()) {
      density.push_back(dummy);
      inFile >> dummy;
    }
    return std::pair<std::string, std::vector<double>>(ss.str(), density);
  } else {
    throw std::runtime_error("Could not open file: " + filename + "\n");
  }
}

BOOST_AUTO_TEST_CASE(gaussianwriter_test) {
  libint2::initialize();

  // 1. Setup Orbitals Object
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gaussianwriter/methane.xyz");
  orbitals.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                           "/gaussianwriter/def2-tzvp.xml");

  orbitals.MOs().eigenvalues() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gaussianwriter/orb_energies.mm");
  orbitals.MOs().eigenvectors() =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gaussianwriter/orb_coeffs.mm");
  orbitals.setQMEnergy(-40.53758914883583);
  orbitals.setNumberOfAlphaElectrons(5);
  orbitals.setNumberOfOccupiedLevels(5);

  // 2. Write object to gaussian checkpoint file (.fchk)
  Logger log;
  GaussianWriter writer(log);
  writer.WriteFile("methane", orbitals);

  // 3.1 Check everything except density matrix
  Index linesToSkip = 2;  // basisset names will be different due to test folder
  std::pair<std::string, std::vector<double>> testFile =
      readFchkFile("methane.fchk", linesToSkip);
  std::pair<std::string, std::vector<double>> refFile = readFchkFile(
      std::string(XTP_TEST_DATA_FOLDER) + "/gaussianwriter/fchkRefMethane.fchk",
      linesToSkip);

  bool filesAreEqual = testFile.first == refFile.first;
  if (!filesAreEqual) {
    std::cout << "GENERATED FILE: " << std::endl;
    std::cout << testFile.first << std::endl;
    std::cout << "REFERENCE FILE: " << std::endl;
    std::cout << refFile.first << std::endl;
  }
  BOOST_CHECK(filesAreEqual);

  // 3.2 Check Density Matrix
  bool densityMatrixIsEqual = true;
  for (Index i = 0; i < static_cast<Index>(testFile.second.size()); ++i) {
    densityMatrixIsEqual =
        densityMatrixIsEqual &&
        (std::abs(testFile.second[i] - refFile.second[i]) < 1e-10);
  }
  if (!densityMatrixIsEqual) {
    std::cout << "GENERATED FILE vs REFERENCE: " << std::endl;
    for (Index i = 0; i < static_cast<Index>(testFile.second.size()); ++i) {
      std::cout << testFile.second[i] << "   " << refFile.second[i] << std::endl;
    }
  }
  BOOST_CHECK(densityMatrixIsEqual);

  libint2::finalize();
}
}