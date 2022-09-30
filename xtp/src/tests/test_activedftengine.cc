

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

#define BOOST_TEST_MODULE dftengine_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(dftengine_test)

archive_file =
    std::string(XTP_TEST_DATA_FOLDER) + "/activedftengine/molecule.orb";

BOOST_AUTO_TEST_CASE(dft_active) {
  libint2::initialize();
  DFTEngine activedft;
  Orbitals orb;
  orb.ReadFromCpt(archive_file) orb.QMAtoms().LoadFromFile(
      std::string(XTP_TEST_DATA_FOLDER) + "/activedftengine/molecule.xyz");
  std::ofstream xml("dftengine.xml");
  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<xtpdft>" << std::endl;
  xml << "<screening_eps>1e-9</screening_eps>\n";
  xml << "<fock_matrix_reset>5</fock_matrix_reset>\n";
  xml << "<convergence>" << std::endl;
  xml << "    <energy>1e-7</energy>" << std::endl;
  xml << "    <method>DIIS</method>" << std::endl;
  xml << "    <DIIS_start>0.002</DIIS_start>" << std::endl;
  xml << "    <ADIIS_start>0.8</ADIIS_start>" << std::endl;
  xml << "    <DIIS_length>20</DIIS_length>" << std::endl;
  xml << "    <levelshift>0.0</levelshift>" << std::endl;
  xml << "    <levelshift_end>0.2</levelshift_end>" << std::endl;
  xml << "    <max_iterations>100</max_iterations>\n";
  xml << "    <error>1e-7</error>\n";
  xml << "    <DIIS_maxout>false</DIIS_maxout>\n";
  xml << "    <mixing>0.7</mixing>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<max_iterations>100</max_iterations>" << std::endl;
  xml << "<dft_in_dft>" << std::endl;
  xml << "    <activeatoms>0</activeatoms>" << std::endl;
  xml << "    <truncate_basis>False</truncate_basis>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");

  Logger log;
  activedft.setLogger(&log);
  activedft.Initialize(prop.get("dftpackage"));
  activedft.EvaluateActiveRegion(orb);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(13);
  mo_eigenvalues << -19.0797217317, -1.0193288277, -0.5207045448, -0.3421671503,
      -0.2737294649, 0.1190570943, 0.2108430632, 0.9536993123, 1.0433048316,
      1.4681715152, 1.5465977240, 1.6726054535, 2.7743814447;
  // Eigen::MatrixXd MOs_coeff_ref =
  //     votca::tools::EigenIO_MatrixMarket::ReadMatrix(
  //         std::string(XTP_TEST_DATA_FOLDER) + "/activedftengine/mo_eigenvectors.mm");
  bool check_eng = MOs_energy_ref.isApprox(orb.getEmbeddedMOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result energy" << std::endl;
    std::cout << orb.getEmbeddedMOs().eigenvalues() << std::endl;
    std::cout << "ref energy" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  libint2::finalize();
}
