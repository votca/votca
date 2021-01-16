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

QMMolecule Water() {
  QMMolecule mol(" ", 1);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/espfit/molecule.xyz");
  return mol;
}

void WriteBasis321G() {
  std::ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << std::endl;
  basisfile << "  <!--Basis set created by xtp_basisset from 3-21G.nwchem at "
               "Thu Sep 15 15:40:33 2016-->"
            << std::endl;
  basisfile << "  <element name=\"H\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "  <element name=\"O\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"3.220370e+02\">" << std::endl;
  basisfile << "        <contractions factor=\"5.923940e-02\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"4.843080e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"3.515000e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.042060e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"7.076580e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"7.402940e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"-4.044530e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"2.445860e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.576200e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"1.221560e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"8.539550e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"3.736840e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();
}

BOOST_AUTO_TEST_CASE(dft_full) {
  libint2::initialize();
  DFTEngine dft;

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

  std::ofstream xml("dftengine.xml");
  xml << "<package>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<use_auxbasisset>false</use_auxbasisset>" << std::endl;
  xml << "<use_ecp>false</use_ecp>" << std::endl;
  xml << "<read_guess>0</read_guess>" << std::endl;
  xml << "<xtpdft>" << std::endl;
  xml << "<use_external_field>false</use_external_field>" << std::endl;
  xml << "<use_external_density>false</use_external_density>" << std::endl;
  xml << "<with_screening choices=\"bool\">true</with_screening>\n";
  xml << "<screening_eps  choices=\"float+\">1e-9</screening_eps>\n";
  xml << "<four_center_method>direct</four_center_method>\n";
  xml << "<fock_matrix_reset>5</fock_matrix_reset>\n";
  xml << "<convergence>" << std::endl;
  xml << "    <energy>1e-7</energy>" << std::endl;
  xml << "    <method>DIIS</method>" << std::endl;
  xml << "    <DIIS_start>0.002</DIIS_start>" << std::endl;
  xml << "    <ADIIS_start>0.8</ADIIS_start>" << std::endl;
  xml << "    <DIIS_length>20</DIIS_length>" << std::endl;
  xml << "    <levelshift>0.0</levelshift>" << std::endl;
  xml << "    <levelshift_end>0.2</levelshift_end>" << std::endl;
  xml << "    <max_iterations choices=\"int+\">100</max_iterations>\n";
  xml << "    <error choices=\"float+\">1e-7</error>\n";
  xml << "    <DIIS_maxout choices=\"bool\">false</DIIS_maxout>\n";
  xml << "    <mixing choices=\"float+\">0.7</mixing>\n";
  xml << "</convergence>" << std::endl;
  xml << "<initial_guess>independent</initial_guess>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<integration_grid_small>0</integration_grid_small>" << std::endl;
  xml << "<max_iterations>200</max_iterations>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</package>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop);
  dft.Evaluate(orb);

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -75.891017293070945, 1e-5);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(13);
  MOs_energy_ref << -19.0739, -1.01904, -0.520731, -0.341996, -0.27356,
      0.118834, 0.210783, 0.953576, 1.04314, 1.46895, 1.54729, 1.67293, 2.77584;
  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  Eigen::MatrixXd MOs_coeff_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/dftengine/MOs_coeff_ref.mm");

  AOBasis basis = orb.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd proj = MOs_coeff_ref.leftCols(5).transpose() *
                         overlap.Matrix() *
                         orb.MOs().eigenvectors().leftCols(5);
  Eigen::VectorXd norms = proj.colwise().norm();
  bool check_coeff = norms.isApproxToConstant(1, 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(density_guess) {
  libint2::initialize();
  DFTEngine dft;

  std::unique_ptr<StaticSite> s =
      std::make_unique<StaticSite>(0, "C", 3 * Eigen::Vector3d::UnitX());
  Vector9d multipoles;
  multipoles << 1.0, 0.5, 1.0, -1.0, 0.1, -0.2, 0.333, 0.1, 0.15;
  s->setMultipole(multipoles, 2);
  std::vector<std::unique_ptr<StaticSite> > multipole_vec;
  multipole_vec.push_back(std::move(s));

  dft.setExternalcharges(&multipole_vec);

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

  std::ofstream xml("dftengine.xml");
  xml << "<package>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<use_auxbasisset>false</use_auxbasisset>" << std::endl;
  xml << "<use_ecp>false</use_ecp>" << std::endl;
  xml << "<read_guess>0</read_guess>" << std::endl;
  xml << "<xtpdft>" << std::endl;
  xml << "<use_external_field>false</use_external_field>" << std::endl;
  xml << "<use_external_density>false</use_external_density>" << std::endl;
  xml << "<with_screening choices=\"bool\">true</with_screening>\n";
  xml << "<screening_eps  choices=\"float+\">1e-9</screening_eps>\n";
  xml << "<four_center_method>direct</four_center_method>\n";
  xml << "<fock_matrix_reset>5</fock_matrix_reset>\n";
  xml << "<convergence>" << std::endl;
  xml << "    <energy>1e-7</energy>" << std::endl;
  xml << "    <method>DIIS</method>" << std::endl;
  xml << "    <DIIS_start>0.002</DIIS_start>" << std::endl;
  xml << "    <ADIIS_start>0.8</ADIIS_start>" << std::endl;
  xml << "    <DIIS_length>20</DIIS_length>" << std::endl;
  xml << "    <levelshift>0.0</levelshift>" << std::endl;
  xml << "    <levelshift_end>0.2</levelshift_end>" << std::endl;
  xml << "    <max_iterations choices=\"int+\">100</max_iterations>\n";
  xml << "    <error choices=\"float+\">1e-7</error>\n";
  xml << "    <DIIS_maxout choices=\"bool\">false</DIIS_maxout>\n";
  xml << "    <mixing choices=\"float+\">0.7</mixing>\n";
  xml << "</convergence>" << std::endl;
  xml << "<initial_guess>atom</initial_guess>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<integration_grid_small>0</integration_grid_small>" << std::endl;
  xml << "<max_iterations>1</max_iterations>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</package>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop);
  dft.Evaluate(orb);

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -75.891684954029387, 1e-5);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(13);
  MOs_energy_ref << -19.3481, -1.30585, -0.789203, -0.59822, -0.555272,
      -0.150066, -0.0346099, 0.687671, 0.766599, 1.17942, 1.28947, 1.41871,
      2.49675;

  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  Eigen::MatrixXd MOs_coeff_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/dftengine/MOs_coeff_ref2.mm");

  AOBasis basis = orb.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd proj = MOs_coeff_ref.leftCols(5).transpose() *
                         overlap.Matrix() *
                         orb.MOs().eigenvectors().leftCols(5);
  Eigen::VectorXd norms = proj.colwise().norm();
  bool check_coeff = norms.isApproxToConstant(1, 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
