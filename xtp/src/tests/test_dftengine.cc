

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

#define BOOST_TEST_MODULE dftengine_test

// Third party includes
#include <boost/test/unit_test.hpp>
#include <sstream>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using votca::Index;

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

  std::ofstream xml("dftengine2.xml");
  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<initial_guess>independent</initial_guess>" << std::endl;
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<max_iterations>200</max_iterations>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine2.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
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

  AOBasis basis = orb.getDftBasis();
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

  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<initial_guess>atom</initial_guess>" << std::endl;
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<max_iterations>1</max_iterations>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
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

  AOBasis basis = orb.getDftBasis();
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

BOOST_AUTO_TEST_CASE(huckel_guess) {
  libint2::initialize();
  DFTEngine dft;

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

  std::ofstream xml("dftengine_huckel.xml");

  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<initial_guess>huckel</initial_guess>" << std::endl;
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<max_iterations>1</max_iterations>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine_huckel.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  dft.Evaluate(orb);

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -75.89101729307095, 1e-5);

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
          std::string(XTP_TEST_DATA_FOLDER) + "/dftengine/MOs_coeff_ref3.mm");

  AOBasis basis = orb.getDftBasis();
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

BOOST_AUTO_TEST_CASE(huckel_dft_guess) {
  libint2::initialize();
  DFTEngine dft;

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

  std::ofstream xml("dftengine_huckel_dft.xml");

  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<initial_guess>huckel_dft</initial_guess>" << std::endl;
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine_huckel_dft.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  dft.Evaluate(orb);

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -75.89101729307095, 1e-5);

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
          std::string(XTP_TEST_DATA_FOLDER) + "/dftengine/MOs_coeff_ref4.mm");

  AOBasis basis = orb.getDftBasis();
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

BOOST_AUTO_TEST_CASE(dft_cation) {
  libint2::initialize();
  DFTEngine dft;

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

  std::ofstream xml("dftengine_cation.xml");

  xml << "<dftpackage>" << std::endl;
  xml << "<spin>2</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>1</charge>" << std::endl;
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<initial_guess>huckel_dft</initial_guess>" << std::endl;
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine_cation.xml");

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  orb.setChargeAndSpin(1, 2);
  dft.Evaluate(orb);

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -75.456593793657646, 1e-5);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(13);
  MOs_energy_ref << -19.6958, -1.58537, -1.03212, -0.899489, -0.873537,
      -0.255014, -0.169427, 0.548747, 0.609376, 0.899063, 1.05047, 1.18105,
      2.23438;

  Eigen::VectorXd MOs_energy_ref_beta = Eigen::VectorXd::Zero(13);
  MOs_energy_ref_beta << -19.668, -1.50742, -1.00909, -0.837772, -0.588854,
      -0.242175, -0.161593, 0.547669, 0.616913, 1.02967, 1.07354, 1.2038,
      2.29571;

  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;

    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  bool check_eng_beta =
      MOs_energy_ref_beta.isApprox(orb.MOs_beta().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng_beta, true);
  if (!check_eng_beta) {
    std::cout << "result eng beta" << std::endl;
    std::cout << orb.MOs_beta().eigenvalues() << std::endl;

    std::cout << "ref eng beta" << std::endl;
    std::cout << MOs_energy_ref_beta << std::endl;
  }

  Eigen::MatrixXd MOs_coeff_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/dftengine/MOs_coeff_cation_alpha.mm");

  Eigen::MatrixXd MOs_coeff_beta_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/dftengine/MOs_coeff_cation_beta.mm");

  AOBasis basis = orb.getDftBasis();
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

  Eigen::MatrixXd proj_beta = MOs_coeff_beta_ref.leftCols(5).transpose() *
                              overlap.Matrix() *
                              orb.MOs_beta().eigenvectors().leftCols(5);
  Eigen::VectorXd norms_beta = proj_beta.colwise().norm();
  bool check_coeff_beta = norms_beta.isApproxToConstant(1, 1e-5);
  BOOST_CHECK_EQUAL(check_coeff_beta, true);
  if (!check_coeff_beta) {
    std::cout << "result coeff beta" << std::endl;
    std::cout << orb.MOs_beta().eigenvectors() << std::endl;
    std::cout << "ref coeff beta" << std::endl;
    std::cout << MOs_coeff_beta_ref << std::endl;
  }

  libint2::finalize();
}

namespace {
void WriteDimerGuessXML(const std::string& filename, const std::string& guess,
                        const std::string& extra) {
  std::ofstream xml(filename);
  xml << "<dftpackage>\n";
  xml << "<spin>1</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>0</charge>\n";
  xml << "<functional>XC_HYB_GGA_XC_PBEH</functional>\n";
  xml << "<basisset>3-21G.xml</basisset>\n";
  xml << "<initial_guess>" << guess << "</initial_guess>\n";
  xml << extra;
  xml << "<xtpdft>\n";
  xml << "<screening_eps>1e-9</screening_eps>\n";
  xml << "<fock_matrix_reset>5</fock_matrix_reset>\n";
  xml << "<convergence>\n";
  xml << "    <energy>1e-8</energy>\n";
  xml << "    <method>DIIS</method>\n";
  xml << "    <DIIS_start>0.002</DIIS_start>\n";
  xml << "    <ADIIS_start>0.8</ADIIS_start>\n";
  xml << "    <DIIS_length>20</DIIS_length>\n";
  xml << "    <levelshift>0.0</levelshift>\n";
  xml << "    <levelshift_end>0.2</levelshift_end>\n";
  xml << "    <max_iterations>100</max_iterations>\n";
  xml << "    <error>1e-7</error>\n";
  xml << "    <DIIS_maxout>false</DIIS_maxout>\n";
  xml << "    <mixing>0.7</mixing>\n";
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>\n";
  xml << "<integration_grid>xcoarse</integration_grid>\n";
  xml << "</xtpdft>\n";
  xml << "</dftpackage>\n";
}

QMMolecule WaterDimer(const Eigen::Vector3d& shift) {
  QMMolecule dimer(" ", 1);
  QMMolecule a = Water();
  QMMolecule b = Water();
  b.Translate(shift);
  Index id = 0;
  for (const QMAtom& at : a) {
    dimer.push_back(QMAtom(id++, at.getElement(), at.getPos()));
  }
  for (const QMAtom& at : b) {
    dimer.push_back(QMAtom(id++, at.getElement(), at.getPos()));
  }
  return dimer;
}

Index SCFIterations(const std::string& log) {
  // the SCF loop logs one "Iteration" line per step
  Index count = 0;
  std::size_t pos = 0;
  while ((pos = log.find(" Iteration ", pos)) != std::string::npos) {
    ++count;
    ++pos;
  }
  return count;
}
}  // namespace

// Closed-shell dimer guess for a restricted run: two converged water monomers
// placed 10 bohr apart. The guess must converge to the same state as a
// standard guess, and in fewer iterations.
BOOST_AUTO_TEST_CASE(dimer_guess_closed_shell) {
  libint2::initialize();
  WriteBasis321G();

  {
    Orbitals mono;
    mono.QMAtoms() = Water();
    WriteDimerGuessXML("dftengine_monomer.xml", "atom", "");
    votca::tools::Property prop;
    prop.LoadFromXML("dftengine_monomer.xml");
    Logger log;
    DFTEngine dft;
    dft.setLogger(&log);
    dft.Initialize(prop.get("dftpackage"));
    BOOST_REQUIRE(dft.Evaluate(mono));
    mono.WriteToCpt("water_monomer.orb");
  }

  const Eigen::Vector3d shift(10.0, 0.0, 0.0);

  auto RunDimer = [&](const std::string& guess, const std::string& extra,
                      Index& iterations) {
    Orbitals orb;
    orb.QMAtoms() = WaterDimer(shift);
    WriteDimerGuessXML("dftengine_dimer.xml", guess, extra);
    votca::tools::Property prop;
    prop.LoadFromXML("dftengine_dimer.xml");
    Logger log;
    log.setReportLevel(votca::Log::info);
    DFTEngine dft;
    dft.setLogger(&log);
    dft.Initialize(prop.get("dftpackage"));
    BOOST_REQUIRE(dft.Evaluate(orb));
    std::stringstream ss;
    ss << log;
    iterations = SCFIterations(ss.str());
    return orb.getDFTTotalEnergy();
  };

  Index it_atom = 0;
  Index it_dimer = 0;
  double e_atom = RunDimer("atom", "", it_atom);
  double e_dimer = RunDimer("dimer_guess",
                            "<dimer_guess_orbA>water_monomer.orb</"
                            "dimer_guess_orbA>\n<dimer_guess_orbB>water_"
                            "monomer.orb</dimer_guess_orbB>\n",
                            it_dimer);
  BOOST_CHECK_SMALL(e_dimer - e_atom, 1e-6);
  BOOST_TEST_MESSAGE("SCF iterations: atom guess "
                     << it_atom << ", dimer guess " << it_dimer);
  BOOST_CHECK_LT(it_dimer, it_atom);

  libint2::finalize();
}

// QM/MM warm start: the converged orbitals handed back in reach the same
// energy in a few iterations, for closed and open shells; orbitals that do not
// fit (here: none at all) fall back to the configured guess.
BOOST_AUTO_TEST_CASE(warm_start_from_previous_orbitals) {
  libint2::initialize();
  WriteBasis321G();

  for (const std::string charge_spin : {"0 1", "1 2"}) {
    std::istringstream cs(charge_spin);
    int charge, spin;
    cs >> charge >> spin;
    WriteDimerGuessXML("dftengine_warm.xml", "atom", "");
    votca::tools::Property prop;
    prop.LoadFromXML("dftengine_warm.xml");
    prop.set("dftpackage.charge", std::to_string(charge));
    prop.set("dftpackage.spin", std::to_string(spin));

    auto Run = [&](Orbitals& orb, bool warm, Index& iterations,
                   std::string& logtext) {
      Logger log;
      log.setReportLevel(votca::Log::info);
      DFTEngine dft;
      dft.setLogger(&log);
      dft.Initialize(prop.get("dftpackage"));
      dft.setWarmStart(warm);
      BOOST_REQUIRE(dft.Evaluate(orb));
      std::stringstream ss;
      ss << log;
      logtext = ss.str();
      iterations = SCFIterations(logtext);
      return orb.getDFTTotalEnergy();
    };

    Orbitals orb;
    orb.QMAtoms() = Water();
    orb.setChargeAndSpin(charge, spin);
    Index it_cold = 0;
    std::string log_cold;
    const double e_cold = Run(orb, false, it_cold, log_cold);

    Index it_warm = 0;
    std::string log_warm;
    const double e_warm = Run(orb, true, it_warm, log_warm);
    BOOST_CHECK(log_warm.find("Starting from the orbitals of the previous") !=
                std::string::npos);
    BOOST_CHECK_SMALL(e_warm - e_cold, 1e-7);
    BOOST_CHECK_LE(it_warm, 3);
    BOOST_TEST_MESSAGE("charge " << charge << ": cold " << it_cold
                                 << " iterations, warm " << it_warm);

    Orbitals fresh;
    fresh.QMAtoms() = Water();
    fresh.setChargeAndSpin(charge, spin);
    Index it_fresh = 0;
    std::string log_fresh;
    const double e_fresh = Run(fresh, true, it_fresh, log_fresh);
    BOOST_CHECK(log_fresh.find("not usable as guess (no MOs)") !=
                std::string::npos);
    BOOST_CHECK(log_fresh.find("Starting from the orbitals of the previous") ==
                std::string::npos);
    BOOST_CHECK_SMALL(e_fresh - e_cold, 1e-7);
    // The iteration counts of the cold and the fallback run are not compared:
    // with several threads, summation order differs from run to run and can
    // change the SCF path (seen for the cation: 10 or 16 iterations).
    BOOST_TEST_MESSAGE("charge " << charge << ": fallback " << it_fresh
                                 << " iterations");
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
