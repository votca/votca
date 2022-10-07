

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

#define BOOST_TEST_MODULE activedftengine_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(activedftengine_test)

void WriteBasisSVP() {
  std::ofstream basisfile("def2-svp.xml");
  basisfile << "<basis name=\"def2-svp\">" << std::endl;
  basisfile << "  <!--Basis set created by xtp_basisset from def2-svp.nw at "
               "Thu Jan  2 17:19:55 2020-->"
            << std::endl;
  basisfile << "  <element name=\"H\">" << std::endl;
  basisfile << "    <shell type=\"S\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"1.301070e+01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"1.968216e-02\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.962257e+00\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"1.379652e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"4.445380e-01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"4.783193e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"S\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"1.219496e-01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"P\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"8.000000e-01\">" << std::endl;
  basisfile << "        <contractions type=\"P\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "  <element name=\"O\">" << std::endl;
  basisfile << "    <shell type=\"S\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"2.266177e+03\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"-5.343181e-03\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"3.408701e+02\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"-3.989004e-02\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"7.736314e+01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"-1.785391e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"2.147964e+01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"-4.642768e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"6.658943e+00\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"-4.430975e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"S\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"8.097598e-01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"S\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"2.553077e-01\">" << std::endl;
  basisfile << "        <contractions type=\"S\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"P\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"1.772150e+01\">" << std::endl;
  basisfile << "        <contractions type=\"P\" factor=\"4.339457e-02\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"3.863551e+00\">" << std::endl;
  basisfile << "        <contractions type=\"P\" factor=\"2.309412e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.048092e+00\">" << std::endl;
  basisfile << "        <contractions type=\"P\" factor=\"5.137531e-01\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"P\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"2.764154e-01\">" << std::endl;
  basisfile << "        <contractions type=\"P\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell type=\"D\" scale=\"1.0\">" << std::endl;
  basisfile << "      <constant decay=\"1.200000e+00\">" << std::endl;
  basisfile << "        <contractions type=\"D\" factor=\"1.000000e+00\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();
}

BOOST_AUTO_TEST_CASE(dft_active) {
  libint2::initialize();
  DFTEngine activedft;
  Orbitals orb;

  std::string archive_file =
      std::string(XTP_TEST_DATA_FOLDER) + "/activedftengine/molecule.orb";
  orb.ReadFromCpt(archive_file);
  WriteBasisSVP();

  std::ofstream xml("dftengine.xml");
  xml << "<dftpackage>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<basisset>def2-svp.xml</basisset>" << std::endl;
  xml << "<initial_guess>atom</initial_guess>" << std::endl;
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
  xml << "    <threshold>0.4</threshold>" << std::endl;
  xml << "    <levelshift>10000</levelshift>" << std::endl;
  xml << "    <truncate_basis>False</truncate_basis>" << std::endl;
  xml << "    <truncation_threshold>1e-4</truncation_threshold>" << std::endl;
  xml << "</dft_in_dft>" << std::endl;
  xml << "</xtpdft>" << std::endl;
  xml << "</dftpackage>" << std::endl;
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine.xml");

  Logger log;
  activedft.setLogger(&log);
  activedft.Initialize(prop.get("dftpackage"));
  activedft.EvaluateActiveRegion(orb);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(24);
  MOs_energy_ref << -19.1916, -1.00827, -0.52987, -0.381083, -0.305608,
      0.0651889, 0.144076, 0.581983, 0.644849, 0.940892, 0.947913, 1.0306,
      1.11042, 1.35566, 1.41993, 1.55998, 1.80324, 2.2424, 2.28424, 2.97774,
      3.02425, 3.20965, 3.5225, 3.84776;

  bool check_eng =
      MOs_energy_ref.isApprox(orb.getEmbeddedMOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result energy" << std::endl;
    std::cout << orb.getEmbeddedMOs().eigenvalues() << std::endl;
    std::cout << "ref energy" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  Eigen::MatrixXd MOs_coeff_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/activedftengine/mo_eigenvectors.mm");
  AOBasis basis = orb.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd proj = MOs_coeff_ref.transpose() * overlap.Matrix() *
                         orb.getEmbeddedMOs().eigenvectors();
  Eigen::VectorXd norms = proj.colwise().norm();
  bool check_coeff = norms.isApproxToConstant(1, 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.getEmbeddedMOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
