/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE dftengine_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/dftengine.h>

#include "votca/xtp/orbitals.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(dftengine_test)

BOOST_AUTO_TEST_CASE(dft_full) {

  DFTEngine dft;

  std::ofstream xyzfile("molecule.xyz");
  xyzfile << "3" << std::endl;
  xyzfile << "Water molecule" << std::endl;
  xyzfile << "O          0.00000        0.00000        0.11779" << std::endl;
  xyzfile << "H          0.00000        0.75545       -0.47116" << std::endl;
  xyzfile << "H          0.00000       -0.75545       -0.47116" << std::endl;

  xyzfile.close();

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
  Orbitals orb;
  orb.QMAtoms().LoadFromFile("molecule.xyz");

  std::ofstream xml("dftengine.xml");
  xml << "<package>" << std::endl;
  xml << "<spin>1</spin>" << std::endl;
  xml << "<name>xtp</name>" << std::endl;
  xml << "<charge>0</charge>" << std::endl;
  xml << "<convergence>" << std::endl;
  xml << "    <energy>1e-7</energy>" << std::endl;
  xml << "    <method>DIIS</method>" << std::endl;
  xml << "    <DIIS_start>0.002</DIIS_start>" << std::endl;
  xml << "    <ADIIS_start>0.8</ADIIS_start>" << std::endl;
  xml << "    <DIIS_length>20</DIIS_length>" << std::endl;
  xml << "    <levelshift>0.0</levelshift>" << std::endl;
  xml << "    <levelshift_end>0.2</levelshift_end>" << std::endl;
  xml << "</convergence>" << std::endl;
  xml << "<initial_guess>independent</initial_guess>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<integration_grid_small>0</integration_grid_small>" << std::endl;
  xml << "<xc_functional>XC_HYB_GGA_XC_PBEH</xc_functional>" << std::endl;
  xml << "<max_iterations>200</max_iterations>" << std::endl;
  xml << "<read_guess>0</read_guess>" << std::endl;
  xml << "<cleanup></cleanup>" << std::endl;
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

  Eigen::MatrixXd MOs_coeff_ref = Eigen::MatrixXd::Zero(13, 13);
  MOs_coeff_ref << 0.982343, 0.225297, 3.39697e-12, 0.104226, -4.02405e-12,
      0.111297, 1.47612e-11, -1.73573e-12, -0.0611901, -1.82245e-12, 0.0515886,
      -6.55285e-12, -0.0865545, 0.101935, -0.213424, 1.90464e-12, -0.0897338,
      3.33978e-12, -0.0568715, -1.26376e-11, 1.12461e-11, 0.0976597,
      1.73492e-11, -0.131099, 1.55304e-11, 1.64108, -0.00360159, 0.116669,
      -5.86782e-13, -0.438404, 2.7592e-11, 0.249979, 3.24481e-11, -1.05194e-11,
      -0.319739, 5.38091e-11, -0.988233, 7.8109e-11, -0.147551, -2.79025e-13,
      3.45211e-12, -0.409602, 3.0704e-12, -5.80283e-12, 2.85189e-11, -0.358104,
      -0.198503, -3.7813e-12, -1.77834e-11, -8.49144e-11, -1.0425, 4.00797e-12,
      1.29673e-14, -1.00926e-12, 7.23423e-12, -2.7662e-11, -0.5242, 1.34007e-11,
      7.15316e-13, 1.60688e-11, -1.50996e-11, -1.02779, -5.12996e-11,
      1.43023e-11, 7.55262e-12, -0.0420683, -0.681791, -1.96952e-11, -0.470854,
      -8.52383e-12, -1.00918, -1.12167e-10, -1.35751e-13, 0.0433163,
      -3.21991e-13, -0.210812, 2.43573e-11, -2.00561, 0.00731426, 0.115649,
      3.63395e-12, -0.484273, 4.21012e-11, 0.461816, 3.83154e-11, -1.34171e-12,
      -0.217276, -5.47088e-11, 1.18317, -9.64031e-11, 0.501366, 2.3837e-12,
      3.59972e-12, -0.348495, 1.18667e-11, -1.82988e-11, 1.0157e-10, -0.745697,
      -0.354101, 1.19913e-11, 2.08924e-11, 1.11025e-10, 1.44212, -3.91477e-12,
      -1.66267e-12, -4.34816e-12, 4.68601e-12, -3.15464e-11, -0.629373,
      2.74539e-11, 1.77388e-11, -2.11844e-11, 1.89289e-11, 0.966973,
      5.04727e-11, -1.68919e-11, -7.27192e-12, 0.00274585, -0.126026, -0.238054,
      0.130359, -1.2134e-11, 0.0855516, 0.0897993, 0.954381, -0.965681,
      1.55138e-11, 0.314424, -0.169014, 0.291888, 0.00750313, -0.0234965,
      -0.19353, 0.118088, 4.80145e-11, 0.828874, 1.11045, -0.776492, 0.534722,
      -4.61076e-11, 0.0817261, -0.524685, 0.3665, 0.00274585, -0.126026,
      0.238054, 0.130359, -9.25038e-13, 0.0855516, -0.0897993, -0.954381,
      -0.965681, -1.21157e-11, 0.314424, 0.169014, 0.291888, 0.00750313,
      -0.0234965, 0.19353, 0.118088, 3.94754e-12, 0.828874, -1.11045, 0.776492,
      0.534722, 1.64278e-11, 0.0817261, 0.524685, 0.3665;

  AOBasis basis = orb.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  Eigen::MatrixXd proj =
      MOs_coeff_ref.transpose() * overlap.Matrix() * orb.MOs().eigenvectors();
  Eigen::VectorXd norms = proj.colwise().norm();
  bool check_coeff = norms.isApproxToConstant(1, 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
