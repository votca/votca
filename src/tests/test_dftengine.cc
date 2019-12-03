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

QMMolecule Water() {
  std::ofstream xyzfile("molecule.xyz");
  xyzfile << "3" << std::endl;
  xyzfile << "Water molecule" << std::endl;
  xyzfile << "O          0.00000        0.00000        0.11779" << std::endl;
  xyzfile << "H          0.00000        0.75545       -0.47116" << std::endl;
  xyzfile << "H          0.00000       -0.75545       -0.47116" << std::endl;

  xyzfile.close();
  QMMolecule mol(" ", 1);
  mol.LoadFromFile("molecule.xyz");
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

  DFTEngine dft;

  WriteBasis321G();

  Orbitals orb;
  orb.QMAtoms() = Water();

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
}

BOOST_AUTO_TEST_CASE(density_guess) {

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
  xml << "<convergence>" << std::endl;
  xml << "    <energy>1e-7</energy>" << std::endl;
  xml << "    <method>DIIS</method>" << std::endl;
  xml << "    <DIIS_start>0.002</DIIS_start>" << std::endl;
  xml << "    <ADIIS_start>0.8</ADIIS_start>" << std::endl;
  xml << "    <DIIS_length>20</DIIS_length>" << std::endl;
  xml << "    <levelshift>0.0</levelshift>" << std::endl;
  xml << "    <levelshift_end>0.2</levelshift_end>" << std::endl;
  xml << "</convergence>" << std::endl;
  xml << "<initial_guess>atom</initial_guess>" << std::endl;
  xml << "<basisset>3-21G.xml</basisset>" << std::endl;
  xml << "<integration_grid>xcoarse</integration_grid>" << std::endl;
  xml << "<integration_grid_small>0</integration_grid_small>" << std::endl;
  xml << "<xc_functional>XC_HYB_GGA_XC_PBEH</xc_functional>" << std::endl;
  xml << "<max_iterations>1</max_iterations>" << std::endl;
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

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -78.611399705809276, 1e-5);

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

  Eigen::MatrixXd MOs_coeff_ref = Eigen::MatrixXd::Zero(13, 13);
  MOs_coeff_ref << 0.982347, -0.223306, 0.011141, 0.109466, 0.0234867, 0.105295,
      0.0232115, -0.0182764, 0.0582594, 0.00059232, 0.0489508, 0.000786996,
      -0.0873035, 0.101932, 0.207993, -0.0126581, -0.0990048, -0.0340845,
      -0.0557366, -0.0168192, 0.021685, -0.085871, -0.0512844, -0.107919,
      0.0169856, 1.64211, -0.00381954, -0.123607, 0.0114843, -0.435627,
      -0.0189638, 0.263385, 0.0440122, -0.0986781, 0.279353, 0.139172,
      -0.982073, -0.0532888, -0.132683, 0.000164598, 0.0116932, 0.41142,
      -0.00374944, -0.0241652, -0.079091, 0.358415, 0.192504, 0.0608463,
      -0.071234, -0.0658821, 1.03301, -0.014661, 0.000466344, 0.0212544,
      0.0156916, -0.0476611, 0.513555, -0.0376232, -0.00115035, -0.0503564,
      0.0669386, -1.01848, -0.113743, -0.0696086, -0.0305003, -0.0420924,
      0.680535, -0.0444622, -0.488698, -0.0832719, -0.967803, -0.194593,
      0.03382, -0.0554861, 0.114266, -0.237959, -0.0390543, -2.00325,
      0.00740548, -0.136959, 0.0167012, -0.471414, -0.0393881, 0.445018,
      0.0765772, -0.0772324, 0.231882, -0.149126, 1.18178, 0.0723819, 0.478636,
      -8.34239e-05, 0.0176715, 0.346204, -0.0199012, -0.00769652, -0.137677,
      0.723167, 0.328144, 0.0914483, 0.0851728, 0.086366, -1.44522, 0.0255437,
      -0.000170684, 0.0465323, 0.0317499, -0.0465383, 0.631769, 0.00435838,
      -0.000600695, 0.0289997, -0.0399768, 0.950868, 0.115385, 0.0724111,
      0.0388676, 0.00276614, 0.125487, 0.239726, 0.138419, -0.0431558,
      0.0732387, -0.079736, -1.19554, 0.648011, 0.0220574, 0.299091, 0.195476,
      0.278534, 0.00753585, 0.0323646, 0.202131, 0.147493, 0.0411125, 1.02626,
      -0.920028, 0.896255, -0.261932, -0.13503, 0.0604202, 0.536989, 0.36072,
      0.00275706, 0.113321, -0.245446, 0.117288, -0.00329392, 0.0765362,
      0.124891, 0.620243, 1.21143, 0.00508303, 0.306816, -0.154032, 0.287174,
      0.00747457, 0.016969, -0.168657, 0.108007, 0.0371845, 0.594021, 1.23285,
      -0.611537, -0.749223, -0.0227042, 0.139562, -0.513734, 0.378741;

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
}

BOOST_AUTO_TEST_SUITE_END()
