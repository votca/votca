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
  xyzfile << " 5" << std::endl;
  xyzfile << " methane" << std::endl;
  xyzfile << " C            .000000     .000000     .000000" << std::endl;
  xyzfile << " H            .629118     .629118     .629118" << std::endl;
  xyzfile << " H           -.629118    -.629118     .629118" << std::endl;
  xyzfile << " H            .629118    -.629118    -.629118" << std::endl;
  xyzfile << " H           -.629118     .629118    -.629118" << std::endl;
  xyzfile.close();

  std::ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << std::endl;
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
  basisfile << "  <element name=\"C\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << std::endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << std::endl;
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

  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), -40.24119548, 1e-5);

  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(17);
  MOs_energy_ref << -10.1468, -0.711178, -0.400278, -0.400278, -0.400278,
      0.162029, 0.207707, 0.207707, 0.207707, 0.726624, 0.726624, 0.726624,
      1.05647, 1.08935, 1.08935, 1.08935, 1.7139;
  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }

  Eigen::MatrixXd MOs_coeff_ref = Eigen::MatrixXd::Zero(17, 17);
  MOs_coeff_ref << 0.985986, -0.212012, -1.41137e-14, 1.87274e-14, 1.29873e-13,
      -0.185237, -7.29025e-14, -1.32229e-13, 6.34381e-12, -8.58718e-13,
      6.38622e-15, 1.39359e-14, 0.122724, 1.00452e-14, -7.25392e-14,
      1.38121e-11, -0.0160637, 0.108446, 0.201993, -1.2858e-14, 3.53637e-14,
      2.61866e-12, 0.0834681, 3.00917e-14, 6.30652e-14, -7.21577e-12,
      5.73411e-12, 5.1607e-14, 3.95273e-14, -0.16, 5.21885e-14, -5.65616e-14,
      -1.54422e-12, -1.99639, 2.53665e-15, 2.23849e-13, 0.0248746, -0.312978,
      -0.217885, -2.66563e-12, 0.00350206, 0.266431, -0.187946, -0.460444,
      -0.055683, -0.643566, 5.48337e-11, -0.00320464, 0.649803, -0.458766,
      -5.08742e-12, 1.77625e-14, 3.88024e-13, -0.280985, 0.132582, -0.222524,
      -2.24466e-12, -0.232532, -0.129714, -0.188215, -0.457016, 0.586622,
      0.27622, 5.4308e-11, 0.56445, -0.321387, -0.459159, -4.93364e-12,
      9.86962e-15, 2.9995e-13, 0.25783, 0.174684, -0.221488, -2.46395e-12,
      0.228557, -0.136053, -0.188609, -0.456524, -0.531093, 0.372575,
      5.44335e-11, -0.560452, -0.327395, -0.459812, -5.01028e-12, -0.0800426,
      0.626197, 7.65165e-14, -1.80841e-13, -1.44075e-11, 2.43668, 9.51464e-13,
      1.71691e-12, -7.55698e-11, -3.13924e-12, -2.10712e-13, -2.76665e-13,
      -0.142423, -1.81012e-13, 4.02015e-13, -3.13604e-11, 3.92764, -3.48211e-13,
      1.13228e-12, 0.0202501, -0.254791, -0.177377, -3.45453e-11, 0.0144232,
      1.0973, -0.774057, 0.771208, 0.0932646, 1.07792, -6.39905e-11, 0.00396605,
      -0.804196, 0.567769, 8.19778e-12, -3.57182e-13, 1.20269e-12, -0.228746,
      0.107934, -0.181154, -3.32464e-11, -0.957683, -0.534227, -0.775161,
      0.765466, -0.982546, -0.462646, -6.33112e-11, -0.698563, 0.397748,
      0.568255, 8.01569e-12, -3.52863e-13, 1.15417e-12, 0.209896, 0.142208,
      -0.18031, -3.39172e-11, 0.94131, -0.560333, -0.776784, 0.764642, 0.88954,
      -0.624034, -6.34903e-11, 0.693615, 0.405184, 0.569062, 8.11001e-12,
      -0.000891932, 0.121179, 0.000710901, -0.00236145, -0.273694, -0.0189428,
      0.000107842, -0.000151336, 0.128651, -0.60076, -6.7597e-05, 0.00228628,
      0.666903, -0.000555113, -0.000714544, 0.96403, -0.149759, 0.0154614,
      0.0188171, 0.00065559, -0.00217772, -0.252399, -0.923329, 0.00143538,
      -0.00201428, 1.71233, -0.0881061, -9.91361e-06, 0.0003353, -0.302447,
      0.000866187, 0.00111496, -1.50425, -0.778387, -0.000891932, 0.121179,
      0.0198603, -0.25647, 0.0935031, -0.0189428, -0.00170333, -0.121231,
      -0.0430248, 0.198111, -0.048626, -0.565071, 0.666903, 0.00503981,
      -0.908646, -0.322014, -0.149759, 0.0154614, 0.0188171, 0.0183151,
      -0.236516, 0.0862281, -0.923329, -0.0226713, -1.61358, -0.572659,
      0.0290545, -0.00713138, -0.082872, -0.302447, -0.00786402, 1.41783,
      0.502465, -0.778387, -0.000891932, 0.121179, 0.212514, 0.146825,
      0.0905238, -0.0189428, -0.104235, 0.0621351, -0.0427231, 0.20154,
      -0.464363, 0.323523, 0.666903, 0.784874, 0.458884, -0.320552, -0.149759,
      0.0154614, 0.0188171, 0.195979, 0.135401, 0.0834807, -0.923329, -1.38737,
      0.827016, -0.568644, 0.0295574, -0.0681024, 0.0474471, -0.302447, -1.2247,
      -0.716033, 0.500182, -0.778387, -0.000891932, 0.121179, -0.233085,
      0.112007, 0.0896668, -0.0189428, 0.105831, 0.0592473, -0.0429026,
      0.201109, 0.513057, 0.239262, 0.666903, -0.789359, 0.450476, -0.321464,
      -0.149759, 0.0154614, 0.0188171, -0.21495, 0.103292, 0.0826904, -0.923329,
      1.4086, 0.78858, -0.571033, 0.0294942, 0.0752436, 0.0350896, -0.302447,
      1.2317, -0.702914, 0.501607, -0.778387;

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
