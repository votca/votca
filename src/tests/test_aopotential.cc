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

#define BOOST_TEST_MODULE aopotential_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/orbitals.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(aopotential_test)

BOOST_AUTO_TEST_CASE(aopotentials_test) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  AOMultipole esp;
  esp.FillPotential(aobasis, orbitals.QMAtoms());
  Eigen::MatrixXd esp_ref = Eigen::MatrixXd::Zero(17, 17);

  esp_ref << -36.3803, -2.68098, 0, 0, 0, -3.82378, 0, 0, 0, -0.366604,
      -1.69872, -0.366604, -1.69872, -0.366604, -1.69872, -0.366604, -1.69872,
      -2.68098, -8.18472, 0, 0, 0, -5.77391, 0, 0, 0, -1.3464, -2.94621,
      -1.3464, -2.94621, -1.3464, -2.94621, -1.3464, -2.94621, 0, 0, -8.70705,
      0, 0, 0, -3.49193, 0, 0, -1.11814, -0.876661, -1.11814, -0.876661,
      1.11814, 0.876661, 1.11814, 0.876661, 0, 0, 0, -8.70705, 0, 0, 0,
      -3.49193, 0, -1.11814, -0.876661, 1.11814, 0.876661, 1.11814, 0.876661,
      -1.11814, -0.876661, 0, 0, 0, 0, -8.70705, 0, 0, 0, -3.49193, -1.11814,
      -0.876661, 1.11814, 0.876661, -1.11814, -0.876661, 1.11814, 0.876661,
      -3.82378, -5.77391, 0, 0, 0, -6.04705, 0, 0, 0, -1.91545, -3.64247,
      -1.91545, -3.64247, -1.91545, -3.64247, -1.91545, -3.64247, 0, 0,
      -3.49193, 0, 0, 0, -4.45574, 0, 0, -1.49469, -1.4276, -1.49469, -1.4276,
      1.49469, 1.4276, 1.49469, 1.4276, 0, 0, 0, -3.49193, 0, 0, 0, -4.45574, 0,
      -1.49469, -1.4276, 1.49469, 1.4276, 1.49469, 1.4276, -1.49469, -1.4276, 0,
      0, 0, 0, -3.49193, 0, 0, 0, -4.45574, -1.49469, -1.4276, 1.49469, 1.4276,
      -1.49469, -1.4276, 1.49469, 1.4276, -0.366604, -1.3464, -1.11814,
      -1.11814, -1.11814, -1.91545, -1.49469, -1.49469, -1.49469, -5.52644,
      -3.23475, -0.0526377, -0.668246, -0.0526377, -0.668246, -0.0526377,
      -0.668246, -1.69872, -2.94621, -0.876661, -0.876661, -0.876661, -3.64247,
      -1.4276, -1.4276, -1.4276, -3.23475, -4.25826, -0.668246, -1.83787,
      -0.668246, -1.83787, -0.668246, -1.83787, -0.366604, -1.3464, -1.11814,
      1.11814, 1.11814, -1.91545, -1.49469, 1.49469, 1.49469, -0.0526377,
      -0.668246, -5.52644, -3.23475, -0.0526377, -0.668246, -0.0526377,
      -0.668246, -1.69872, -2.94621, -0.876661, 0.876661, 0.876661, -3.64247,
      -1.4276, 1.4276, 1.4276, -0.668246, -1.83787, -3.23475, -4.25826,
      -0.668246, -1.83787, -0.668246, -1.83787, -0.366604, -1.3464, 1.11814,
      1.11814, -1.11814, -1.91545, 1.49469, 1.49469, -1.49469, -0.0526377,
      -0.668246, -0.0526377, -0.668246, -5.52644, -3.23475, -0.0526377,
      -0.668246, -1.69872, -2.94621, 0.876661, 0.876661, -0.876661, -3.64247,
      1.4276, 1.4276, -1.4276, -0.668246, -1.83787, -0.668246, -1.83787,
      -3.23475, -4.25826, -0.668246, -1.83787, -0.366604, -1.3464, 1.11814,
      -1.11814, 1.11814, -1.91545, 1.49469, -1.49469, 1.49469, -0.0526377,
      -0.668246, -0.0526377, -0.668246, -0.0526377, -0.668246, -5.52644,
      -3.23475, -1.69872, -2.94621, 0.876661, -0.876661, 0.876661, -3.64247,
      1.4276, -1.4276, 1.4276, -0.668246, -1.83787, -0.668246, -1.83787,
      -0.668246, -1.83787, -3.23475, -4.25826;

  bool check_esp = esp.Matrix().isApprox(esp_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_esp, 1);
  if (!check_esp) {
    std::cout << "esp Ref" << endl;
    std::cout << esp_ref << endl;
    std::cout << "esp" << endl;
    std::cout << esp.Matrix() << endl;
  }

  ofstream ecpfile("ecp.xml");
  ecpfile << "<pseudopotential name=\"ECP_STUTTGART\">" << endl;
  ecpfile << "  <element name=\"C\" lmax=\"3\" ncore=\"2\">" << endl;
  ecpfile << "    <shell type=\"F\"><constant power=\"2\" decay=\"1.0\" "
             "contraction=\"0.0\"></constant></shell>"
          << endl;
  ecpfile << "    <shell type=\"S\"><constant power=\"2\" decay=\"6.40105200\" "
             "contraction=\"33.12163800\"></constant></shell>"
          << endl;
  ecpfile << "    <shell type=\"P\"><constant power=\"2\" decay=\"7.30774700\" "
             "contraction=\"-1.98625700\"></constant></shell>"
          << endl;
  ecpfile << "    <shell type=\"D\"><constant power=\"2\" decay=\"5.96179600\" "
             "contraction=\"-9.45431800\"></constant></shell>"
          << endl;
  ecpfile << "  </element>" << endl;
  ecpfile << "</pseudopotential>" << endl;
  ecpfile.close();
  BasisSet ecps;
  ecps.LoadPseudopotentialSet("ecp.xml");
  AOBasis ecpbasis;
  ecpbasis.ECPFill(ecps, orbitals.QMAtoms());

  AOECP ecp;
  ecp.FillPotential(aobasis, ecpbasis);
  Eigen::MatrixXd ecp_ref = Eigen::MatrixXd::Zero(17, 17);
  ecp_ref << 21.6188, 1.34835, 0, 0, 0, 2.29744, 0, 0, 0, 0.209711, 1.01592,
      0.209711, 1.01592, 0.209711, 1.01592, 0.209711, 1.01592, 1.34835,
      0.702249, 0, 0, 0, 0.4993, 0, 0, 0, 0.0564639, 0.225665, 0.0564639,
      0.225665, 0.0564639, 0.225665, 0.0564639, 0.225665, 0, 0, -0.0737545, 0,
      0, 0, -0.00882987, 0, 0, -0.00178626, -0.00193605, -0.00178626,
      -0.00193605, 0.00178626, 0.00193605, 0.00178626, 0.00193605, 0, 0, 0,
      -0.0737545, 0, 0, 0, -0.00882987, 0, -0.00178626, -0.00193605, 0.00178626,
      0.00193605, 0.00178626, 0.00193605, -0.00178626, -0.00193605, 0, 0, 0, 0,
      -0.0737545, 0, 0, 0, -0.00882987, -0.00178626, -0.00193605, 0.00178626,
      0.00193605, -0.00178626, -0.00193605, 0.00178626, 0.00193605, 2.29744,
      0.4993, 0, 0, 0, 0.458665, 0, 0, 0, 0.0477375, 0.20545, 0.0477375,
      0.20545, 0.0477375, 0.20545, 0.0477375, 0.20545, 0, 0, -0.00882987, 0, 0,
      0, -0.0011596, 0, 0, -0.000240513, -0.000255319, -0.000240513,
      -0.000255319, 0.000240513, 0.000255319, 0.000240513, 0.000255319, 0, 0, 0,
      -0.00882987, 0, 0, 0, -0.0011596, 0, -0.000240513, -0.000255319,
      0.000240513, 0.000255319, 0.000240513, 0.000255319, -0.000240513,
      -0.000255319, 0, 0, 0, 0, -0.00882987, 0, 0, 0, -0.0011596, -0.000240513,
      -0.000255319, 0.000240513, 0.000255319, -0.000240513, -0.000255319,
      0.000240513, 0.000255319, 0.209711, 0.0564639, -0.00178626, -0.00178626,
      -0.00178626, 0.0477375, -0.000240513, -0.000240513, -0.000240513,
      0.00468574, 0.0212243, 0.0052935, 0.0215396, 0.0052935, 0.0215396,
      0.0052935, 0.0215396, 1.01592, 0.225665, -0.00193605, -0.00193605,
      -0.00193605, 0.20545, -0.000255319, -0.000255319, -0.000255319, 0.0212243,
      0.0918741, 0.0215396, 0.0921252, 0.0215396, 0.0921252, 0.0215396,
      0.0921252, 0.209711, 0.0564639, -0.00178626, 0.00178626, 0.00178626,
      0.0477375, -0.000240513, 0.000240513, 0.000240513, 0.0052935, 0.0215396,
      0.00468574, 0.0212243, 0.0052935, 0.0215396, 0.0052935, 0.0215396,
      1.01592, 0.225665, -0.00193605, 0.00193605, 0.00193605, 0.20545,
      -0.000255319, 0.000255319, 0.000255319, 0.0215396, 0.0921252, 0.0212243,
      0.0918741, 0.0215396, 0.0921252, 0.0215396, 0.0921252, 0.209711,
      0.0564639, 0.00178626, 0.00178626, -0.00178626, 0.0477375, 0.000240513,
      0.000240513, -0.000240513, 0.0052935, 0.0215396, 0.0052935, 0.0215396,
      0.00468574, 0.0212243, 0.0052935, 0.0215396, 1.01592, 0.225665,
      0.00193605, 0.00193605, -0.00193605, 0.20545, 0.000255319, 0.000255319,
      -0.000255319, 0.0215396, 0.0921252, 0.0215396, 0.0921252, 0.0212243,
      0.0918741, 0.0215396, 0.0921252, 0.209711, 0.0564639, 0.00178626,
      -0.00178626, 0.00178626, 0.0477375, 0.000240513, -0.000240513,
      0.000240513, 0.0052935, 0.0215396, 0.0052935, 0.0215396, 0.0052935,
      0.0215396, 0.00468574, 0.0212243, 1.01592, 0.225665, 0.00193605,
      -0.00193605, 0.00193605, 0.20545, 0.000255319, -0.000255319, 0.000255319,
      0.0215396, 0.0921252, 0.0215396, 0.0921252, 0.0215396, 0.0921252,
      0.0212243, 0.0918741;

  bool check_ecp = ecp.Matrix().isApprox(ecp_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_ecp, 1);
  if (!check_ecp) {
    std::cout << "ecp Ref" << endl;
    std::cout << ecp_ref << endl;
    std::cout << "ecp" << endl;
    std::cout << ecp.Matrix() << endl;
  }

  ofstream mpsfile("polarsite.mps");
  mpsfile << "! One Site" << endl;
  mpsfile << "! N=1 " << endl;
  mpsfile << "Units angstrom" << endl;
  mpsfile << "  C +0 0 3 Rank 1" << endl;
  mpsfile << "+0" << endl;
  mpsfile << "10 0 0" << endl;
  mpsfile << "0 0 0 0 0" << endl;
  mpsfile
      << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 +1.9445387 "
      << endl;
  mpsfile.close();
  StaticSegment seg("", 0);
  seg.LoadFromFile("polarsite.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites;
  for (const StaticSite& site : seg) {
    externalsites.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
  }
  AOMultipole dip;
  dip.FillPotential(aobasis, externalsites);

  Eigen::MatrixXd dip_ref = Eigen::MatrixXd::Zero(17, 17);
  dip_ref << 0.31114997753, 0.059568868026, 0.0090978711864, 0, 0,
      0.056104697636, 0.0013498178976, 0, 0, 0.0061933281198, 0.025459181656,
      0.0061933281198, 0.025459181656, 0.0056130806569, 0.024860733171,
      0.0056130806569, 0.024860733171, 0.059568868026, 0.31114997753,
      0.066842196408, 0, 0, 0.2368963398, 0.042609848798, 0, 0, 0.073695437658,
      0.13588175121, 0.073695437658, 0.13588175121, 0.047134059517,
      0.11392427425, 0.047134059517, 0.11392427425, 0.0090978711864,
      0.066842196408, 0.32666220712, 0, 0, 0.065599720802, 0.17980265473, 0, 0,
      0.075224551189, 0.083547839473, 0.075224551189, 0.083547839473,
      -0.035337083996, -0.0087046753266, -0.035337083996, -0.0087046753266, 0,
      0, 0, 0.30339386273, 0, 0, 0, 0.15697695635, 0, 0.061346561258,
      0.043126415371, -0.061346561258, -0.043126415371, -0.040820728861,
      -0.037257483146, 0.040820728861, 0.037257483146, 0, 0, 0, 0,
      0.30339386273, 0, 0, 0, 0.15697695635, 0.061346561258, 0.043126415371,
      -0.061346561258, -0.043126415371, 0.040820728861, 0.037257483146,
      -0.040820728861, -0.037257483146, 0.056104697636, 0.2368963398,
      0.065599720802, 0, 0, 0.31114557005, 0.12399363244, 0, 0, 0.13566608278,
      0.24813270723, 0.13566608278, 0.24813270723, 0.072222745084,
      0.16730389611, 0.072222745084, 0.16730389611, 0.0013498178976,
      0.042609848798, 0.17980265473, 0, 0, 0.12399363244, 0.38517412728, 0, 0,
      0.13766652824, 0.2355220863, 0.13766652824, 0.2355220863, -0.05335528609,
      -0.024078433641, -0.05335528609, -0.024078433641, 0, 0, 0, 0.15697695635,
      0, 0, 0, 0.27407784996, 0, 0.10959962132, 0.10744652669, -0.10959962132,
      -0.10744652669, -0.060016250362, -0.076591108649, 0.060016250362,
      0.076591108649, 0, 0, 0, 0, 0.15697695635, 0, 0, 0, 0.27407784996,
      0.10959962132, 0.10744652669, -0.10959962132, -0.10744652669,
      0.060016250362, 0.076591108649, -0.060016250362, -0.076591108649,
      0.0061933281198, 0.073695437658, 0.075224551189, 0.061346561258,
      0.061346561258, 0.13566608278, 0.13766652824, 0.10959962132,
      0.10959962132, 0.40885081105, 0.26407630555, 0.0038749978449,
      0.053453882722, 0.002270661607, 0.043208755867, 0.002270661607,
      0.043208755867, 0.025459181656, 0.13588175121, 0.083547839473,
      0.043126415371, 0.043126415371, 0.24813270723, 0.2355220863,
      0.10744652669, 0.10744652669, 0.26407630555, 0.40853020343,
      0.053453882722, 0.17647932404, 0.026292260419, 0.10354578683,
      0.026292260419, 0.10354578683, 0.0061933281198, 0.073695437658,
      0.075224551189, -0.061346561258, -0.061346561258, 0.13566608278,
      0.13766652824, -0.10959962132, -0.10959962132, 0.0038749978449,
      0.053453882722, 0.40885081105, 0.26407630555, 0.002270661607,
      0.043208755867, 0.002270661607, 0.043208755867, 0.025459181656,
      0.13588175121, 0.083547839473, -0.043126415371, -0.043126415371,
      0.24813270723, 0.2355220863, -0.10744652669, -0.10744652669,
      0.053453882722, 0.17647932404, 0.26407630555, 0.40853020343,
      0.026292260419, 0.10354578683, 0.026292260419, 0.10354578683,
      0.0056130806569, 0.047134059517, -0.035337083996, -0.040820728861,
      0.040820728861, 0.072222745084, -0.05335528609, -0.060016250362,
      0.060016250362, 0.002270661607, 0.026292260419, 0.002270661607,
      0.026292260419, 0.19479973371, 0.12582094161, 0.001654409346,
      0.023957628139, 0.024860733171, 0.11392427425, -0.0087046753266,
      -0.037257483146, 0.037257483146, 0.16730389611, -0.024078433641,
      -0.076591108649, 0.076591108649, 0.043208755867, 0.10354578683,
      0.043208755867, 0.10354578683, 0.12582094161, 0.19479972247,
      0.023957628139, 0.075477416664, 0.0056130806569, 0.047134059517,
      -0.035337083996, 0.040820728861, -0.040820728861, 0.072222745084,
      -0.05335528609, 0.060016250362, -0.060016250362, 0.002270661607,
      0.026292260419, 0.002270661607, 0.026292260419, 0.001654409346,
      0.023957628139, 0.19479973371, 0.12582094161, 0.024860733171,
      0.11392427425, -0.0087046753266, 0.037257483146, -0.037257483146,
      0.16730389611, -0.024078433641, 0.076591108649, -0.076591108649,
      0.043208755867, 0.10354578683, 0.043208755867, 0.10354578683,
      0.023957628139, 0.075477416664, 0.12582094161, 0.19479972247;
  bool dip_check = dip_ref.isApprox(dip.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(dip_check, 1);
  if (!dip_check) {
    std::cout << "dip Ref" << endl;
    std::cout << dip_ref << endl;
    std::cout << "Dip" << endl;
    std::cout << dip.Matrix() << endl;
  }

  ofstream mpsfile2("polarsite2.mps");
  mpsfile2 << "! One Site" << endl;
  mpsfile2 << "! N=1 " << endl;
  mpsfile2 << "Units angstrom" << endl;
  mpsfile2 << "  C +0 0 3 Rank 2" << endl;
  mpsfile2 << "+0" << endl;
  mpsfile2 << "0 0 0" << endl;
  mpsfile2 << "100 0 0 0 0" << endl;
  mpsfile2
      << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 +1.9445387 "
      << endl;
  mpsfile2.close();
  StaticSegment seg2("", 0);
  seg2.LoadFromFile("polarsite2.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites2;
  for (const StaticSite& site : seg2) {
    externalsites2.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
  }

  AOMultipole quad;
  quad.FillPotential(aobasis, externalsites2);

  Eigen::MatrixXd quad_ref = Eigen::MatrixXd::Zero(17, 17);

  quad_ref << -0.411624, -0.0788044, -0.0180535, 0, 0, -0.0742216, -0.00267854,
      0, 0, -0.00838389, -0.0338781, -0.00838389, -0.0338781, -0.00723469,
      -0.0326907, -0.00723469, -0.0326907, -0.0788044, -0.411624, -0.13264, 0,
      0, -0.313393, -0.0845537, 0, 0, -0.104821, -0.186907, -0.104821,
      -0.186907, -0.0540895, -0.143531, -0.0540895, -0.143531, -0.0180535,
      -0.13264, -0.452667, 0, 0, -0.130174, -0.257994, 0, 0, -0.114539,
      -0.140379, -0.114539, -0.140379, 0.037363, -0.00805777, 0.037363,
      -0.00805777, 0, 0, 0, -0.391103, 0, 0, 0, -0.197601, 0, -0.0838419,
      -0.0558751, 0.0838419, 0.0558751, 0.0463037, 0.0451165, -0.0463037,
      -0.0451165, 0, 0, 0, 0, -0.391103, 0, 0, 0, -0.197601, -0.0838419,
      -0.0558751, 0.0838419, 0.0558751, -0.0463037, -0.0451165, 0.0463037,
      0.0451165, -0.0742216, -0.313393, -0.130174, 0, 0, -0.411571, -0.245827,
      0, 0, -0.190962, -0.351334, -0.190962, -0.351334, -0.0766038, -0.195938,
      -0.0766038, -0.195938, -0.00267854, -0.0845537, -0.257994, 0, 0,
      -0.245827, -0.60644, 0, 0, -0.206616, -0.406573, -0.206616, -0.406573,
      0.0523077, -0.00658091, 0.0523077, -0.00658091, 0, 0, 0, -0.197601, 0, 0,
      0, -0.313544, 0, -0.144221, -0.125379, 0.144221, 0.125379, 0.0615623,
      0.0796695, -0.0615623, -0.0796695, 0, 0, 0, 0, -0.197601, 0, 0, 0,
      -0.313544, -0.144221, -0.125379, 0.144221, 0.125379, -0.0615623,
      -0.0796695, 0.0615623, 0.0796695, -0.00838389, -0.104821, -0.114539,
      -0.0838419, -0.0838419, -0.190962, -0.206616, -0.144221, -0.144221,
      -0.557673, -0.3602, -0.00648275, -0.0819761, -0.00281427, -0.0584669,
      -0.00281427, -0.0584669, -0.0338781, -0.186907, -0.140379, -0.0558751,
      -0.0558751, -0.351334, -0.406573, -0.125379, -0.125379, -0.3602,
      -0.554921, -0.0819761, -0.292576, -0.0285113, -0.128308, -0.0285113,
      -0.128308, -0.00838389, -0.104821, -0.114539, 0.0838419, 0.0838419,
      -0.190962, -0.206616, 0.144221, 0.144221, -0.00648275, -0.0819761,
      -0.557673, -0.3602, -0.00281427, -0.0584669, -0.00281427, -0.0584669,
      -0.0338781, -0.186907, -0.140379, 0.0558751, 0.0558751, -0.351334,
      -0.406573, 0.125379, 0.125379, -0.0819761, -0.292576, -0.3602, -0.554921,
      -0.0285113, -0.128308, -0.0285113, -0.128308, -0.00723469, -0.0540895,
      0.037363, 0.0463037, -0.0463037, -0.0766038, 0.0523077, 0.0615623,
      -0.0615623, -0.00281427, -0.0285113, -0.00281427, -0.0285113, -0.194913,
      -0.125894, -0.00180873, -0.0252281, -0.0326907, -0.143531, -0.00805777,
      0.0451165, -0.0451165, -0.195938, -0.00658091, 0.0796695, -0.0796695,
      -0.0584669, -0.128308, -0.0584669, -0.128308, -0.125894, -0.194913,
      -0.0252281, -0.0825406, -0.00723469, -0.0540895, 0.037363, -0.0463037,
      0.0463037, -0.0766038, 0.0523077, -0.0615623, 0.0615623, -0.00281427,
      -0.0285113, -0.00281427, -0.0285113, -0.00180873, -0.0252281, -0.194913,
      -0.125894, -0.0326907, -0.143531, -0.00805777, -0.0451165, 0.0451165,
      -0.195938, -0.00658091, -0.0796695, 0.0796695, -0.0584669, -0.128308,
      -0.0584669, -0.128308, -0.0252281, -0.0825406, -0.125894, -0.194913;

  bool quad_check = quad_ref.isApprox(quad.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(quad_check, 1);
  if (!quad_check) {
    std::cout << "Quad Ref" << endl;
    std::cout << quad_ref << endl;
    std::cout << "Quad" << endl;
    std::cout << quad.Matrix() << endl;
  }
}

BOOST_AUTO_TEST_CASE(aomultipole_comparison) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  {
    ofstream mpsfile("polarsite.mps");
    mpsfile << "! One Site" << endl;
    mpsfile << "! N=1 " << endl;
    mpsfile << "Units angstrom" << endl;
    mpsfile << "  C +0 0 15 Rank 1" << endl;
    mpsfile << "+0" << endl;
    mpsfile << "10 10 10" << endl;
    mpsfile << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 "
               "+1.9445387 "
            << endl;
    mpsfile.close();
    StaticSegment seg("", 0);
    seg.LoadFromFile("polarsite.mps");

    std::vector<std::unique_ptr<StaticSite> > externalsites;
    for (const StaticSite& site : seg) {
      externalsites.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site)));
    }
    AOMultipole dip;
    dip.FillPotential(aobasis, externalsites);

    Eigen::MatrixXd dip_ref = Eigen::MatrixXd::Zero(17, 17);

    double a = 0.1;                            // this is in a0
    double mag_d = seg[0].getDipole().norm();  // this is in e * a0
    const Eigen::Vector3d dir_d = seg[0].getDipole().normalized();
    const Eigen::Vector3d A = seg[0].getPos() + 0.5 * a * dir_d;
    const Eigen::Vector3d B = seg[0].getPos() - 0.5 * a * dir_d;
    double qA = mag_d / a;
    double qB = -qA;
    StaticSite site1 = StaticSite(0, "", A);
    site1.setCharge(qA);
    StaticSite site2 = StaticSite(1, "", B);
    site2.setCharge(qB);
    std::vector<std::unique_ptr<StaticSite> > externalsites_mono;
    externalsites_mono.push_back(
        std::unique_ptr<StaticSite>(new StaticSite(site1)));
    externalsites_mono.push_back(
        std::unique_ptr<StaticSite>(new StaticSite(site2)));
    AOMultipole mono2;
    mono2.FillPotential(aobasis, externalsites_mono);

    bool dip_check = mono2.Matrix().isApprox(dip.Matrix(), 1e-4);
    BOOST_CHECK_EQUAL(dip_check, 1);
    if (!dip_check) {
      std::cout << "mono2 Ref" << endl;
      std::cout << mono2.Matrix() << endl;
      std::cout << "Dip" << endl;
      std::cout << dip.Matrix() << endl;
    }
  }

  {
    ofstream mpsfile2("polarsite2.mps");
    mpsfile2 << "! One Site" << endl;
    mpsfile2 << "! N=1 " << endl;
    mpsfile2 << "Units angstrom" << endl;
    mpsfile2 << "  C +0 0 15 Rank 2" << endl;
    mpsfile2 << "+0" << endl;
    mpsfile2 << "0 0 0" << endl;
    mpsfile2 << "100 100 100 100 100" << endl;
    mpsfile2 << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 "
                "+1.9445387 "
             << endl;
    mpsfile2.close();
    StaticSegment seg2("", 0);
    seg2.LoadFromFile("polarsite2.mps");

    std::vector<std::unique_ptr<StaticSite> > externalsites2;
    for (const StaticSite& site : seg2) {
      externalsites2.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site)));
    }
    AOMultipole quad;
    quad.FillPotential(aobasis, externalsites2);

    std::vector<std::unique_ptr<StaticSite> > externalsites_mono6;
    const Eigen::Matrix3d components = seg2[0].CalculateCartesianMultipole();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(components);
    double a = 2 * 0.01;
    for (int i = 0; i < 3; i++) {
      double q = es.eigenvalues()[i] / (a * a);
      const Eigen::Vector3d vec1 =
          seg2[0].getPos() + 0.5 * a * es.eigenvectors().col(i);
      const Eigen::Vector3d vec2 =
          seg2[0].getPos() - 0.5 * a * es.eigenvectors().col(i);
      StaticSite site1 = StaticSite(0, "", vec1);
      site1.setCharge(q);
      StaticSite site2 = StaticSite(1, "", vec2);
      site2.setCharge(q);
      externalsites_mono6.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site1)));
      externalsites_mono6.push_back(
          std::unique_ptr<StaticSite>(new StaticSite(site2)));
    }

    AOMultipole mono6;
    mono6.FillPotential(aobasis, externalsites_mono6);

    bool quad_check = mono6.Matrix().isApprox(quad.Matrix(), 1e-4);
    BOOST_CHECK_EQUAL(quad_check, 1);
    if (!quad_check) {
      std::cout << "mono6 Ref" << endl;
      std::cout << mono6.Matrix() << endl;
      std::cout << "Quad" << endl;
      std::cout << quad.Matrix() << endl;
      std::cout << "diff" << endl;
      std::cout << mono6.Matrix().cwiseQuotient(quad.Matrix()) << endl;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
