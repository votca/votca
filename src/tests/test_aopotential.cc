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
  basis.Load("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

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
  ECPBasisSet ecps;
  ecps.Load("ecp.xml");
  ECPAOBasis ecpbasis;
  ecpbasis.Fill(ecps, orbitals.QMAtoms());
  AOECP ecp;
  OPENMP::setMaxThreads(1);
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

  AOPlanewave planewave;
  std::vector<Eigen::Vector3d> kpoints = {
      {1, 1, 1}, {2, 1, 1}, {-1, -1, -1}, {-2, -1, -1}};
  Eigen::MatrixXcd planewave_ref = Eigen::MatrixXcd::Zero(17, 17);
  planewave_ref.real() << 3.75219, 0.593498, 0, 0, 0, 0.608505, 0, 0, 0,
      0.0567792, 0.27062, 0.0632057, 0.272261, 0.064459, 0.272572, 0.0632057,
      0.272261, 0.593498, 1.49476, 0, 0, 0, 0.832298, 0, 0, 0, -0.238135,
      0.280867, 0.219652, 0.41493, 0.336769, 0.433781, 0.219652, 0.41493, 0, 0,
      1.73441, -0.561535, -0.782038, 0, 0.40241, -0.358572, -0.473183, -0.29074,
      -0.133159, 0.34088, 0.301536, -0.26042, -0.125687, -0.100203, -0.0636456,
      0, 0, -0.561535, 1.73441, -0.782038, 0, -0.358572, 0.40241, -0.473183,
      -0.29074, -0.133159, -0.100203, -0.0636456, -0.26042, -0.125687, 0.34088,
      0.301536, 0, 0, -0.782038, -0.782038, 1.0729, 0, -0.473183, -0.473183,
      0.0585781, -0.332886, -0.235043, -0.0390516, -0.00714159, 0.378793,
      0.250569, -0.0390516, -0.00714159, 0.608505, 0.832298, 0, 0, 0, 0.338232,
      0, 0, 0, -0.442088, -0.0451742, 0.121971, 0.165715, 0.361259, 0.18082,
      0.121971, 0.165715, 0, 0, 0.40241, -0.358572, -0.473183, 0, -0.0935013,
      -0.431733, -0.487193, -0.371122, -0.258685, 0.277042, 0.22887, -0.200825,
      0.0253968, 0.0711417, 0.0603183, 0, 0, -0.358572, 0.40241, -0.473183, 0,
      -0.431733, -0.0935013, -0.487193, -0.371122, -0.258685, 0.0711417,
      0.0603183, -0.200825, 0.0253968, 0.277042, 0.22887, 0, 0, -0.473183,
      -0.473183, 0.0585781, 0, -0.487193, -0.487193, -0.25988, -0.327681,
      -0.280982, 0.133092, 0.0875579, 0.425108, 0.209312, 0.133092, 0.0875579,
      0.0567792, -0.238135, -0.29074, -0.29074, -0.332886, -0.442088, -0.371122,
      -0.371122, -0.327681, -1.20105, -0.569133, 0.00595096, -0.153839,
      -0.000835787, -0.14224, 0.00595096, -0.153839, 0.27062, 0.280867,
      -0.133159, -0.133159, -0.235043, -0.0451742, -0.258685, -0.258685,
      -0.280982, -0.569133, -0.233828, 0.125534, 0.0385781, 0.145059, 0.025619,
      0.125534, 0.0385781, 0.0632057, 0.219652, 0.34088, -0.100203, -0.0390516,
      0.121971, 0.277042, 0.0711417, 0.133092, 0.00595096, 0.125534, -0.177345,
      0.00224974, 0.0060307, 0.0196205, -0.000862638, 0.0011456, 0.272261,
      0.41493, 0.301536, -0.0636456, -0.00714159, 0.165715, 0.22887, 0.0603183,
      0.0875579, -0.153839, 0.0385781, 0.00224974, 0.0721698, 0.0956966,
      0.0385781, 0.0011456, 0.025619, 0.064459, 0.336769, -0.26042, -0.26042,
      0.378793, 0.361259, -0.200825, -0.200825, 0.425108, -0.000835787,
      0.145059, 0.0060307, 0.0956966, 1.45719, 0.565884, 0.0060307, 0.0956966,
      0.272572, 0.433781, -0.125687, -0.125687, 0.250569, 0.18082, 0.0253968,
      0.0253968, 0.209312, -0.14224, 0.025619, 0.0196205, 0.0385781, 0.565884,
      0.129591, 0.0196205, 0.0385781, 0.0632057, 0.219652, -0.100203, 0.34088,
      -0.0390516, 0.121971, 0.0711417, 0.277042, 0.133092, 0.00595096, 0.125534,
      -0.000862638, 0.0011456, 0.0060307, 0.0196205, -0.177345, 0.00224974,
      0.272261, 0.41493, -0.0636456, 0.301536, -0.00714159, 0.165715, 0.0603183,
      0.22887, 0.0875579, -0.153839, 0.0385781, 0.0011456, 0.025619, 0.0956966,
      0.0385781, 0.00224974, 0.0721698;

  planewave.FillPotential(aobasis, kpoints);
  bool planewave_check = planewave_ref.isApprox(planewave.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(planewave_check, 1);
  if (!planewave_check) {
    std::cout << "planewave Ref real" << endl;
    std::cout << planewave_ref.real() << endl;
    std::cout << "planewave real" << endl;
    std::cout << planewave.Matrix().real() << endl;
    std::cout << "planewave Ref imag" << endl;
    std::cout << planewave_ref.imag() << endl;
    std::cout << "planewave imag" << endl;
    std::cout << planewave.Matrix().imag() << endl;
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
  basis.Load("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

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

BOOST_AUTO_TEST_CASE(large_l_test) {
  std::ofstream xyzfile("C2.xyz");
  xyzfile << " 2" << std::endl;
  xyzfile << " C2" << std::endl;
  xyzfile << " C            .000000     .000000     .000000" << std::endl;
  xyzfile << " C            1.000000     .000000     .000000" << std::endl;
  xyzfile.close();

  QMMolecule mol("C", 0);
  mol.LoadFromFile("C2.xyz");

  ofstream basisfile("G.xml");
  basisfile << "<basis name=\"G\">" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"G\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"G\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();

  BasisSet basisset;
  basisset.Load("G.xml");
  AOBasis dftbasis;
  dftbasis.Fill(basisset, mol);

  int dftbasissize = 18;

  ofstream mpsfile2("polarsite_all.mps");
  mpsfile2 << "! One Site" << endl;
  mpsfile2 << "! N=1 " << endl;
  mpsfile2 << "Units angstrom" << endl;
  mpsfile2 << "  C +0 0 15 Rank 2" << endl;
  mpsfile2 << "+1" << endl;
  mpsfile2 << "1 2 3" << endl;
  mpsfile2 << "100 100 100 100 100" << endl;
  mpsfile2 << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 "
              "+1.9445387 "
           << endl;
  mpsfile2.close();
  StaticSegment seg2("", 0);
  seg2.LoadFromFile("polarsite2.mps");

  std::vector<std::unique_ptr<StaticSite> > externalsites2;
  for (const StaticSite& site : seg2) {
    externalsites2.push_back(std::unique_ptr<StaticSite>(new StaticSite(site)));
  }

  AOMultipole esp;
  esp.FillPotential(dftbasis, externalsites2);
  Eigen::MatrixXd esp_ref = Eigen::MatrixXd::Zero(dftbasissize, dftbasissize);
  esp_ref << -0.00329622, 5.90075e-07, 5.90075e-07, 6.2501e-07, 6.2501e-07,
      -4.25931e-22, 3.08906e-23, 2.49204e-23, -1.6263e-19, -2.28713e-05,
      -1.0835e-09, -3.65012e-08, 2.90862e-07, 2.9385e-05, 3.78966e-10,
      -1.26539e-07, -5.00105e-07, -2.12666e-05, 5.90075e-07, -0.00329527,
      -4.66135e-07, 1.18679e-06, -1.18679e-06, 3.69095e-07, -3.69095e-07,
      4.07203e-23, 4.23516e-22, -1.1282e-09, -3.16016e-06, -1.87545e-07,
      5.22697e-08, 2.1684e-09, 7.87395e-06, 7.53975e-08, -3.17261e-07,
      -3.50329e-10, 5.90075e-07, -4.66135e-07, -0.0032962, 1.18679e-06,
      1.18679e-06, 3.69095e-07, 3.69095e-07, 6.47532e-24, 4.23516e-22,
      4.95528e-08, 1.86943e-07, 4.48941e-05, -3.1373e-09, -2.69242e-07,
      -4.65838e-07, -3.65944e-05, 2.39705e-08, 9.97351e-07, 6.2501e-07,
      1.18679e-06, 1.18679e-06, -0.00329428, 3.22744e-23, 1.74232e-06,
      -1.74232e-06, 2.45327e-07, -2.45327e-07, -2.91598e-07, -4.99958e-08,
      -3.04765e-09, 1.26168e-05, 2.9921e-07, 6.50238e-08, -7.883e-10,
      -2.43652e-05, 2.64349e-08, 6.2501e-07, -1.18679e-06, 1.18679e-06,
      3.22744e-23, -0.00329428, 1.74232e-06, 1.74232e-06, 2.45327e-07,
      2.45327e-07, 2.93851e-05, 2.1181e-09, 2.53129e-07, -2.97891e-07,
      -4.00233e-05, -3.37862e-09, -3.59518e-08, 6.05178e-07, 3.69522e-05,
      -4.25931e-22, 3.69095e-07, 3.69095e-07, 1.74232e-06, 1.74232e-06,
      -0.00329186, -7.86115e-23, 1.84059e-06, -1.84059e-06, 4.69653e-10,
      7.87399e-06, 4.67315e-07, -7.04073e-08, -3.48705e-09, -2.10232e-05,
      -2.325e-07, 8.65301e-07, 3.90158e-09, 3.08906e-23, -3.69095e-07,
      3.69095e-07, -1.74232e-06, 1.74232e-06, -7.86115e-23, -0.00329186,
      1.84059e-06, 1.84059e-06, 1.17066e-07, -7.46429e-08, -3.65947e-05,
      -7.54333e-10, 4.87653e-08, 2.30366e-07, 3.56938e-05, -1.18238e-08,
      -1.10373e-06, 2.49204e-23, 4.07203e-23, 6.47532e-24, 2.45327e-07,
      2.45327e-07, 1.84059e-06, 1.84059e-06, -0.00328848, -1.05879e-22,
      5.01174e-07, 3.143e-07, 2.39367e-08, -2.43654e-05, -6.07642e-07,
      -8.5722e-07, -1.19134e-08, 6.79029e-05, 2.06212e-07, -1.6263e-19,
      4.23516e-22, 4.23516e-22, -2.45327e-07, 2.45327e-07, -1.84059e-06,
      1.84059e-06, -1.05879e-22, -0.00328848, -2.1267e-05, -2.99221e-10,
      -9.88047e-07, -2.85798e-08, 3.69526e-05, 3.73364e-09, 1.09343e-06,
      -1.99306e-07, -6.99141e-05, -2.28713e-05, -1.1282e-09, 4.95528e-08,
      -2.91598e-07, 2.93851e-05, 4.69653e-10, 1.17066e-07, 5.01174e-07,
      -2.1267e-05, -0.00288848, 5.31864e-07, 6.79875e-07, 8.0963e-07,
      8.41929e-07, -1.71405e-10, -1.77116e-10, -5.07417e-12, -5.18613e-12,
      -1.0835e-09, -3.16016e-06, 1.86943e-07, -4.99958e-08, 2.1181e-09,
      7.87399e-06, -7.46429e-08, 3.143e-07, -2.99221e-10, 5.31864e-07,
      -0.00288747, -6.03853e-07, 1.36732e-06, -1.0697e-06, 4.97096e-07,
      -4.7805e-07, -1.31435e-10, 1.27205e-10, -3.65012e-08, -1.87545e-07,
      4.48941e-05, -3.04765e-09, 2.53129e-07, 4.67315e-07, -3.65947e-05,
      2.39367e-08, -9.88047e-07, 6.79875e-07, -6.03853e-07, -0.00288873,
      1.06979e-06, 1.36742e-06, 4.78061e-07, 4.97108e-07, -1.27208e-10,
      -1.31439e-10, 2.90862e-07, 5.22697e-08, -3.1373e-09, 1.26168e-05,
      -2.97891e-07, -7.04073e-08, -7.54333e-10, -2.43654e-05, -2.85798e-08,
      8.0963e-07, 1.36732e-06, 1.06979e-06, -0.00288695, -6.51611e-12,
      2.00737e-06, -1.57058e-06, 3.30307e-07, -3.17677e-07, 2.9385e-05,
      2.1684e-09, -2.69242e-07, 2.9921e-07, -4.00233e-05, -3.48705e-09,
      4.87653e-08, -6.07642e-07, 3.69526e-05, 8.41929e-07, -1.0697e-06,
      1.36742e-06, -6.51611e-12, -0.00288695, 1.57058e-06, 2.00737e-06,
      3.17677e-07, 3.30307e-07, 3.78966e-10, 7.87395e-06, -4.65838e-07,
      6.50238e-08, -3.37862e-09, -2.10232e-05, 2.30366e-07, -8.5722e-07,
      3.73364e-09, -1.71405e-10, 4.97096e-07, 4.78061e-07, 2.00737e-06,
      1.57058e-06, -0.00288504, -5.18336e-17, 2.12047e-06, -1.65932e-06,
      -1.26539e-07, 7.53975e-08, -3.65944e-05, -7.883e-10, -3.59518e-08,
      -2.325e-07, 3.56938e-05, -1.19134e-08, 1.09343e-06, -1.77116e-10,
      -4.7805e-07, 4.97108e-07, -1.57058e-06, 2.00737e-06, -5.18336e-17,
      -0.00288504, 1.65932e-06, 2.12047e-06, -5.00105e-07, -3.17261e-07,
      2.39705e-08, -2.43652e-05, 6.05178e-07, 8.65301e-07, -1.18238e-08,
      6.79029e-05, -1.99306e-07, -5.07417e-12, -1.31435e-10, -1.27208e-10,
      3.30307e-07, 3.17677e-07, 2.12047e-06, 1.65932e-06, -0.00288238,
      -8.47033e-22, -2.12666e-05, -3.50329e-10, 9.97351e-07, 2.64349e-08,
      3.69522e-05, 3.90158e-09, -1.10373e-06, 2.06212e-07, -6.99141e-05,
      -5.18613e-12, 1.27205e-10, -1.31439e-10, -3.17677e-07, 3.30307e-07,
      -1.65932e-06, 2.12047e-06, -8.47033e-22, -0.00288238;

  bool check_esp = esp.Matrix().isApprox(esp_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_esp, 1);
  if (!check_esp) {
    cout << "esp ref" << endl;
    cout << esp_ref << endl;
    cout << "esp result" << endl;
    cout << esp.Matrix() << endl;
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
  ecpfile << "    <shell type=\"D\"><constant power=\"1\" decay=\"5.96179600\" "
             "contraction=\"-9.45431800\"></constant></shell>"
          << endl;
  ecpfile << "  </element>" << endl;
  ecpfile << "</pseudopotential>" << endl;
  ecpfile.close();
  ECPBasisSet ecps;
  ecps.Load("ecp.xml");
  ECPAOBasis ecpbasis;
  ecpbasis.Fill(ecps, mol);

  AOECP ecp;
  ecp.FillPotential(dftbasis, ecpbasis);

  Eigen::MatrixXd ecp_ref = Eigen::MatrixXd::Zero(dftbasissize, dftbasissize);
  // apparently ecps are zero, sounds weird

  bool check_ecp = ecp.Matrix().isApprox(ecp_ref, 0.00001);

  BOOST_CHECK_EQUAL(check_ecp, 1);
  if (!check_ecp) {
    cout << "ref" << endl;
    cout << ecp_ref << endl;
    cout << "ecp result" << endl;
    cout << ecp.Matrix() << endl;
  }

  AOPlanewave planewave;
  std::vector<Eigen::Vector3d> kpoints = {
      {1, 1, 1}, {2, 1, 1}, {-1, -1, -1}, {-2, -1, -1}};
  Eigen::MatrixXcd planewave_ref =
      Eigen::MatrixXcd::Zero(dftbasissize, dftbasissize);
  planewave_ref.real() << 2.78686, -0.0725054, -0.10791, 0.422843, 0.202853,
      -0.0250416, 0.000303226, 0.010992, -0.00526033, 0.00392097, -4.96257e-05,
      -0.000536107, -0.00483455, -0.00509484, 0.000137166, -0.00105514,
      0.00856805, 0.00392341, -0.0725054, 2.92595, -0.313275, -0.213136,
      0.150911, 0.115488, -0.266102, 0.000805192, 0.0186204, 9.60884e-05,
      0.000503997, 0.00299424, 0.000493185, -0.00012728, -0.00126233,
      -0.00123335, -0.00304358, 9.45352e-05, -0.10791, -0.313275, 2.62469,
      -0.13664, -0.213264, 0.242433, 0.126715, -0.0178129, -0.000128433,
      0.000544835, -0.00299301, -0.00756059, 0.000205998, -0.00277786,
      0.00743251, 0.0062102, -0.000533323, 0.0103551, 0.422843, -0.213136,
      -0.13664, 2.74363, 0.0132207, -0.302104, 0.204886, 0.0817605, -0.173207,
      0.00483672, -0.000492351, -8.25166e-05, -0.00208217, -0.00501665,
      0.000700739, 7.01576e-05, 0.00413476, -0.00024633, 0.202853, 0.150911,
      -0.213264, 0.0132207, 2.73111, -0.205297, -0.301633, 0.173395, 0.0811361,
      -0.00509547, 3.85713e-05, 0.00276634, 0.00501296, 0.00698738,
      -0.000107299, -0.000656976, -0.0103686, -0.00668892, -0.0250416, 0.115488,
      0.242433, -0.302104, -0.205297, 2.66405, -0.000118649, -0.301713,
      0.207342, -0.000161553, -0.00126233, -0.00743581, -0.000702956,
      0.000253862, 0.00336668, 0.00372282, 0.00834485, -0.000359769,
      0.000303226, -0.266102, 0.126715, 0.204886, -0.301633, -0.000118649,
      2.66483, -0.207347, -0.301744, 0.00104784, 0.00123199, 0.00620541,
      -3.88906e-05, 0.000666773, -0.00371921, -0.0059922, 9.25436e-05,
      -0.011252, 0.010992, 0.000805192, -0.0178129, 0.0817605, 0.173395,
      -0.301713, -0.207347, 2.54589, -6.76574e-06, -0.00857285, 0.00304192,
      -0.000424275, 0.00413956, 0.0103775, -0.00834049, 0.000381058, -0.0114706,
      -0.00368562, -0.00526033, 0.0186204, -0.000128433, -0.173207, 0.0811361,
      0.207342, -0.301744, -6.76574e-06, 2.54587, 0.00392389, -7.63477e-05,
      -0.0103458, 0.000253077, -0.00668838, 0.000195053, 0.0112433, 0.00366679,
      0.0122079, 0.00392097, 9.60884e-05, 0.000544835, 0.00483672, -0.00509547,
      -0.000161553, 0.00104784, -0.00857285, 0.00392389, -1.48692, 0.0400761,
      0.0685193, -0.265065, -0.162967, 0.018105, 0.00176914, -0.00883063,
      0.00321974, -4.96257e-05, 0.000503997, -0.00299301, -0.000492351,
      3.85713e-05, -0.00126233, 0.00123199, 0.00304192, -7.63477e-05, 0.0400761,
      -1.59791, 0.196603, 0.135373, -0.084255, -0.0938385, 0.168324,
      0.000889671, -0.0134697, -0.000536107, 0.00299424, -0.00756059,
      -8.25166e-05, 0.00276634, -0.00743581, 0.00620541, -0.000424275,
      -0.0103458, 0.0685193, 0.196603, -1.35589, 0.0739181, 0.134347, -0.149309,
      -0.10074, 0.0127738, 0.00159265, -0.00483455, 0.000493185, 0.000205998,
      -0.00208217, 0.00501296, -0.000702956, -3.88906e-05, 0.00413956,
      0.000253077, -0.265065, 0.135373, 0.0739181, -1.44864, -0.0106212,
      0.189636, -0.111534, -0.0656841, 0.107445, -0.00509484, -0.00012728,
      -0.00277786, -0.00501665, 0.00698738, 0.000253862, 0.000666773, 0.0103775,
      -0.00668838, -0.162967, -0.084255, 0.134347, -0.0106212, -1.44094,
      0.111887, 0.189281, -0.107619, -0.0651825, 0.000137166, -0.00126233,
      0.00743251, 0.000700739, -0.000107299, 0.00336668, -0.00371921,
      -0.00834049, 0.000195053, 0.018105, -0.0938385, -0.149309, 0.189636,
      0.111887, -1.38511, 0.000109377, 0.187057, -0.111243, -0.00105514,
      -0.00123335, 0.0062102, 7.01576e-05, -0.000656976, 0.00372282, -0.0059922,
      0.000381058, 0.0112433, 0.00176914, 0.168324, -0.10074, -0.111534,
      0.189281, 0.000109377, -1.38574, 0.111246, 0.187082, 0.00856805,
      -0.00304358, -0.000533323, 0.00413476, -0.0103686, 0.00834485,
      9.25436e-05, -0.0114706, 0.00366679, -0.00883063, 0.000889671, 0.0127738,
      -0.0656841, -0.107619, 0.187057, 0.111246, -1.29299, 5.43541e-06,
      0.00392341, 9.45352e-05, 0.0103551, -0.00024633, -0.00668892,
      -0.000359769, -0.011252, -0.00368562, 0.0122079, 0.00321974, -0.0134697,
      0.00159265, 0.107445, -0.0651825, -0.111243, 0.187082, 5.43541e-06,
      -1.29297;

  planewave.FillPotential(dftbasis, kpoints);
  bool planewave_check = planewave_ref.isApprox(planewave.Matrix(), 1e-4);
  BOOST_CHECK_EQUAL(planewave_check, 1);
  if (!planewave_check) {
    std::cout << "planewave Ref real" << endl;
    std::cout << planewave_ref.real() << endl;
    std::cout << "planewave real" << endl;
    std::cout << planewave.Matrix().real() << endl;
    std::cout << "planewave Ref imag" << endl;
    std::cout << planewave_ref.imag() << endl;
    std::cout << "planewave imag" << endl;
    std::cout << planewave.Matrix().imag() << endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
