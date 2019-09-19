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

#define BOOST_TEST_MODULE ppm_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/ppm.h>
#include <votca/xtp/threecenter.h>

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(ppm_test)

BOOST_AUTO_TEST_CASE(ppm_full) {

  std::ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  std::ofstream basisfile("3-21G.xml");
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

  Orbitals orb;
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);

  AOKinetic kinetic;
  kinetic.Fill(aobasis);

  AOMultipole esp;
  esp.FillPotential(aobasis, orbitals.QMAtoms());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(kinetic.Matrix() +
                                                    esp.Matrix());

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, es.eigenvectors());

  RPA rpa = RPA(Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(es.eigenvalues());

  PPM ppm;
  ppm.PPM_construct_parameters(rpa);

  Eigen::VectorXd ppm_freq = Eigen::VectorXd::Zero(17);
  ppm_freq << 19.4503, 12.2429, 10.2167, 10.2167, 10.2167, 17.0832, 12.7117,
      12.7117, 12.7117, 10.9455, 10.9455, 10.9455, 11.1861, 9.60843, 9.60843,
      9.60843, 9.61518;
  Eigen::VectorXd ppm_w = Eigen::VectorXd::Zero(17);
  ppm_w << 3.59422e-05, 0.00121795, 0.00343632, 0.00343632, 0.00343632,
      0.00394659, 0.012066, 0.012066, 0.012066, 0.0241118, 0.0241118, 0.0241118,
      0.0286786, 0.11036, 0.11036, 0.11036, 0.191732;

  bool f_check = ppm_freq.isApprox(ppm.getPpm_freq(), 0.0001);

  bool w_check = ppm_w.isApprox(ppm.getPpm_weight(), 0.0001);

  if (!f_check) {
    cout << "ppm_freq" << endl;
    cout << ppm.getPpm_freq() << endl;
    cout << "ppm_freq_ref" << endl;
    cout << ppm_freq << endl;
  }

  if (!w_check) {
    cout << "ppm_w" << endl;
    cout << ppm.getPpm_weight() << endl;
    cout << "ppm_w_ref" << endl;
    cout << ppm_w << endl;
  }

  BOOST_CHECK_EQUAL(f_check, 1);
  BOOST_CHECK_EQUAL(w_check, 1);
}

BOOST_AUTO_TEST_SUITE_END()
