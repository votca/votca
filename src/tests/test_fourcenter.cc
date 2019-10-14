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

#define BOOST_TEST_MODULE fourcenter_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/fourcenter.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(fourcenter_test)

QMMolecule Methane() {
  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile("molecule.xyz");
  return mol;
}

BOOST_AUTO_TEST_CASE(small_l_test) {

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

  QMMolecule mol = Methane();
  BasisSet basis;
  basis.Load("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  FCMatrix fcenter;
  Eigen::Tensor<double, 4> block(
      aobasis.getShell(0).getNumFunc(), aobasis.getShell(4).getNumFunc(),
      aobasis.getShell(2).getNumFunc(), aobasis.getShell(1).getNumFunc());
  block.setZero();
  // S, S,SP,SP
  fcenter.FillFourCenterRepBlock(block, aobasis.getShell(0),
                                 aobasis.getShell(4), aobasis.getShell(2),
                                 aobasis.getShell(1));

  Eigen::Map<Eigen::VectorXd> mapped_result(block.data(), block.size());
  Eigen::VectorXd ref = Eigen::VectorXd::Zero(block.size());
  ref << 0.0569081, 0.000525033, 0.000525033, 0.000525033, 0.00112254,
      0.0329851, 1.75149e-05, 1.75149e-05, 0.00112254, 1.75149e-05, 0.0329851,
      1.75149e-05, 0.00112254, 1.75149e-05, 1.75149e-05, 0.0329851;

  Eigen::TensorMap<Eigen::Tensor<double, 4> > ref_block(
      ref.data(), aobasis.getShell(0).getNumFunc(),
      aobasis.getShell(4).getNumFunc(), aobasis.getShell(2).getNumFunc(),
      aobasis.getShell(1).getNumFunc());

  bool check = mapped_result.isApprox(ref, 0.0001);
  BOOST_CHECK_EQUAL(check, 1);
  if (!check) {
    cout << "ref" << endl;
    cout << ref_block << endl;
    cout << "result" << endl;
    cout << block << endl;
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

  FCMatrix fcenter;
  Eigen::Tensor<double, 4> block(
      dftbasis.getShell(0).getNumFunc(), dftbasis.getShell(1).getNumFunc(),
      dftbasis.getShell(0).getNumFunc(), dftbasis.getShell(1).getNumFunc());
  block.setZero();
  fcenter.FillFourCenterRepBlock(block, dftbasis.getShell(0),
                                 dftbasis.getShell(1), dftbasis.getShell(0),
                                 dftbasis.getShell(1));
  // we only check the first and last 600 values because this gets silly quite
  // quickly
  Eigen::Map<Eigen::VectorXd> mapped_result(block.data(), block.size());
  Eigen::VectorXd ref_head = Eigen::VectorXd::Zero(600);
  ref_head << 0.000154879, 0, 0, 0, -0.000207289, 0, 0, 0, 0.00018609, 0,
      1.46214e-05, 0, 0, 0, -3.66305e-05, 0, 0, 0, 0, 0, -0.00024582, 0, 0, 0,
      0.000202586, 0, 0, 0, 0, 0, -7.2164e-05, 0, 0, 0, 0.000151013, 0,
      -0.000207289, 0, 0, 0, 0.000288865, 0, 0, 0, -0.000300126, 0,
      -3.66305e-05, 0, 0, 0, 9.68421e-05, 0, 0, 0, 0, 0, 0.000202586, 0, 0, 0,
      -0.000190589, 0, 0, 0, 0, 0, 0.000151013, 0, 0, 0, -0.00040428, 0,
      0.00018609, 0, 0, 0, -0.000300126, 0, 0, 0, 0.000486627, 0, 7.06567e-07,
      0, 0, 0, -1.63569e-06, 0, 0, 0, 4.97533e-07, 0, 0, 0, -6.84264e-07, 0, 0,
      0, 3.99569e-07, 0, 0, 0, 1.54412e-06, 0, 0, 0, -5.28594e-06, 0, 0, 0,
      1.81524e-06, 0, 0, 0, -8.90228e-07, 0, 0, 0, -9.65282e-07, 0, 0, 0,
      2.31786e-06, 0, 0, 0, -1.14551e-06, 0, 0, 0, 1.63645e-06, 0, 0, 0,
      -1.14625e-06, 0, 0, 0, -8.0714e-07, 0, 0, 0, 3.26105e-06, 0, 0, 0,
      -5.4175e-06, 0, 0, 0, 2.96428e-06, 0, 0, 0, 7.41209e-07, 0, 0, 0,
      -2.09356e-06, 0, 0, 0, 0, 0, 1.76361e-05, 0, 0, 0, -1.34772e-05, 0, 0, 0,
      0, 0, -1.99285e-06, 0, 0, 0, 6.69284e-06, 0, -1.26266e-05, 0, 0, 0,
      2.11377e-05, 0, 0, 0, -3.46685e-05, 0, 2.54723e-06, 0, 0, 0, -6.05119e-06,
      0, 0, 0, 0, 0, -2.76912e-05, 0, 0, 0, 2.20562e-05, 0, 0, 0, 0, 0,
      4.6536e-06, 0, 0, 0, -1.65484e-05, 0, 9.40844e-06, 0, 0, 0, -1.65093e-05,
      0, 0, 0, 3.04757e-05, 0, -7.59235e-06, 0, 0, 0, 1.90044e-05, 0, 0, 0, 0,
      0, 3.98164e-05, 0, 0, 0, -3.55039e-05, 0, 0, 0, 0, 0, 2.17269e-05, 0, 0,
      0, -4.60101e-05, 0, 0, 0, -8.72453e-06, 0, 0, 0, 4.51535e-06, 0, 0, 0,
      8.88534e-06, 0, 0, 0, -2.22662e-05, 0, 0, 0, -2.13718e-05, 0, 0, 0,
      2.48627e-05, 0, 0, 0, -9.34983e-06, 0, 0, 0, -2.5474e-05, 0, 0, 0,
      5.81377e-05, 0, 0, 0, 2.18666e-05, 0, 0, 0, -1.24967e-05, 0, 0, 0,
      -4.69753e-06, 0, 0, 0, 1.29828e-05, 0, 0, 0, 4.5322e-05, 0, 0, 0,
      -5.68227e-05, 0, 0, 0, 3.37806e-05, 0, 0, 0, 1.03893e-05, 0, 0, 0,
      -3.65769e-05, 0, -0.000207289, 0, 0, 0, 0.000278892, 0, 0, 0,
      -0.000255849, 0, -1.90244e-05, 0, 0, 0, 4.77843e-05, 0, 0, 0, 0, 0,
      0.00032722, 0, 0, 0, -0.000270709, 0, 0, 0, 0, 0, 9.38385e-05, 0, 0, 0,
      -0.000198009, 0, 0.000278696, 0, 0, 0, -0.000389881, 0, 0, 0, 0.000410996,
      0, 4.7755e-05, 0, 0, 0, -0.000126478, 0, 0, 0, 0, 0, -0.000270475, 0, 0,
      0, 0.000254326, 0, 0, 0, 0, 0, -0.000197671, 0, 0, 0, 0.000531983, 0,
      -0.000255013, 0, 0, 0, 0.000410047, 0, 0, 0, -0.000660739, 0,
      -1.63569e-06, 0, 0, 0, 3.93349e-06, 0, 0, 0, -1.14551e-06, 0, 0, 0,
      1.56139e-06, 0, 0, 0, -9.64519e-07, 0, 0, 0, -3.54162e-06, 0, 0, 0,
      1.15576e-05, 0, 0, 0, -4.20943e-06, 0, 0, 0, 2.24771e-06, 0, 0, 0,
      2.22163e-06, 0, 0, 0, -5.47348e-06, 0, 0, 0, 2.74742e-06, 0, 0, 0,
      -3.82355e-06, 0, 0, 0, 2.60145e-06, 0, 0, 0, 1.98797e-06, 0, 0, 0,
      -7.23079e-06, 0, 0, 0, 1.1799e-05, 0, 0, 0, -6.57267e-06, 0, 0, 0,
      -1.77545e-06, 0, 0, 0, 4.8898e-06, 0, 0, 0, 0, 0, -1.34772e-05, 0, 0, 0,
      1.15948e-05, 0, 0, 0, 0, 0, 1.02595e-06, 0, 0, 0, -3.0199e-06, 0,
      9.40844e-06, 0, 0, 0, -1.49465e-05, 0, 0, 0, 2.20935e-05, 0, -1.38396e-06,
      0, 0, 0, 3.59972e-06, 0, 0, 0, 0, 0, 2.02644e-05, 0, 0, 0, -1.75648e-05,
      0, 0, 0, 0, 0, -2.62081e-06, 0, 0, 0, 7.14787e-06, 0, -8.05938e-06, 0, 0,
      0, 1.27386e-05, 0, 0, 0, -1.80383e-05, 0, 3.51684e-06, 0, 0, 0,
      -8.74311e-06, 0, 0, 0, 0, 0, -2.5891e-05, 0, 0, 0, 2.216e-05, 0, 0, 0, 0,
      0, -4.60101e-05, 0, 0, 0, 0.000100836, 0, 0, 0, 1.81445e-05, 0, 0, 0,
      -9.6702e-06, 0, 0, 0, -1.8603e-05, 0, 0, 0, 4.67064e-05, 0, 0, 0,
      4.5322e-05, 0, 0, 0, -5.38802e-05, 0;

  Eigen::VectorXd ref_tail = Eigen::VectorXd::Zero(600);
  ref_tail << 0, 1.91895e-06, 0, 0, 0, -3.35984e-06, 0, 0, 0, -1.12417e-06, 0,
      0, 0, 6.45611e-06, 0, 0, 0, -8.40705e-06, 0, 0, 0, 5.98931e-06, 0, 0, 0,
      1.75656e-06, 0, 0, 0, -4.92594e-06, 0, 0, 0, 0, 0, 3.98164e-05, 0, 0, 0,
      -2.5891e-05, 0, 0, 0, 0, 0, -5.35609e-06, 0, 0, 0, 2.60192e-05, 0,
      -3.46685e-05, 0, 0, 0, 7.21411e-05, 0, 0, 0, -0.000172952, 0, 6.04413e-06,
      0, 0, 0, -1.36934e-05, 0, 0, 0, 0, 0, -7.8967e-05, 0, 0, 0, 5.98723e-05,
      0, 0, 0, 0, 0, 1.20808e-05, 0, 0, 0, -7.11702e-05, 0, 2.20935e-05, 0, 0,
      0, -5.47172e-05, 0, 0, 0, 0.000166337, 0, -2.73454e-05, 0, 0, 0,
      7.43369e-05, 0, 0, 0, 0, 0, 0.000179694, 0, 0, 0, -0.000171822, 0, 0, 0,
      0, 0, 1.03893e-05, 0, 0, 0, -2.58539e-05, 0, 0, 0, -3.34526e-06, 0, 0, 0,
      2.32969e-06, 0, 0, 0, 3.75158e-06, 0, 0, 0, -9.81597e-06, 0, 0, 0,
      -9.34983e-06, 0, 0, 0, 1.24723e-05, 0, 0, 0, -1.05839e-05, 0, 0, 0,
      -1.43844e-05, 0, 0, 0, 3.64727e-05, 0, 0, 0, 9.08306e-06, 0, 0, 0,
      -6.39519e-06, 0, 0, 0, -2.8213e-06, 0, 0, 0, 7.47498e-06, 0, 0, 0,
      2.43554e-05, 0, 0, 0, -3.30882e-05, 0, 0, 0, 3.00837e-05, 0, 0, 0,
      1.40516e-05, 0, 0, 0, -3.8288e-05, 0, -0.000300126, 0, 0, 0, 0.000410047,
      0, 0, 0, -0.000400827, 0, -2.50483e-05, 0, 0, 0, 6.334e-05, 0, 0, 0, 0, 0,
      0.000472519, 0, 0, 0, -0.000395247, 0, 0, 0, 0, 0, 0.000122209, 0, 0, 0,
      -0.000263552, 0, 0.000410996, 0, 0, 0, -0.0005815, 0, 0, 0, 0.000639999,
      0, 6.34912e-05, 0, 0, 0, -0.000168793, 0, 0, 0, 0, 0, -0.000396745, 0, 0,
      0, 0.000370954, 0, 0, 0, 0, 0, -0.000264962, 0, 0, 0, 0.000723024, 0,
      -0.000405204, 0, 0, 0, 0.000645062, 0, 0, 0, -0.0010162, 0, -2.09356e-06,
      0, 0, 0, 4.8898e-06, 0, 0, 0, -1.14625e-06, 0, 0, 0, 2.17032e-06, 0, 0, 0,
      -3.35984e-06, 0, 0, 0, -4.46628e-06, 0, 0, 0, 2.45479e-05, 0, 0, 0,
      -5.74841e-06, 0, 0, 0, 3.41706e-06, 0, 0, 0, 3.44954e-06, 0, 0, 0,
      -8.61514e-06, 0, 0, 0, 2.60145e-06, 0, 0, 0, -5.42124e-06, 0, 0, 0,
      9.74276e-06, 0, 0, 0, 2.91801e-06, 0, 0, 0, -1.8978e-05, 0, 0, 0,
      2.50108e-05, 0, 0, 0, -1.7715e-05, 0, 0, 0, -4.92594e-06, 0, 0, 0,
      1.41673e-05, 0, 0, 0, 0, 0, -3.55039e-05, 0, 0, 0, 2.216e-05, 0, 0, 0, 0,
      0, 5.09326e-06, 0, 0, 0, -2.52448e-05, 0, 3.04757e-05, 0, 0, 0,
      -6.59293e-05, 0, 0, 0, 0.000166337, 0, -5.60962e-06, 0, 0, 0, 1.24753e-05,
      0, 0, 0, 0, 0, 7.25763e-05, 0, 0, 0, -5.45388e-05, 0, 0, 0, 0, 0,
      -1.1181e-05, 0, 0, 0, 6.90232e-05, 0, -1.80383e-05, 0, 0, 0, 4.88879e-05,
      0, 0, 0, -0.00016137, 0, 2.61236e-05, 0, 0, 0, -7.13757e-05, 0, 0, 0, 0,
      0, -0.000171822, 0, 0, 0, 0.000166651, 0, 0, 0, 0, 0, -3.65769e-05, 0, 0,
      0, 8.65917e-05, 0, 0, 0, 1.60594e-05, 0, 0, 0, -1.03721e-05, 0, 0, 0,
      -1.75306e-05, 0, 0, 0, 4.68748e-05, 0, 0, 0, 3.37806e-05, 0, 0, 0,
      -4.31343e-05, 0, 0, 0, 3.00837e-05, 0, 0, 0, 4.79355e-05, 0, 0, 0,
      -0.000120333, 0, 0, 0, -4.32202e-05, 0, 0, 0, 2.94858e-05, 0, 0, 0,
      1.20523e-05, 0, 0, 0, -3.39703e-05, 0, 0, 0, -8.09325e-05, 0, 0, 0,
      0.000109618, 0, 0, 0, -9.58214e-05, 0, 0, 0, -3.8288e-05, 0, 0, 0,
      0.000118792, 0, 0.000486627, 0, 0, 0, -0.000660739, 0, 0, 0, 0.000627084,
      0, 4.48388e-05, 0, 0, 0, -0.000113808, 0, 0, 0, 0, 0, -0.00083548, 0, 0,
      0, 0.00070618, 0, 0, 0, 0, 0, -0.000205877, 0, 0, 0, 0.000437748, 0,
      -0.000660739, 0, 0, 0, 0.000933188, 0, 0, 0, -0.0010162, 0, -0.000113808,
      0, 0, 0, 0.000304671, 0, 0, 0, 0, 0, 0.00070618, 0, 0, 0, -0.000670916, 0,
      0, 0, 0, 0, 0.000437748, 0, 0, 0, -0.00121331, 0, 0.000627084, 0, 0, 0,
      -0.0010162, 0, 0, 0, 0.00167288;
  bool check_head = mapped_result.head<600>().isApprox(ref_head, 0.0001);
  BOOST_CHECK_EQUAL(check_head, 1);
  if (!check_head) {
    cout << "ref" << endl;
    cout << ref_head.transpose() << endl;
    cout << "result" << endl;
    cout << mapped_result.head<600>().transpose() << endl;
  }

  bool check_tail = mapped_result.tail<600>().isApprox(ref_tail, 0.0001);
  BOOST_CHECK_EQUAL(check_tail, 1);
  if (!check_tail) {
    cout << "ref" << endl;
    cout << ref_tail.transpose() << endl;
    cout << "result" << endl;
    cout << mapped_result.tail<600>().transpose() << endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
