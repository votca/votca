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

#define BOOST_TEST_MODULE rpa_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>

using namespace std;
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_h2p) {

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

  Eigen::VectorXd eigenvals = Eigen::VectorXd::Zero(17);
  eigenvals << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
      0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418,
      3.4418, 4.81084, 17.1838;

  Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Zero(17, 17);
  eigenvectors << 0.0185815, 2.9133e-17, 8.49354e-17, -0.00312916, 0.0420075,
      1.11356e-16, 1.85886e-17, -0.0334732, 0.0485113, -8.71556e-18,
      -3.79994e-17, -0.0346485, -0.0248392, -3.32286e-22, -4.62643e-17,
      0.0144472, 0.996183, 0.166534, 2.80578e-16, 7.85515e-16, -0.0299557,
      0.409156, 1.10045e-15, 1.67666e-16, -0.336895, 0.608646, -1.00343e-16,
      -5.04622e-16, -0.437568, -0.299751, -2.77682e-17, -6.06812e-16, 0.176594,
      -0.0866677, 0.0010572, -0.0210402, -0.0345975, 0.035778, 0.0611836,
      0.0374747, -0.154443, 0.0892921, 0.0842611, -0.35309, -0.0759572,
      0.278374, -0.409082, -0.64367, 0.308248, -0.261525, -0.000315534,
      -0.0010572, -0.0404824, 0.000922593, -0.035778, -0.0611836, -0.115015,
      -0.109676, -0.0892921, -0.0842611, -0.242326, 0.267806, -0.278374,
      0.409082, -0.0548848, 0.711558, 0.261525, 0.000315534, 0.0010572,
      -0.0194422, 0.0355201, 0.035778, 0.0611836, -0.152489, 0.0447677,
      0.0892921, 0.0842611, 0.110764, 0.343764, 0.278374, -0.409082, 0.588785,
      0.403311, -0.261525, -0.000315534, -0.823783, -9.8891e-16, -3.34692e-15,
      0.103497, -0.277613, -5.51463e-16, 4.95594e-17, 0.163544, 0.121215,
      -7.23985e-17, -2.63149e-16, -0.259891, -0.284396, -1.74149e-16,
      -6.65818e-16, 0.208987, 0.00782842, -0.0333718, 0.22696, 0.373203,
      -0.337332, -0.251625, -0.144131, 0.594004, -0.329076, -0.0456626, 0.18588,
      0.0399869, -0.0631275, -0.0704844, -0.231899, 0.111054, -0.189161,
      0.000129868, 0.0333718, 0.436683, -0.009952, 0.337332, 0.251625, 0.442357,
      0.421824, 0.329076, 0.0456626, 0.12757, -0.140984, 0.0631275, 0.0704844,
      -0.0197737, 0.256357, 0.189161, -0.000129868, -0.0333718, 0.209723,
      -0.383155, -0.337332, -0.251625, 0.586489, -0.172181, -0.329076,
      -0.0456626, -0.0583106, -0.180971, -0.0631275, -0.0704844, 0.212125,
      0.145303, -0.189161, 0.000129868, -0.00177478, 0.0553645, -0.00126176,
      -0.0164247, 0.23154, -0.262519, -0.250334, -0.0135392, -0.429472, 0.45567,
      -0.503583, -0.223493, -0.211802, -0.020461, 0.265268, 0.0023362,
      -0.00241145, 0.294363, -0.686239, 0.0156394, 0.204055, -0.360136,
      0.267096, 0.254698, 0.074687, -0.0228668, 0.132236, -0.14614, -0.174986,
      -0.185046, -0.0109958, 0.142556, 0.0661743, 0.0022999, -0.00177478,
      -0.0265895, 0.0485779, -0.0164247, 0.23154, 0.348055, -0.102182,
      -0.0135392, -0.429472, 0.208281, 0.646413, -0.223493, -0.211802,
      -0.219498, -0.150354, 0.0023362, -0.00241145, 0.294363, 0.329576,
      -0.60212, 0.204055, -0.360136, -0.354123, 0.103963, 0.074687, -0.0228668,
      0.0604434, 0.18759, -0.174986, -0.185046, -0.117959, -0.0808008,
      0.0661743, 0.0022999, -0.00177478, -0.028775, -0.0473162, -0.0164247,
      0.23154, -0.0855356, 0.352515, -0.0135392, -0.429472, -0.663951, -0.14283,
      -0.223493, -0.211802, 0.239959, -0.114914, 0.0023362, -0.00241145,
      0.294363, 0.356664, 0.586481, 0.204055, -0.360136, 0.0870267, -0.358661,
      0.074687, -0.0228668, -0.192679, -0.0414494, -0.174986, -0.185046,
      0.128955, -0.0617554, 0.0661743, 0.0022999, 0.00741062, -3.87173e-16,
      -4.31863e-16, 0.0468488, -0.0476991, 7.27357e-16, 1.23654e-15, -0.43422,
      -0.159247, -4.34945e-17, 1.2743e-16, 0.503528, -0.228856, -7.97629e-17,
      -2.53026e-16, 0.689669, -0.00301027, 0.173046, 7.91486e-15, 8.39419e-15,
      -0.717804, -0.0195249, -1.10754e-15, -1.66789e-15, 0.551371, 0.0684292,
      4.15572e-17, -1.84233e-16, 0.0105378, -0.148396, -1.63792e-16,
      -4.6499e-16, 0.351571, 0.00210309;

  Logger log;
  TCMatrix_gwbse Mmn(log);
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, eigenvectors);

  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(eigenvals);
  rpa.configure(4, 0, 16);

  Eigen::VectorXd rpa_omega_ref = Eigen::VectorXd(60);
  rpa_omega_ref << 0.104192, 0.104192, 0.187814, 0.559693, 0.559693, 0.572575,
      0.577988, 0.577989, 0.579088, 0.618403, 0.618403, 0.67005, 0.678538,
      0.678538, 0.722771, 1.10797, 1.41413, 1.41413, 1.58866, 1.60381, 1.60381,
      1.64709, 1.87331, 1.87331, 1.88646, 1.8926, 1.89268, 1.89268, 1.933,
      1.933, 2.01832, 2.40974, 2.42192, 2.42192, 2.46371, 2.85829, 2.8853,
      2.8853, 2.90367, 2.90367, 2.92541, 2.94702, 3.3382, 3.3382, 3.35102,
      3.3566, 3.3566, 3.35835, 3.39617, 3.39617, 4.22882, 4.71607, 4.72233,
      4.72233, 4.76567, 16.5917, 17.0793, 17.093, 17.093, 17.1377;

  Eigen::VectorXd rpa_XpY_diag_ref = Eigen::VectorXd(60);
  rpa_XpY_diag_ref << 0.00898067, -0.00898068, 0.00228621, -7.40059e-10,
      -0.000447591, -7.21786e-10, 2.86571e-09, 1.16917e-08, 5.18349e-09,
      0.000124051, 1.89797e-10, -1.549e-05, -0.00513947, -0.00513956,
      -7.71176e-08, -2.09605e-08, 0.00412762, 0.00412738, 6.253e-09,
      3.93273e-05, 2.09406e-05, -0.000532391, 3.70347e-05, 2.0857e-06,
      -6.13653e-10, 0.000723575, 3.41096e-05, 0.00942648, 0.0497694, 0.0497693,
      -1.28726e-07, 3.68417e-09, -7.00924e-06, 7.01603e-06, -1.36252e-09,
      -1.08353e-09, 0.000549202, 0.000549203, -6.22315e-09, -3.83321e-08,
      -2.4032e-08, -2.19679e-09, 5.76794e-10, -1.46987e-08, -1.27183e-08,
      -0.138355, -6.49432e-09, -2.71175e-05, 2.69024e-05, -2.69024e-05,
      0.000155952, 0.000311876, -0.0002498, -0.000249804, 0.000816289,
      1.05559e-06, 1.84996e-12, 5.76753e-06, -9.90455e-11, -0.000800706;

  RPA::rpa_eigensolution sol = rpa.Diagonalize_H2p();

  bool check_rpa_eigenvalues = rpa_omega_ref.isApprox(sol.omega, 0.0001);
  if (!check_rpa_eigenvalues) {
    cout << "rpa_omega" << endl;
    cout << sol.omega << endl;
    cout << "rpa_omega_ref" << endl;
    cout << rpa_omega_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_rpa_eigenvalues, 1);

  Eigen::VectorXd rpa_XpY_diag = sol.XpY.diagonal();

  bool check_rpa_XpY_diag =
      rpa_XpY_diag_ref.cwiseAbs().isApprox(rpa_XpY_diag.cwiseAbs(), 0.0001);
  if (!check_rpa_XpY_diag) {
    cout << "rpa_XpY_diag" << endl;
    cout << rpa_XpY_diag.transpose() << endl;
    cout << "rpa_XpY_diag_ref" << endl;
    cout << rpa_XpY_diag_ref.transpose() << endl;
  }
  BOOST_CHECK_EQUAL(check_rpa_XpY_diag, 1);
}

BOOST_AUTO_TEST_SUITE_END()