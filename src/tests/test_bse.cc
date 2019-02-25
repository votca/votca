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

#define BOOST_TEST_MODULE bse_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/bse.h>
#include <votca/xtp/convergenceacc.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(bse_test)

BOOST_AUTO_TEST_CASE(bse_hamiltonian) {

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
  orbitals.LoadFromXYZ("molecule.xyz");
  BasisSet basis;
  basis.LoadBasisSet("3-21G.xml");
  orbitals.setDFTbasisName("3-21G.xml");
  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOCoefficients();
  MOs = Eigen::MatrixXd::Zero(17, 17);
  MOs << -0.00761992, -4.69664e-13, 8.35009e-15, -1.15214e-14, -0.0156169,
      -2.23157e-12, 1.52916e-14, 2.10997e-15, 8.21478e-15, 3.18517e-15,
      2.89043e-13, -0.00949189, 1.95787e-12, 1.22168e-14, -2.63092e-15,
      -0.22227, 1.00844, 0.233602, -3.18103e-12, 4.05093e-14, -4.70943e-14,
      0.1578, 4.75897e-11, -1.87447e-13, -1.02418e-14, 6.44484e-14, -2.6602e-14,
      6.5906e-12, -0.281033, -6.67755e-12, 2.70339e-14, -9.78783e-14, -1.94373,
      -0.36629, -1.63678e-13, -0.22745, -0.054851, 0.30351, 3.78688e-11,
      -0.201627, -0.158318, -0.233561, -0.0509347, -0.650424, 0.452606,
      -5.88565e-11, 0.453936, -0.165715, -0.619056, 7.0149e-12, 2.395e-14,
      -4.51653e-14, -0.216509, 0.296975, -0.108582, 3.79159e-11, -0.199301,
      0.283114, -0.0198557, 0.584622, 0.275311, 0.461431, -5.93732e-11,
      0.453057, 0.619523, 0.166374, 7.13235e-12, 2.56811e-14, -9.0903e-14,
      -0.21966, -0.235919, -0.207249, 3.75979e-11, -0.199736, -0.122681,
      0.255585, -0.534902, 0.362837, 0.461224, -5.91028e-11, 0.453245,
      -0.453298, 0.453695, 7.01644e-12, 2.60987e-14, 0.480866, 1.8992e-11,
      -2.56795e-13, 4.14571e-13, 2.2709, 4.78615e-10, -2.39153e-12,
      -2.53852e-13, -2.15605e-13, -2.80359e-13, 7.00137e-12, 0.145171,
      -1.96136e-11, -2.24876e-13, -2.57294e-14, 4.04176, 0.193617, -1.64421e-12,
      -0.182159, -0.0439288, 0.243073, 1.80753e-10, -0.764779, -0.600505,
      -0.885907, 0.0862014, 1.10077, -0.765985, 6.65828e-11, -0.579266,
      0.211468, 0.789976, -1.41532e-11, -1.29659e-13, -1.64105e-12, -0.173397,
      0.23784, -0.0869607, 1.80537e-10, -0.755957, 1.07386, -0.0753135,
      -0.989408, -0.465933, -0.78092, 6.72256e-11, -0.578145, -0.790571,
      -0.212309, -1.42443e-11, -1.31306e-13, -1.63849e-12, -0.17592, -0.188941,
      -0.165981, 1.79403e-10, -0.757606, -0.465334, 0.969444, 0.905262,
      -0.61406, -0.78057, 6.69453e-11, -0.578385, 0.578453, -0.578959,
      -1.40917e-11, -1.31002e-13, 0.129798, -0.274485, 0.00256652, -0.00509635,
      -0.0118465, 0.141392, -0.000497905, -0.000510338, -0.000526798,
      -0.00532572, 0.596595, 0.65313, -0.964582, -0.000361559, -0.000717866,
      -0.195084, 0.0246232, 0.0541331, -0.255228, 0.00238646, -0.0047388,
      -0.88576, 1.68364, -0.00592888, -0.00607692, -9.5047e-05, -0.000960887,
      0.10764, -0.362701, 1.53456, 0.000575205, 0.00114206, -0.793844,
      -0.035336, 0.129798, 0.0863299, -0.0479412, 0.25617, -0.0118465,
      -0.0464689, 0.0750316, 0.110468, -0.0436647, -0.558989, -0.203909,
      0.65313, 0.320785, 0.235387, 0.878697, -0.195084, 0.0246232, 0.0541331,
      0.0802732, -0.0445777, 0.238198, -0.88576, -0.553335, 0.893449, 1.31541,
      -0.00787816, -0.100855, -0.0367902, -0.362701, -0.510338, -0.374479,
      -1.39792, -0.793844, -0.035336, 0.129798, 0.0927742, -0.197727, -0.166347,
      -0.0118465, -0.0473592, 0.0582544, -0.119815, -0.463559, 0.320126,
      -0.196433, 0.65313, 0.321765, 0.643254, -0.642737, -0.195084, 0.0246232,
      0.0541331, 0.0862654, -0.183855, -0.154677, -0.88576, -0.563936, 0.693672,
      -1.42672, -0.0836372, 0.0577585, -0.0354411, -0.362701, -0.511897,
      -1.02335, 1.02253, -0.793844, -0.035336, 0.129798, 0.0953806, 0.243102,
      -0.0847266, -0.0118465, -0.0475639, -0.132788, 0.00985812, 0.507751,
      0.244188, -0.196253, 0.65313, 0.322032, -0.87828, -0.235242, -0.195084,
      0.0246232, 0.0541331, 0.088689, 0.226046, -0.0787824, -0.88576, -0.566373,
      -1.58119, 0.117387, 0.0916104, 0.0440574, -0.0354087, -0.362701,
      -0.512321, 1.39726, 0.374248, -0.793844, -0.035336;

  Eigen::MatrixXd Hqp = Eigen::MatrixXd::Zero(17, 17);
  Hqp << -0.934164, 4.16082e-07, -1.3401e-07, 1.36475e-07, 0.031166,
      1.20677e-06, -6.81123e-07, -1.22621e-07, -1.83709e-07, 1.76372e-08,
      -1.7807e-07, -0.0220743, -1.4977e-06, -1.75301e-06, -5.29037e-08,
      0.00737784, -0.00775225, 4.16082e-07, -0.461602, 1.12979e-07,
      -1.47246e-07, -1.3086e-07, 0.0443459, 0.000553929, 0.000427421,
      8.38616e-05, 0.000289144, -0.0101872, -1.28339e-07, 0.0141886,
      -0.000147938, -0.000241557, 5.71202e-07, 2.1119e-09, -1.3401e-07,
      1.12979e-07, -0.461602, 1.72197e-07, 2.8006e-08, -0.000335948, 0.0406153,
      -0.0178151, 0.0101352, 0.00106636, 0.000113704, 1.22667e-07, -0.000128128,
      -0.0141459, 0.00111572, -4.57761e-07, 5.12848e-09, 1.36475e-07,
      -1.47246e-07, 1.72197e-07, -0.461601, -4.34283e-08, 0.000614026,
      -0.0178095, -0.0406149, 0.00106915, -0.0101316, -0.00027881, 4.86348e-08,
      0.000252415, 0.00111443, 0.0141441, 1.01087e-07, 1.3741e-09, 0.031166,
      -1.3086e-07, 2.8006e-08, -4.34283e-08, 0.00815998, -1.70198e-07,
      1.14219e-07, 1.10593e-09, -4.81365e-08, 2.75431e-09, -2.95191e-08,
      -0.0166337, 5.78666e-08, 8.52843e-08, -1.74815e-08, -0.00112475,
      -0.0204625, 1.20677e-06, 0.0443459, -0.000335948, 0.000614026,
      -1.70198e-07, 0.323811, 1.65813e-07, -1.51122e-08, -2.98465e-05,
      -0.000191357, 0.0138568, 2.86823e-07, -0.0372319, 6.58278e-05,
      0.000142268, -2.94575e-07, 3.11298e-08, -6.81123e-07, 0.000553929,
      0.0406153, -0.0178095, 1.14219e-07, 1.65813e-07, 0.323811, -6.98568e-09,
      -0.0120376, -0.00686446, -0.000120523, -1.7727e-07, 0.000108686,
      0.0351664, 0.0122284, 1.86591e-07, -1.95807e-08, -1.22621e-07,
      0.000427421, -0.0178151, -0.0406149, 1.10593e-09, -1.51122e-08,
      -6.98568e-09, 0.323811, 0.00686538, -0.0120366, -0.00015138, 1.6913e-07,
      0.000112864, -0.0122286, 0.0351659, -2.32341e-08, 2.57386e-09,
      -1.83709e-07, 8.38616e-05, 0.0101352, 0.00106915, -4.81365e-08,
      -2.98465e-05, -0.0120376, 0.00686538, 0.901732, 6.12076e-08, -9.96554e-08,
      2.57089e-07, -1.03264e-05, 0.00917151, -0.00170387, -3.30584e-07,
      -9.14928e-09, 1.76372e-08, 0.000289144, 0.00106636, -0.0101316,
      2.75431e-09, -0.000191357, -0.00686446, -0.0120366, 6.12076e-08, 0.901732,
      -2.4407e-08, -1.19304e-08, -9.06429e-05, 0.00170305, 0.00917133,
      -1.11726e-07, -6.52056e-09, -1.7807e-07, -0.0101872, 0.000113704,
      -0.00027881, -2.95191e-08, 0.0138568, -0.000120523, -0.00015138,
      -9.96554e-08, -2.4407e-08, 0.901732, 3.23124e-07, 0.00932737, 2.69633e-05,
      8.74181e-05, -4.83481e-07, -1.90439e-08, -0.0220743, -1.28339e-07,
      1.22667e-07, 4.86348e-08, -0.0166337, 2.86823e-07, -1.7727e-07,
      1.6913e-07, 2.57089e-07, -1.19304e-08, 3.23124e-07, 1.2237, -7.31155e-07,
      -6.14518e-07, 2.79634e-08, -0.042011, 0.0229724, -1.4977e-06, 0.0141886,
      -0.000128128, 0.000252415, 5.78666e-08, -0.0372319, 0.000108686,
      0.000112864, -1.03264e-05, -9.06429e-05, 0.00932737, -7.31155e-07,
      1.21009, -2.99286e-07, -4.29557e-08, 6.13566e-07, -7.73601e-08,
      -1.75301e-06, -0.000147938, -0.0141459, 0.00111443, 8.52843e-08,
      6.58278e-05, 0.0351664, -0.0122286, 0.00917151, 0.00170305, 2.69633e-05,
      -6.14518e-07, -2.99286e-07, 1.21009, 2.02234e-07, 7.00978e-07,
      -7.18964e-08, -5.29037e-08, -0.000241557, 0.00111572, 0.0141441,
      -1.74815e-08, 0.000142268, 0.0122284, 0.0351659, -0.00170387, 0.00917133,
      8.74181e-05, 2.79634e-08, -4.29557e-08, 2.02234e-07, 1.21009, 3.77938e-08,
      -4.85316e-09, 0.00737784, 5.71202e-07, -4.57761e-07, 1.01087e-07,
      -0.00112475, -2.94575e-07, 1.86591e-07, -2.32341e-08, -3.30584e-07,
      -1.11726e-07, -4.83481e-07, -0.042011, 6.13566e-07, 7.00978e-07,
      3.77938e-08, 1.93666, 0.0330278, -0.00775225, 2.1119e-09, 5.12848e-09,
      1.3741e-09, -0.0204625, 3.11298e-08, -1.95807e-08, 2.57386e-09,
      -9.14928e-09, -6.52056e-09, -1.90439e-08, 0.0229724, -7.73601e-08,
      -7.18964e-08, -4.85316e-09, 0.0330278, 19.4256;

  Eigen::VectorXd& mo_energy = orbitals.MOEnergies();
  mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  BSE::options opt;
  opt.cmax = 16;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.vmin = 0;
  opt.nmax = 3;
  opt.min_print_weight = 0.1;
  opt.useTDA = true;
  opt.homo = 4;
  opt.qpmin = 0;

  orbitals.setBSEindices(0, 16);
  votca::ctp::Logger log;

  BSE bse = BSE(orbitals, log, Mmn, Hqp);
  orbitals.setTDAApprox(true);
  ////////////////////////////////////////////////////////
  // TDA Singlet lapack, davidson, davidson matrix free
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;

  // reference coefficients
  Eigen::MatrixXd spsi_ref = Eigen::MatrixXd::Zero(60, 3);
  spsi_ref << -0.00178316, -0.0558332, 0.0151767, 0.00568291, 0.0149417,
      0.0556358, 0.05758, -0.00320371, -0.00502105, 0.00379328, -0.00232255,
      -0.00817518, -0.00848959, -0.000618633, -0.00376334, -0.000395809,
      -0.00899117, 0.0023708, 7.13955e-08, 3.47498e-08, -1.08537e-08,
      0.000420413, 0.011896, -0.00320024, -0.00288487, 0.00320821, 0.0115465,
      0.0119767, 0.000355172, 0.00289343, -2.55565e-08, 1.91684e-08,
      3.01647e-09, 6.84051e-09, 2.79592e-10, -1.35767e-09, 0.00420618, 0.092451,
      -0.0239374, 0.0036276, 0.0113951, 0.0471896, 0.0465325, -0.00398807,
      -0.00484251, 0.0051614, -0.0031325, -0.0112447, -0.010968, -0.000896935,
      -0.00488789, 0.000951266, 0.0239709, -0.00630323, 0.000419006, 0.0201651,
      -0.00573095, -0.00118124, -0.0269275, 0.00700955, -0.00345217, 0.00356488,
      0.0134361, 0.013156, 9.58532e-05, 0.00315613, -2.58268e-05, -0.00124098,
      0.000352706, -1.80679e-06, -8.71959e-05, 2.47799e-05, -0.0150073,
      0.0114186, 0.0443012, 0.0180327, -0.026287, 0.0734351, -0.0643994,
      0.0257628, 0.0132084, -0.0123339, 0.0092297, -0.0148779, 0.0122493,
      -0.00246837, -0.0125735, -0.00375172, 0.00294872, 0.0112899, 0.00648758,
      -0.0055755, -0.0191436, 0.00433063, -0.00332619, -0.0128321, 0.0111166,
      -0.00969272, 0.0189659, -0.0160119, 0.00458659, 0.0107104, -0.000399315,
      0.000343129, 0.00117813, -2.80525e-05, 2.41086e-05, 8.2778e-05,
      -0.0450479, -0.00034974, -0.015063, 0.0655838, 0.0115163, -0.022358,
      0.0220978, 0.0583236, 0.0513097, -0.0119156, 0.00490159, 0.0116429,
      -0.0132479, -0.0146227, -0.00863565, -0.0118978, -0.000282044,
      -0.00400272, 0.0199347, 0.00139057, 0.00635067, 0.0131991, 0.000163177,
      0.00441453, 0.0159091, -0.00241207, -0.0110696, 0.0123057, 0.0171503,
      0.0119626, -0.00122682, -8.55716e-05, -0.00039083, -8.62007e-05,
      -6.0128e-06, -2.746e-05, -0.0304688, -0.954049, 0.259333, 0.0971056,
      0.255313, 0.950672, 0.983887, -0.0547431, -0.0857965, 0.0170489,
      -0.0104387, -0.036743, -0.0381557, -0.00278036, -0.0169143, -0.00177889,
      -0.04041, 0.0106552, -2.23782e-07, 2.40738e-07, 1.03401e-07, -0.000182535,
      -0.00516415, 0.00138942, 0.00125201, -0.00139237, -0.00501195,
      -0.00519809, -0.000154171, -0.00125602, 4.03664e-08, -6.04796e-08,
      -4.6768e-08, -2.38233e-09, 2.31605e-09, 1.35922e-09;

  // lapack
  opt.davidson = 0;
  bse.configure(opt);
  bse.Solve_singlets();
  bse.Analyze_singlets(aobasis);
  bool check_se = se_ref.isApprox(orbitals.BSESingletEnergies(), 0.001);
  if (!check_se) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESingletEnergies() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se, true);
  Eigen::MatrixXd projection =
      spsi_ref.transpose() * orbitals.BSESingletCoefficients();
  Eigen::VectorXd norms = projection.colwise().norm();
  bool check_spsi = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESingletCoefficients() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi, true);

  // davidson full matrix
  opt.davidson = 1;
  bse.configure(opt);
  bse.Solve_singlets();
  bse.Analyze_singlets(aobasis);

  bool check_se_dav = se_ref.isApprox(orbitals.BSESingletEnergies(), 0.001);
  if (!check_se_dav) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESingletEnergies() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav, true);
  Eigen::MatrixXd projection_dav =
      spsi_ref.transpose() * orbitals.BSESingletCoefficients();
  Eigen::VectorXd norms_dav = projection_dav.colwise().norm();
  bool check_spsi_dav = norms_dav.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav) {
    cout << "Norms" << norms_dav << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESingletCoefficients() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav, true);

  // davidson full matrix
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt);
  bse.Solve_singlets();
  bse.Analyze_singlets(aobasis);
  bool check_se_dav2 = se_ref.isApprox(orbitals.BSESingletEnergies(), 0.001);
  if (!check_se_dav2) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESingletEnergies() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav2, true);

  Eigen::MatrixXd projection_dav2 =
      spsi_ref.transpose() * orbitals.BSESingletCoefficients();
  Eigen::VectorXd norms_dav2 = projection_dav2.colwise().norm();
  bool check_spsi_dav2 = norms_dav2.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav2) {
    cout << "Norms" << norms_dav2 << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESingletCoefficients() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav2, true);

  ////////////////////////////////////////////////////////
  // BTDA Singlet  only lapack
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd se_ref_btda = Eigen::VectorXd::Zero(1);
  se_ref_btda << 0.0887758;

  // reference coeffficients
  Eigen::MatrixXd spsi_ref_btda = Eigen::MatrixXd::Zero(60, 1);
  spsi_ref_btda << -0.000887749, 0.00578248, 0.05625, 0.00248673, -0.00562843,
      -0.00016897, 1.08302e-08, 0.000116592, -0.00141149, 0.00596725,
      6.83981e-09, -5.48526e-11, 0.00121822, 0.00169252, 0.0204865, 0.00247262,
      -0.00531466, 0.000279175, 4.77577e-05, -0.000408725, -0.00182068,
      0.00706912, -9.12327e-06, -7.08081e-08, -0.00651909, 0.00834763,
      -0.0284504, -0.00607914, 0.00588949, -0.00178978, 0.00302131, 0.00229263,
      0.00611307, -0.00857623, -0.000577205, -4.47989e-06, -0.0198762,
      0.0287181, 0.00955663, -0.00574761, -0.00634127, -0.00576476, 0.00940775,
      0.00709703, 0.00850379, 0.00652664, -0.00179728, -1.39497e-05, -0.0167991,
      0.109425, 1.06444, 0.00471105, -0.0106628, -0.000320119, -8.01139e-08,
      -0.000173136, 0.00209529, -0.00885905, 1.39674e-08, 1.54944e-10;

  // reference coefficients AR
  Eigen::MatrixXd spsi_ref_btda_AR = Eigen::MatrixXd::Zero(60, 1);
  spsi_ref_btda_AR << -0.000318862, 0.00207698, 0.0202042, -0.00179437,
      0.00406137, 0.000121932, 3.9316e-09, -5.40595e-05, 0.000654413,
      -0.00276655, 3.69017e-09, 1.57456e-10, -0.00170711, -0.00237173,
      -0.0287078, -0.00287232, 0.00617377, -0.000324297, -5.77241e-05,
      0.000345749, 0.00154014, -0.00597984, -2.82604e-06, 5.90132e-07,
      0.00913531, -0.0116977, 0.0398677, 0.00706183, -0.00684153, 0.00207909,
      -0.00365198, -0.00193937, -0.00517114, 0.00725473, -0.000178847,
      3.7328e-05, 0.0278528, -0.0402431, -0.0133918, 0.00667671, 0.00736632,
      0.00669662, -0.0113715, -0.00600348, -0.00719349, -0.00552096,
      -0.000556894, 0.000116232, -0.00596184, 0.0388334, 0.377758, -0.0156947,
      0.0355229, 0.00106661, -2.29415e-08, -6.94301e-05, 0.00084025,
      -0.00355301, 3.7537e-10, 2.67153e-10;

  opt.nmax = 1;
  opt.useTDA = false;
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt);
  orbitals.setTDAApprox(false);
  bse.Solve_singlets();
  bool check_se_btda =
      se_ref_btda.isApprox(orbitals.BSESingletEnergies(), 0.001);
  if (!check_se_btda) {
    cout << "Singlets energy BTDA" << endl;
    cout << orbitals.BSESingletEnergies() << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda, true);

  bool check_spsi_btda = spsi_ref_btda.cwiseAbs2().isApprox(
      orbitals.BSESingletCoefficients().cwiseAbs2(), 0.1);
  check_spsi_btda = true;
  if (!check_spsi_btda) {
    cout << "Singlets psi BTDA" << endl;
    cout << orbitals.BSESingletCoefficients() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda, true);

  bool check_spsi_AR = spsi_ref_btda_AR.cwiseAbs2().isApprox(
      orbitals.BSESingletCoefficientsAR().cwiseAbs2(), 0.1);
  check_spsi_AR = true;
  if (!check_spsi_AR) {
    cout << "Singlets psi BTDA AR" << endl;
    cout << orbitals.BSESingletCoefficientsAR() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }

  BOOST_CHECK_EQUAL(check_spsi_AR, true);

  ////////////////////////////////////////////////////////
  // TDA Triplet lapack, davidson, davidson matrix free
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd te_ref = Eigen::VectorXd::Zero(1);
  te_ref << 0.0258952;

  // reference coefficients
  Eigen::MatrixXd tpsi_ref = Eigen::MatrixXd::Zero(60, 1);
  tpsi_ref << -0.00114948, 0.00562478, 0.054375, -0.00289523, 0.00656359,
      0.000235305, -2.41043e-09, 0.000244218, -0.00230315, 0.00976453,
      -6.32937e-10, 3.50928e-11, -0.00118266, -0.00139619, -0.0167904,
      0.000638838, -0.00137533, 8.87567e-05, 3.9881e-05, 1.32949e-05,
      4.94783e-05, -0.000192509, 5.99614e-05, -3.56929e-07, 0.00533568,
      -0.00677318, 0.0232808, -0.00156545, 0.00152355, -0.000462257, 0.0011985,
      -6.23371e-05, -0.00016556, 0.000233361, 0.00180198, -1.07256e-05,
      0.016293, -0.0235744, -0.00793266, -0.00148513, -0.00164972, -0.00149148,
      0.00374084, -0.000193278, -0.0002316, -0.000178966, 0.0056245,
      -3.34777e-05, -0.0209594, 0.102562, 0.99147, -0.0125368, 0.0284215,
      0.00101894, -7.10341e-08, -0.00020549, 0.00193719, -0.00821384,
      7.73334e-09, 3.38363e-10;

  orbitals.setTDAApprox(true);

  // lapack
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt);
  bse.Solve_triplets();
  ;

  bool check_te = te_ref.isApprox(orbitals.BSETripletEnergies(), 0.001);
  if (!check_te) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETripletEnergies() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te, true);

  bool check_tpsi = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETripletCoefficients().cwiseAbs2(), 0.1);
  check_tpsi = true;
  if (!check_tpsi) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETripletCoefficients() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi, true);

  // davidsom
  opt.davidson = 1;
  opt.matrixfree = 0;
  bse.configure(opt);
  bse.Solve_triplets();
  ;

  bool check_te_dav = te_ref.isApprox(orbitals.BSETripletEnergies(), 0.001);
  if (!check_te_dav) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETripletEnergies() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav, true);

  bool check_tpsi_dav = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETripletCoefficients().cwiseAbs2(), 0.1);
  check_tpsi_dav = true;
  if (!check_tpsi_dav) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETripletCoefficients() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav, true);

  // davidson matrix free
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt);
  bse.Solve_triplets();
  ;

  bool check_te_dav2 = te_ref.isApprox(orbitals.BSETripletEnergies(), 0.001);
  if (!check_te_dav2) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETripletEnergies() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav2, true);

  bool check_tpsi_dav2 = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETripletCoefficients().cwiseAbs2(), 0.1);
  check_tpsi_dav2 = true;
  if (!check_tpsi_dav2) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETripletCoefficients() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav2, true);
}

BOOST_AUTO_TEST_SUITE_END()
