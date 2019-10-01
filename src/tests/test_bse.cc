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
#include <votca/xtp/qmfragment.h>

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
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  BasisSet basis;
  basis.Load("3-21G.xml");
  orbitals.setDFTbasisName("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
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

  Eigen::VectorXd& mo_energy = orbitals.MOs().eigenvalues();
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
  Logger log;

  BSE bse = BSE(log, Mmn, Hqp);
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
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_singlets(orbitals);
  bool check_se = se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se, true);
  Eigen::MatrixXd projection =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms = projection.colwise().norm();
  bool check_spsi = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi, true);

  // davidson full matrix
  opt.davidson = 1;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_singlets(orbitals);

  std::vector<QMFragment<BSE_Population> > singlets;
  bse.Analyze_singlets(singlets, orbitals);

  bool check_se_dav =
      se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_dav) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav, true);
  Eigen::MatrixXd projection_dav =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_dav = projection_dav.colwise().norm();
  bool check_spsi_dav = norms_dav.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav) {
    cout << "Norms" << norms_dav << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav, true);

  // davidson matrix free
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_singlets(orbitals);
  bool check_se_dav2 =
      se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_dav2) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav2, true);

  Eigen::MatrixXd projection_dav2 =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_dav2 = projection_dav2.colwise().norm();
  bool check_spsi_dav2 = norms_dav2.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav2) {
    cout << "Norms" << norms_dav2 << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav2, true);

  ////////////////////////////////////////////////////////
  // BTDA Singlet  lanczos and lapack
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd se_ref_btda = Eigen::VectorXd::Zero(3);
  se_ref_btda << 0.0887758, 0.0887758, 0.0887758;

  // reference coeffficients
  Eigen::MatrixXd spsi_ref_btda = Eigen::MatrixXd::Zero(60, 3);
  spsi_ref_btda <<    0.00228773,   -0.0572672,    0.0185901,
                     -0.00606479,    0.0182899,    0.0570883,
                      -0.0599032,  -0.00403878,  -0.00506982,
                     -0.00239507,  -0.00174279,  -0.00511526,
                        0.005396, -0.000467089,   -0.0023674,
                     0.000293773,  -0.00562893,   0.00178028,
                    -1.18994e-08,   3.6883e-08,  1.39133e-09,
                     -0.00038655,   0.00890738,  -0.00286588,
                      0.00216846,   0.00287559,   0.00864496,
                     -0.00910256,  0.000306734,   0.00218106,
                    -6.67006e-09,  1.50464e-08, -4.23043e-09,
                     5.73874e-11,  1.02039e-09,  3.42402e-10,
                     -0.00300256,    0.0559254,   -0.0174141,
                     -0.00229461,   0.00836668,    0.0285808,
                       -0.028529,  -0.00280458,  -0.00285455,
                     -0.00278921,  -0.00201767,  -0.00602331,
                      0.00596683, -0.000570189,  -0.00262914,
                     -0.00061392,    0.0128447,  -0.00405267,
                    -0.000330257,     0.011532,  -0.00389023,
                     0.000903715,   -0.0173665,   0.00543237,
                      0.00223576,   0.00277158,    0.0086669,
                     -0.00860769,  9.71086e-05,   0.00205762,
                     5.28565e-05,  -0.00184564,  0.000622629,
                     6.78691e-07, -2.36991e-05,  7.99472e-06,
                      0.00913939,   0.00838224,    0.0268121,
                      -0.0109663,   -0.0140137,    0.0457817,
                       0.0393868,      0.01648,   0.00733715,
                       0.0066705,   0.00467073,  -0.00831923,
                      -0.0066442,  -0.00172131,  -0.00676294,
                      0.00202369,   0.00191325,   0.00605252,
                     -0.00373624,  -0.00379982,   -0.0109468,
                     -0.00281175,  -0.00260087,  -0.00827985,
                     -0.00723077,  -0.00578823,    0.0126979,
                       0.0104477,   0.00341532,   0.00685506,
                     0.000597973,  0.000608139,   0.00175196,
                     7.67794e-06,  7.80883e-06,  2.24965e-05,
                       0.0276451, -0.000451369,  -0.00921566,
                      -0.0402838,   0.00608812,    -0.014051,
                      -0.0138464,    0.0371613,    0.0297107,
                      0.00645305,   0.00301651,   0.00619365,
                      0.00726364,   -0.0081113,  -0.00430282,
                      0.00646942, -0.000209234,  -0.00216462,
                      -0.0115828,  0.000896897,   0.00364201,
                     -0.00863663,  0.000180708,   0.00287774,
                       -0.010391,  -0.00200537,  -0.00715858,
                     -0.00813919,    0.0115232,   0.00727228,
                      0.00185377, -0.000143538, -0.000582886,
                     2.38027e-05, -1.84326e-06, -7.48457e-06,
                       0.0391594,    -0.980256,     0.318212,
                       -0.103812,      0.31307,     0.977192,
                        -1.02537,   -0.0691326,   -0.0867812,
                     -0.00843227,  -0.00613583,   -0.0180091,
                       0.0189974,  -0.00164438,  -0.00833481,
                      0.00103433,   -0.0198172,   0.00626757,
                     7.10218e-08,  2.05723e-07,  5.20066e-08,
                     0.000176091,  -0.00405684,   0.00130541,
                    -0.000987275,  -0.00130929,  -0.00393686,
                      0.00414458, -0.000139691, -0.000993321,
                    -9.14836e-09, -7.65102e-08, -3.81041e-08,
                    -1.73462e-10,  2.76218e-09,  6.98481e-10;

  // // reference coefficients AR
  Eigen::MatrixXd spsi_ref_btda_AR = Eigen::MatrixXd::Zero(60, 3);
  spsi_ref_btda_AR <<  0.000584144,   -0.0146225,   0.00474677,
                       -0.00154856,   0.00466997,    0.0145766,
                        -0.0152954,  -0.00103125,  -0.00129451,
                         0.0017559,   0.00127769,   0.00375014,
                       -0.00395598,  0.000342451,   0.00173561,
                      -0.000215381,   0.00412671,  -0.00130515,
                      -5.02586e-09,  4.52834e-08,  1.34149e-08,
                       0.000156138,  -0.00359781,   0.00115758,
                      -0.000875892,  -0.00116148,  -0.00349181,
                        0.00367662, -0.000123883, -0.000880951,
                      -5.20546e-09,  2.56561e-08,   6.1721e-09,
                      -1.30377e-10, -3.42705e-10,  -1.4581e-10,
                        0.00225763,   -0.0420501,    0.0130936,
                        0.00172531,  -0.00629088,   -0.0214897,
                         0.0214509,   0.00210875,   0.00214631,
                        0.00443924,   0.00321126,   0.00958646,
                       -0.00949663,  0.000907487,   0.00418444,
                       0.000977068,    -0.020443,   0.00645007,
                       0.000411481,   -0.0143682,   0.00484706,
                      -0.000752112,    0.0144531,  -0.00452103,
                       -0.00186068,  -0.00230664,   -0.0072129,
                        0.00716363,  -8.0822e-05,  -0.00171243,
                       4.85758e-05,  -0.00169622,  0.000572232,
                      -4.64045e-06,  0.000162032, -5.46606e-05,
                       -0.00687191,   -0.0063025,   -0.0201599,
                        0.00824556,    0.0105368,    -0.034423,
                        -0.0296148,   -0.0123913,  -0.00551674,
                        -0.0106166,  -0.00743383,    0.0132406,
                         0.0105747,   0.00273958,    0.0107636,
                       -0.00322083,  -0.00304507,  -0.00963294,
                        0.00465513,   0.00473433,     0.013639,
                        0.00234006,   0.00216454,   0.00689084,
                        0.00601773,   0.00481718,   -0.0105676,
                       -0.00869493,  -0.00284236,  -0.00570504,
                       0.000549541,  0.000558892,    0.0016101,
                      -5.24965e-05, -5.33896e-05, -0.000153809,
                        -0.0207863,  0.000339395,    0.0069292,
                         0.0302893,  -0.00457762,    0.0105648,
                          0.010411,   -0.0279414,   -0.0223393,
                        -0.0102705,  -0.00480099,  -0.00985759,
                        -0.0115606,    0.0129097,   0.00684819,
                        -0.0102966,  0.000333027,   0.00344514,
                         0.0144314,  -0.00111747,  -0.00453772,
                        0.00718776, -0.000150389,  -0.00239498,
                        0.00864786,   0.00166895,   0.00595764,
                        0.00677374,     -0.00959,  -0.00605226,
                        0.00170364, -0.000131921, -0.000535686,
                      -0.000162745,   1.2602e-05,  5.11726e-05,
                        0.00995914,    -0.249299,    0.0809272,
                        -0.0264016,      0.07962,     0.248518,
                         -0.260775,   -0.0175818,   -0.0220701,
                         0.0109283,   0.00795199,    0.0233398,
                        -0.0246207,   0.00213124,    0.0108019,
                       -0.00134058,    0.0256833,  -0.00812285,
                       1.49164e-08, -4.69903e-09,  5.23125e-09,
                       6.14995e-07, -1.40366e-05,  4.50767e-06,
                      -3.26276e-06, -4.40871e-06, -1.35086e-05,
                       1.40497e-05, -4.71623e-07, -3.42445e-06,
                       1.99993e-09, -3.51811e-08, -4.56933e-09,
                      -2.19573e-10,  6.10086e-10, -2.35986e-11;


  opt.nmax = 3;
  opt.useTDA = false;
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  orbitals.setTDAApprox(false);
  bse.Solve_singlets(orbitals);

  bool check_se_btda =
      se_ref_btda.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_btda) {
    cout << "Singlets energy BTDA" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda, true);

  projection = spsi_ref_btda.transpose() * orbitals.BSESinglets().eigenvectors();
  norms = projection.colwise().norm();
  bool check_spsi_btda = norms.isApproxToConstant(1, 1e-5);

  //check_spsi_btda = true;
  if (!check_spsi_btda) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA (Lapack)" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda, true);

  projection = spsi_ref_btda_AR.transpose() * orbitals.BSESinglets().eigenvectors2();
  norms = projection.colwise().norm();
  bool check_spsi_btda_AR = norms.isApproxToConstant(1, 1e-5);
  
  //check_spsi_AR = true;
  if (!check_spsi_btda_AR) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA AR (Lapack)" << endl;
    cout << orbitals.BSESinglets().eigenvectors2() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_AR, true);



  // Davidson matrix free
  opt.matrixfree = 1;
  opt.davidson = 1;
  opt.nmax = 3; 

  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_singlets(orbitals);
  
  bool check_se_btda_mf =
      se_ref_btda.isApprox(orbitals.BSESinglets().eigenvalues().head(opt.nmax), 0.001);
  if (!check_se_btda_mf) {
    cout << "Singlets energy BTDA (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvalues().head(opt.nmax) << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda_mf, true);
  
  projection = spsi_ref_btda.transpose() * orbitals.BSESinglets().eigenvectors();
  norms = projection.colwise().norm();
  bool check_spsi_btda_mf = norms.isApproxToConstant(1, 1e-5);

  if (!check_spsi_btda_mf) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_mf, true);

  projection = spsi_ref_btda_AR.transpose() * orbitals.BSESinglets().eigenvectors2();
  norms = projection.colwise().norm();
  bool check_spsi_btda_AR_mf = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi_btda_AR_mf) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA AR (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvectors2() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_AR_mf, true);

  ////////////////////////////////////////////////////////
  // TDA Triplet lapack, davidson, davidson matrix free
  ////////////////////////////////////////////////////////

  // reference energy
  opt.nmax = 1;
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
  opt.useTDA = true;

  // lapack
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_triplets(orbitals);
  std::vector<QMFragment<BSE_Population> > triplets;
  bse.Analyze_triplets(triplets, orbitals);

  bool check_te = te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te, true);

  bool check_tpsi = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  check_tpsi = true;
  if (!check_tpsi) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi, true);

  // davidson
  opt.davidson = 1;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_triplets(orbitals);

  bool check_te_dav =
      te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te_dav) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav, true);

  bool check_tpsi_dav = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  if (!check_tpsi_dav) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav, true);

  // davidson matrix free
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt, orbitals.MOs().eigenvalues());
  bse.Solve_triplets(orbitals);

  bool check_te_dav2 =
      te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te_dav2) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav2, true);

  bool check_tpsi_dav2 = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  if (!check_tpsi_dav2) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav2, true);
}

BOOST_AUTO_TEST_SUITE_END()
