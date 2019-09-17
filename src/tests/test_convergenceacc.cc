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

#define BOOST_TEST_MODULE convergenceacc_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/convergenceacc.h>
#include <votca/xtp/orbitals.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(convergenceacc_test)

BOOST_AUTO_TEST_CASE(levelshift_test) {

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
  AOOverlap overlap;
  overlap.Fill(aobasis);

  AOKinetic kinetic;
  kinetic.Fill(aobasis);
  AOMultipole esp;
  esp.FillPotential(aobasis, orbitals.QMAtoms());
  Eigen::MatrixXd H = kinetic.Matrix() + esp.Matrix();
  ConvergenceAcc d;
  Orbitals orb;
  int occlevels = 5;
  Logger log;
  ConvergenceAcc::options opt;
  opt.mode = ConvergenceAcc::KSmode::closed;
  opt.levelshift = 0.1;
  opt.histlength = 10;
  opt.numberofelectrons = occlevels * 2;
  d.setLogger(&log);
  d.Configure(opt);
  d.setOverlap(overlap, 1e-8);
  orb.MOs() = d.SolveFockmatrix(H);

  Eigen::VectorXd mo_eng_ref = Eigen::VectorXd::Zero(17);
  mo_eng_ref << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094, -3.72889,
      -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048, -2.63214,
      -2.08206, -2.08206, -2.08206, -2.03268;

  bool eng_comp = mo_eng_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  if (!eng_comp) {
    std::cout << "result energies" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref energies" << std::endl;
    std::cout << mo_eng_ref << std::endl;
  }

  Eigen::MatrixXd mo_ref = Eigen::MatrixXd::Zero(17, 17);
  mo_ref << -0.996559, -0.223082, 4.81443e-15, 2.21045e-15, -6.16146e-17,
      -3.16966e-16, 5.46703e-18, -1.09681e-15, -0.0301914, 6.45993e-16,
      1.05377e-16, 3.41154e-16, -0.102052, -5.79826e-16, 9.38593e-16,
      -4.69346e-15, -0.111923, -0.0445146, 0.88316, -1.94941e-14, -8.67388e-15,
      -7.26679e-16, 1.16326e-14, -3.35886e-15, 2.37877e-14, 0.866126,
      3.2068e-15, 3.80914e-15, 3.24563e-15, -0.938329, -6.4404e-15, 1.10811e-14,
      -5.5056e-14, -1.28767, 8.15798e-17, 2.30849e-14, 1.04169, 0.117804,
      0.0951759, 0.179467, 0.147031, 0.39183, -1.02927e-14, 0.32699, -0.223689,
      -0.130009, 1.0375e-15, -0.0940179, 0.126956, 0.0122904, 1.41709e-15,
      4.60157e-17, -7.1203e-15, 0.143338, -0.980459, -0.355251, 0.41611,
      -0.10826, -0.149964, 2.41546e-16, 0.12214, -0.0512447, 0.39537,
      1.1054e-15, -0.0996828, -0.0636092, -0.105478, 5.10746e-15, -5.25872e-18,
      4.8424e-15, 0.0488925, 0.364515, -0.9863, 0.0447336, 0.417155, -0.177023,
      5.76117e-15, -0.228081, -0.348136, 0.0253377, -1.05286e-15, 0.079576,
      0.0703157, -0.117608, 5.31327e-15, 0.0472716, 0.235837, -3.58018e-15,
      -1.68354e-15, 2.3989e-15, -9.86879e-15, 4.52519e-15, -1.6106e-14,
      -0.599523, -1.31237e-14, -8.63443e-15, -8.61196e-15, 1.8333, 2.05275e-14,
      -3.9562e-14, 1.89874e-13, 4.24316, -2.74184e-16, -1.53939e-15, -0.162416,
      -0.0183675, -0.0148395, -0.151162, -0.123842, -0.330032, 1.10084e-15,
      -1.45092, 0.992556, 0.576875, -3.82954e-15, 0.604373, -0.816111,
      -0.0790061, -8.89474e-15, -2.24862e-16, 3.23655e-15, -0.0223487, 0.152869,
      0.0553894, -0.350483, 0.0911859, 0.126313, -5.48468e-15, -0.541961,
      0.227383, -1.75434, -3.89443e-15, 0.640788, 0.408897, 0.67804,
      -3.17156e-14, -2.81346e-17, -1.09423e-15, -0.00762313, -0.0568338,
      0.15378, -0.0376785, -0.351364, 0.149104, -4.94425e-15, 1.01204, 1.54475,
      -0.112429, 8.50653e-15, -0.511536, -0.452008, 0.756019, -3.3547e-14,
      -0.00106227, 0.0237672, 0.00345981, -0.00139675, -0.00349474, -0.597906,
      -0.425733, -0.0605479, -0.343823, 0.162103, -0.45692, 0.21318, -0.600309,
      0.310843, -0.36406, 0.574148, 0.0554949, -0.00928842, -0.0414346,
      0.0619999, -0.0250297, -0.0626259, 0.00227746, 0.00162164, 0.00023063,
      -0.0301047, 0.273177, -0.770004, 0.359253, 0.0095153, -0.8783, 1.02867,
      -1.62228, -1.24406, -0.00106227, 0.0237672, 0.00238182, 0.00205737,
      0.00402848, 0.262742, 0.151145, -0.671213, -0.343823, 0.317484, 0.12884,
      -0.40386, -0.600309, 0.201313, -0.327527, -0.641099, 0.0554949,
      -0.00928842, -0.0414346, 0.0426822, 0.0368682, 0.0721904, -0.0010008,
      -0.000575719, 0.00255669, -0.0301047, 0.535026, 0.217123, -0.680588,
      0.0095153, -0.568818, 0.925441, 1.81145, -1.24406, -0.00106227, 0.0237672,
      -0.00318563, 0.0034409, -0.00203628, 0.514364, -0.353326, 0.391148,
      -0.343823, -0.496623, -0.0536813, -0.176018, -0.600309, -0.744328,
      -0.01898, 0.0665156, 0.0554949, -0.00928842, -0.0414346, -0.0570866,
      0.0616609, -0.0364902, -0.00195924, 0.00134584, -0.0014899, -0.0301047,
      -0.836913, -0.0904642, -0.296627, 0.0095153, 2.10313, 0.0536287,
      -0.187943, -1.24406, -0.00106227, 0.0237672, -0.002656, -0.00410152,
      0.00150255, -0.1792, 0.627913, 0.340613, -0.343823, 0.0170366, 0.38176,
      0.366698, -0.600309, 0.232172, 0.710567, 0.000435528, 0.0554949,
      -0.00928842, -0.0414346, -0.0475955, -0.0734994, 0.0269257, 0.000682583,
      -0.00239176, -0.00129742, -0.0301047, 0.0287103, 0.643346, 0.617962,
      0.0095153, -0.656011, -2.00774, -0.0012306, -1.24406;

  Eigen::MatrixXd proj =
      mo_ref.transpose() * overlap.Matrix() * orb.MOs().eigenvectors();
  Eigen::VectorXd norms = proj.colwise().norm();
  bool mo_comp = norms.isApproxToConstant(1, 1e-5);

  if (!mo_comp) {
    std::cout << "result mos" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref mos" << std::endl;
    std::cout << mo_ref << std::endl;
    std::cout << "norms" << std::endl;
    std::cout << norms << std::endl;
    std::cout << "proj" << std::endl;
    std::cout << proj << std::endl;
  }

  for (unsigned i = occlevels; i < 17; i++) {
    orb.MOs().eigenvalues()(i) += opt.levelshift;
  }

  d.Levelshift(H, orb.MOs().eigenvectors());
  votca::tools::EigenSystem result = d.SolveFockmatrix(H);

  bool check_level =
      result.eigenvalues().isApprox(orb.MOs().eigenvalues(), 0.00001);
  BOOST_CHECK_EQUAL(check_level, 1);
}

BOOST_AUTO_TEST_SUITE_END()
