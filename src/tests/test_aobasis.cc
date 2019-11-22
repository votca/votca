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

#define BOOST_TEST_MODULE aobasis_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/convergenceacc.h>
#include <votca/xtp/orbitals.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(aobasis_test)

BOOST_AUTO_TEST_CASE(FillNormBasis_test) {
  std::ofstream basisfile("notnormalized.xml");
  basisfile << "<basis name=\"def2-TZVP\">" << std::endl;
  basisfile << "  <element name=\"Al\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
  basisfile << "      <constant decay=\"1.570000e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"3.330000e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"D\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell> " << std::endl;
  basisfile << "  </element> " << std::endl;
  basisfile << "</basis> " << std::endl;
  basisfile.close();

  std::ofstream xyzfile("Al.xyz");
  xyzfile << " 1" << std::endl;
  xyzfile << " Al" << std::endl;
  xyzfile << " Al            .000000     .000000     .000000" << std::endl;
  xyzfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("Al.xyz");
  BasisSet basis;
  basis.Load("notnormalized.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  const AOShell& shell = aobasis.getShell(0);
  std::vector<double> ref_results = {0.1831079647, 0.9155398233};
  Index i = 0;
  bool check_norm = true;
  for (const AOGaussianPrimitive& gaussian : shell) {
    if (std::abs(ref_results[i] - gaussian.getContraction()[2]) > 1e-7) {
      check_norm = false;
      break;
    }
    i++;
  }

  i = 0;
  if (!check_norm) {
    for (const AOGaussianPrimitive& gaussian : shell) {
      std::cout << "Ref:" << ref_results[i]
                << " result:" << gaussian.getContraction()[2] << std::endl;
      i++;
    }
  }
  BOOST_CHECK_EQUAL(check_norm, 1);
}

BOOST_AUTO_TEST_CASE(ReorderMos_test) {

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
  basisfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
  basisfile << "      <constant decay=\"1.570000e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell> " << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  BasisSet basis;
  basis.Load("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  orbitals.MOs().eigenvectors() =
      Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  orbitals.MOs().eigenvectors() << 0.996559, -0.223082, -5.74064e-15,
      1.70342e-15, 2.96852e-16, -3.57518e-16, -4.15482e-16, -1.07176e-17,
      -0.0301914, 1.01049e-16, 1.89625e-16, -1.73748e-16, 2.38718e-16,
      1.62375e-17, 0.102052, 9.62277e-16, -6.0957e-16, -1.06426e-15, 0.111923,
      -1.5113e-15, -2.10033e-15, -4.58257e-15, 0.0445146, 0.88316, 2.08298e-14,
      -4.52658e-15, 1.11277e-15, 3.26653e-15, 1.00375e-14, 2.77629e-15,
      0.866126, 5.6158e-15, -5.62939e-15, -7.52681e-16, 2.37035e-15,
      -2.50255e-17, 0.938329, 1.04884e-14, -6.90948e-15, -1.21867e-14, 1.28767,
      -1.95522e-14, -2.43333e-14, -5.31215e-14, 1.17833e-16, -2.62332e-16,
      0.231688, 0.60693, 0.827952, -0.216204, -0.0399226, 0.286201,
      -3.69205e-15, 0.106469, -0.42034, 0.21456, -3.06369e-16, 3.43687e-16,
      -1.20522e-15, 0.00988083, -0.10169, -0.149274, 2.71906e-15, -0.0154314,
      0.0340686, 0.0918581, 1.53082e-16, 2.78929e-14, -1.00703, 0.299218,
      0.0624584, 0.0937755, 0.328402, 0.11665, -5.17067e-15, 0.440979,
      0.0102622, -0.198718, -3.59947e-16, 4.92066e-16, 1.59132e-17, 0.173706,
      -0.0356156, 0.0357605, -5.96263e-16, 0.0939227, -0.0213146, 0.0236835,
      4.20694e-17, 6.11932e-16, 0.199382, 0.80601, -0.646639, -0.273334,
      0.144247, -0.186362, -9.99121e-16, -0.168103, -0.239303, -0.385398,
      1.34626e-15, 4.28838e-16, 1.93238e-17, -0.049494, -0.145299, 0.0957059,
      -5.74112e-16, -0.0278764, -0.0906737, 0.0289463, -0.0472716, 0.235837,
      1.11456e-14, -6.49632e-15, -4.93927e-15, 2.93059e-16, -9.29852e-15,
      -4.62943e-15, -0.599523, -1.21396e-14, 5.36322e-15, 2.68772e-15,
      -5.69707e-15, 2.52476e-16, -1.8333, -3.11059e-14, 1.50032e-14,
      2.96739e-14, -4.24316, 6.46996e-14, 8.35461e-14, 1.79304e-13,
      -2.05526e-16, -1.86284e-15, -0.0360124, -0.0943383, -0.128693, 0.0723654,
      0.0133625, -0.0957943, 1.59924e-14, -0.348785, 1.377, -0.702883,
      2.73335e-16, 4.12143e-17, 9.70315e-15, -0.0748836, 0.770673, 1.1313,
      -1.56785e-14, 0.0762958, -0.168441, -0.454163, -4.61587e-16, -7.91469e-15,
      0.156528, -0.046509, -0.00970823, -0.0313875, -0.109919, -0.0390437,
      1.03899e-14, -1.44462, -0.0336184, 0.650987, 2.38336e-15, -3.75347e-15,
      1.19166e-15, -1.31646, 0.269919, -0.271017, 1.65778e-16, -0.464371,
      0.105383, -0.117096, -1.12124e-16, -1.56626e-15, -0.030991, -0.125282,
      0.10051, 0.0914873, -0.0482808, 0.062377, 4.17567e-15, 0.550694, 0.783942,
      1.26254, -4.93981e-15, -2.05029e-15, -1.57029e-15, 0.375099, 1.10117,
      -0.725324, 1.34493e-15, 0.137826, 0.448307, -0.143116, 5.73096e-18,
      1.02238e-16, -2.48464e-16, -1.10921e-16, -1.52385e-17, -2.91919e-16,
      -4.23503e-16, 5.26081e-17, -1.95459e-16, -1.93453e-15, -5.85128e-16,
      -3.01357e-15, -0.998193, 0.060094, 1.59178e-15, -4.71758e-16,
      -3.52931e-16, -7.53139e-16, 4.14142e-16, -1.88079e-15, -1.622e-16,
      3.54387e-17, 1.23411e-17, 1.96963e-16, 0.00308147, 0.0124569, -0.00999387,
      0.242938, -0.128206, 0.165638, -1.35965e-15, -0.139405, -0.198451,
      -0.319605, 1.5009e-15, -1.26406e-16, 6.70813e-16, 0.20618, 0.605281,
      -0.398688, -6.86196e-15, -0.133052, -0.432778, 0.138158, -2.40587e-17,
      2.86328e-16, -0.0155638, 0.00462444, 0.0009653, -0.0833475, -0.291883,
      -0.103678, 3.92148e-15, 0.365698, 0.00851032, -0.164794, -6.99698e-16,
      -1.95342e-16, -2.30698e-16, -0.723617, 0.148366, -0.14897, 1.34703e-14,
      0.448285, -0.101733, 0.113039, -4.00384e-18, -2.56438e-17, 0.00358075,
      0.00938016, 0.0127961, 0.192161, 0.0354831, -0.254375, -2.70117e-15,
      0.0882931, -0.348582, 0.177931, -8.4513e-16, 1.17851e-15, 5.49385e-15,
      -0.0411612, 0.423615, 0.621838, 2.45073e-14, -0.0736529, 0.162606,
      0.438431, 9.53106e-18, 7.87376e-17, 4.50653e-16, -1.1446e-16, 3.98938e-16,
      3.51619e-16, -9.23785e-18, -8.69969e-16, -2.47306e-16, 1.14344e-15,
      -1.91458e-15, -8.51659e-16, -0.060094, -0.998193, 4.38673e-16,
      9.70656e-16, 3.30285e-16, 1.02556e-15, -1.00514e-16, 3.88181e-16,
      1.69959e-18, -1.19517e-15, 0.00106227, 0.0237672, -0.001459, 0.00433715,
      0.000617507, 0.461796, -0.504929, -0.252611, -0.343823, 0.0850202,
      -0.145542, -0.0828265, 1.94519e-15, 1.78636e-15, 0.600309, 0.0703459,
      -0.148256, -0.00934183, -0.0554949, -0.271103, 0.417352, -0.773903,
      0.00928842, -0.0414346, -0.0288797, 0.0858503, 0.012223, 0.047698,
      -0.0521531, -0.0260917, -0.0301047, 0.287304, -0.49182, -0.279891,
      3.08213e-16, 2.95652e-16, -0.0095153, 0.627532, -1.32254, -0.0833354,
      1.24406, 0.484301, -0.745562, 1.38251, 0.00106227, 0.0237672, 0.0026328,
      -0.00126226, 0.00357714, 0.0427602, 0.598097, -0.415299, -0.343823,
      -0.0372958, -0.0428745, 0.179002, -3.86028e-16, 1.11668e-15, 0.600309,
      -0.0599788, 0.0415619, -0.147278, -0.0554949, 0.43641, -0.782307,
      -0.210113, 0.00928842, -0.0414346, 0.052114, -0.0249854, 0.0708065,
      0.00441662, 0.0617762, -0.0428954, -0.0301047, -0.126031, -0.144883,
      0.604892, 2.31849e-15, -4.33324e-15, -0.0095153, -0.535051, 0.37076,
      -1.31382, 1.24406, -0.779607, 1.39752, 0.375349, 0.00106227, 0.0237672,
      0.00246913, -0.000253666, -0.00389357, 0.176084, 0.168299, 0.687525,
      -0.343823, -0.160372, 0.0382744, -0.0899273, 9.65416e-17, -6.69665e-16,
      0.600309, -0.122276, -0.00419351, 0.109758, -0.0554949, 0.569725,
      0.553977, 0.46382, 0.00928842, -0.0414346, 0.0488744, -0.00502111,
      -0.0770701, 0.0181874, 0.0173832, 0.0710131, -0.0301047, -0.541935,
      0.129339, -0.303886, 4.5522e-15, -1.4397e-16, -0.0095153, -1.09078,
      -0.0374089, 0.979111, 1.24406, -1.01776, -0.989631, -0.828573, 0.00106227,
      0.0237672, -0.00364293, -0.00282122, -0.000301075, -0.680641, -0.261466,
      -0.0196149, -0.343823, 0.112648, 0.150142, -0.00624844, 2.57639e-15,
      -1.04514e-15, 0.600309, 0.111909, 0.110888, 0.0468622, -0.0554949,
      -0.735032, -0.189022, 0.520197, 0.00928842, -0.0414346, -0.0721087,
      -0.0558438, -0.00595952, -0.0703021, -0.0270063, -0.00202599, -0.0301047,
      0.380663, 0.507365, -0.021115, -4.23727e-15, 3.19439e-15, -0.0095153,
      0.998299, 0.989191, 0.418042, 1.24406, 1.31307, 0.337671, -0.929285;

  Eigen::MatrixXd ref = orbitals.MOs().eigenvectors();
  aobasis.ReorderMOs(orbitals.MOs().eigenvectors(), "xtp", "gaussian");
  aobasis.ReorderMOs(orbitals.MOs().eigenvectors(), "gaussian", "xtp");

  bool check_reorder_gaus = ref.isApprox(orbitals.MOs().eigenvectors(), 1e-7);
  if (!check_reorder_gaus) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "reordered" << std::endl;
    std::cout << orbitals.MOs().eigenvectors() << std::endl;
  }

  BOOST_CHECK_EQUAL(check_reorder_gaus, 1);

  aobasis.ReorderMOs(orbitals.MOs().eigenvectors(), "xtp", "nwchem");
  aobasis.ReorderMOs(orbitals.MOs().eigenvectors(), "nwchem", "xtp");

  bool check_reorder_nwchem = ref.isApprox(orbitals.MOs().eigenvectors(), 1e-7);

  BOOST_CHECK_EQUAL(check_reorder_nwchem, 1);
}

BOOST_AUTO_TEST_SUITE_END()
