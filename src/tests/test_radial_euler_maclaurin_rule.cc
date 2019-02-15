/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE aomatrix_test
#include "votca/xtp/orbitals.h"
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/xtp/radial_euler_maclaurin_rule.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(radial_euler_maclaurin_test)

BOOST_AUTO_TEST_CASE(setup_test) {

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
  AOBasis aobasis;
  aobasis.AOBasisFill(basis, orbitals.QMAtoms());

  EulerMaclaurinGrid radialgrid;
  auto grid = radialgrid.CalculateAtomicRadialGrids(aobasis, orbitals.QMAtoms(),
                                                    "medium");

  auto Cgrid = grid.at("C");
  auto Hgrid = grid.at("H");

  Eigen::VectorXd C_radius_ref = Eigen::VectorXd::Zero(49);
  C_radius_ref << 2.80419e-05, 0.000224341, 0.00075721, 0.00179513, 0.00350698,
      0.00606227, 0.00963155, 0.0143869, 0.0205023, 0.0281546, 0.0375238,
      0.0487943, 0.0621556, 0.0778038, 0.0959422, 0.116783, 0.14055, 0.167478,
      0.197817, 0.231835, 0.269818, 0.312078, 0.358952, 0.410809, 0.468057,
      0.531147, 0.600581, 0.676925, 0.760816, 0.852983, 0.954259, 1.06561,
      1.18816, 1.32325, 1.47244, 1.63766, 1.82121, 2.02599, 2.25563, 2.51479,
      2.80962, 3.14844, 3.54292, 4.01013, 4.57655, 5.28651, 6.22317, 7.57314,
      9.93197;

  Eigen::VectorXd C_weight_ref = Eigen::VectorXd::Zero(49);
  C_weight_ref << 6.61523e-14, 1.69369e-11, 4.34206e-10, 4.33973e-09,
      2.58921e-08, 1.11494e-07, 3.8345e-07, 1.11898e-06, 2.88109e-06,
      6.72222e-06, 1.44868e-05, 2.92464e-05, 5.59082e-05, 0.000102053,
      0.000179067, 0.000303666, 0.000499919, 0.00080193, 0.00125739, 0.00193227,
      0.00291702, 0.00433476, 0.00635227, 0.00919452, 0.0131643, 0.0186686,
      0.0262548, 0.0366601, 0.05088, 0.0702642, 0.096653, 0.132572, 0.181518,
      0.24838, 0.340075, 0.466532, 0.64224, 0.888765, 1.23897, 1.74433, 2.48828,
      3.61166, 5.36479, 8.22247, 13.1661, 22.479, 42.4806, 96.4423, 338.808;

  Eigen::VectorXd H_radius_ref = Eigen::VectorXd::Zero(49);
  H_radius_ref << 3.36503e-05, 0.00026921, 0.000908652, 0.00215416, 0.00420837,
      0.00727472, 0.0115579, 0.0172643, 0.0246028, 0.0337855, 0.0450285,
      0.0585531, 0.0745868, 0.0933646, 0.115131, 0.14014, 0.16866, 0.200973,
      0.23738, 0.278202, 0.323782, 0.374494, 0.430742, 0.492971, 0.561669,
      0.637376, 0.720697, 0.81231, 0.912979, 1.02358, 1.14511, 1.27873, 1.4258,
      1.5879, 1.76693, 1.96519, 2.18545, 2.43119, 2.70675, 3.01774, 3.37154,
      3.77813, 4.25151, 4.81216, 5.49186, 6.34382, 7.4678, 9.08777, 11.9184;

  Eigen::VectorXd H_weight_ref = Eigen::VectorXd::Zero(49);
  H_weight_ref << 1.14311e-13, 2.9267e-11, 7.50308e-10, 7.49906e-09,
      4.47416e-08, 1.92661e-07, 6.62602e-07, 1.9336e-06, 4.97853e-06,
      1.1616e-05, 2.50333e-05, 5.05378e-05, 9.66094e-05, 0.000176347,
      0.000309428, 0.000524735, 0.00086386, 0.00138573, 0.00217277, 0.00333897,
      0.0050406, 0.00749047, 0.0109767, 0.0158881, 0.0227479, 0.0322593,
      0.0453684, 0.0633487, 0.0879206, 0.121417, 0.167016, 0.229085, 0.313664,
      0.429201, 0.58765, 0.806167, 1.10979, 1.53579, 2.14094, 3.01421, 4.29975,
      6.24096, 9.27036, 14.2084, 22.7511, 38.8437, 73.4066, 166.652, 585.46;

  BOOST_CHECK_EQUAL(Cgrid.radius.size(), C_radius_ref.size());
  BOOST_CHECK_EQUAL(Cgrid.weight.size(), C_weight_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.radius.size(), H_radius_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.weight.size(), H_weight_ref.size());

  bool Cradius = C_radius_ref.isApprox(Cgrid.radius, 0.001);
  bool Cweight = C_weight_ref.isApprox(Cgrid.weight, 0.001);
  BOOST_CHECK_EQUAL(Cradius, true);
  BOOST_CHECK_EQUAL(Cweight, true);

  bool Hradius = H_radius_ref.isApprox(Hgrid.radius, 0.001);
  bool Hweight = H_weight_ref.isApprox(Hgrid.weight, 0.001);
  BOOST_CHECK_EQUAL(Hradius, true);
  BOOST_CHECK_EQUAL(Hweight, true);
}

BOOST_AUTO_TEST_SUITE_END()
