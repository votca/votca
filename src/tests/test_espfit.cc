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

#define BOOST_TEST_MODULE espfit_test
#include <boost/test/unit_test.hpp>
#include <votca/ctp/logger.h>
#include <votca/xtp/espfit.h>

#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(espfit_test)

BOOST_AUTO_TEST_CASE(esp_charges) {

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

  Eigen::MatrixXd dmat = Eigen::MatrixXd::Zero(17, 17);
  dmat << 0.00157507, 0.0337454, 4.48905e-16, -5.93152e-16, 7.87133e-17,
      0.030876, 2.51254e-16, -1.49094e-16, 5.77899e-17, 0.00415998, -0.00445632,
      0.00415998, -0.00445632, 0.00415998, -0.00445632, 0.00415998, -0.00445632,
      0.0337454, 0.722983, 2.66427e-15, -4.44783e-15, 3.45846e-16, 0.661507,
      4.39854e-15, -2.02475e-15, 1.04832e-15, 0.0891262, -0.095475, 0.0891262,
      -0.095475, 0.0891262, -0.095475, 0.0891262, -0.095475, 4.48905e-16,
      2.66427e-15, 1.52199, 2.88658e-15, 2.09034e-15, -7.94212e-15, 0.215492,
      2.8727e-15, -1.40513e-15, 0.141933, -0.0402359, 0.141933, -0.0402359,
      -0.141933, 0.0402359, -0.141933, 0.0402359, -5.93152e-16, -4.44783e-15,
      2.88658e-15, 1.52199, -2.31759e-15, 9.21105e-15, -2.22045e-15, 0.215492,
      1.6263e-15, 0.141933, -0.0402359, -0.141933, 0.0402359, -0.141933,
      0.0402359, 0.141933, -0.0402359, 7.87133e-17, 3.45846e-16, 2.09034e-15,
      -2.31759e-15, 1.52199, 2.98902e-15, -2.04958e-15, 4.79738e-15, 0.215492,
      0.141933, -0.0402359, -0.141933, 0.0402359, 0.141933, -0.0402359,
      -0.141933, 0.0402359, 0.030876, 0.661507, -7.94212e-15, 9.21105e-15,
      2.98902e-15, 0.605259, 2.55488e-15, 2.7779e-17, 1.33759e-15, 0.0815477,
      -0.0873567, 0.0815477, -0.0873567, 0.0815477, -0.0873567, 0.0815477,
      -0.0873567, 2.51254e-16, 4.39854e-15, 0.215492, -2.22045e-15,
      -2.04958e-15, 2.55488e-15, 0.0305108, 3.29597e-17, -5.29036e-16,
      0.0200958, -0.00569686, 0.0200958, -0.00569686, -0.0200958, 0.00569686,
      -0.0200958, 0.00569686, -1.49094e-16, -2.02475e-15, 2.8727e-15, 0.215492,
      4.79738e-15, 2.7779e-17, 3.29597e-17, 0.0305108, 9.55941e-16, 0.0200958,
      -0.00569686, -0.0200958, 0.00569686, -0.0200958, 0.00569686, 0.0200958,
      -0.00569686, 5.77899e-17, 1.04832e-15, -1.40513e-15, 1.6263e-15, 0.215492,
      1.33759e-15, -5.29036e-16, 9.55941e-16, 0.0305108, 0.0200958, -0.00569686,
      -0.0200958, 0.00569686, 0.0200958, -0.00569686, -0.0200958, 0.00569686,
      0.00415998, 0.0891262, 0.141933, 0.141933, 0.141933, 0.0815477, 0.0200958,
      0.0200958, 0.0200958, 0.0506951, -0.0230264, -0.00224894, -0.00801753,
      -0.00224894, -0.00801753, -0.00224894, -0.00801753, -0.00445632,
      -0.095475, -0.0402359, -0.0402359, -0.0402359, -0.0873567, -0.00569686,
      -0.00569686, -0.00569686, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.00801753, 0.0115445, 0.00415998, 0.0891262,
      0.141933, -0.141933, -0.141933, 0.0815477, 0.0200958, -0.0200958,
      -0.0200958, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, -0.00445632, -0.095475, -0.0402359,
      0.0402359, 0.0402359, -0.0873567, -0.00569686, 0.00569686, 0.00569686,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, 0.00415998, 0.0891262, -0.141933, -0.141933,
      0.141933, 0.0815477, -0.0200958, -0.0200958, 0.0200958, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00224894,
      -0.00801753, -0.00445632, -0.095475, 0.0402359, 0.0402359, -0.0402359,
      -0.0873567, 0.00569686, 0.00569686, -0.00569686, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992, -0.00801753, 0.0115445,
      0.00415998, 0.0891262, -0.141933, 0.141933, -0.141933, 0.0815477,
      -0.0200958, 0.0200958, -0.0200958, -0.00224894, -0.00801753, -0.00224894,
      -0.00801753, -0.00224894, -0.00801753, 0.0506951, -0.0230264, -0.00445632,
      -0.095475, 0.0402359, -0.0402359, 0.0402359, -0.0873567, 0.00569686,
      -0.00569686, 0.00569686, -0.00801753, 0.0115445, -0.00801753, 0.0115445,
      -0.00801753, 0.0115445, -0.0230264, 0.0157992;

  votca::ctp::Logger log;

  Espfit esp = Espfit(&log);
  esp.setUseSVD(1e-8);
  esp.Fit2Density(orbitals.QMAtoms(), dmat, aobasis, "medium");
  Eigen::VectorXd pcharges = Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  int index = 0;
  for (QMAtom* atom : orbitals.QMAtoms()) {
    pcharges(index) = atom->getPartialcharge();
    index++;
    atom->setPartialcharge(0.0);
  }

  Eigen::VectorXd p_ref = Eigen::VectorXd::Zero(5);
  p_ref << -1.52561, 0.881377, 0.881411, 0.881411, 0.881411;

  bool check_esp_num = p_ref.isApprox(pcharges, 0.01);
  if (!check_esp_num) {
    cout << "ref" << endl;
    cout << p_ref << endl;
    cout << "calc" << endl;
    cout << pcharges << endl;
  }
  BOOST_CHECK_EQUAL(check_esp_num, 1);


  std::vector<std::pair<int, int> > pairconstraint;
  std::pair<int, int> p1;
  p1.first = 1;
  p1.second = 2;
  pairconstraint.push_back(p1);
  std::pair<int, int> p2;
  p2.first = 3;
  p2.second = 4;
  pairconstraint.push_back(p2);
  Espfit esp2 = Espfit(&log);
  esp2.setUseSVD(1e-8);
  esp2.setPairConstraint(pairconstraint);
  esp2.Fit2Density(orbitals.QMAtoms(), dmat, aobasis, "medium");
  Eigen::VectorXd pcharges_equal =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  index = 0;
  for (const QMAtom* atom : orbitals.QMAtoms()) {
    pcharges_equal(index) = atom->getPartialcharge();
    index++;
  }

  bool check_p1 = (std::abs(pcharges_equal(1) - pcharges_equal(2)) < 1e-6);
  bool check_p2 = (std::abs(pcharges_equal(3) - pcharges_equal(4)) < 1e-6);
  BOOST_CHECK_EQUAL(check_p1 && check_p2, 1);

  std::vector<Espfit::region> regionconstraint;
  Espfit::region reg;
  reg.atomindices = {1, 2, 3};
  reg.charge = 1.0;
  regionconstraint.push_back(reg);
  Espfit esp3 = Espfit(&log);
  esp3.setRegionConstraint(regionconstraint);
  esp3.setUseSVD(1e-8);
  esp3.Fit2Density(orbitals.QMAtoms(), dmat, aobasis, "medium");
  Eigen::VectorXd pcharges_reg =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  index = 0;
  for (const QMAtom* atom : orbitals.QMAtoms()) {
    pcharges_reg(index) = atom->getPartialcharge();
    index++;
  }
  bool check_reg = (std::abs(pcharges_reg.segment(1, 3).sum() - 1.0) < 1e-6);
  if (!check_reg) {
    std::cout << "All charges " << pcharges_reg << std::endl;
    std::cout << "Sum of charges 1,2,3 should equal 1:"
              << pcharges_reg.segment(1, 3).sum() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_reg, 1);
}


BOOST_AUTO_TEST_CASE(analytic_vs_numeric) {

ofstream xyzfile("molecule.xyz");
  xyzfile << " 1" << endl;
  xyzfile << " carbon" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
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

  Eigen::MatrixXd dmat = 0.01*Eigen::MatrixXd::Ones(9, 9);
  Eigen::VectorXd diag=Eigen::VectorXd::Zero(9);
  diag<<1,2,3,4,5,8,1,2,2;
  dmat.diagonal()=diag;
  votca::ctp::Logger log;

  Espfit esp = Espfit(&log);
  esp.setUseSVD(1e-8);


  esp.Fit2Density_analytic(orbitals.QMAtoms(), dmat, aobasis);
  Eigen::VectorXd pcharges_anal =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  int index = 0;
  for (const QMAtom* atom : orbitals.QMAtoms()) {
    pcharges_anal(index) = atom->getPartialcharge();
    index++;
  }

  esp.Fit2Density(orbitals.QMAtoms(), dmat, aobasis, "medium");
  Eigen::VectorXd pcharges = Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
   index = 0;
  for (QMAtom* atom : orbitals.QMAtoms()) {
    pcharges(index) = atom->getPartialcharge();
    index++;
    atom->setPartialcharge(0.0);
  }


  bool check_esp_ana = pcharges.isApprox(pcharges_anal, 0.01);
  if (!check_esp_ana) {
    cout << "numeric" << endl;
    cout << pcharges << endl;
    cout << "analytic" << endl;
    cout << pcharges_anal << endl;
  }
  BOOST_CHECK_EQUAL(check_esp_ana, 1);
}

BOOST_AUTO_TEST_SUITE_END()
