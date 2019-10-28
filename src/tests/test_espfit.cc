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

#define BOOST_TEST_MODULE espfit_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/espfit.h>
#include <votca/xtp/logger.h>

#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(espfit_test)

BOOST_AUTO_TEST_CASE(esp_charges) {
  std::ofstream xyzfile("molecule.xyz");
  xyzfile << "3" << std::endl;
  xyzfile << "Water molecule" << std::endl;
  xyzfile << "O          0.00000        0.00000        0.11779" << std::endl;
  xyzfile << "H          0.00000        0.75545       -0.47116" << std::endl;
  xyzfile << "H          0.00000       -0.75545       -0.47116" << std::endl;

  xyzfile.close();

  std::ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << std::endl;
  basisfile << "  <!--Basis set created by xtp_basisset from 3-21G.nwchem at "
               "Thu Sep 15 15:40:33 2016-->"
            << std::endl;
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
  basisfile << "  <element name=\"O\">" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << std::endl;
  basisfile << "      <constant decay=\"3.220370e+02\">" << std::endl;
  basisfile << "        <contractions factor=\"5.923940e-02\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"4.843080e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"3.515000e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.042060e+01\">" << std::endl;
  basisfile << "        <contractions factor=\"7.076580e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"7.402940e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"-4.044530e-01\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"2.445860e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "      <constant decay=\"1.576200e+00\">" << std::endl;
  basisfile << "        <contractions factor=\"1.221560e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"8.539550e-01\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << std::endl;
  basisfile << "      <constant decay=\"3.736840e-01\">" << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << std::endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << std::endl;
  basisfile << "      </constant>" << std::endl;
  basisfile << "    </shell>" << std::endl;
  basisfile << "  </element>" << std::endl;
  basisfile << "</basis>" << std::endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  orbitals.setDFTbasisName("3-21G.xml");
  orbitals.setBasisSetSize(13);
  orbitals.setNumberOfOccupiedLevels(5);

  Eigen::MatrixXd MOs = Eigen::MatrixXd::Zero(13, 13);
  MOs << 0.982343, 0.225297, 3.39697e-12, 0.104226, -4.02405e-12, 0.111297,
      1.47612e-11, -1.73573e-12, -0.0611901, -1.82245e-12, 0.0515886,
      -6.55285e-12, -0.0865545, 0.101935, -0.213424, 1.90464e-12, -0.0897338,
      3.33978e-12, -0.0568715, -1.26376e-11, 1.12461e-11, 0.0976597,
      1.73492e-11, -0.131099, 1.55304e-11, 1.64108, -0.00360159, 0.116669,
      -5.86782e-13, -0.438404, 2.7592e-11, 0.249979, 3.24481e-11, -1.05194e-11,
      -0.319739, 5.38091e-11, -0.988233, 7.8109e-11, -0.147551, -2.79025e-13,
      3.45211e-12, -0.409602, 3.0704e-12, -5.80283e-12, 2.85189e-11, -0.358104,
      -0.198503, -3.7813e-12, -1.77834e-11, -8.49144e-11, -1.0425, 4.00797e-12,
      1.29673e-14, -1.00926e-12, 7.23423e-12, -2.7662e-11, -0.5242, 1.34007e-11,
      7.15316e-13, 1.60688e-11, -1.50996e-11, -1.02779, -5.12996e-11,
      1.43023e-11, 7.55262e-12, -0.0420683, -0.681791, -1.96952e-11, -0.470854,
      -8.52383e-12, -1.00918, -1.12167e-10, -1.35751e-13, 0.0433163,
      -3.21991e-13, -0.210812, 2.43573e-11, -2.00561, 0.00731426, 0.115649,
      3.63395e-12, -0.484273, 4.21012e-11, 0.461816, 3.83154e-11, -1.34171e-12,
      -0.217276, -5.47088e-11, 1.18317, -9.64031e-11, 0.501366, 2.3837e-12,
      3.59972e-12, -0.348495, 1.18667e-11, -1.82988e-11, 1.0157e-10, -0.745697,
      -0.354101, 1.19913e-11, 2.08924e-11, 1.11025e-10, 1.44212, -3.91477e-12,
      -1.66267e-12, -4.34816e-12, 4.68601e-12, -3.15464e-11, -0.629373,
      2.74539e-11, 1.77388e-11, -2.11844e-11, 1.89289e-11, 0.966973,
      5.04727e-11, -1.68919e-11, -7.27192e-12, 0.00274585, -0.126026, -0.238054,
      0.130359, -1.2134e-11, 0.0855516, 0.0897993, 0.954381, -0.965681,
      1.55138e-11, 0.314424, -0.169014, 0.291888, 0.00750313, -0.0234965,
      -0.19353, 0.118088, 4.80145e-11, 0.828874, 1.11045, -0.776492, 0.534722,
      -4.61076e-11, 0.0817261, -0.524685, 0.3665, 0.00274585, -0.126026,
      0.238054, 0.130359, -9.25038e-13, 0.0855516, -0.0897993, -0.954381,
      -0.965681, -1.21157e-11, 0.314424, 0.169014, 0.291888, 0.00750313,
      -0.0234965, 0.19353, 0.118088, 3.94754e-12, 0.828874, -1.11045, 0.776492,
      0.534722, 1.64278e-11, 0.0817261, 0.524685, 0.3665;

  orbitals.MOs().eigenvectors() = MOs;
  orbitals.MOs().eigenvalues() = Eigen::VectorXd::Ones(13);
  QMState gs = QMState("n");
  Logger log;
  Espfit esp = Espfit(log);
  esp.setUseSVD(1e-8);
  StaticSegment result = esp.Fit2Density(orbitals, gs, "xcoarse");
  Eigen::VectorXd pcharges = Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  Index index = 0;
  for (const auto& site : result) {
    pcharges(index) = site.getCharge();
    index++;
  }
  Eigen::VectorXd p_ref = Eigen::VectorXd::Zero(3);
  p_ref << -0.827774, 0.413985, 0.41379;

  bool check_esp_num = p_ref.isApprox(pcharges, 0.01);
  if (!check_esp_num) {
    std::cout << "ref" << std::endl;
    std::cout << p_ref << std::endl;
    std::cout << "calc" << std::endl;
    std::cout << pcharges << std::endl;
  }
  BOOST_CHECK_EQUAL(check_esp_num, 1);

  std::vector<std::pair<Index, Index> > pairconstraint;
  std::pair<Index, Index> p1;
  p1.first = 1;
  p1.second = 2;
  pairconstraint.push_back(p1);
  Espfit esp2 = Espfit(log);
  esp2.setUseSVD(1e-8);
  esp2.setPairConstraint(pairconstraint);
  StaticSegment result2 = esp2.Fit2Density(orbitals, gs, "xcoarse");
  Eigen::VectorXd pcharges_equal = Eigen::VectorXd::Zero(result2.size());
  index = 0;
  for (const auto& site : result2) {
    pcharges_equal(index) = site.getCharge();
    index++;
  }

  bool check_p1 = (std::abs(pcharges_equal(1) - pcharges_equal(2)) < 1e-6);
  BOOST_CHECK_EQUAL(check_p1, 1);

  std::vector<QMFragment<double> > regionconstraint;

  std::string indeces = "1...2";
  QMFragment<double> reg = QMFragment<double>(0, indeces);
  reg.value() = 1.0;
  regionconstraint.push_back(reg);
  Espfit esp3 = Espfit(log);
  esp3.setRegionConstraint(regionconstraint);
  esp3.setUseSVD(1e-8);
  StaticSegment result3 = esp3.Fit2Density(orbitals, gs, "xcoarse");
  Eigen::VectorXd pcharges_reg =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  index = 0;

  for (const auto& site : result3) {
    pcharges_reg(index) = site.getCharge();
    index++;
  }

  bool check_reg = (std::abs(pcharges_reg.segment(1, 2).sum() - 1.0) < 1e-6);
  if (!check_reg) {
    std::cout << "All charges " << pcharges_reg << std::endl;
    std::cout << "Sum of charges 1,2,3 should equal 1:"
              << pcharges_reg.segment(1, 2).sum() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_reg, 1);
}

BOOST_AUTO_TEST_SUITE_END()
