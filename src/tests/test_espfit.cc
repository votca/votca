/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/espfit.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include <libint2/initialize.h>
using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(espfit_test)

BOOST_AUTO_TEST_CASE(esp_charges) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/espfit/molecule.xyz");
  orbitals.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                           "/espfit/3-21G.xml");
  orbitals.setBasisSetSize(13);
  orbitals.setNumberOfOccupiedLevels(5);

  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/espfit/MOs.mm");
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

  std::string indeces = "1:2";
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

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
