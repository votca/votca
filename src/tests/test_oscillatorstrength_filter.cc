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

#define BOOST_TEST_MODULE oscillatorstrength_filter_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/filterfactory.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(oscillatorstrength_filter_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {
  FilterFactory::RegisterAll();

  libint2::initialize();
  std::unique_ptr<StateFilter_base> osc_f =
      std::unique_ptr<StateFilter_base>(Filter().Create("oscillatorstrength"));

  std::ofstream opt("oscillatorstrength.xml");
  opt << "<oscillatorstrength>0.0045</oscillatorstrength>" << std::endl;
  opt.close();
  votca::tools::Property prop;
  prop.LoadFromXML("oscillatorstrength.xml");
  osc_f->Initialize(prop.get("oscillatorstrength"));

  Orbitals A;
  A.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                    "/oscillatorstrength_filter/3-21G.xml");
  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                           "/oscillatorstrength_filter/molecule.xyz");
  A.setBasisSetSize(17);
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/oscillatorstrength_filter/MOs_A.mm");

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/oscillatorstrength_filter/spsi_ref.mm");

  A.BSESinglets().eigenvectors() = spsi_ref;

  // reference energy
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;

  A.BSESinglets().eigenvalues() = se_ref;
  A.CalcCoupledTransition_Dipoles();

  BOOST_CHECK_EQUAL(osc_f->NeedsInitialState(), false);

  std::vector<votca::Index> ref = {0};
  std::vector<votca::Index> results =
      osc_f->CalcIndeces(A, QMStateType::Singlet);

  BOOST_CHECK_EQUAL(results.size(), ref.size());

  for (votca::Index i = 0; i < votca::Index(ref.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i], results[i]);
  }

  BOOST_REQUIRE_THROW(osc_f->CalcIndeces(A, QMStateType::Gstate),
                      std::runtime_error);
  BOOST_REQUIRE_THROW(osc_f->CalcIndeces(A, QMStateType::Triplet),
                      std::runtime_error);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
