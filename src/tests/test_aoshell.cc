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

#define BOOST_TEST_MODULE aoshell_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(aoshell_test)

BOOST_AUTO_TEST_CASE(EvalAOspace) {

  QMMolecule mol = QMMolecule("", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aoshell/Al.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aoshell/largeshell.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  const AOShell& shell = aobasis.getShell(0);

  Eigen::Vector3d gridpos = Eigen::Vector3d::Ones();
  Eigen::VectorXd aoval = Eigen::VectorXd::Zero(aobasis.AOBasisSize());
  Eigen::MatrixX3d aograd = Eigen::MatrixX3d::Zero(aobasis.AOBasisSize(), 3);
  Eigen::Block<Eigen::MatrixX3d> grad_block =
      aograd.block(0, 0, aobasis.AOBasisSize(), 3);
  Eigen::VectorBlock<Eigen::VectorXd> ao_block =
      aoval.segment(0, aobasis.AOBasisSize());

  shell.EvalAOspace(ao_block, grad_block, gridpos);

  Eigen::VectorXd aoval_2 = Eigen::VectorXd::Zero(aobasis.AOBasisSize());
  Eigen::VectorBlock<Eigen::VectorXd> ao_block_2 =
      aoval_2.segment(0, aobasis.AOBasisSize());
  shell.EvalAOspace(ao_block_2, gridpos);

  Eigen::VectorXd aoval_ref = Eigen::VectorXd::Zero(aobasis.AOBasisSize());
  aoval_ref << 0.0680316, 0.0895832, 0.0895832, 0.0895832, 0, 0.126153,
      0.126153, 0.126153, 0, -0.102376, 0.0626925, 0.0626925, 0.198251, 0,
      0.0809357, -0.0809357, -0.122215, -0.0552111, -0.0552111, 0.156161, 0,
      0.146075, -0.146075, 0, -0.103291;

  Eigen::MatrixX3d aograd_ref =
      Eigen::MatrixX3d::Zero(aobasis.AOBasisSize(), 3);

  aograd_ref << -0.057521967096, -0.057521967096, -0.057521967096,
      -0.091846116024, -0.091846116024, -0.0022629076774, -0.091846116024,
      -0.0022629076774, -0.091846116024, -0.0022629076774, -0.091846116024,
      -0.091846116024, -0.072834589679, -0.072834589679, 0.14566917936,
      -0.16812154688, -0.041968337007, -0.041968337007, -0.041968337007,
      -0.16812154688, -0.041968337007, -0.041968337007, -0.041968337007,
      -0.16812154688, 0.12615320987, -0.12615320987, 0, 0.027261217239,
      0.027261217239, 0.18082590591, -0.17342532206, -0.11073280044,
      0.14003728606, -0.11073280044, -0.17342532206, 0.14003728606,
      -0.15191670048, -0.15191670048, -0.15191670048, 0.19825116059,
      -0.19825116059, 0, 0.099851661525, -0.14295543066, -0.14295543066,
      0.14295543066, -0.099851661525, 0.14295543066, 0.16861461687,
      0.16861461687, -0.0059782842903, -0.042137230049, -0.097348353145,
      0.28912950853, -0.097348353145, -0.042137230049, 0.28912950853,
      -0.27121951095, -0.27121951095, 0.11918208443, 0.15616063815,
      -0.15616063815, 0, 0.11148463165, -0.32674007231, -0.18066517099,
      0.32674007231, -0.11148463165, 0.18066517099, 0.20658110657,
      -0.20658110657, 0, 0.024459014248, 0.024459014248, 0.23104012081;

  bool ao_check = aoval_ref.isApprox(aoval, 1e-5);
  if (!ao_check) {
    std::cout << "ref" << std::endl;
    std::cout << aoval_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << aoval << std::endl;
  }
  BOOST_CHECK_EQUAL(ao_check, 1);
  bool aograd_check = aograd_ref.isApprox(aograd, 1e-5);
  if (!aograd_check) {
    std::cout << "ref" << std::endl;
    std::cout << aograd_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << aograd << std::endl;
  }

  BOOST_CHECK_EQUAL(aograd_check, 1);

  bool ao1vsao2_check = aoval_2.isApprox(aoval, 1e-5);
  if (!ao1vsao2_check) {
    std::cout << "ref" << std::endl;
    std::cout << aoval << std::endl;
    std::cout << "result" << std::endl;
    std::cout << aoval_2 << std::endl;
  }

  BOOST_CHECK_EQUAL(aograd_check, 1);
}

BOOST_AUTO_TEST_SUITE_END()
