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

#define BOOST_TEST_MODULE bse_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bseoperator_btda.h"
#include "votca/xtp/matrixfreeoperator.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(bse_test_operatorbtda)

class HermitianBlockOperator : public MatrixFreeOperator {
 public:
  HermitianBlockOperator() = default;

  void attach_matrix(const Eigen::MatrixXd &mat);
  Eigen::MatrixXd matmul(const Eigen::MatrixXd &input) const {
    return _mat * input;
  }

  Eigen::VectorXd diagonal() const { return _mat.diagonal(); }

 private:
  Eigen::MatrixXd _mat;
};

void HermitianBlockOperator::attach_matrix(const Eigen::MatrixXd &mat) {
  _mat = mat;
}

BOOST_AUTO_TEST_CASE(bse_test_operatorbtda) {

  votca::Index size = 7;

  Eigen::MatrixXd r1 = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd symm1 = r1 + r1.transpose();

  Eigen::MatrixXd r2 = Eigen::MatrixXd::Random(size, size);
  Eigen::MatrixXd symm2 = r2 + r2.transpose();

  HermitianBlockOperator Rop;
  Rop.set_size(size);
  Rop.attach_matrix(symm1);

  HermitianBlockOperator Cop;
  Cop.set_size(size);
  Cop.attach_matrix(symm2);

  // create Hamiltonian operator
  HamiltonianOperator<HermitianBlockOperator, HermitianBlockOperator> Hop(Rop,
                                                                          Cop);

  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(2 * size, 2 * size);

  Eigen::MatrixXd H=Hop*identity;

  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(2 * size, 2 * size);
  ref.topLeftCorner(size,size)=symm1;
  ref.topRightCorner(size,size)=symm2;
  ref.bottomLeftCorner(size,size)=-symm2;
  ref.bottomRightCorner(size,size)=-symm1;

  bool check = ref.isApprox(H, 1E-9);
  BOOST_CHECK_EQUAL(check, true);
  if (!check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << H << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
