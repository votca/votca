/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#define BOOST_TEST_MODULE matrix_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(matrix_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  vec v1(1, 2, 3);
  vec v2(4, 5, 6);
  vec v3(7, 8, 9);
  matrix mat2(v1, v2, v3);
  matrix mat3(mat2);
  matrix mat4 = mat2;
  std::cout << "mat2" << mat2 << std::endl;
  std::cout << "mat4" << mat4 << std::endl;
}

BOOST_AUTO_TEST_CASE(eigen_test) {
  vec v1(1, 2, 3);
  vec v2(4, 5, 6);
  vec v3(7, 8, 9);
  matrix mat(v1, v2, v3);
  matrix mat2(mat.ToEigenMatrix());

  BOOST_CHECK(mat.isClose(mat2, 0.01));
  std::cout << "mat" << mat << std::endl;
  std::cout << "mat eigen" << mat.ToEigenMatrix() << std::endl;
}

BOOST_AUTO_TEST_CASE(overloadoperator_test) {

  vec v1(1, 1, 1);
  matrix mat1(v1, v1, v1);
  matrix mat2(v1, v1, v1);

  { BOOST_CHECK(mat1.isClose(mat2, 0.01)); }

  vec v2(0, 0, 0);
  matrix mat3(v2, v2, v2);
  matrix mat4(v2, v2, v2);

  {
    mat3 += mat1;
    BOOST_CHECK(mat3.isClose(mat1, 0.01));
    mat3 -= mat1;
    BOOST_CHECK(mat3.isClose(mat4, 0.01));
  }

  vec v3(-1, -1, -1);
  matrix mat5(v3, v3, v3);
  {
    mat1 /= (-1);
    BOOST_CHECK(mat1.isClose(mat5, 0.01));
    mat1 *= (-1);
    BOOST_CHECK(mat1.isClose(mat2, 0.01));
  }
}

BOOST_AUTO_TEST_SUITE_END()
