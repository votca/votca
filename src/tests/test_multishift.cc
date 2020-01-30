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

#define BOOST_TEST_MODULE multishift_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <votca/xtp/multishift.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <votca/xtp/orbitals.h>
using namespace votca::xtp;
using namespace std; 

BOOST_AUTO_TEST_SUITE(multishift_test)

BOOST_AUTO_TEST_CASE(multishift_test){
    
    Multishift multi;
    
    multi.setMatrixSize(30);
    
    Eigen::MatrixXcd I;

    I.setIdentity(30, 30);
    
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(30, 30);
    
    std::cout << "A=" << std::endl << A << std::endl;

    Eigen::VectorXcd b = Eigen::VectorXcd::Random(30);

    std::cout << "b=" << std::endl << b << std::endl;
    
    std::cout << "Starting BiCG" << std::endl << std::endl;

    Multishift::MultiShiftResult result = multi.ComplexBiCG(A, b);
    
    std::cout << "BiCG complete" << std::endl << std::endl;
    
    BOOST_CHECK_EQUAL(result.converged, true);
    
    std::cout << "x=" << std::endl << result._x << std::endl;
    std::cout << "diff x= " <<std::endl<<A.colPivHouseholderQr().solve(b)-result._x;
    
    Eigen::VectorXcd res = (A)*result._x - b;

    std::cout << "res=" << std::endl << res << std::endl;

    BOOST_CHECK_EQUAL(res.norm() < 1e-13, true);
    
    std::complex<double> w(10, 0);

    Eigen::VectorXcd x_w = multi.DoMultishift(A, b, w, result);

    std::cout << "x_w=" << std::endl << x_w << std::endl;

    Eigen::VectorXcd res_w = (A + w * I) * x_w - b;

    std::cout << "res_w=" << std::endl << res_w << std::endl;

    BOOST_CHECK_EQUAL(res_w.norm() < 1e-13, true);
    

}BOOST_AUTO_TEST_SUITE_END()