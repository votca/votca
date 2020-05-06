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

#define BOOST_TEST_MODULE rspace_test
#include "../src/RSpace.hpp"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(rspace_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Kokkos::initialize();
  {
    RSpace<double> a;
    a.init_params(0.3, 0.7, 1.2);

    // testing particle positions
    std::vector<double> xyz{0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.50,
                            0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.50, 0.50,
                            0.00, 0.50, 0.00, 0.50, 0.50, 0.50, 0.50, 0.50};

    // testing charges
    std::vector<double> q{-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};

    // testing dipole momenta
    std::vector<double> d{-0.5, 0.3,  -0.1, -0.4, 0.4,  -0.2, -0.3, 0.5,
                          -0.3, -0.2, 0.4,  -0.4, -0.1, 0.3,  -0.5, -0.0,
                          0.2,  -0.4, 0.1,  0.1,  -0.3, 0.2,  0.0,  -0.2};

    std::vector<double> Q{
        -0.5, 0.3,  -0.1, -0.1, 0.5,  -0.3, 0.1,  0.1,  -0.3, -0.4, 0.4,  -0.2,
        0.0,  0.4,  -0.4, 0.2,  0.0,  -0.2, -0.3, 0.5,  -0.3, 0.1,  0.3,  -0.5,
        0.3,  -0.1, -0.1, -0.2, 0.4,  -0.4, 0.2,  0.2,  -0.4, 0.4,  -0.2, -0.0,
        -0.1, 0.3,  -0.5, 0.3,  0.1,  -0.3, 0.5,  -0.3, 0.1,  0.0,  0.2,  -0.4,
        0.4,  0.0,  -0.2, 0.4,  -0.4, 0.2,  0.1,  0.1,  -0.3, 0.5,  -0.1, -0.1,
        0.3,  -0.5, 0.3,  0.2,  0.0,  -0.2, 0.4,  -0.2, 0.0,  0.2,  -0.4, 0.4};

    a.compute(xyz, q, d, Q);
    double result = 0.5;
    BOOST_CHECK_CLOSE(result, 0.5, 1e-5);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
