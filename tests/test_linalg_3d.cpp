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

#define BOOST_TEST_MODULE linalg_test
#include "../src/kokkos_linalg.hpp"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(linalg_test)

BOOST_AUTO_TEST_CASE(cross) {
  Kokkos::initialize();
  {
    Kokkos::View<double[3]> a("a");
    a(0) = -1;
    a(1) = 7;
    a(2) = 4;
    Kokkos::View<double[3]> b("b");
    b(0) = -5;
    b(1) = 8;
    b(2) = 4;

    auto result = kokkos_linalg_3d::cross(a, b);
    BOOST_CHECK_CLOSE(result[0], -4, 1e-9);
    BOOST_CHECK_CLOSE(result[1], -16, 1e-9);
    BOOST_CHECK_CLOSE(result[2], 27, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(dot) {
  Kokkos::initialize();
  {
    Kokkos::View<double[3]> a("a");
    a(0) = -1;
    a(1) = 7;
    a(2) = 4;
    Kokkos::View<double[3]> b("b");
    b(0) = -5;
    b(1) = 8;
    b(2) = 4;

    auto result = kokkos_linalg_3d::dot(a, b);
    BOOST_CHECK_CLOSE(result, 77, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(gemv) {
  Kokkos::initialize();
  {

    Kokkos::View<double[3]> a("a");
    a(0) = 1;
    a(1) = 2;
    a(2) = 3;
    Kokkos::View<double[9]> A("A");
    A(0) = 5;
    A(1) = 1;
    A(2) = 3;
    A(3) = 1;
    A(4) = 1;
    A(5) = 1;
    A(6) = 1;
    A(7) = 2;
    A(8) = 1;

    auto result = kokkos_linalg_3d::gemv(A, a);
    BOOST_CHECK_CLOSE(result[0], 16, 1e-9);
    BOOST_CHECK_CLOSE(result[1], 6, 1e-9);
    BOOST_CHECK_CLOSE(result[2], 8, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(norm) {
  Kokkos::initialize();
  {
    Kokkos::View<double[3]> a("a");
    a(0) = -1;
    a(1) = 7;
    a(2) = 4;
    auto result = kokkos_linalg_3d::norm(a);
    BOOST_CHECK_CLOSE(result, std::sqrt(66), 1e-5);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(trace) {
  Kokkos::initialize();
  {
    Kokkos::View<double[9]> A("A");
    A(0) = 5;
    A(1) = 1;
    A(2) = 3;
    A(3) = 1;
    A(4) = 1;
    A(5) = 1;
    A(6) = 1;
    A(7) = 2;
    A(8) = 1;

    auto result = kokkos_linalg_3d::trace(A);
    BOOST_CHECK_CLOSE(result, 7, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(cross_matrix_product) {
  Kokkos::initialize();
  {
    Kokkos::View<double[9]> A("A");
    A(0) = 5;
    A(1) = 1;
    A(2) = 3;
    A(3) = 1;
    A(4) = 1;
    A(5) = 1;
    A(6) = 1;
    A(7) = 2;
    A(8) = 1;

    Kokkos::View<double[9]> B("B");
    B(0) = -1;
    B(1) = 1;
    B(2) = 3;
    B(3) = 1;
    B(4) = 1;
    B(5) = 1;
    B(6) = 1;
    B(7) = 2;
    B(8) = 1;
    auto result = kokkos_linalg_3d::cross_matrix_product(A, B);
    BOOST_CHECK_CLOSE(result[0], 0, 1e-9);
    BOOST_CHECK_CLOSE(result[1], -6, 1e-9);
    BOOST_CHECK_CLOSE(result[2], 6, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(scale_3d) {
  Kokkos::initialize();
  {
    Kokkos::View<double[3]> v("v");
    v(0) = -1.2;
    v(1) = 2.4;
    v(2) = 1000.0;
    double s = 0.25;

    auto result = kokkos_linalg_3d::scale_3d(s, v);

    BOOST_CHECK_CLOSE(result[0], -0.3, 1e-9);
    BOOST_CHECK_CLOSE(result[1], 0.6, 1e-9);
    BOOST_CHECK_CLOSE(result[2], 250.0, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_CASE(dualbase_3d) {
  Kokkos::initialize();
  {
    Kokkos::View<double[3]> a("a");
    a(0) = 1.0;
    a(1) = 1.0;
    a(2) = 1.0;

    Kokkos::View<double[3]> b("b");
    b(0) = -1.0;
    b(1) = 1.0;
    b(2) = 1.0;

    Kokkos::View<double[3]> c("c");
    c(0) = 1.0;
    c(1) = -1.0;
    c(2) = 1.0;

    auto result = kokkos_linalg_3d::dualbase_3d(a, b, c);

    BOOST_CHECK_CLOSE(result[0][0], 0.5, 1e-9);
    BOOST_CHECK_CLOSE(result[0][1], 0.5, 1e-9);
    BOOST_CHECK_CLOSE(result[0][2], 0.0, 1e-9);

    BOOST_CHECK_CLOSE(result[1][0], -0.5, 1e-9);
    BOOST_CHECK_CLOSE(result[1][1], 0.0, 1e-9);
    BOOST_CHECK_CLOSE(result[1][2], 0.5, 1e-9);

    BOOST_CHECK_CLOSE(result[2][0], 0.0, 1e-9);
    BOOST_CHECK_CLOSE(result[2][1], -0.5, 1e-9);
    BOOST_CHECK_CLOSE(result[2][2], 0.5, 1e-9);
  }
  Kokkos::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
