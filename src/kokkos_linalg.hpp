
#ifndef KOKKOS_LINALG_CLASS_INCLUDED
#define KOKKOS_LINALG_CLASS_INCLUDED

#include "Kokkos_Core.hpp"
#include <cmath>

namespace kokkos_linalg_3d {

template <class Vector1, class Vector2>
Kokkos::View<typename Vector1::value_type[3]> cross(const Vector1& x,
                                                    const Vector2& y) {
  Kokkos::View<typename Vector1::value_type[3]> result("crossproduct");
  assert(x.extent(0) == y.extent(0) && "Dimensions of vectors do not match!");
  assert(x.extent(0) == 3 && "Only three dimensions allowed");
  result[0] = x(1) * y(2) - x(2) * y(1);
  result[1] = x(2) * y(0) - x(0) * y(2);
  result[2] = x(0) * y(1) - x(1) * y(0);
  return result;
}

template <class Vector1, class Vector2>
typename Vector1::value_type dot(const Vector1& x, const Vector2& y) {
  assert(x.extent(0) == y.extent(0) && "Dimensions of arrays do not match!");
  typename Vector1::value_type result = 0.0;
  for (auto i = 0; i < x.extent(0); i++) {
    result += x(i) * y(i);
  }
  return result;
}

template <class Vector>
typename Vector::value_type norm(const Vector& x) {
  return std::sqrt(dot(x, x));
}

template <class Matrix, class Vector>
Kokkos::View<typename Vector::value_type[3]> gemv(const Matrix& A,
                                                  const Vector& x) {
  Kokkos::View<typename Vector::value_type[3]> result("matmul");
  result(0) = A(0) * x(0) + A(1) * x(1) + A(2) * x(2);
  result(1) = A(3) * x(0) + A(4) * x(1) + A(5) * x(2);
  result(2) = A(6) * x(0) + A(7) * x(1) + A(8) * x(2);
  return result;
}

// eqn 21 in Smith Point Multipoles in Ewald Summation(Revisited)
// generalized cross product
template <class Matrix1, class Matrix2>
Kokkos::View<typename Matrix1::value_type[3]> cross_matrix_product(
    const Matrix1& A, const Matrix2& B) {
  assert(A.extent(0) == B.extent(0) && "Dimensions of matrices do not match!");
  assert(A.extent(0) == 9 && "Only 3x3 matrices allowed");
  Kokkos::View<typename Matrix1::value_type[3]> result("cross_matrix_product");

  result(0) = A(1) * B(2) + A(4) * B(5) + A(7) * B(8) - A(2) * B(1) -
              A(5) * B(4) - A(8) * B(7);
  result(1) = A(2) * B(0) + A(5) * B(3) + A(8) * B(6) - A(0) * B(2) -
              A(3) * B(5) - A(6) * B(8);
  result(2) = A(0) * B(1) + A(3) * B(4) + A(6) * B(7) - A(1) * B(0) -
              A(4) * B(3) - A(7) * B(6);
  return result;
}

template <class Matrix>
typename Matrix::value_type trace(const Matrix& A) {
  assert(A.extent(0) == 9 && "Only 3x3 matrices allowed");
  return A(0) + A(4) + A(8);
}
}  // namespace kokkos_linalg_3d

#endif