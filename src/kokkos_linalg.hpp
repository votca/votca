
#ifndef KOKKOS_LINALG_CLASS_INCLUDED
#define KOKKOS_LINALG_CLASS_INCLUDED

#include "Kokkos_Core.hpp"
#include <cmath>

namespace kokkos_linalg_3d {

template <class T, class Vector1, class Vector2>
Kokkos::View<T[3]> cross(const Vector1& x, const Vector2& y) {
  Kokkos::View<T[3]> result("crossproduct");
  static_assert(x.extend(0) == y.extend(0),
                "Dimensions of vectors do not match!");
  static_assert(x.extend(0) == 3, "Only three dimensions allowed");
  result[0] = x(1) * y(2) - x(2) * y(1);
  result[1] = x(2) * y(0) - x(0) * y(2);
  result[2] = x(0) * y(1) - x(1) * y(0);
  return result;
}

template <class T, class Vector1, class Vector2>
T dot(const Vector1& x, const Vector2& y) {
  static_assert(x.extend(0) == y.extend(0),
                "Dimensions of arrays do not match!");
  T result = 0.0;
  for (auto i = 0; i < x.extend(0); i++) {
    result += x(i) * y(i);
  }
  return result;
}

template <class T, class Vector>
T norm(const Vector& x) {
  return std::sqrt(dot(x, x));
}

template <class T, class Matrix, class Vector>
Kokkos::View<T[3]> gemv(const Matrix& A, const Vector& x) {
  Kokkos::View<T[3]> result("matmul");
  result(0) = A(0) * x(0) + A(1) * x(1) + A(2) * x(2);
  result(1) = A(3) * x(0) + A(4) * x(1) + A(5) * x(2);
  result(2) = A(6) * x(0) + A(7) * x(1) + A(8) * x(2);
  return result;
}

// eqn 21 in Smith Point Multipoles in Ewald Summation(Revisited)
// generalized cross product
template <class T, class Matrix1, class Matrix2>
Kokkos::View<T[3]> cross_matrix_product(const Matrix1& A, const Matrix2& B) {
  static_assert(A.extend(0) == B.extend(0),
                "Dimensions of matrices do not match!");
  static_assert(A.extend(0) == 9, "Only 3x3 matrices allowed");
  Kokkos::View<T[3]> result("cross_matrix_product");

  result(0) = A(1) * B(2) + A(4) * B(5) + A(7) * B(8) - A(2) * B(1) -
              A(5) * B(4) - A(8) * B(7);
  result(1) = A(2) * B(0) + A(5) * B(3) + A(8) * B(6) - A(0) * B(2) -
              A(3) * B(5) - A(6) * B(8);
  result(2) = A(0) * B(1) + A(3) * B(4) + A(6) * B(7) - A(1) * B(0) -
              A(4) * B(3) - A(7) * B(6);
  return result;
}

template <class T, class Matrix>
T trace(const Matrix& A) {
  static_assert(A.extend(0) == 9, "Only 3x3 matrices allowed");
  return A(0) + A(4) + A(8);
}
}  // namespace kokkos_linalg_3d

#endif