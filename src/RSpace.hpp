#ifndef RSPACE_CLASS_INCLUDED
#define RSPACE_CLASS_INCLUDED

#include "kokkos_linalg.hpp"
#include <vector>

/**
 *
 * Class to compute the r-space part of an Ewald sum for
 * quadrupole, dipole and/or monopole systems
 *
 */

template <class T>
class RSpace {
 public:
  /**
   *
   * Constructor for the class (base-version)
   *
   */
  RSpace() {}
  /**
   *
   * Constructor for the class
   *
   * @param alpha     Ewald sum splitting parameter
   * @param r_max     floating point cut-off radius for the r-space
   *                  part
   * @param _l         system length
   */
  RSpace(T alpha, T r_max, T l) { init_params(alpha, r_max, l); }
  /*
   *
   * Initialization of the internal parameters
   *
   * @param alpha     Ewald sum splitting parameter
   * @param r_max     floating point cut-off radius for the r-space
   *                  part
   * @param _l         system length
   */
  void init_params(T, T, T);
  /**
   *
   * Computation of the rspace part of the Ewald sum, fills the
   * private variables of the class which can be accessed by the
   * corresponding getters after the computation
   *
   * @param xyz   particle positions (xyz(0), xyz(1), ...)
   * @param q     particle charges(q,q,q,...)
   * @param d     particle dipole moments(d(0,0),d(0,1),
   *              d(1,0),d(1,1),d(2,0),d(2,1),...)
   * @param q     particle quadrupole moments(
   * q(0,0),...,q(0,8),q(1,0),...,q(1,8),...)
   */
  void compute(const std::vector<T>&, const std::vector<T>&,
               const std::vector<T>&, const std::vector<T>&);

 private:
  T alpha;  ///< splitting parameter of the Ewald sum
  T r_max;  ///< floating-point version of the r-space cutoff
  T l;      // boxsize

  T pot_energy;             ///< potential energy computed for the r-space part
  Kokkos::View<T * [3]> f;  ///< forces computed for the k-space part
  Kokkos::View<T * [3]> tqe;  ///< torque computed for the k-space part

  bool is_monopole;    ///< compute the monopole contributions
  bool is_dipole;      ///< compute the dipole contributions
  bool is_quadrupole;  ///< compute the quadrupole contribution

  template <class Vector>
  Kokkos::View<T[3]> cross(const Vector& x, const Vector& y) const;
};

template <class T>
void RSpace<T>::init_params(T _alpha, T _r_max, T _l) {
  alpha = _alpha;
  r_max = _r_max;
  l = _l;
}

template <class T>
void RSpace<T>::compute(const std::vector<T>& _xyz, const std::vector<T>& _q,
                        const std::vector<T>& _d, const std::vector<T>& _Q) {

  // get number of particles
  size_t N = _xyz.size() / 3;

  // store C++ vectors to Kokkos Views
  Kokkos::View<T * [3]> xyz("positions", N);
  Kokkos::View<T*> q("charges", N);
  Kokkos::View<T * [3]> d("dipole moments", N);
  Kokkos::View<T * [9]> Q("quadrupole moments", N);

  Kokkos::parallel_for(
      N, KOKKOS_LAMBDA(const int n) {
        q(n) = _q.at(n);
        for (int i = 0; i < 3; ++i) {
          xyz(n, i) = _xyz.at(3 * n + i);
          d(n, i) = _d.at(3 * n + i);
        }
        for (int i = 0; i < 9; ++i) {
          Q(n, i) = _Q.at(9 * n + i);
        }
      });

  // create Kokkos views of sufficient size to store the factors
  f = Kokkos::View<T * [3]>("r-space force contribution", N);
  tqe = Kokkos::View<T * [3]>("r-space torque contribution", N);

  Kokkos::parallel_for(
      N, KOKKOS_LAMBDA(const int n) {
        for (int i = 0; i < 3; ++i) {
          f(n, i) = 0.0;
          tqe(n, i) = 0.0;
        }
      });

  Kokkos::parallel_for(
      N, KOKKOS_LAMBDA(const int i) {
        for (int j = i + 1; j < N; ++j) {
          Kokkos::View<T[3]> dR("distance vector");
          for (int d = 0; d < 3; d++) {
            dR(d) = xyz(i, d) - xyz(j, d);
            dR(d) -= std::round(dR[d] / l) * l;
          }
          T dist2 = dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2];
          T inv_dist = 1.0 / std::sqrt(inv_dist);

          // eqn 69-71 in Smith Point Multipoles in Ewald Summation(Revisited)
          T Bn[6];
          Bn[0] = inv_dist * std::erfc(alpha / inv_dist);
          T expfactor = std::exp(-alpha * alpha * dist2);
          for (int b = 1; b < 6; b++) {
            Bn[b] = inv_dist * inv_dist * (2 * b - 1) * Bn[b - 1] +
                    std::pow(2 * alpha * alpha, b) / (alpha * std::sqrt(M_PI)) *
                        expfactor;
          }

          auto Di = Kokkos::subview(d, i, Kokkos::ALL());
          auto Dj = Kokkos::subview(d, j, Kokkos::ALL());
          auto Qi = Kokkos::subview(Q, i, Kokkos::ALL());
          auto Qj = Kokkos::subview(Q, j, Kokkos::ALL());

          Kokkos::View<T[3]> dixdj = cross(Di, Dj);
          Kokkos::View<T[3]> dixR = cross(Di, dR);
          Kokkos::View<T[3]> djxR = cross(Dj, dR);

          Kokkos::View<T[3]> QIR = gemv(Qi, dR);
          Kokkos::View<T[3]> QJR = gemv(Qj, dR);
          Kokkos::View<T[3]> QIQJR = gemv(Qi, QJR);
          Kokkos::View<T[3]> QJQIR = gemv(Qj, QIR);
          Kokkos::View<T[3]> QIXQJ = cross_matrix_product(Qi, Qj);

          Kokkos::View<T[3]> RxQIR = cross(dR, QIR);
          Kokkos::View<T[3]> RxQJR = cross(dR, QJR);

          Kokkos::View<T[3]> RxQIJR = cross(dR, QIQJR);
          Kokkos::View<T[3]> RxQJIR = cross(dR, QJQIR);

          Kokkos::View<T[3]> QJRXQIR = cross(QJR, QIR);

          Kokkos::View<T[3]> QIDJ = gemv(Qi, Dj);
          Kokkos::View<T[3]> QJDI = gemv(Qj, Di);

          Kokkos::View<T[3]> DIXQJR = cross(Di, QJR);
          Kokkos::View<T[3]> DJXQIR = cross(Dj, QIR);

          Kokkos::View<T[3]> RXQIDJ = cross(dR, QIDJ);

          Kokkos::View<T[3]> RXQJDI = cross(dR, QJDI);

          T QII = trace(Qi);
          T QJJ = trace(Qj);

          T DD = dot(Di, Dj);
          T DIR = dot(Di, dR);
          T DJR = dot(Dj, dR);
          T QIRR = dot(QIR, dR);
          T QJRR = dot(QJR, dR);
          T QIRDJ = dot(QIR, dR);
        }
      });
}

#endif
