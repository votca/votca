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

  /**
   *
   * Getter for the r-space contribution to forces
   *
   * @return std::vector for type T containing the forces
   */
  std::vector<T> get_forces();

  /**
   *
   * Getter for the r-space contribution to torque
   *
   * @return std::vector for type T containing the torque
   */
  std::vector<T> get_torque();

  /**
   *
   * Getter for the r-space contribution to energy
   *
   * @return type T containing the energy
   */
  T get_energy();

  /**
   *
   * Getter for the r-space contribution to forces
   *
   * @return type T containing the virial
   */
  T get_virial();

 private:
  T alpha;  ///< splitting parameter of the Ewald sum
  T r_max;  ///< floating-point version of the r-space cutoff
  T l;      // boxsize

  T pot_energy = 0.0;       ///< potential energy computed for the r-space part
  T virial = 0.0;           ///< virial computed for the r-space part
  Kokkos::View<T * [3]> f;  ///< forces computed for the r-space part
  Kokkos::View<T * [3]> tqe;  ///< torque computed for the r-space part

  bool is_monopole;    ///< compute the monopole contributions
  bool is_dipole;      ///< compute the dipole contributions
  bool is_quadrupole;  ///< compute the quadrupole contribution
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

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
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

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    for (int i = 0; i < 3; ++i) {
      f(n, i) = 0.0;
      tqe(n, i) = 0.0;
    }
  });

  pot_energy = 0.0;
  virial = 0.0;
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int i) {
    T thread_Virial = 0.0;
    T thread_Potential = 0.0;
    for (int j = i + 1; j < N; ++j) {
      std::array<T, 3> dR;
      for (int d = 0; d < 3; d++) {
        dR[d] = xyz(i, d) - xyz(j, d);
        dR[d] -= std::round(dR[d] / l) * l;
      }
      T dist2 = dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2];
      T dist = std::sqrt(dist2);
      T inv_dist = 1.0 / dist;

      // eqn 69-71 in Smith Point Multipoles in Ewald Summation(Revisited)
      std::array<T, 6> Bn;
      Bn[0] = inv_dist * std::erfc(alpha / inv_dist);
      T expfactor = std::exp(-alpha * alpha * dist2);
      for (int b = 1; b < 6; b++) {
        Bn[b] = inv_dist * inv_dist *
                ((T)(2 * b - 1) * Bn[b - 1] + std::pow(2 * alpha * alpha, b) /
                                                  (alpha * std::sqrt(M_PI)) *
                                                  expfactor);
      }

      auto Di = Kokkos::subview(d, i, Kokkos::ALL());
      auto Dj = Kokkos::subview(d, j, Kokkos::ALL());
      auto Qi = Kokkos::subview(Q, i, Kokkos::ALL());
      auto Qj = Kokkos::subview(Q, j, Kokkos::ALL());

      std::array<T, 3> dixdj = kokkos_linalg_3d::cross(Di, Dj);
      std::array<T, 3> dixR = kokkos_linalg_3d::cross(Di, dR);
      std::array<T, 3> djxR = kokkos_linalg_3d::cross(Dj, dR);

      std::array<T, 3> QIR = kokkos_linalg_3d::gemv(Qi, dR);
      std::array<T, 3> QJR = kokkos_linalg_3d::gemv(Qj, dR);
      std::array<T, 3> QIQJR = kokkos_linalg_3d::gemv(Qi, QJR);
      std::array<T, 3> QJQIR = kokkos_linalg_3d::gemv(Qj, QIR);
      std::array<T, 3> QIXQJ = kokkos_linalg_3d::cross_matrix_product(Qi, Qj);

      std::array<T, 3> RxQIR = kokkos_linalg_3d::cross(dR, QIR);
      std::array<T, 3> RxQJR = kokkos_linalg_3d::cross(dR, QJR);

      std::array<T, 3> RxQIJR = kokkos_linalg_3d::cross(dR, QIQJR);
      std::array<T, 3> RxQJIR = kokkos_linalg_3d::cross(dR, QJQIR);

      std::array<T, 3> QJRXQIR = kokkos_linalg_3d::cross(QJR, QIR);

      std::array<T, 3> QIDJ = kokkos_linalg_3d::gemv(Qi, Dj);
      std::array<T, 3> QJDI = kokkos_linalg_3d::gemv(Qj, Di);

      std::array<T, 3> DIXQJR = kokkos_linalg_3d::cross(Di, QJR);
      std::array<T, 3> DJXQIR = kokkos_linalg_3d::cross(Dj, QIR);

      std::array<T, 3> RXQIDJ = kokkos_linalg_3d::cross(dR, QIDJ);

      std::array<T, 3> RXQJDI = kokkos_linalg_3d::cross(dR, QJDI);

      T QII = kokkos_linalg_3d::trace(Qi);
      T QJJ = kokkos_linalg_3d::trace(Qj);  // sc1

      T DD = kokkos_linalg_3d::dot(Di, Dj);        // sc2
      T DIR = kokkos_linalg_3d::dot(Di, dR);       // sc3
      T DJR = kokkos_linalg_3d::dot(Dj, dR);       // sc4
      T QIRR = kokkos_linalg_3d::dot(QIR, dR);     // sc5
      T QJRR = kokkos_linalg_3d::dot(QJR, dR);     // sc6
      T QIRDJ = kokkos_linalg_3d::dot(QIR, Dj);    // sc7
      T QJRDI = kokkos_linalg_3d::dot(QJR, Di);    // sc8
      T QJRQIR = kokkos_linalg_3d::dot(QJR, QIR);  // sc9
      T QMIQMJ = kokkos_linalg_3d::dot(Qi, Qj);    // sc10

      // eqn 38-42
      std::array<T, 9> GL;
      GL[0] = q(i) * q(j);
      GL[1] = q(j) * DIR - q(i) * DJR;
      GL[2] = q(i) * QJRR + q(j) * QIRR - DIR * DJR;
      GL[3] = DIR * QJRR - DJR * QIRR;
      GL[4] = QIRR * QJRR;
      GL[5] = -(QJJ * QIRR + QJRR * QII + 4 * QJRQIR);
      GL[6] = DD - q(j) * QII - q(i) * QJJ;
      GL[7] = DJR * QII - QJJ * DIR + 2 * (QIRDJ - QJRDI);
      GL[8] = 2 * QMIQMJ + QII * QJJ;

      thread_Potential += Bn[0] * GL[0] + Bn[1] * (GL[1] + GL[6]) +
                          Bn[2] * (GL[2] + GL[7] + GL[8]) +
                          Bn[3] * (GL[3] + GL[5]) + Bn[4] * GL[4];
#ifndef DEBUG_OUTPUT_ENABLED
      std::cout << "r-space DEBUG: " << i << " <> " << j << " " << pot_energy
                << " + " << thread_Potential << " ( " << (Bn[0] * GL[0])
                << " ) "
                << " " << Bn[0] << " "
                << " " << Bn[1] << " "
                << " " << Bn[2] << " "
                << " " << Bn[3] << " "
                << " " << Bn[4] << " "
                << " | "
                << " " << QII << " "
                << " " << QJJ << " "
                << " " << GL[0] << " "
                << " " << GL[1] << " "
                << " " << GL[2] << " "
                << " " << GL[3] << " "
                << " " << GL[4] << " "
                << " " << GL[5] << " "
                << " " << GL[6] << " "
                << " " << GL[7] << " " << std::endl;
#endif
      thread_Virial += Bn[0] * GL[0] + Bn[1] * (2 * GL[1] + 3 * GL[6]) +
                       Bn[2] * (3 * GL[2] + 4 * GL[7] + 5 * GL[8]) +
                       Bn[3] * (4 * GL[3] + 5 * GL[5]) + Bn[4] * 5 * GL[4];
    }
    Kokkos::atomic_add(&pot_energy, thread_Potential);
    Kokkos::atomic_add(&virial, thread_Virial);
  });
}

/**
 *
 * Getter for the r-space contribution to forces
 *
 * @return std::vector for type T containing the forces
 */
template <class T>
std::vector<T> RSpace<T>::get_forces() {
  // get number of particles
  int N = f.extent(0);

  std::vector<T> k_forces(3 * N);

  for (int n = 0; n < N; ++n) {
    for (int d = 0; d < 3; ++d) {
      k_forces.at(3 * n + d) = f(n, d);
    }
  }

  return k_forces;
}

/**
 *
 * Getter for the r-space contribution to torque
 *
 * @return std::vector for type T containing the torque
 */
template <class T>
std::vector<T> RSpace<T>::get_torque() {
  // get number of particles
  int N = tqe.extent(0);

  std::vector<T> k_torque(3 * N);

  for (int n = 0; n < N; ++n) {
    for (int d = 0; d < 3; ++d) {
      k_torque.at(3 * n + d) = tqe(n, d);
    }
  }

  return k_torque;
}

/**
 *
 * Getter for the r-space contribution to energy
 *
 * @return type T containing the energy
 */
template <class T>
T RSpace<T>::get_energy() {
  return pot_energy;
}

/**
 *
 * Getter for the r-space contribution to forces
 *
 * @return type T containing the virial
 */
template <class T>
T RSpace<T>::get_virial() {
  return virial;
}

#endif
