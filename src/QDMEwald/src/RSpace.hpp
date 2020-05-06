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

  namespace kl3 = kokkos_linalg_3d;

  // get number of particles
  size_t N = _xyz.size() / 3;

  // store C++ vectors to Kokkos Views
  Kokkos::View<T * [3]> xyz("positions", N);
  Kokkos::View<T*> q("charges", N);
  Kokkos::View<T * [3]> d("dipole moments", N);
  Kokkos::View<T * [9]> Q("quadrupole moments", N);

  Kokkos::View<T*> pe("potential energy", N);

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    pe(n) = 0.0;
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

    std::array<T, 3> force_i = {0, 0, 0};
    std::array<T, 3> torque_i = {0, 0, 0};
    for (int j = i + 1; j < N; ++j) {
      std::array<T, 3> dR;
      for (int d = 0; d < 3; d++) {
        dR[d] = xyz(i, d) - xyz(j, d);
        dR[d] -= std::round(dR[d] / l) * l;
      }
      T dist2 = dR[0] * dR[0] + dR[1] * dR[1] + dR[2] * dR[2];
      //      if (dist2 > r_max * r_max) continue;
      T inv_dist = 1.0 / std::sqrt(dist2);

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

      std::array<T, 3> dixdj = kl3::cross(Di, Dj);
      std::array<T, 3> dixR = kl3::cross(Di, dR);
      std::array<T, 3> djxR = kl3::cross(Dj, dR);

      std::array<T, 3> QIR = kl3::gemv(Qi, dR);
      std::array<T, 3> QJR = kl3::gemv(Qj, dR);
      std::array<T, 3> QIQJR = kl3::gemv(Qi, QJR);
      std::array<T, 3> QJQIR = kl3::gemv(Qj, QIR);
      std::array<T, 3> QIXQJ = kl3::cross_matrix_product(Qi, Qj);

      std::array<T, 3> RxQIR = kl3::cross(dR, QIR);
      std::array<T, 3> RxQJR = kl3::cross(dR, QJR);

      std::array<T, 3> RxQIJR = kl3::cross(dR, QIQJR);
      std::array<T, 3> RxQJIR = kl3::cross(dR, QJQIR);

      std::array<T, 3> QJRXQIR = kl3::cross(QJR, QIR);

      std::array<T, 3> QIDJ = kl3::gemv(Qi, Dj);
      std::array<T, 3> QJDI = kl3::gemv(Qj, Di);

      std::array<T, 3> DIXQJR = kl3::cross(Di, QJR);
      std::array<T, 3> DJXQIR = kl3::cross(Dj, QIR);

      std::array<T, 3> RXQIDJ = kl3::cross(dR, QIDJ);

      std::array<T, 3> RXQJDI = kl3::cross(dR, QJDI);

      T QII = kl3::trace(Qi);
      T QJJ = kl3::trace(Qj);  // sc1

      T DD = kl3::dot(Di, Dj);        // sc2
      T DIR = kl3::dot(Di, dR);       // sc3
      T DJR = kl3::dot(Dj, dR);       // sc4
      T QIRR = kl3::dot(QIR, dR);     // sc5
      T QJRR = kl3::dot(QJR, dR);     // sc6
      T QIRDJ = kl3::dot(QIR, Dj);    // sc7
      T QJRDI = kl3::dot(QJR, Di);    // sc8
      T QJRQIR = kl3::dot(QJR, QIR);  // sc9
      T QMIQMJ = kl3::dot(Qi, Qj);    // sc10

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

      thread_Potential = Bn[0] * GL[0] + Bn[1] * (GL[1] + GL[6]) +
                         Bn[2] * (GL[2] + GL[7] + GL[8]) +
                         Bn[3] * (GL[3] + GL[5]) + Bn[4] * GL[4];

#ifdef DEBUG_OUTPUT_ENABLED
      std::cout << "r-space DEBUG: " << i << " <> " << j << " " << pot_energy
                << " { " << xyz(i, 0) << ", " << xyz(i, 1) << ", " << xyz(i, 2)
                << " } "
                << " { " << xyz(j, 0) << ", " << xyz(j, 1) << ", " << xyz(j, 2)
                << " } "
                << " + " << thread_Potential << " " << dist2 << " "
                << " " << (Bn[0] * GL[0]) << " "
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

      Kokkos::atomic_add(&pe(i), 0.5 * thread_Potential);
      Kokkos::atomic_add(&pe(j), 0.5 * thread_Potential);
      Kokkos::atomic_add(&pot_energy, thread_Potential);

      thread_Virial += Bn[0] * GL[0] + Bn[1] * (2 * GL[1] + 3 * GL[6]) +
                       Bn[2] * (3 * GL[2] + 4 * GL[7] + 5 * GL[8]) +
                       Bn[3] * (4 * GL[3] + 5 * GL[5]) + Bn[4] * 5 * GL[4];

      std::array<T, 7> F_M_coeff;

      F_M_coeff[0] = Bn[1] * GL[0] + Bn[2] * (GL[1] + GL[6]) +
                     Bn[3] * (GL[2] + GL[7] + GL[8]) + Bn[4] * (GL[3] + GL[5]) +
                     Bn[5] * GL[4];
      F_M_coeff[1] = -q(j) * Bn[1] + (DJR + QJJ) * Bn[2] - QJRR * Bn[3];
      F_M_coeff[2] = q(i) * Bn[1] + (DIR - QII) * Bn[2] + QIRR * Bn[3];
      F_M_coeff[3] = 2 * Bn[2];
      F_M_coeff[4] = 2 * (-q(j) * Bn[2] + (QJJ + DJR) * Bn[3] - QJRR * Bn[4]);
      F_M_coeff[5] = 2 * (-q(i) * Bn[2] + (QII - DIR) * Bn[3] - QIRR * Bn[4]);
      F_M_coeff[6] = 4 * Bn[3];

      std::array<T, 3> thread_force = {0, 0, 0};
      kl3::add_to(thread_force, kl3::scale_3d(F_M_coeff[0], dR));
      kl3::add_to(thread_force, kl3::scale_3d(F_M_coeff[1], Di));
      kl3::add_to(thread_force, kl3::scale_3d(F_M_coeff[2], Dj));
      kl3::add_to(thread_force,
                  kl3::scale_3d(F_M_coeff[3], kl3::subtract(QJDI, QIDJ)));
      kl3::add_to(thread_force, kl3::scale_3d(F_M_coeff[4], QIR));
      kl3::add_to(thread_force, kl3::scale_3d(F_M_coeff[5], QJR));
      kl3::add_to(thread_force,
                  kl3::scale_3d(F_M_coeff[6], kl3::add(QIQJR, QJQIR)));

      std::array<T, 3> thread_M_i = {0, 0, 0};

      kl3::subtract_from(thread_M_i, kl3::scale_3d(Bn[1], dixdj));
      kl3::add_to(thread_M_i, kl3::scale_3d(F_M_coeff[1], dixR));
      std::array<T, 3> temp_sum = kl3::add(DIXQJR, DJXQIR);
      kl3::add_to(temp_sum, RXQIDJ);
      kl3::subtract_from(temp_sum, kl3::scale_3d(2, QIXQJ));
      kl3::add_to(thread_M_i, kl3::scale_3d(F_M_coeff[3], temp_sum));
      kl3::subtract_from(thread_M_i, kl3::scale_3d(F_M_coeff[4], RxQIR));
      kl3::subtract_from(
          thread_M_i, kl3::scale_3d(F_M_coeff[6], kl3::add(RxQIJR, QJRXQIR)));

      std::array<T, 3> thread_M_j = {0, 0, 0};

      kl3::add_to(thread_M_j, kl3::scale_3d(Bn[1], dixdj));
      kl3::add_to(thread_M_j, kl3::scale_3d(F_M_coeff[2], djxR));
      temp_sum = kl3::add(DIXQJR, DJXQIR);
      kl3::add_to(temp_sum, RXQJDI);
      kl3::subtract_from(temp_sum, kl3::scale_3d(2, QIXQJ));
      kl3::subtract_from(thread_M_j, kl3::scale_3d(F_M_coeff[3], temp_sum));
      kl3::subtract_from(thread_M_j, kl3::scale_3d(F_M_coeff[5], RxQJR));
      kl3::subtract_from(
          thread_M_j,
          kl3::scale_3d(F_M_coeff[6], kl3::subtract(RxQIJR, QJRXQIR)));

      kl3::subtract_from(force_i, thread_force);
      kl3::add_to(torque_i, thread_M_i);

      for (int k = 0; k < 3; k++) {
        Kokkos::atomic_add(&f(j, k), thread_force[k]);
        Kokkos::atomic_add(&tqe(j, k), thread_M_j[k]);
      }
    }

    for (int k = 0; k < 3; k++) {
      Kokkos::atomic_add(&f(i, k), force_i[k]);
      Kokkos::atomic_add(&tqe(i, k), torque_i[k]);
    }

    Kokkos::atomic_add(&virial, thread_Virial);
  });

  // self energy correction
  T U_self = 0.0;
  Kokkos::parallel_reduce(N,
                          KOKKOS_LAMBDA(const int i, T& self_tmp) {
                            T tmp = -q(i) * q(i) * alpha / std::sqrt(M_PI);
                            pe(i) += tmp;
                            self_tmp += tmp;
                          },
                          U_self);

  pot_energy += U_self;
#ifdef QDMEWALD_DEBUG_ENABLED
  std::cout << "DEBUG: r-space energy particle 0:   " << pe(0) << " "
            << pot_energy << std::endl;
  std::cout << "DEBUG: r-space energy particle N/2: " << pe(N / 2) << " "
            << pot_energy << std::endl;
  std::cout << "DEBUG: r-space energy particle N-1: " << pe(N - 1) << " "
            << pot_energy << std::endl;
#endif
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
