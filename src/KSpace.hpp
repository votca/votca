#ifndef KSPACE_CLASS_INCLUDED
#define KSPACE_CLASS_INCLUDED

#include "Kokkos_Core.hpp"
#include <cmath>
#include <vector>

/**
 *
 * Class to compute the k-space part of an Ewald sum for
 * quadrupole, dipole and/or monopole systems
 *
 */

template <class T>
class KSpace {
 public:
  /**
   *
   * Constructor for the class (base-version)
   *
   */
  KSpace() {}
  /**
   *
   * Constructor for the class
   *
   * @param _alpha     Ewald sum splitting parameter
   * @param _k_max     floating point cut-off radius for the k-space
   *                  part
   * @param _l         system length
   */
  KSpace(const T _alpha, const T _k_max, const T _l) {
    init_params(_alpha, _k_max, _l);
  }
  /*
   *
   * Initialization and computation of internal parameters
   *
   * @param _alpha     Ewald sum splitting parameter
   * @param _k_max     floating point cut-off radius for the k-space
   *                  part
   * @param _l         system length
   */
  void init_params(const T, const T, const T);
  /**
   *
   * Computation of the kspace part of the Ewald sum, fills the
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
  T alpha;        ///< splitting parameter of the Ewald sum
  T gamma;        ///< transformed splitting parameter of the Ewald sum
  int offset;     ///< offset for the navigation in the k-based multi-dim arrays
  int k_sq_int;   ///< square value of the k-space limit
  int k_max_int;  ///< k-space image limit

  T l;                      ///< system length
  T rcl;                    ///< wave length of cell
  T pot_energy;             ///< potential energy computed for the k-space part
  T virial;                 ///< virial computed for the k-space part
  Kokkos::View<T * [3]> f;  ///< forces computed for the k-space part
  Kokkos::View<T * [3]> tqe;        ///< torque computed for the k-space part
  Kokkos::View<T * [15]> vec_comp;  ///< vector components for the computation
                                    ///< of forces
  Kokkos::View<T[6]> vec_sums;      ///< vector summation for the computation of
                                    ///< pot. energy

  bool is_monopole;    ///< compute the monopole contributions
  bool is_dipole;      ///< compute the dipole contributions
  bool is_quadrupole;  ///< compute the quadrupole contributions

  Kokkos::View<T***> cos_fac;  ///< cosine based exponential factors
  Kokkos::View<T***> sin_fac;  ///< sine based exponential factors

  Kokkos::View<T*> ak;  ///< AK coefficients

  /**
   *
   * Compute the exponential factors for each particle
   *
   * @param xyz   std::vector of type T containing the positions of the
   *              particles to be used in the computation of the exponential
   *              factors
   * @param l     system length
   */
  void compute_exponentials(Kokkos::View<T * [3]>);

  /**
   *
   * Computation of the AK coefficients for the matrix multiplications
   *
   * @param l     system length
   */
  void compute_ak();

  /**
   *
   * Computation of the vector components for the given k-cell and storing the
   * results to compute the forces
   *
   * @param x     x coordinate of the k-cell
   * @param y     y coordinate of the k-cell
   * @param z     z coordinate of the k-cell
   * @param q     array of charges
   * @param d     array of dipole moments
   * @param Q     array of quadrupole moments
   */
  void compute_vector_components(const int, const int, const int,
                                 const Kokkos::View<T*>&,
                                 const Kokkos::View<T * [3]>&,
                                 const Kokkos::View<T * [9]>&);

  /**
   *
   * Computation of the vector sums for the given cell and storing the
   * results to compute the energy
   *
   */
  void compute_vector_sums();
};

/*
 *
 * Initialization and computation of internal parameters
 *
 * @param alpha     Ewald sum splitting parameter
 * @param k_max     floating point cut-off radius for the k-space
 *                  part
 */
template <class T>
void KSpace<T>::init_params(const T _alpha, const T _k_max, const T _l) {
  alpha = _alpha;
  gamma = -0.25 / (alpha * alpha);
  k_sq_int = (int)_k_max;
  k_max_int = (int)std::floor(sqrt(_k_max));
  // transformed length
  T rcl = 2.0 * M_PI / _l;
  l = _l;
}

/*
 *
 * Computation of the AK coefficients for the matrix multiplications
 * (currently only implemented for host-space)
 */
template <class T>
void KSpace<T>::compute_ak() {
  T expf = 0.0;

  if (alpha > (T)1e-12) {
    expf = std::exp(gamma * rcl * rcl);
  }

  ak = Kokkos::View<T*>("AK coefficients", k_sq_int);

  Kokkos::parallel_for(k_sq_int, KOKKOS_LAMBDA(const int k) {
    T rksq = (T)k * rcl * rcl;
    T eksq = std::pow(expf, (T)k);
    ak(k) = eksq / rksq;
  });
}

/*
 *
 * Computation of the exponential factors for the k-space part
 * (currently only implemented for host-space)
 *
 * @param xyz       particles positons
 */
template <class T>
void KSpace<T>::compute_exponentials(Kokkos::View<T * [3]> xyz) {
  // get number of particles
  size_t N = xyz.size() / 3.0;
  // create Kokkos views of sufficient size to store the factors
  cos_fac = Kokkos::View<T***>("cosine exponential factors", N,
                               2 * (k_max_int) + 1, 3);
  sin_fac =
      Kokkos::View<T***>("sine exponential factors", N, 2 * (k_max_int) + 1, 3);

  // to deal with negative k-values
  offset = k_max_int;

  // initialize the first factors (k == 0)
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    for (int d = 0; d < 3; ++d) {
      cos_fac(n, 0 + offset, d) = 1.0;
      sin_fac(n, 0 + offset, d) = 0.0;
    }
  });

  // compute the exponential factors (k == 1 / k == -1)
  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    for (int d = 0; d < 3; ++d) {
      cos_fac(n, 1 + offset, d) = std::cos(rcl * xyz(3 * n, d));
      sin_fac(n, 1 + offset, d) = std::sin(rcl * xyz(3 * n, d));
      cos_fac(n, -1 + offset, d) = cos_fac(n, 1 + offset, d);
      sin_fac(n, -1 + offset, d) = -sin_fac(n, 1 + offset, d);
    }
  });

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    for (int k = 2; k <= k_max_int; ++k) {
      for (int d = 0; d < 3; ++d) {
        cos_fac(n, k + offset, d) =
            cos_fac(n, k - 1 + offset, d) * cos_fac(n, 1 + offset, d) -
            sin_fac(n, k - 1 + offset, d) * sin_fac(n, 1 + offset, d);
        sin_fac(n, k + offset, d) =
            sin_fac(n, k - 1 + offset, d) * cos_fac(n, 1 + offset, d) -
            cos_fac(n, k - 1 + offset, d) * sin_fac(n, 1 + offset, d);
        cos_fac(n, -k + offset, d) = cos_fac(n, k + offset, d);
        sin_fac(n, -k + offset, d) = -sin_fac(n, k + offset, d);
      }
    }
  });
}

/**
 *
 * Computation of the vector components for the given k-cell and storing the
 * results to compute the forces
 *
 * @param x     x coordinate of the k-cell
 * @param y     y coordinate of the k-cell
 * @param z     z coordinate of the k-cell
 * @param q     array of charges
 * @param d     array of dipole moments
 * @param Q     array of quadrupole moments
 */
template <class T>
void KSpace<T>::compute_vector_components(const int x, const int y, const int z,
                                          const Kokkos::View<T*>& q,
                                          const Kokkos::View<T * [3]>& d,
                                          const Kokkos::View<T * [9]>& Q) {
  size_t N = q.extent(0);

  // create vector component view
  // [ CKR SKR DK[x] DK[y] DK[z] QK[x] QK[y] QK[z] CKC CKS DKC DKS QKC QKS ]
  // [  0   1   2     3     4      5     6     7    8   9   10  11  12  13 ]
  vec_comp = Kokkos::View<T * [15]>("vector components", N);

  T rx = rcl * (T)x;
  T ry = rcl * (T)y;
  T rz = rcl * (T)z;

  Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n) {
    T cxy = cos_fac(n, x, 1) * cos_fac(n, y, 2) -
            sin_fac(n, x, 1) * sin_fac(n, y, 2);
    T sxy = sin_fac(n, x, 1) * cos_fac(n, y, 2) -
            sin_fac(n, y, 2) * cos_fac(n, x, 1);
    // scalar product dipole moment and k-vector
    T dk = rx * d(n, 0) + ry * d(n, 1) + rz * d(n, 2);
    // tensor product quadrupole moment and k-vector
    T Qk = rx * (rx * q(n, 0) + ry * q(n, 3) + rz * q(n, 6)) +
           ry * (rx * q(n, 1) + ry * q(n, 4) + rz * q(n, 7)) +
           rz * (rx * q(n, 2) + ry * q(n, 5) + rz * q(n, 8));
    // cosine based vector component (monopole)
    vec_comp(n, 0) = cxy * cos_fac(n, z, 3) - sxy * sin_fac(n, z, 3);
    // sine based vector component (monopole)
    vec_comp(n, 1) = sxy * cos_fac(n, z, 3) + cxy * sin_fac(n, z, 3);
    // vector component(s) (dipole)
    vec_comp(n, 2) = rz * d(n, 1) - ry * d(n, 2);
    vec_comp(n, 3) = rx * d(n, 2) - ry * d(n, 0);
    vec_comp(n, 4) = ry * d(n, 0) - ry * d(n, 1);
    // vector component(s) (quadrupole)
    vec_comp(n, 5) = rz * (rx * q(n, 1) + ry * q(n, 4) + rz * q(n, 7)) -
                     ry * (rx * q(n, 2) + ry * q(n, 5) + rz * q(n, 8));
    vec_comp(n, 6) = rx * (rx * q(n, 2) + ry * q(n, 5) + rz * q(n, 8)) -
                     rz * (rx * q(n, 0) + ry * q(n, 3) + rz * q(n, 6));
    vec_comp(n, 7) = ry * (rx * q(n, 0) + ry * q(n, 3) + rz * q(n, 6)) -
                     rx * (rx * q(n, 1) + ry * q(n, 4) + rz * q(n, 7));
    // cosine component (monopole)
    vec_comp(n, 8) = q(n) * vec_comp(n, 0);
    // sine component (monopole)
    vec_comp(n, 9) = q(n) * vec_comp(n, 1);
    // cosine component (dipole)
    vec_comp(n, 10) = dk * vec_comp(n, 0);
    // sine component (dipole)
    vec_comp(n, 11) = dk * vec_comp(n, 1);
    // cosine component (quadrupole)
    vec_comp(n, 12) = Qk * vec_comp(n, 0);
    // sine component (quadrupole)
    vec_comp(n, 13) = Qk * vec_comp(n, 1);
  });
}

/**
 *
 * Computation of the vector components for the given k-cell and storing the
 * results to compute the forces
 *
 */
template <class T>
void KSpace<T>::compute_vector_sums() {}

/**
 *
 * Computation of the kspace part of the Ewald sum, fills the
 * private variables of the class which can be accessed by the
 * corresponding getters after the computation
 *
 * @param _xyz   particle positions (xyz(0), xyz(1), ...)
 * @param _q     particle charges(q,q,q,...)
 * @param _d     particle dipole moments(d(0,0),d(0,1),
 *              d(1,0),d(1,1),d(2,0),d(2,1),...)
 * @param _Q     particle quadrupole moments(
 * Q(0,0),...,Q(0,8),Q(1,0),...,Q(1,8),...)
 */
template <class T>
void KSpace<T>::compute(const std::vector<T>& _xyz, const std::vector<T>& _q,
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

  // compute exponential factors
  compute_exponentials(xyz);

  // compute AK coefficients
  compute_ak();

  // loop over the k-images to compute the necessary
  // parameters (energy, virial, torque, forces)
  int miny = 0;
  int minz = 0;
  double rx, ry, rz;
  for (int ix = 0; ix <= k_max_int; ++ix) {
    rx = rcl * (T)ix;
    for (int iy = miny; iy <= k_max_int; ++iy) {
      ry = rcl * (T)iy;
      for (int iz = minz; iz <= k_max_int; ++iz) {
        rz = rcl * (T)iz;
        // -1 to conform to C++ style array indexing
        // opposed to original Fortran indexing
        int kk = ix * ix + iy * iy + iz * iz - 1;
        if (kk > k_sq_int) continue;
        // compute exp (ikr) and the corresponding scalar and vector
        // products for the different multipoles
        compute_vector_components(ix, iy, iz, q, d, Q);
        compute_vector_sums();
      }
    }
  }
}

#endif
