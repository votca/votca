#ifndef QDMEWALD_INCLUDED
#define QDMEWALD_INCLUDED

#include "KSpace.hpp"
#include "RSpace.hpp"

template <class T>
class QDMEwald {
 public:
  /**
   *
   * Constructor for the class (base-version)
   *
   */
  QDMEwald() {}
  /**
   *
   * Constructor for the class
   *
   * @param alpha     Ewald sum splitting parameter
   * @param k_max     floating point cut-off radius for the k-space
   *                  part
   * @param r_max     floating point cut-off radius for the r-space
   *                  part
   */
  QDMEwald(T alpha, T k_max, T r_max)
      : kspace(KSpace<T>(alpha, k_max)), rspace(RSpace<T>(alpha, r_max)) {}

  /**
   *
   * Computation of the Ewald sum for quadrupoles, dipoles and monopoles
   *
   * @param l         system size (cubic system)
   * @param xyz       particle positions
   * @param q         particle charges
   * @param d         particle dipole moments
   * @param Q         particle quadrupole moments
   */
  void compute(const T, const std::vector<T>&, const std::vector<T>&,
               const std::vector<T>&, const std::vector<T>&);

 private:
  KSpace<T> kspace;
  RSpace<T> rspace;
};

/**
 *
 * Computation of the Ewald sum for quadrupoles, dipoles and monopoles
 *
 * @param l         system size (cubic system)
 * @param xyz       particle positions
 * @param q         particle charges
 * @param d         particle dipole moments
 * @param Q         particle quadrupole moments
 */
template <class T>
void QDMEwald<T>::compute(const T l, const std::vector<T>& xyz,
                          const std::vector<T>& q, const std::vector<T>& d,
                          const std::vector<T>& Q) {
  // compute k-space contribution
  kspace.compute(l, xyz, q, d, Q);
}
#endif
