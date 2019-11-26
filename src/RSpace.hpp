#ifndef RSPACE_CLASS_INCLUDED
#define RSPACE_CLASS_INCLUDED

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
   */
  RSpace(T alpha, T r_max) { init_params(alpha, r_max); }
  /*
   *
   * Initialization of the internal parameters
   *
   * @param alpha     Ewald sum splitting parameter
   * @param r_max     floating point cut-off radius for the r-space
   *                  part
   */
  void init_params(T, T);
  /**
   *
   * Computation of the kspace part of the Ewald sum, fills the
   * private variables of the class which can be accessed by the
   * corresponding getters after the computation
   *
   */
  void compute();

 private:
  T alpha;  ///< splitting parameter of the Ewald sum
  T r_max;  ///< floating-point version of the r-space cutoff

  T pot_energy;        ///< potential energy computed for the r-space part
  std::vector<T> f;    ///< forces computed for the r-space part
  std::vector<T> tqe;  ///< torque computed for the r-space part

  bool is_monopole;    ///< compute the monopole contributions
  bool is_dipole;      ///< compute the dipole contributions
  bool is_quadrupole;  ///< compute the quadrupole contributions
};

template <class T>
void RSpace<T>::init_params(T _alpha, T _r_max) {
  alpha = _alpha;
  r_max = _r_max;
}

#endif
