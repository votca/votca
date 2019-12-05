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
   * @param l         system size (cubic system)
   */
  QDMEwald(T alpha, T k_max, T r_max, T l)
      : kspace(KSpace<T>(alpha, k_max, l)),
        rspace(RSpace<T>(alpha, r_max, l)) {}

  /**
   *
   * Computation of the Ewald sum for quadrupoles, dipoles and monopoles
   *
   * @param xyz       particle positions
   * @param q         particle charges
   * @param d         particle dipole moments
   * @param Q         particle quadrupole moments
   */
  void compute(const std::vector<T>&, const std::vector<T>&,
               const std::vector<T>&, const std::vector<T>&);

  /**
   *
   * Getter for the k-space contribution to energy
   *
   * @return type T containing the energy
   */

  T get_energy();

 private:
  KSpace<T> kspace;
  RSpace<T> rspace;

  std::vector<T> force;
  std::vector<T> torque;

  T total_energy;
  T total_virial;
};

/**
 *
 * Computation of the Ewald sum for quadrupoles, dipoles and monopoles
 *
 * @param xyz       particle positions
 * @param q         particle charges
 * @param d         particle dipole moments
 * @param Q         particle quadrupole moments
 */
template <class T>
void QDMEwald<T>::compute(const std::vector<T>& xyz, const std::vector<T>& q,
                          const std::vector<T>& d, const std::vector<T>& Q) {

  total_energy = 0.0;
  total_virial = 0.0;

  // compute k-space contribution
  kspace.compute(xyz, q, d, Q);
  // compute r-space contribution
  rspace.compute(xyz, q, d, Q);

  total_energy += kspace.get_energy();
  std::cout << "KSpace energy: " << kspace.get_energy() << std::endl;
  total_energy += rspace.get_energy();
  std::cout << "RSpace energy: " << rspace.get_energy() << std::endl;

  force = kspace.get_forces();
  torque = kspace.get_torque();

  std::vector<T> total_force(3, (T)0);

  std::vector<T> r_force = rspace.get_forces();
  std::vector<T> r_torque = rspace.get_torque();

  for (int i = 0; i < force.size(); ++i) {
    force.at(i) += r_force.at(i);
    torque.at(i) += r_torque.at(i);
    total_force.at(i % 3) += r_force.at(i);
  }

  std::cout << "total forces: " << total_force.at(0) << " " << total_force.at(1)
            << " " << total_force.at(2) << std::endl;
}

/**
 *
 * Getter for the k-space contribution to energy
 *
 * @return type T containing the energy
 */
template <class T>
T QDMEwald<T>::get_energy() {
  return total_energy;
}

#endif
