/*
 *            Copyright 2019 Forschungszentrum Juelich GmbH
 *            Copyright 2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_QDMEWALD_PRIVATE_H
#define VOTCA_XTP_QDMEWALD_PRIVATE_H

// Local private VOTCA includes
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

  /**
   *
   * Getter for the k-space contribution to energy
   *
   * @return std::vector<T> containing the forces
   */
  std::vector<T> get_forces();

  /**
   *
   * Getter for the k-space contribution to energy
   *
   * @return std::vector<T> containing the torque
   */
  std::vector<T> get_torque();

  /**
   *
   * Getter for the k-space contribution to energy
   *
   * @return std::vector<T> containing the sum over all forces
   */
  std::vector<T> get_total_force();

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
  total_energy += rspace.get_energy();
#ifndef QDMEWALD_DEBUG_ENABLED
  std::cout << "KSpace energy: " << kspace.get_energy() << std::endl;
  std::cout << "RSpace energy: " << rspace.get_energy() << std::endl;
#endif

  force = kspace.get_forces();
  torque = kspace.get_torque();

  std::vector<T> r_force = rspace.get_forces();
  std::vector<T> r_torque = rspace.get_torque();

  for (size_t i = 0; i < force.size(); ++i) {
    force.at(i) += r_force.at(i);
    torque.at(i) += r_torque.at(i);
  }
}

/**
 *
 * Getter for the total energy
 *
 * @return type T containing the energy
 */
template <class T>
T QDMEwald<T>::get_energy() {
  return total_energy;
}

/**
 *
 * Getter for the k-space contribution to energy
 *
 * @return std::vector<T> containing the forces
 */
template <class T>
std::vector<T> QDMEwald<T>::get_forces() {
  return force;
}

/**
 *
 * Getter for the k-space contribution to energy
 *
 * @return std::vector<T> containing the torque
 */
template <class T>
std::vector<T> QDMEwald<T>::get_torque() {
  return torque;
}

/**
 *
 * Getter for the k-space contribution to energy
 *
 * @return std::vector<T> containing the sum over all forces
 */
template <class T>
std::vector<T> QDMEwald<T>::get_total_force() {
  std::vector<T> total_force(3, (T)0);
  for (size_t i = 0; i < force.size(); ++i) {
#ifdef QDMEWALD_DEBUG_ENABLED
    if (i % 3 == 0) std::cout << "idx " << i / 3 << ": ";
    std::cout << force.at(i) << " ";
#endif
    total_force.at(i % 3) += force.at(i);
#ifdef QDMEWALD_DEBUG_ENABLED
    if (i % 3 == 2)
      std::cout << " => " << total_force.at(0) << " " << total_force.at(1)
                << " " << total_force.at(2) << std::endl;
#endif
  }
  return total_force;
}

#endif // VOTCA_XTP_QDMEWALD_PRIVATE_H

