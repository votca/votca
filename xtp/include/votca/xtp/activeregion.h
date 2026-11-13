/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_ACTIVEREGION_H
#define VOTCA_XTP_ACTIVEREGION_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "orbitals.h"

namespace votca {
namespace xtp {

/**
 * \brief Checks and grows the active region of a projection-based embedding
 * from the frontier orbitals of the full reference calculation.
 *
 * Atoms are indexed by their position in Orbitals::QMAtoms(); the
 * dft_in_dft.activeatoms indices are these positions.
 *
 * Without a truncated basis the embedded virtual space is the full one, so
 * only the occupied frontier orbitals decide whether an excitation fits into
 * the active region (its holes must lie there). With a truncated basis the
 * virtual frontier orbitals count as well.
 */
class ActiveRegion {
 public:
  struct Options {
    bool grow = false;
    // frontier orbitals: occupied levels within this window below the HOMO
    // (and, with include_virtuals, virtual levels within it above the
    // LUMO), in Hartree
    double window = 0.1 / 27.211386245988;
    // per-atom Mulliken population of a frontier orbital above which a
    // neighbouring atom is added
    double threshold = 0.0025;
    bool include_virtuals = false;
  };

  struct Report {
    std::vector<Index> atoms;        // final active atoms, sorted
    std::vector<Index> added;        // atoms added by the growth
    std::vector<Index> orbitals;     // frontier orbitals (MO indices)
    Eigen::VectorXd energies;        // their energies
    Eigen::VectorXd coverage_input;  // population on the input region
    Eigen::VectorXd coverage_final;  // population on the final region
    // atoms above the threshold in molecules without active atoms (no
    // bond path to the region; add them by hand if they belong to it)
    std::vector<Index> unreachable;
    // inactive heavy atoms bonded to two or more active atoms: leaving them
    // out cuts several bonds (e.g. a ring atom); grow adds them
    std::vector<Index> enclosed;
    // inactive hydrogens bonded to an active atom; grow adds them
    std::vector<Index> cut_hydrogens;
    std::string Format() const;
  };

  /// Bonded neighbours of each atom (covalent radii, factor 1.2).
  static std::vector<std::vector<Index>> BondGraph(const QMMolecule& mol);

  /// Mulliken population of each atom (rows) in each orbital (columns)
  /// of `mos` (AO coefficients as columns).
  static Eigen::MatrixXd AtomPopulations(const Orbitals& orb,
                                         const Eigen::MatrixXd& mos);

  /**
   * \brief Throws unless the active atoms of every molecule (connected
   * component of the bond graph) form a single connected piece: a region
   * that splits a molecule cuts it in more places than intended. Several
   * molecules, each partly or fully active, are fine.
   */
  static void CheckConnected(const QMMolecule& mol,
                             const std::vector<Index>& active);

  /**
   * \brief Checks the input region, and with options.grow adds every atom
   * whose population in a frontier orbital exceeds the threshold, in the
   * molecules that have active atoms, with the shortest bond path that
   * connects it to the region (so the region stays connected); then heavy
   * atoms bonded to two or more active atoms (so that no atom is cut out of
   * the region by several bonds), and all hydrogens bonded to active atoms.
   * `full` must hold the MOs of the full reference calculation.
   */
  static Report Evaluate(const Orbitals& full, const std::vector<Index>& active,
                         const Options& options);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ACTIVEREGION_H
