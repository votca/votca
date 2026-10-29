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
#ifndef VOTCA_XTP_LATTICEMODEL_H
#define VOTCA_XTP_LATTICEMODEL_H

// Standard includes
#include <array>
#include <cstdint>
#include <random>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "qmstate.h"
#include "topology.h"

namespace votca {
namespace xtp {

/**
 * \brief Surrogate lattice model of a hopping system, written as an ordinary
 * XTP topology (segments with one pseudo atom each, a periodic box and a
 * neighbour list with couplings), so that every transport method in XTP can be
 * run and cross-checked on it.
 *
 * Follows the setup of Chen, van der Hoorn, Baumeier, arXiv:2405.15836:
 * simple cubic lattice with spacing a, sites at (i+1/2)a, periodic in all
 * three directions (a time-of-flight calculation cuts one of them), pairs for
 * |r| < cutoff, |J|^2 = J0^2 exp(-(|r|-a)/decay), one reorganization energy
 * for all pairs, and Gaussian site energies without or with spatial
 * correlation (moving-average construction of Baumeier et al., PRB 86,
 * 184202 (2012); eq. 2.5 of the paper excludes the site itself from the
 * averages, see Options::include_self).
 *
 * All quantities are in atomic units (bohr, Hartree).
 */
class LatticeModel {
 public:
  enum class Disorder { None, Uncorrelated, Correlated };

  struct Options {
    std::array<Index, 3> dims = {8, 8, 8};
    double spacing = 1.0;       // a
    double cutoff = 1.5;        // pairs with |r| < cutoff
    double J2 = 0.0;            // |J|^2 at |r| = a
    double decay_length = 1.0;  // of |J|^2
    double reorganization = 0.0;
    double mean_energy = 0.0;
    double sigma = 0.0;
    Disorder disorder = Disorder::None;
    double eta_a = 0.2;  // weight of the uncorrelated part
    double eta_b = 0.4;  // weight of the short-range moving average
    Index l_b = 9;       // number of neighbours in the short-range average
    Index l_c = 280;     // number of neighbours in the long-range average
    /// the l_b and l_c sites of the averages include the site itself
    /// (Baumeier et al., PRB 86, 184202 (2012)); false: the l nearest other
    /// sites (eq. 2.5 of Chen et al., arXiv:2405.15836)
    bool include_self = true;
    std::uint64_t seed = 1;
    std::vector<QMStateType> states = {QMStateType::Electron,
                                       QMStateType::Hole};
    std::string segment_name = "site";
  };

  explicit LatticeModel(const Options& opt);

  /// Fill an empty topology with segments, box and neighbour list.
  void BuildTopology(Topology& top) const;

  /// Site energies in segment order
  const std::vector<double>& SiteEnergies() const { return energies_; }

  /// Positions in segment order
  const std::vector<Eigen::Vector3d>& Positions() const { return positions_; }

  Index NumberOfSites() const { return Index(positions_.size()); }

  static Disorder DisorderFromString(const std::string& name);

  /// Deterministic standard normal numbers (Box-Muller on mt19937_64).
  /// std::normal_distribution differs between standard libraries, this
  /// does not, so a seed gives the same lattice on every platform.
  class NormalGenerator {
   public:
    explicit NormalGenerator(std::uint64_t seed) : engine_(seed) {}
    double operator()();

   private:
    double Uniform();  // in (0,1]
    std::mt19937_64 engine_;
    bool has_spare_ = false;
    double spare_ = 0.0;
  };

 private:
  Options opt_;
  Eigen::Matrix3d box_;
  std::vector<Eigen::Vector3d> positions_;
  std::vector<double> energies_;

  Eigen::Vector3d ShortestConnection(const Eigen::Vector3d& r1,
                                     const Eigen::Vector3d& r2) const;
  void GenerateEnergies();
  /// indices of the l nearest other sites of every site, ties broken by index
  std::vector<std::vector<Index>> NearestNeighbours(Index l) const;
  /// lattice offsets (in units of a, reduced to the box) of all other sites
  /// within radius of a site, with their distances; by translation symmetry
  /// the same for every site, so neighbour searches cost O(n) instead of
  /// O(n^2)
  std::vector<std::pair<double, std::array<Index, 3>>> Offsets(
      double radius) const;
  /// site index of lattice coordinates (taken modulo the dimensions)
  Index SiteIndex(Index ix, Index iy, Index iz) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LATTICEMODEL_H
