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
#ifndef VOTCA_XTP_STOCHASTICNETWORK_H
#define VOTCA_XTP_STOCHASTICNETWORK_H

// Standard includes
#include <array>
#include <cstdint>
#include <functional>
#include <random>
#include <string>
#include <vector>

// Local VOTCA includes
#include "eigen.h"
#include "qmstate.h"
#include "topology.h"

namespace votca {
namespace tools {
class Property;
}
}  // namespace votca

namespace votca {
namespace xtp {

/**
 * \brief Platform independent random numbers for the stochastic models:
 * mt19937_64 with self-written transforms (the std:: distributions differ
 * between standard libraries).
 */
class ModelRandom {
 public:
  /// independent streams for the same seed (points, energies, pairs)
  ModelRandom(std::uint64_t seed, std::uint64_t stream);
  /// in (0,1]
  double Uniform() {
    return (double(engine_() >> 11) + 1.0) * (1.0 / 9007199254740992.0);
  }
  double Normal();
  /// Poisson distributed count
  Index Poisson(double mean);
  /// integer in [0, n)
  Index Integer(Index n);

 private:
  std::mt19937_64 engine_;
  bool has_spare_ = false;
  double spare_ = 0.0;
};

/**
 * \brief Periodic orthorhombic cell list for neighbour searches.
 */
class PeriodicCells {
 public:
  /// box lengths and the largest search radius used with Within()
  PeriodicCells(const Eigen::Vector3d& box, double radius,
                const std::vector<Eigen::Vector3d>& points);
  /// shortest periodic connection r2 - r1
  Eigen::Vector3d Connection(const Eigen::Vector3d& r1,
                             const Eigen::Vector3d& r2) const;
  /// calls f(j, connection, distance) for all points within radius of p
  /// (radius <= the radius given at construction)
  void Within(const Eigen::Vector3d& p, double radius,
              const std::function<void(Index, const Eigen::Vector3d&, double)>&
                  f) const;
  /// indices of the l points nearest to point i (including i itself),
  /// ties broken by index
  std::vector<Index> Nearest(Index i, Index l) const;
  /// indices of the l points nearest to any position p
  std::vector<Index> NearestTo(const Eigen::Vector3d& p, Index l) const;
  Index NumberOfPoints() const { return Index(points_.size()); }

 private:
  Eigen::Vector3d box_;
  std::array<Index, 3> n_;
  Eigen::Vector3d cell_;
  const std::vector<Eigen::Vector3d>& points_;
  std::vector<Index> start_;
  std::vector<Index> members_;
  Index CellIndex(const Eigen::Vector3d& p, std::array<Index, 3>& c) const;
  template <class F>
  void Shell(const std::array<Index, 3>& c, Index s, F&& f) const;
};

/**
 * \brief Stochastic network model of a molecular charge transport network
 * (Baumeier, Stenzel, Poelking, Andrienko, Schmidt, Phys. Rev. B 86, 184202
 * (2012); Kordt et al., J. Chem. Theory Comput. 10, 2508 (2014)): a random
 * weighted graph with the statistics of a microscopic (atomistic) system,
 * in a periodic box of any size.
 *
 * - Sites: iterative dominance-competition point process. Poisson points
 *   get radii R = r_h + X with X Gamma distributed (by inversion); a point
 * inside the sphere of a point with a larger radius is removed (all points of
 * an iteration count); later iterations fill the remaining space without
 *   violating the spheres of accepted points, until the number of sites is
 *   reached.
 * - Site energies per carrier: moving average of three Gaussian fields,
 *   E_i = mean + sigma (sqrt(w_a) A_i + sqrt(w_b/l_b) sum_{N_b(i)} B_j +
 *   sqrt(w_c/l_c) sum_{N_c(i)} C_j), N_l(i) the l sites nearest to i,
 *   i itself included, w_c = 1 - w_a - w_b. A second carrier can be
 *   correlated with the first: z_2 = rho z_1 + sqrt(1 - rho^2) w with w a
 *   moving average with the second carrier's parameters.
 * - Pairs: two sites at distance r are connected with probability f(r)
 *   (table, linear interpolation).
 * - Couplings: log10(J^2 / eV^2) = m(r) + s(r) Q(u), with m and s^2 linear
 *   in r and Q the quantile function of the standardized residuals
 *   (table; a standard normal if empty), u = Phi(g) for a standard normal
 *   g; g of two carriers correlated with coefficient rho_J (Gaussian
 *   copula).
 * - Reorganization and internal energies per carrier: constants.
 * - Traps (optional): templates of trap sites or clusters from a
 *   microscopic system (TrapFinder) are inserted with a density per site:
 *   floor(density N) + Bernoulli(rest) traps of each template, each at a
 *   random free site and its nearest free sites, whose energies (all
 *   carriers) are set to the bulk mean plus the template offsets.
 *
 * Options are in atomic units (bohr, Hartree) except the coupling model,
 * which is in log10(J^2 / eV^2) per nm as in the papers.
 */
class StochasticNetwork {
 public:
  struct PointProcess {
    double hard_core = 0.0;        // r_h
    double radius_mean = 0.0;      // mean of X
    double radius_variance = 0.0;  // variance of X (0: X = mean)
    /// candidates per iteration in units of the number of sites
    double candidates = 3.0;
    Index max_iterations = 200;
  };

  struct EnergyModel {
    double mean = 0.0;
    double sigma = 0.0;  // 0: all sites at mean
    double w_a = 1.0;
    double w_b = 0.0;
    Index l_b = 1;
    Index l_c = 1;
  };

  struct CouplingModel {
    // log10(J^2/eV^2) mean = mean_offset + mean_slope r[nm]
    double mean_offset = -6.0;
    double mean_slope = 0.0;
    // variance = var_offset + var_slope r[nm] (clipped at 0)
    double var_offset = 1.0;
    double var_slope = 0.0;
    /// quantiles of the standardized residual at u = (k + 0.5) / K,
    /// k = 0..K-1 (empty: standard normal)
    std::vector<double> residual_quantiles;
  };

  struct Carrier {
    QMStateType state;
    EnergyModel energy;
    CouplingModel coupling;
    double U_xX_nN = 0.0;  // internal energy
    double U_nX_nN = 0.0;  // reorganization parts
    double U_xN_xX = 0.0;
  };

  /// a trap template: energies of its sites relative to the mean of the
  /// (Gaussian) bulk, for every carrier, taken from a microscopic system
  struct Trap {
    Index carrier = 0;     // the carrier it traps (index into carriers)
    double density = 0.0;  // traps of this kind per site
    /// offsets[k][c]: site k (deepest first), carrier c, Hartree
    std::vector<std::vector<double>> offsets;
    double release_rate = 0.0;  // microscopic release rate (1/s, information)
  };

  struct Options {
    Eigen::Vector3d box = Eigen::Vector3d::Zero();  // orthorhombic, bohr
    Index sites = 0;
    PointProcess points;
    /// connection probability table: f at distances r (ascending, bohr)
    std::vector<double> connection_r;
    std::vector<double> connection_p;
    std::vector<Carrier> carriers;
    /// correlation of the site energies / coupling residuals of the second
    /// carrier with the first
    double energy_correlation = 0.0;
    double coupling_correlation = 0.0;
    std::uint64_t seed = 1;
    std::string segment_name = "site";
    /// trap templates, inserted after the Gaussian energies
    std::vector<Trap> traps;
    /// multiplies all trap densities (0: no traps)
    double trap_factor = 1.0;
  };

  explicit StochasticNetwork(const Options& opt);

  Index NumberOfSites() const { return Index(positions_.size()); }
  const std::vector<Eigen::Vector3d>& Positions() const { return positions_; }
  /// radii of the point process (r_h + X)
  const std::vector<double>& Radii() const { return radii_; }
  /// inserted traps: template index and sites (deepest first)
  struct InsertedTrap {
    Index type;
    std::vector<Index> sites;
  };
  const std::vector<InsertedTrap>& Traps() const { return inserted_; }
  /// site energies of carrier c (index into Options::carriers)
  const std::vector<double>& SiteEnergies(Index c) const {
    return energies_[c];
  }
  struct Pair {
    Index i;
    Index j;
    Eigen::Vector3d R;       // shortest connection from i to j
    std::vector<double> J2;  // per carrier, Hartree^2
  };
  const std::vector<Pair>& Pairs() const { return pairs_; }
  /// iterations the point process needed
  Index Iterations() const { return iterations_; }

  /// fill an empty topology (segments with one pseudo atom, box, pairs)
  void BuildTopology(Topology& top) const;

  /// options from a parameter file written by networkfit (units nm, eV);
  /// sites and box are left as in the file and can be overwritten
  static Options FromProperty(const tools::Property& p);
  /// the parameter file (root element "stochasticnetwork")
  static std::string ToXML(const Options& opt);

  /// standardized moving average z_i (unit variance) of the standard normal
  /// fields A, B, C (one value per site); E_i = mean + sigma z_i
  static std::vector<double> MovingAverage(const PeriodicCells& cells,
                                           const EnergyModel& m, Index sites,
                                           const std::vector<double>& A,
                                           const std::vector<double>& Bf,
                                           const std::vector<double>& C);
  /// one realization of the point process (as used by the constructor);
  /// returns the number of iterations, throws if the box is saturated
  static Index SamplePoints(const PointProcess& pp, const Eigen::Vector3d& box,
                            Index sites, ModelRandom& random,
                            std::vector<Eigen::Vector3d>& positions,
                            std::vector<double>& radii);
  /// connection probability at distance r (bohr)
  double Connection(double r) const;
  /// standardized residual for a standard normal g
  static double Residual(const CouplingModel& m, double g);
  /// log10(J^2/eV^2) for distance r (bohr) and standard normal g
  static double Log10J2(const CouplingModel& m, double r, double g);

 private:
  Options opt_;
  std::vector<Eigen::Vector3d> positions_;
  std::vector<double> radii_;
  std::vector<std::vector<double>> energies_;
  std::vector<Pair> pairs_;
  Index iterations_ = 0;
  std::vector<InsertedTrap> inserted_;

  void Check() const;
  void GenerateEnergies(ModelRandom& random);
  void GeneratePairs(ModelRandom& random);
  void InsertTraps(ModelRandom& random);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STOCHASTICNETWORK_H
