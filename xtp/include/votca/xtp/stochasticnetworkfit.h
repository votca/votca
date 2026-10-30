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
#ifndef VOTCA_XTP_STOCHASTICNETWORKFIT_H
#define VOTCA_XTP_STOCHASTICNETWORKFIT_H

// Standard includes
#include <cstdint>
#include <string>
#include <vector>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "eigen.h"
#include "qmstate.h"
#include "stochasticnetwork.h"
#include "topology.h"

namespace votca {
namespace xtp {

/**
 * \brief Fits the parameters of a StochasticNetwork to a microscopic
 * system (a state file with site energies and couplings).
 *
 * - Point process: hard core, mean and variance of the radii by a coarse
 *   grid and Nelder-Mead on the squared differences of the nearest
 *   neighbour distance distribution D(r), the spherical contact
 *   distribution H(r) and (weight 1/4) the pair correlation function g(r),
 *   averaged over a few realizations with fixed random numbers.
 * - Site energies (per carrier): mean, standard deviation and the spatial
 *   correlation kappa(r) = mean of z_i z_j over pairs at distance r. The
 *   moving average gives
 *   kappa(r) = w_b O_{l_b}(r) + w_c O_{l_c}(r), with O_l(r) the mean
 *   overlap |N_l(i) cap N_l(j)| / l of the neighbour sets of two sites at
 *   distance r, evaluated on the microscopic sites. For every (l_b, l_c) of
 *   the grids the weights (w_b, w_c >= 0, w_b + w_c <= 1) are a weighted
 *   least squares fit; the best pair is taken. The microscopic kappa is
 *   estimated with the sample mean and variance, which biases it by the mean
 *   correlation cbar = (w_a + w_b l_b + w_c l_c) / N over all pairs; the
 *   model is compared in the same form, (kappa - cbar) / (1 - cbar) (fixed
 *   point iteration), so that large model networks have the unbiased
 *   kappa; the standard deviation is corrected in the same way, sigma =
 *   s / sqrt(1 - cbar). Correlations over distances comparable to the box
 *   are not determined by the microscopic system; the fit reports when l_c
 *   is at the end of its grid. The site energies are the full site energies
 *   (electrostatic, polarization and internal).
 * - Connections: f(r) = (number of pairs) / (number of site pairs) in
 *   distance bins.
 * - Couplings (per carrier): y = log10(J^2/eV^2) linear in r, the squared
 *   residuals linear in r, quantiles of the standardized residuals.
 * - Two carriers: correlation of the site energies (Pearson) and of the
 *   coupling residuals (normal scores, Gaussian copula). The moving average
 *   of the second carrier is fitted to (kappa_2 - rho^2 kappa_1) /
 *   (1 - rho^2), so that the model reproduces kappa_2.
 * - Reorganization energies: means over the segments.
 * - Traps (Settings::traps, e.g. from TrapFinder): their sites are left out
 *   of the site energy statistics of all carriers (mean, sigma, kappa,
 *   carrier correlation), and each becomes a StochasticNetwork::Trap
 *   template (energies of its sites relative to the bulk mean, all
 *   carriers) with the density 1/N of the microscopic system.
 */
class StochasticNetworkFit {
 public:
  /// a trap of the microscopic system (from TrapFinder): its sites are
  /// left out of the bulk energy statistics and become a template
  struct TrapInput {
    Index carrier = 0;
    std::vector<Index> sites;   // deepest first
    double release_rate = 0.0;  // 1/s
  };

  struct Settings {
    std::vector<QMStateType> carriers;
    std::vector<TrapInput> traps;
    double bin = 0.05 * tools::conv::nm2bohr;
    /// largest distance for kappa(r) and O_l(r); 0: just below half the box
    double correlation_cutoff = 0.0;
    std::vector<Index> l_b;  // grids
    std::vector<Index> l_c;
    Index quantiles = 100;  // 0: Gaussian residuals
    /// false: the point process is not fitted but taken from points
    bool fit_points = true;
    StochasticNetwork::PointProcess points;
    Index realizations = 4;
    Index max_evaluations = 150;
    std::uint64_t seed = 1;
    double candidates = 3.0;
    Index max_iterations = 200;
  };

  /// microscopic system
  struct Data {
    Eigen::Vector3d box = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> positions;  // wrapped into the box
    std::vector<QMStateType> carriers;
    std::vector<std::vector<double>> energies;  // per carrier, Hartree
    std::vector<double> U_xX_nN;                // per carrier, means
    std::vector<double> U_nX_nN;
    std::vector<double> U_xN_xX;
    struct Pair {
      Index i;
      Index j;
      double dist;
      std::vector<double> J2;  // per carrier, Hartree^2
    };
    std::vector<Pair> pairs;
  };
  static Data FromTopology(const Topology& top,
                           const std::vector<QMStateType>& carriers);
  /// the same data from a model network (validation)
  static Data FromNetwork(const StochasticNetwork& net,
                          const std::vector<QMStateType>& carriers,
                          const Eigen::Vector3d& box);

  /// a function on a grid
  struct Curve {
    std::vector<double> r;  // bohr
    std::vector<double> value;
    std::vector<double> count;  // weights (pairs per bin), may be empty
  };

  struct EnergyFit {
    double mean = 0.0;
    double sigma = 0.0;        // sample standard deviation
    double sigma_model = 0.0;  // sigma / sqrt(1 - cbar) (see kappa)
    double cbar = 0.0;         // mean correlation over all pairs
    Index l_b = 1;
    Index l_c = 1;
    double w_a = 1.0;
    double w_b = 0.0;
    double w_c = 0.0;
    double rms = 0.0;  // rms deviation of the kappa fit
    Curve kappa;       // microscopic (sample mean and variance)
    /// the model on the microscopic sites, in the same (sample) form
    Curve kappa_model;
    bool fitted = false;  // false: sigma zero, no correlation fit
    /// rms of the best fit for each l_c (r: l_c)
    Curve l_c_scan;
    bool l_b_at_edge = false;  // the best l_b / l_c is the end of its grid
    bool l_c_at_edge = false;
  };

  struct CouplingFit {
    StochasticNetwork::CouplingModel model;
    Index used = 0;
    Index excluded = 0;  // J^2 <= 0
    double skewness = 0.0;
    double excess_kurtosis = 0.0;
  };

  struct PointFit {
    StochasticNetwork::PointProcess process;
    double objective = 0.0;
    Index evaluations = 0;
    Curve D, H, g;                    // microscopic
    Curve D_model, H_model, g_model;  // model, mean over the realizations
  };

  struct Result {
    StochasticNetwork::Options options;
    PointFit points;
    std::vector<EnergyFit> energies;
    std::vector<CouplingFit> couplings;
    Curve connection;
  };

  StochasticNetworkFit(const Data& data, const Settings& settings);
  Result Fit() const;

  // statistics, also used for the validation (all distances bohr)
  /// nearest neighbour distances
  static std::vector<double> NearestDistances(
      const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box);
  /// distances of `probes` random positions to the nearest point
  static std::vector<double> ContactDistances(
      const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box,
      Index probes, std::uint64_t seed);
  /// empirical distribution function on a grid
  static Curve Cdf(const std::vector<double>& x, const std::vector<double>& r);
  /// pair correlation function in bins of width bin up to rmax
  static Curve PairCorrelation(const std::vector<Eigen::Vector3d>& points,
                               const Eigen::Vector3d& box, double bin,
                               double rmax);
  /// kappa(r) of standardized values z in bins up to rmax
  /// (pairs with a site marked in skip are left out)
  static Curve Correlation(const std::vector<Eigen::Vector3d>& points,
                           const Eigen::Vector3d& box,
                           const std::vector<double>& z, double bin,
                           double rmax, const std::vector<char>& skip = {});
  /// inverse of the standard normal distribution function
  static double NormalQuantile(double u);
  /// overlap curves O_l(r) for the given l
  static std::vector<Curve> Overlaps(const std::vector<Eigen::Vector3d>& points,
                                     const Eigen::Vector3d& box,
                                     const std::vector<Index>& l, double bin,
                                     double rmax);
  /// weights w_b, w_c >= 0, w_b + w_c <= 1 minimizing
  /// sum count (w_b Ob + w_c Oc - kappa)^2; returns the weighted mean square
  static double FitWeights(const Curve& kappa, const Curve& Ob, const Curve& Oc,
                           double& w_b, double& w_c);
  /// pairs per site in distance bins (connected pairs of the data)
  static Curve PairsPerSite(const Data& data, double bin, double rmax);
  /// mean of log10(J^2/eV^2) in distance bins
  static Curve MeanLog10J2(const Data& data, Index carrier, double bin,
                           double rmax);

 private:
  const Data& data_;
  Settings settings_;

  PointFit FitPoints() const;
  Curve FitConnection() const;
  CouplingFit FitCoupling(Index carrier,
                          std::vector<double>& standardized) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STOCHASTICNETWORKFIT_H
