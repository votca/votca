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

// Standard includes
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

// Third party includes
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/citations.h"
#include "votca/xtp/stochasticnetwork.h"

namespace votca {
namespace xtp {

// ---------------------------------------------------------------- random

ModelRandom::ModelRandom(std::uint64_t seed, std::uint64_t stream) {
  // std::seed_seq and mt19937_64 are fully specified by the standard
  std::seed_seq seq{
      std::uint32_t(seed & 0xffffffffu), std::uint32_t(seed >> 32),
      std::uint32_t(stream & 0xffffffffu), std::uint32_t(stream >> 32)};
  engine_.seed(seq);
}

double ModelRandom::Normal() {
  if (has_spare_) {
    has_spare_ = false;
    return spare_;
  }
  double r = std::sqrt(-2.0 * std::log(Uniform()));
  double phi = 2.0 * tools::conv::Pi * Uniform();
  spare_ = r * std::sin(phi);
  has_spare_ = true;
  return r * std::cos(phi);
}

Index ModelRandom::Poisson(double mean) {
  if (!(mean >= 0.0) || !std::isfinite(mean)) {
    throw std::runtime_error("ModelRandom::Poisson: mean must be >= 0");
  }
  // sum of Poisson variables with means <= 500 (Knuth's product of uniforms)
  Index n = 0;
  while (mean > 0.0) {
    double m = std::min(mean, 500.0);
    mean -= m;
    const double limit = std::exp(-m);
    double p = Uniform();
    while (p > limit) {
      ++n;
      p *= Uniform();
    }
  }
  return n;
}

Index ModelRandom::Integer(Index n) {
  if (n <= 0) {
    throw std::runtime_error("ModelRandom::Integer: n must be positive");
  }
  const std::uint64_t un = std::uint64_t(n);
  const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max() -
                              std::numeric_limits<std::uint64_t>::max() % un;
  std::uint64_t x;
  do {
    x = engine_();
  } while (x >= limit);
  return Index(x % un);
}

// ---------------------------------------------------------------- cells

PeriodicCells::PeriodicCells(const Eigen::Vector3d& box, double radius,
                             const std::vector<Eigen::Vector3d>& points)
    : box_(box), points_(points) {
  if (!(box.minCoeff() > 0.0)) {
    throw std::runtime_error("PeriodicCells: box lengths must be positive");
  }
  if (!(radius > 0.0) || !(2.0 * radius < box.minCoeff())) {
    throw std::runtime_error(
        "PeriodicCells: the search radius must be positive and smaller than "
        "half the box");
  }
  for (Index d = 0; d < 3; ++d) {
    n_[d] = std::max(Index(1), Index(std::floor(box[d] / radius)));
  }
  // not many more cells than points (larger cells are fine for Within)
  while (double(n_[0]) * double(n_[1]) * double(n_[2]) >
         2.0 * double(points.size()) + 27.0) {
    Index d = Index(std::max_element(n_.begin(), n_.end()) - n_.begin());
    n_[d] = std::max(Index(1), n_[d] / 2);
  }
  for (Index d = 0; d < 3; ++d) {
    cell_[d] = box[d] / double(n_[d]);
  }
  Index ncells = n_[0] * n_[1] * n_[2];
  std::vector<Index> cell_of(points.size());
  start_.assign(ncells + 1, 0);
  std::array<Index, 3> c;
  for (Index i = 0; i < Index(points.size()); ++i) {
    cell_of[i] = CellIndex(points[i], c);
    ++start_[cell_of[i] + 1];
  }
  for (Index k = 0; k < ncells; ++k) {
    start_[k + 1] += start_[k];
  }
  members_.resize(points.size());
  std::vector<Index> fill(start_.begin(), start_.end() - 1);
  for (Index i = 0; i < Index(points.size()); ++i) {
    members_[fill[cell_of[i]]++] = i;
  }
}

Index PeriodicCells::CellIndex(const Eigen::Vector3d& p,
                               std::array<Index, 3>& c) const {
  for (Index d = 0; d < 3; ++d) {
    double x = p[d] - box_[d] * std::floor(p[d] / box_[d]);
    Index k = Index(std::floor(x / cell_[d]));
    c[d] = std::min(std::max(k, Index(0)), n_[d] - 1);
  }
  return (c[0] * n_[1] + c[1]) * n_[2] + c[2];
}

template <class F>
void PeriodicCells::Shell(const std::array<Index, 3>& c, Index s, F&& f) const {
  // offsets -lo..hi cover every cell exactly once also for few cells
  std::array<Index, 3> lo, hi;
  for (Index d = 0; d < 3; ++d) {
    lo[d] = std::min(s, (n_[d] - 1) / 2);
    hi[d] = std::min(s, n_[d] - 1 - (n_[d] - 1) / 2);
  }
  for (Index o0 = -lo[0]; o0 <= hi[0]; ++o0) {
    for (Index o1 = -lo[1]; o1 <= hi[1]; ++o1) {
      for (Index o2 = -lo[2]; o2 <= hi[2]; ++o2) {
        if (std::max({std::abs(o0), std::abs(o1), std::abs(o2)}) != s) {
          continue;
        }
        Index k0 = (c[0] + o0 + n_[0]) % n_[0];
        Index k1 = (c[1] + o1 + n_[1]) % n_[1];
        Index k2 = (c[2] + o2 + n_[2]) % n_[2];
        f((k0 * n_[1] + k1) * n_[2] + k2);
      }
    }
  }
}

Eigen::Vector3d PeriodicCells::Connection(const Eigen::Vector3d& r1,
                                          const Eigen::Vector3d& r2) const {
  Eigen::Vector3d d = r2 - r1;
  for (Index k = 0; k < 3; ++k) {
    // points in the box: one image shift at most
    if (d[k] > 0.5 * box_[k]) {
      d[k] -= box_[k];
    } else if (d[k] < -0.5 * box_[k]) {
      d[k] += box_[k];
    }
    if (std::abs(d[k]) > 0.5 * box_[k]) {
      d[k] -= box_[k] * std::round(d[k] / box_[k]);
    }
  }
  return d;
}

void PeriodicCells::Within(
    const Eigen::Vector3d& p, double radius,
    const std::function<void(Index, const Eigen::Vector3d&, double)>& f) const {
  if (radius > cell_.minCoeff() * (1.0 + 1e-12)) {
    throw std::runtime_error(
        "PeriodicCells::Within: radius larger than the cell list radius");
  }
  std::array<Index, 3> c;
  CellIndex(p, c);
  for (Index s = 0; s <= 1; ++s) {
    Shell(c, s, [&](Index cell) {
      for (Index k = start_[cell]; k < start_[cell + 1]; ++k) {
        Index j = members_[k];
        Eigen::Vector3d conn = Connection(p, points_[j]);
        double dist = conn.norm();
        if (dist < radius) {
          f(j, conn, dist);
        }
      }
    });
  }
}

std::vector<Index> PeriodicCells::Nearest(Index i, Index l) const {
  return NearestTo(points_[i], l);
}

std::vector<Index> PeriodicCells::NearestTo(const Eigen::Vector3d& p,
                                            Index l) const {
  const Index n = Index(points_.size());
  if (l > n) {
    throw std::runtime_error("PeriodicCells::Nearest: more neighbours (" +
                             std::to_string(l) + ") than points (" +
                             std::to_string(n) + ")");
  }
  if (l <= 0) {
    return {};
  }
  std::array<Index, 3> c;
  CellIndex(p, c);
  std::vector<std::pair<double, Index>> cand;
  const double inf = std::numeric_limits<double>::infinity();
  for (Index s = 0;; ++s) {
    Shell(c, s, [&](Index cell) {
      for (Index k = start_[cell]; k < start_[cell + 1]; ++k) {
        Index j = members_[k];
        cand.emplace_back(Connection(p, points_[j]).norm(), j);
      }
    });
    // every point outside the shells 0..s is farther than `guaranteed`
    double guaranteed = inf;
    bool all = true;
    for (Index d = 0; d < 3; ++d) {
      Index full = std::max((n_[d] - 1) / 2, n_[d] - 1 - (n_[d] - 1) / 2);
      if (s < full) {
        all = false;
        guaranteed = std::min(guaranteed, double(s) * cell_[d]);
      }
    }
    if (all) {
      break;
    }
    if (Index(cand.size()) >= l &&
        std::count_if(cand.begin(), cand.end(),
                      [guaranteed](const std::pair<double, Index>& x) {
                        return x.first <= guaranteed;
                      }) >= l) {
      break;
    }
  }
  std::partial_sort(cand.begin(), cand.begin() + l, cand.end());
  std::vector<Index> result(l);
  for (Index k = 0; k < l; ++k) {
    result[k] = cand[k].second;
  }
  return result;
}

// ------------------------------------------------------------ the model

StochasticNetwork::StochasticNetwork(const Options& opt) : opt_(opt) {
  Check();
  citations::Cite("Baumeier:2012");
  citations::Cite("Kordt:2014");
  if (!opt_.traps.empty() && opt_.trap_factor > 0.0) {
    citations::Cite("Chen:2025b");
  }
  ModelRandom points(opt_.seed, 0);
  iterations_ = SamplePoints(opt_.points, opt_.box, opt_.sites, points,
                             positions_, radii_);
  ModelRandom energies(opt_.seed, 1);
  GenerateEnergies(energies);
  ModelRandom traps(opt_.seed, 3);
  InsertTraps(traps);
  ModelRandom pairs(opt_.seed, 2);
  GeneratePairs(pairs);
}

void StochasticNetwork::Check() const {
  auto fail = [](const std::string& m) {
    throw std::runtime_error("StochasticNetwork: " + m);
  };
  if (!(opt_.box.minCoeff() > 0.0)) {
    fail("box lengths must be positive");
  }
  if (opt_.sites <= 0) {
    fail("need at least one site");
  }
  const PointProcess& pp = opt_.points;
  if (!(pp.hard_core >= 0.0) || !(pp.radius_mean >= 0.0) ||
      !(pp.radius_variance >= 0.0) || !(pp.candidates > 0.0) ||
      pp.max_iterations <= 0) {
    fail(
        "point process: hard core, radius mean and variance must be >= 0, "
        "candidates and iterations > 0");
  }
  if (pp.radius_variance > 0.0 && pp.radius_mean == 0.0) {
    fail("point process: a radius variance needs a positive radius mean");
  }
  if (opt_.connection_r.empty() ||
      opt_.connection_r.size() != opt_.connection_p.size()) {
    fail("connection table: r and p need the same, nonzero length");
  }
  for (std::size_t k = 0; k < opt_.connection_r.size(); ++k) {
    if (!(opt_.connection_p[k] >= 0.0 && opt_.connection_p[k] <= 1.0)) {
      fail("connection probabilities must be in [0, 1]");
    }
    if (!(opt_.connection_r[k] >= 0.0) ||
        (k > 0 && !(opt_.connection_r[k] > opt_.connection_r[k - 1]))) {
      fail("connection distances must be >= 0 and ascending");
    }
  }
  if (!(2.0 * opt_.connection_r.back() < opt_.box.minCoeff())) {
    fail("the largest connection distance must be smaller than half the box");
  }
  if (opt_.carriers.empty()) {
    fail("need at least one carrier");
  }
  for (std::size_t c = 0; c < opt_.carriers.size(); ++c) {
    const Carrier& car = opt_.carriers[c];
    if (!car.state.isKMCState()) {
      fail(car.state.ToString() + " is not a hopping state (e, h, s, t)");
    }
    for (std::size_t d = 0; d < c; ++d) {
      if (opt_.carriers[d].state == car.state) {
        fail("carrier " + car.state.ToString() + " given twice");
      }
    }
    const EnergyModel& e = car.energy;
    double w_c = 1.0 - e.w_a - e.w_b;
    if (!(e.sigma >= 0.0) || !(e.w_a >= 0.0) || !(e.w_b >= 0.0) ||
        w_c < -1e-12) {
      fail("energy model: sigma, w_a, w_b >= 0 and w_a + w_b <= 1");
    }
    if (e.l_b <= 0 || e.l_c <= 0 || (e.w_b > 0.0 && e.l_b > opt_.sites) ||
        (w_c > 1e-12 && e.l_c > opt_.sites)) {
      fail("energy model: l_b and l_c must be in [1, number of sites]");
    }
    const CouplingModel& m = car.coupling;
    for (std::size_t k = 1; k < m.residual_quantiles.size(); ++k) {
      if (m.residual_quantiles[k] < m.residual_quantiles[k - 1]) {
        fail("coupling model: residual quantiles must be ascending");
      }
    }
  }
  if (!(opt_.trap_factor >= 0.0)) {
    fail("the trap factor must be >= 0");
  }
  for (const Trap& t : opt_.traps) {
    if (t.carrier < 0 || t.carrier >= Index(opt_.carriers.size()) ||
        !(t.density >= 0.0) || t.offsets.empty()) {
      fail("traps need a carrier, a density >= 0 and at least one site");
    }
    for (const auto& o : t.offsets) {
      if (o.size() != opt_.carriers.size()) {
        fail("trap offsets need one value per carrier");
      }
    }
  }
  if (!(std::abs(opt_.energy_correlation) <= 1.0) ||
      !(std::abs(opt_.coupling_correlation) <= 1.0)) {
    fail("correlations must be in [-1, 1]");
  }
}

Index StochasticNetwork::SamplePoints(const PointProcess& pp,
                                      const Eigen::Vector3d& box, Index n,
                                      ModelRandom& random,
                                      std::vector<Eigen::Vector3d>& positions,
                                      std::vector<double>& radii) {
  // radii by inversion of the Gamma distribution function: continuous in the
  // parameters for the same random numbers (common random numbers in fits);
  // tabulated quantiles, exact in the tails
  const bool gamma = pp.radius_variance > 0.0;
  const double shape =
      gamma ? pp.radius_mean * pp.radius_mean / pp.radius_variance : 0.0;
  const double scale = gamma ? pp.radius_variance / pp.radius_mean : 0.0;
  constexpr Index M = 512;
  std::vector<double> table;
  if (gamma) {
    table.resize(M + 1);
    for (Index k = 1; k < M; ++k) {
      table[k] = boost::math::gamma_p_inv(shape, double(k) / double(M));
    }
  }
  auto radius = [&]() {
    if (!gamma) {
      return pp.hard_core + pp.radius_mean;
    }
    double u = random.Uniform();
    double t = u * double(M);
    Index k = Index(t);
    double x;
    if (k < 1 || k >= M - 1) {
      x = boost::math::gamma_p_inv(shape, std::min(u, 1.0 - 1e-16));
    } else {
      x = table[k] + (t - double(k)) * (table[k + 1] - table[k]);
    }
    return pp.hard_core + x * scale;
  };
  const double half_box = 0.5 * box.minCoeff();
  positions.clear();
  radii.clear();
  positions.reserve(n);
  radii.reserve(n);
  double max_accepted = 0.0;
  std::vector<Eigen::Vector3d> cand;
  std::vector<double> cand_r;
  std::vector<char> keep;
  std::vector<Index> survivors;
  for (Index it = 1; it <= pp.max_iterations; ++it) {
    Index m = random.Poisson(pp.candidates * double(n));
    cand.resize(m);
    cand_r.resize(m);
    double max_cand = 0.0;
    for (Index k = 0; k < m; ++k) {
      for (Index d = 0; d < 3; ++d) {
        cand[k][d] = (1.0 - random.Uniform()) * box[d];
      }
      cand_r[k] = radius();
      max_cand = std::max(max_cand, cand_r[k]);
    }
    double search = std::max(max_cand, max_accepted);
    if (!(search < half_box)) {
      throw std::runtime_error("StochasticNetwork: point radii (" +
                               std::to_string(search * tools::conv::bohr2nm) +
                               " nm) must be smaller than half the box");
    }
    // a zero search radius means no exclusion
    search = std::max(search, 1e-9 * half_box);
    keep.assign(m, 1);
    // against accepted points: distance >= max(R_candidate, R_accepted)
    if (!positions.empty()) {
      PeriodicCells acc(box, search, positions);
      for (Index k = 0; k < m; ++k) {
        acc.Within(cand[k], std::max(cand_r[k], max_accepted),
                   [&](Index j, const Eigen::Vector3d&, double dist) {
                     if (dist < std::max(cand_r[k], radii[j])) {
                       keep[k] = 0;
                     }
                   });
      }
    }
    // dominance among all candidates of the iteration: inside the sphere of
    // a larger point (equal radii: the earlier point wins)
    {
      PeriodicCells cells(box, std::max(max_cand, 1e-9 * half_box), cand);
      for (Index k = 0; k < m; ++k) {
        if (!keep[k]) {
          continue;  // removed already; it still dominates others
        }
        cells.Within(cand[k], std::max(max_cand, 1e-9 * half_box),
                     [&](Index j, const Eigen::Vector3d&, double dist) {
                       if (j == k) {
                         return;
                       }
                       bool larger = cand_r[j] > cand_r[k] ||
                                     (j < k && cand_r[j] == cand_r[k]);
                       if (larger && dist < cand_r[j]) {
                         keep[k] = 0;
                       }
                     });
      }
    }
    survivors.clear();
    for (Index k = 0; k < m; ++k) {
      if (keep[k]) {
        survivors.push_back(k);
      }
    }
    Index need = n - Index(positions.size());
    if (Index(survivors.size()) > need) {
      // random subset (partial Fisher-Yates)
      for (Index k = 0; k < need; ++k) {
        Index r = k + random.Integer(Index(survivors.size()) - k);
        std::swap(survivors[k], survivors[r]);
      }
      survivors.resize(need);
    }
    for (Index k : survivors) {
      positions.push_back(cand[k]);
      radii.push_back(cand_r[k]);
      max_accepted = std::max(max_accepted, cand_r[k]);
    }
    if (Index(positions.size()) == n) {
      return it;
    }
  }
  throw std::runtime_error(
      "StochasticNetwork: the point process placed only " +
      std::to_string(positions.size()) + " of " + std::to_string(n) +
      " sites in " + std::to_string(pp.max_iterations) +
      " iterations; the box is too small for the radii (increase the box or "
      "max_iterations)");
}

std::vector<double> StochasticNetwork::MovingAverage(
    const PeriodicCells& cells, const EnergyModel& m, Index sites,
    const std::vector<double>& A, const std::vector<double>& Bf,
    const std::vector<double>& C) {
  if (Index(A.size()) != sites || Index(Bf.size()) != sites ||
      Index(C.size()) != sites) {
    throw std::runtime_error(
        "StochasticNetwork::MovingAverage: need one field value per site");
  }
  const double w_c = std::max(0.0, 1.0 - m.w_a - m.w_b);
  const double fa = std::sqrt(m.w_a);
  const double fb = std::sqrt(m.w_b / double(m.l_b));
  const double fc = std::sqrt(w_c / double(m.l_c));
  const Index lb = m.w_b > 0.0 ? m.l_b : 0;
  const Index lc = w_c > 0.0 ? m.l_c : 0;
  const Index lmax = std::max(lb, lc);
  std::vector<double> z(sites);
#pragma omp parallel for schedule(dynamic, 256)
  for (Index i = 0; i < sites; ++i) {
    std::vector<Index> nb = cells.Nearest(i, lmax);
    double b = 0.0;
    for (Index k = 0; k < lb; ++k) {
      b += Bf[nb[k]];
    }
    double c = 0.0;
    for (Index k = 0; k < lc; ++k) {
      c += C[nb[k]];
    }
    z[i] = fa * A[i] + fb * b + fc * c;
  }
  return z;
}

void StochasticNetwork::GenerateEnergies(ModelRandom& random) {
  const Index n = NumberOfSites();
  // cells of about eight points for the nearest-neighbour searches
  double volume = opt_.box.prod();
  double radius =
      std::min(0.49 * opt_.box.minCoeff(), 2.0 * std::cbrt(volume / double(n)));
  PeriodicCells cells(opt_.box, radius, positions_);
  energies_.clear();
  std::vector<double> first;
  std::vector<double> A(n), B(n), C(n);
  for (std::size_t c = 0; c < opt_.carriers.size(); ++c) {
    const EnergyModel& m = opt_.carriers[c].energy;
    // the fields are always drawn, so that the energies of one carrier do
    // not depend on the model of another
    for (std::vector<double>* f : {&A, &B, &C}) {
      for (double& x : *f) {
        x = random.Normal();
      }
    }
    std::vector<double> z = MovingAverage(cells, m, n, A, B, C);
    if (c == 0) {
      first = z;
    } else if (c == 1 && opt_.energy_correlation != 0.0) {
      double rho = opt_.energy_correlation;
      double s = std::sqrt(std::max(0.0, 1.0 - rho * rho));
      for (Index i = 0; i < n; ++i) {
        z[i] = rho * first[i] + s * z[i];
      }
    }
    for (double& x : z) {
      x = m.mean + m.sigma * x;
    }
    energies_.push_back(std::move(z));
  }
}

double StochasticNetwork::Connection(double r) const {
  const std::vector<double>& rs = opt_.connection_r;
  const std::vector<double>& ps = opt_.connection_p;
  if (r <= rs.front()) {
    return ps.front();
  }
  if (r >= rs.back()) {
    return 0.0;
  }
  std::size_t k = std::size_t(std::upper_bound(rs.begin(), rs.end(), r) -
                              rs.begin());  // rs[k-1] < r < rs[k]
  double t = (r - rs[k - 1]) / (rs[k] - rs[k - 1]);
  return ps[k - 1] + t * (ps[k] - ps[k - 1]);
}

namespace {
double NormalCdf(double g) { return 0.5 * std::erfc(-g / std::sqrt(2.0)); }
double NormalQuantile(double u) {
  return -std::sqrt(2.0) * boost::math::erfc_inv(2.0 * u);
}
}  // namespace

double StochasticNetwork::Residual(const CouplingModel& m, double g) {
  const std::vector<double>& q = m.residual_quantiles;
  const Index K = Index(q.size());
  if (K == 0) {
    return g;
  }
  if (K == 1) {
    return q[0] + g;
  }
  // q[k] belongs to the normal quantile g_k = Phi^-1((k + 1/2) / K);
  // linear interpolation in g, beyond the table like a normal distribution
  auto node = [K](Index k) {
    return NormalQuantile((double(k) + 0.5) / double(K));
  };
  double t = NormalCdf(g) * double(K) - 0.5;
  if (t <= 0.0) {
    return q[0] + (g - node(0));
  }
  if (t >= double(K - 1)) {
    return q[K - 1] + (g - node(K - 1));
  }
  Index k = std::min(Index(std::floor(t)), K - 2);
  double g0 = node(k);
  double g1 = node(k + 1);
  double f = std::min(1.0, std::max(0.0, (g - g0) / (g1 - g0)));
  return q[k] + f * (q[k + 1] - q[k]);
}

double StochasticNetwork::Log10J2(const CouplingModel& m, double r, double g) {
  double rnm = r * tools::conv::bohr2nm;
  double mean = m.mean_offset + m.mean_slope * rnm;
  double var = std::max(0.0, m.var_offset + m.var_slope * rnm);
  return mean + std::sqrt(var) * Residual(m, g);
}

void StochasticNetwork::InsertTraps(ModelRandom& random) {
  inserted_.clear();
  if (opt_.traps.empty() || opt_.trap_factor == 0.0) {
    return;
  }
  const Index n = NumberOfSites();
  double volume = opt_.box.prod();
  double radius =
      std::min(0.49 * opt_.box.minCoeff(), 2.0 * std::cbrt(volume / double(n)));
  PeriodicCells cells(opt_.box, radius, positions_);
  std::vector<char> used(n, 0);
  for (std::size_t t = 0; t < opt_.traps.size(); ++t) {
    const Trap& trap = opt_.traps[t];
    double expected = trap.density * opt_.trap_factor * double(n);
    Index count = Index(std::floor(expected));
    if (random.Uniform() <= expected - double(count)) {
      ++count;
    }
    const Index size = Index(trap.offsets.size());
    for (Index k = 0; k < count; ++k) {
      // a random free site whose nearest sites are free as well
      std::vector<Index> sites;
      for (Index attempt = 0; attempt < 1000 && sites.empty(); ++attempt) {
        Index centre = random.Integer(n);
        std::vector<Index> nb = cells.Nearest(centre, size);
        if (std::none_of(nb.begin(), nb.end(),
                         [&used](Index i) { return used[i] != 0; })) {
          sites = nb;
        }
      }
      if (sites.empty()) {
        throw std::runtime_error(
            "StochasticNetwork: no free place for a trap; the trap density "
            "is too high");
      }
      for (Index q = 0; q < size; ++q) {
        used[sites[q]] = 1;
        for (std::size_t c = 0; c < opt_.carriers.size(); ++c) {
          energies_[c][sites[q]] =
              opt_.carriers[c].energy.mean + trap.offsets[q][c];
        }
      }
      inserted_.push_back({Index(t), sites});
    }
  }
}

void StochasticNetwork::GeneratePairs(ModelRandom& random) {
  const Index n = NumberOfSites();
  const double rmax = opt_.connection_r.back();
  pairs_.clear();
  if (!(rmax > 0.0)) {
    return;
  }
  PeriodicCells cells(opt_.box, rmax, positions_);
  const double ev2 = tools::conv::ev2hrt * tools::conv::ev2hrt;
  const std::size_t ncar = opt_.carriers.size();
  const double rho = opt_.coupling_correlation;
  const double s = std::sqrt(std::max(0.0, 1.0 - rho * rho));
  struct Partner {
    Index j;
    Eigen::Vector3d R;
    double dist;
  };
  std::vector<Partner> partners;
  std::vector<double> g(ncar);
  for (Index i = 0; i < n; ++i) {
    partners.clear();
    cells.Within(positions_[i], rmax,
                 [&](Index j, const Eigen::Vector3d& conn, double dist) {
                   if (j > i) {
                     partners.push_back({j, conn, dist});
                   }
                 });
    std::sort(partners.begin(), partners.end(),
              [](const Partner& a, const Partner& b) { return a.j < b.j; });
    for (const Partner& p : partners) {
      double f = Connection(p.dist);
      if (f <= 0.0 || (f < 1.0 && random.Uniform() > f)) {
        continue;
      }
      for (std::size_t c = 0; c < ncar; ++c) {
        g[c] = random.Normal();
      }
      if (ncar > 1) {
        g[1] = rho * g[0] + s * g[1];
      }
      Pair pair{i, p.j, p.R, std::vector<double>(ncar)};
      for (std::size_t c = 0; c < ncar; ++c) {
        pair.J2[c] =
            std::pow(10.0, Log10J2(opt_.carriers[c].coupling, p.dist, g[c])) *
            ev2;
      }
      pairs_.push_back(std::move(pair));
    }
  }
}

void StochasticNetwork::BuildTopology(Topology& top) const {
  if (!top.Segments().empty() || top.NBList().size() > 0) {
    throw std::runtime_error("StochasticNetwork: topology must be empty");
  }
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() = opt_.box;
  top.setBox(box);
  top.setStep(0);
  top.setTime(0.0);
  for (Index i = 0; i < NumberOfSites(); ++i) {
    Segment& seg = top.AddSegment(opt_.segment_name);
    seg.push_back(Atom(0, "C", positions_[i]));
    for (std::size_t c = 0; c < opt_.carriers.size(); ++c) {
      const Carrier& car = opt_.carriers[c];
      seg.setEMpoles(car.state, energies_[c][i]);
      seg.setU_xX_nN(car.U_xX_nN, car.state);
      seg.setU_nX_nN(car.U_nX_nN, car.state);
      seg.setU_xN_xX(car.U_xN_xX, car.state);
    }
  }
  const std::vector<Segment>& segs = top.Segments();
  QMNBList& nblist = top.NBList();
  for (const Pair& p : pairs_) {
    QMPair& pair = nblist.Add(segs[p.i], segs[p.j], p.R);
    for (std::size_t c = 0; c < opt_.carriers.size(); ++c) {
      pair.setJeff2(p.J2[c], opt_.carriers[c].state);
    }
  }
}

// ------------------------------------------------------- parameter file

StochasticNetwork::Options StochasticNetwork::FromProperty(
    const tools::Property& p) {
  const double nm = tools::conv::nm2bohr;
  const double ev = tools::conv::ev2hrt;
  Options opt;
  opt.sites = p.ifExistsReturnElseReturnDefault<Index>("sites", 0);
  if (p.exists("box")) {
    std::vector<double> box = p.get("box").as<std::vector<double>>();
    if (box.size() != 3) {
      throw std::runtime_error("stochasticnetwork: box needs three lengths");
    }
    opt.box = Eigen::Vector3d(box[0], box[1], box[2]) * nm;
  }
  opt.seed = p.ifExistsReturnElseReturnDefault<std::uint64_t>("seed", 1);
  opt.segment_name =
      p.ifExistsReturnElseReturnDefault<std::string>("segmentname", "site");
  const tools::Property& pp = p.get("points");
  opt.points.hard_core = pp.get("hardcore").as<double>() * nm;
  opt.points.radius_mean = pp.get("radiusmean").as<double>() * nm;
  opt.points.radius_variance = pp.get("radiusvariance").as<double>() * nm * nm;
  opt.points.candidates =
      pp.ifExistsReturnElseReturnDefault<double>("candidates", 3.0);
  opt.points.max_iterations =
      pp.ifExistsReturnElseReturnDefault<Index>("maxiterations", 200);
  opt.connection_r = p.get("connection.r").as<std::vector<double>>();
  for (double& r : opt.connection_r) {
    r *= nm;
  }
  opt.connection_p = p.get("connection.p").as<std::vector<double>>();
  opt.energy_correlation =
      p.ifExistsReturnElseReturnDefault<double>("energycorrelation", 0.0);
  opt.coupling_correlation =
      p.ifExistsReturnElseReturnDefault<double>("couplingcorrelation", 0.0);
  for (const tools::Property* c : p.Select("carrier")) {
    Carrier car;
    car.state = QMStateType(c->get("state").as<std::string>());
    const tools::Property& e = c->get("energy");
    car.energy.mean = e.get("mean").as<double>() * ev;
    car.energy.sigma = e.get("sigma").as<double>() * ev;
    car.energy.w_a = e.get("w_a").as<double>();
    car.energy.w_b = e.get("w_b").as<double>();
    car.energy.l_b = e.get("l_b").as<Index>();
    car.energy.l_c = e.get("l_c").as<Index>();
    const tools::Property& j = c->get("coupling");
    car.coupling.mean_offset = j.get("meanoffset").as<double>();
    car.coupling.mean_slope = j.get("meanslope").as<double>();
    car.coupling.var_offset = j.get("varoffset").as<double>();
    car.coupling.var_slope = j.get("varslope").as<double>();
    if (j.exists("quantiles")) {
      car.coupling.residual_quantiles =
          j.get("quantiles").as<std::vector<double>>();
    }
    car.U_xX_nN =
        c->ifExistsReturnElseReturnDefault<double>("U_xX_nN", 0.0) * ev;
    car.U_nX_nN =
        c->ifExistsReturnElseReturnDefault<double>("U_nX_nN", 0.0) * ev;
    car.U_xN_xX =
        c->ifExistsReturnElseReturnDefault<double>("U_xN_xX", 0.0) * ev;
    opt.carriers.push_back(car);
  }
  for (const tools::Property* t : p.Select("trap")) {
    Trap trap;
    QMStateType state(t->get("carrier").as<std::string>());
    trap.carrier = -1;
    for (std::size_t c = 0; c < opt.carriers.size(); ++c) {
      if (opt.carriers[c].state == state) {
        trap.carrier = Index(c);
      }
    }
    if (trap.carrier < 0) {
      throw std::runtime_error("stochasticnetwork: trap for carrier " +
                               state.ToString() + " without a carrier block");
    }
    trap.density = t->get("density").as<double>();
    trap.release_rate =
        t->ifExistsReturnElseReturnDefault<double>("releaserate", 0.0);
    for (const tools::Property* site : t->Select("site")) {
      std::vector<double> o = site->as<std::vector<double>>();
      for (double& x : o) {
        x *= ev;
      }
      trap.offsets.push_back(o);
    }
    opt.traps.push_back(trap);
  }
  return opt;
}

std::string StochasticNetwork::ToXML(const Options& opt) {
  const double nm = tools::conv::bohr2nm;
  const double ev = tools::conv::hrt2ev;
  std::ostringstream out;
  out << std::setprecision(10);
  auto list = [&out](const std::vector<double>& v, double f) {
    for (std::size_t k = 0; k < v.size(); ++k) {
      out << (k ? " " : "") << v[k] * f;
    }
  };
  out << "<?xml version=\"1.0\"?>\n";
  out << "<!-- stochastic network model; lengths nm, energies eV, couplings "
         "log10(J^2/eV^2) with r in nm -->\n";
  out << "<stochasticnetwork>\n";
  out << "  <sites>" << opt.sites << "</sites>\n";
  out << "  <box>" << opt.box[0] * nm << " " << opt.box[1] * nm << " "
      << opt.box[2] * nm << "</box>\n";
  out << "  <seed>" << opt.seed << "</seed>\n";
  out << "  <segmentname>" << opt.segment_name << "</segmentname>\n";
  out << "  <points>\n";
  out << "    <hardcore>" << opt.points.hard_core * nm << "</hardcore>\n";
  out << "    <radiusmean>" << opt.points.radius_mean * nm << "</radiusmean>\n";
  out << "    <radiusvariance>" << opt.points.radius_variance * nm * nm
      << "</radiusvariance>\n";
  out << "    <candidates>" << opt.points.candidates << "</candidates>\n";
  out << "    <maxiterations>" << opt.points.max_iterations
      << "</maxiterations>\n";
  out << "  </points>\n";
  out << "  <connection>\n    <r>";
  list(opt.connection_r, nm);
  out << "</r>\n    <p>";
  list(opt.connection_p, 1.0);
  out << "</p>\n  </connection>\n";
  out << "  <energycorrelation>" << opt.energy_correlation
      << "</energycorrelation>\n";
  out << "  <couplingcorrelation>" << opt.coupling_correlation
      << "</couplingcorrelation>\n";
  for (const Carrier& c : opt.carriers) {
    out << "  <carrier>\n";
    out << "    <state>" << c.state.ToString() << "</state>\n";
    out << "    <energy>\n";
    out << "      <mean>" << c.energy.mean * ev << "</mean>\n";
    out << "      <sigma>" << c.energy.sigma * ev << "</sigma>\n";
    out << "      <w_a>" << c.energy.w_a << "</w_a>\n";
    out << "      <w_b>" << c.energy.w_b << "</w_b>\n";
    out << "      <l_b>" << c.energy.l_b << "</l_b>\n";
    out << "      <l_c>" << c.energy.l_c << "</l_c>\n";
    out << "    </energy>\n";
    out << "    <coupling>\n";
    out << "      <meanoffset>" << c.coupling.mean_offset << "</meanoffset>\n";
    out << "      <meanslope>" << c.coupling.mean_slope << "</meanslope>\n";
    out << "      <varoffset>" << c.coupling.var_offset << "</varoffset>\n";
    out << "      <varslope>" << c.coupling.var_slope << "</varslope>\n";
    if (!c.coupling.residual_quantiles.empty()) {
      out << "      <quantiles>";
      list(c.coupling.residual_quantiles, 1.0);
      out << "</quantiles>\n";
    }
    out << "    </coupling>\n";
    out << "    <U_xX_nN>" << c.U_xX_nN * ev << "</U_xX_nN>\n";
    out << "    <U_nX_nN>" << c.U_nX_nN * ev << "</U_nX_nN>\n";
    out << "    <U_xN_xX>" << c.U_xN_xX * ev << "</U_xN_xX>\n";
    out << "  </carrier>\n";
  }
  for (const Trap& t : opt.traps) {
    out << "  <!-- trap: site energies relative to the bulk mean (eV) for the "
           "carriers in the order of the carrier blocks, deepest site first "
           "-->\n";
    out << "  <trap>\n";
    out << "    <carrier>" << opt.carriers[t.carrier].state.ToString()
        << "</carrier>\n";
    out << "    <density>" << t.density << "</density>\n";
    out << "    <releaserate>" << t.release_rate << "</releaserate>\n";
    for (const auto& o : t.offsets) {
      out << "    <site>";
      list(o, ev);
      out << "</site>\n";
    }
    out << "  </trap>\n";
  }
  out << "</stochasticnetwork>\n";
  return out.str();
}

}  // namespace xtp
}  // namespace votca
