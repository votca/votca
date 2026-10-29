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
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

// Third party includes
#include <boost/math/special_functions/erf.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/stochasticnetworkfit.h"

namespace votca {
namespace xtp {

namespace {

double HalfBox(const Eigen::Vector3d& box) { return 0.4999 * box.minCoeff(); }

double Log10J2eV(double J2) {
  return std::log10(J2 * tools::conv::hrt2ev * tools::conv::hrt2ev);
}

double Mean(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0) / double(x.size());
}

double Variance(const std::vector<double>& x, double mean) {
  double v = 0.0;
  for (double y : x) {
    v += (y - mean) * (y - mean);
  }
  return v / double(x.size());
}

/// least squares y = a + b x
std::array<double, 2> Line(const std::vector<double>& x,
                           const std::vector<double>& y) {
  double mx = Mean(x);
  double my = Mean(y);
  double sxy = 0.0;
  double sxx = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    sxy += (x[k] - mx) * (y[k] - my);
    sxx += (x[k] - mx) * (x[k] - mx);
  }
  double b = sxx > 0.0 ? sxy / sxx : 0.0;
  return {my - b * mx, b};
}

StochasticNetworkFit::Curve Bins(double bin, double rmax) {
  StochasticNetworkFit::Curve c;
  Index n = Index(std::ceil(rmax / bin));
  for (Index k = 0; k < n; ++k) {
    c.r.push_back((double(k) + 0.5) * bin);
  }
  c.value.assign(n, 0.0);
  c.count.assign(n, 0.0);
  return c;
}

}  // namespace

// ------------------------------------------------------------------ data

StochasticNetworkFit::Data StochasticNetworkFit::FromTopology(
    const Topology& top, const std::vector<QMStateType>& carriers) {
  Data d;
  const Eigen::Matrix3d& box = top.getBox();
  Eigen::Matrix3d off = box;
  off.diagonal().setZero();
  if (off.cwiseAbs().maxCoeff() > 1e-9 * box.diagonal().maxCoeff()) {
    throw std::runtime_error(
        "networkfit: the box of the state file must be orthorhombic");
  }
  d.box = box.diagonal();
  d.carriers = carriers;
  const std::vector<Segment>& segs = top.Segments();
  const Index n = Index(segs.size());
  if (n < 2) {
    throw std::runtime_error("networkfit: need at least two segments");
  }
  d.positions.resize(n);
  d.energies.assign(carriers.size(), std::vector<double>(n));
  d.U_xX_nN.assign(carriers.size(), 0.0);
  d.U_nX_nN.assign(carriers.size(), 0.0);
  d.U_xN_xX.assign(carriers.size(), 0.0);
  for (Index i = 0; i < n; ++i) {
    const Segment& seg = segs[i];
    if (seg.getId() != i) {
      throw std::runtime_error(
          "networkfit: segment ids must be consecutive and start at 0");
    }
    Eigen::Vector3d p = seg.getPos();
    for (Index k = 0; k < 3; ++k) {
      p[k] -= d.box[k] * std::floor(p[k] / d.box[k]);
      if (p[k] >= d.box[k]) {
        p[k] = 0.0;
      }
    }
    d.positions[i] = p;
    for (std::size_t c = 0; c < carriers.size(); ++c) {
      d.energies[c][i] = seg.getSiteEnergy(carriers[c]);
      d.U_xX_nN[c] += seg.getU_xX_nN(carriers[c]) / double(n);
      d.U_nX_nN[c] += seg.getU_nX_nN(carriers[c]) / double(n);
      d.U_xN_xX[c] += seg.getU_xN_xX(carriers[c]) / double(n);
    }
  }
  for (const QMPair* pair : top.NBList()) {
    Data::Pair p;
    p.i = pair->Seg1()->getId();
    p.j = pair->Seg2()->getId();
    p.dist = pair->R().norm();
    for (QMStateType s : carriers) {
      p.J2.push_back(pair->getJeff2(s));
    }
    d.pairs.push_back(std::move(p));
  }
  if (d.pairs.empty()) {
    throw std::runtime_error("networkfit: the state file has no pairs");
  }
  return d;
}

StochasticNetworkFit::Data StochasticNetworkFit::FromNetwork(
    const StochasticNetwork& net, const std::vector<QMStateType>& carriers,
    const Eigen::Vector3d& box) {
  Data d;
  d.box = box;
  d.carriers = carriers;
  d.positions = net.Positions();
  for (std::size_t c = 0; c < carriers.size(); ++c) {
    d.energies.push_back(net.SiteEnergies(Index(c)));
  }
  d.U_xX_nN.assign(carriers.size(), 0.0);
  d.U_nX_nN.assign(carriers.size(), 0.0);
  d.U_xN_xX.assign(carriers.size(), 0.0);
  for (const StochasticNetwork::Pair& p : net.Pairs()) {
    d.pairs.push_back({p.i, p.j, p.R.norm(), p.J2});
  }
  return d;
}

// ------------------------------------------------------------ statistics

std::vector<double> StochasticNetworkFit::NearestDistances(
    const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box) {
  double radius = std::min(HalfBox(box),
                           2.0 * std::cbrt(box.prod() / double(points.size())));
  PeriodicCells cells(box, radius, points);
  std::vector<double> d(points.size());
  for (Index i = 0; i < Index(points.size()); ++i) {
    std::vector<Index> nb = cells.Nearest(i, 2);
    Index j = nb[0] == i ? nb[1] : nb[0];
    d[i] = cells.Connection(points[i], points[j]).norm();
  }
  return d;
}

std::vector<double> StochasticNetworkFit::ContactDistances(
    const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box,
    Index probes, std::uint64_t seed) {
  double radius = std::min(HalfBox(box),
                           2.0 * std::cbrt(box.prod() / double(points.size())));
  PeriodicCells cells(box, radius, points);
  ModelRandom random(seed, 7);
  std::vector<double> d(probes);
  for (Index k = 0; k < probes; ++k) {
    Eigen::Vector3d p;
    for (Index c = 0; c < 3; ++c) {
      p[c] = (1.0 - random.Uniform()) * box[c];
    }
    Index j = cells.NearestTo(p, 1)[0];
    d[k] = cells.Connection(p, points[j]).norm();
  }
  return d;
}

StochasticNetworkFit::Curve StochasticNetworkFit::Cdf(
    const std::vector<double>& x, const std::vector<double>& r) {
  std::vector<double> s = x;
  std::sort(s.begin(), s.end());
  Curve c;
  c.r = r;
  for (double y : r) {
    c.value.push_back(
        double(std::upper_bound(s.begin(), s.end(), y) - s.begin()) /
        double(s.size()));
  }
  return c;
}

StochasticNetworkFit::Curve StochasticNetworkFit::PairCorrelation(
    const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box,
    double bin, double rmax) {
  Curve c = Bins(bin, rmax);
  const Index nb = Index(c.r.size());
  PeriodicCells cells(box, rmax, points);
  for (Index i = 0; i < Index(points.size()); ++i) {
    cells.Within(points[i], rmax,
                 [&](Index j, const Eigen::Vector3d&, double dist) {
                   Index b = Index(dist / bin);
                   if (j > i && b < nb) {
                     c.count[b] += 1.0;
                   }
                 });
  }
  const double n = double(points.size());
  const double density = n / box.prod();
  for (Index b = 0; b < nb; ++b) {
    double r0 = double(b) * bin;
    double r1 = r0 + bin;
    double shell = 4.0 / 3.0 * tools::conv::Pi * (r1 * r1 * r1 - r0 * r0 * r0);
    c.value[b] = c.count[b] / (shell * density * 0.5 * n);
  }
  return c;
}

StochasticNetworkFit::Curve StochasticNetworkFit::Correlation(
    const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box,
    const std::vector<double>& z, double bin, double rmax,
    const std::vector<char>& skip) {
  Curve c = Bins(bin, rmax);
  const Index nb = Index(c.r.size());
  PeriodicCells cells(box, rmax, points);
  auto skipped = [&skip](Index i) { return !skip.empty() && skip[i] != 0; };
  for (Index i = 0; i < Index(points.size()); ++i) {
    if (skipped(i)) {
      continue;
    }
    cells.Within(points[i], rmax,
                 [&](Index j, const Eigen::Vector3d&, double dist) {
                   Index b = Index(dist / bin);
                   if (j > i && b < nb && !skipped(j)) {
                     c.count[b] += 1.0;
                     c.value[b] += z[i] * z[j];
                   }
                 });
  }
  for (Index b = 0; b < nb; ++b) {
    if (c.count[b] > 0.0) {
      c.value[b] /= c.count[b];
    }
  }
  return c;
}

std::vector<StochasticNetworkFit::Curve> StochasticNetworkFit::Overlaps(
    const std::vector<Eigen::Vector3d>& points, const Eigen::Vector3d& box,
    const std::vector<Index>& ls, double bin, double rmax) {
  const Index n = Index(points.size());
  const Index lmax = *std::max_element(ls.begin(), ls.end());
  if (lmax > n) {
    throw std::runtime_error("networkfit: l larger than the number of sites");
  }
  PeriodicCells cells(box, rmax, points);
  std::vector<Index> nearest(n * lmax);
#pragma omp parallel for schedule(dynamic, 64)
  for (Index i = 0; i < n; ++i) {
    std::vector<Index> nb = cells.Nearest(i, lmax);
    std::copy(nb.begin(), nb.end(), nearest.begin() + i * lmax);
  }
  Curve bins = Bins(bin, rmax);
  const Index nbins = Index(bins.r.size());
  std::vector<Curve> result(ls.size(), bins);
  std::vector<Index> rank(n, -1);
  std::vector<double> hist(lmax);
  for (Index i = 0; i < n; ++i) {
    for (Index k = 0; k < lmax; ++k) {
      rank[nearest[i * lmax + k]] = k;
    }
    cells.Within(points[i], rmax,
                 [&](Index j, const Eigen::Vector3d&, double dist) {
                   Index b = Index(dist / bin);
                   if (j <= i || b >= nbins) {
                     return;
                   }
                   // a common neighbour with ranks ri, rj belongs to both
                   // sets of size l if max(ri, rj) < l
                   std::fill(hist.begin(), hist.end(), 0.0);
                   for (Index k = 0; k < lmax; ++k) {
                     Index ri = rank[nearest[j * lmax + k]];
                     if (ri >= 0) {
                       hist[std::max(ri, k)] += 1.0;
                     }
                   }
                   std::partial_sum(hist.begin(), hist.end(), hist.begin());
                   for (std::size_t q = 0; q < ls.size(); ++q) {
                     Curve& c = result[q];
                     c.count[b] += 1.0;
                     c.value[b] += hist[ls[q] - 1] / double(ls[q]);
                   }
                 });
    for (Index k = 0; k < lmax; ++k) {
      rank[nearest[i * lmax + k]] = -1;
    }
  }
  for (Curve& c : result) {
    for (Index b = 0; b < nbins; ++b) {
      if (c.count[b] > 0.0) {
        c.value[b] /= c.count[b];
      }
    }
  }
  return result;
}

double StochasticNetworkFit::FitWeights(const Curve& kappa, const Curve& Ob,
                                        const Curve& Oc, double& w_b,
                                        double& w_c) {
  // normal equations of the weighted problem
  double bb = 0.0, bc = 0.0, cc = 0.0, bk = 0.0, ck = 0.0, kk = 0.0, n = 0.0;
  for (std::size_t q = 0; q < kappa.r.size(); ++q) {
    double w = kappa.count.empty() ? 1.0 : kappa.count[q];
    if (w <= 0.0) {
      continue;
    }
    double b = Ob.value[q];
    double c = Oc.value[q];
    double k = kappa.value[q];
    bb += w * b * b;
    bc += w * b * c;
    cc += w * c * c;
    bk += w * b * k;
    ck += w * c * k;
    kk += w * k * k;
    n += w;
  }
  if (n <= 0.0) {
    throw std::runtime_error("networkfit: no pairs for the correlation fit");
  }
  auto cost = [&](double x, double y) {
    return (x * x * bb + 2.0 * x * y * bc + y * y * cc - 2.0 * x * bk -
            2.0 * y * ck + kk) /
           n;
  };
  auto clamp01 = [](double x) { return std::min(1.0, std::max(0.0, x)); };
  std::vector<std::array<double, 2>> cand;
  cand.push_back({0.0, 0.0});
  double det = bb * cc - bc * bc;
  if (det > 1e-14 * bb * cc) {
    cand.push_back({(bk * cc - ck * bc) / det, (ck * bb - bk * bc) / det});
  }
  if (bb > 0.0) {
    cand.push_back({clamp01(bk / bb), 0.0});
  }
  if (cc > 0.0) {
    cand.push_back({0.0, clamp01(ck / cc)});
  }
  // on w_b + w_c = 1: w_b = t, w_c = 1 - t
  double a2 = bb - 2.0 * bc + cc;
  if (a2 > 0.0) {
    double t = clamp01((bk - ck - bc + cc) / a2);
    cand.push_back({t, 1.0 - t});
  }
  double best = std::numeric_limits<double>::infinity();
  for (const auto& x : cand) {
    if (x[0] < 0.0 || x[1] < 0.0 || x[0] + x[1] > 1.0 + 1e-12) {
      continue;
    }
    double f = cost(x[0], x[1]);
    if (f < best) {
      best = f;
      w_b = x[0];
      w_c = x[1];
    }
  }
  return std::max(0.0, best);
}

double StochasticNetworkFit::NormalQuantile(double u) {
  return -std::sqrt(2.0) * boost::math::erfc_inv(2.0 * u);
}

StochasticNetworkFit::Curve StochasticNetworkFit::PairsPerSite(const Data& data,
                                                               double bin,
                                                               double rmax) {
  Curve c = Bins(bin, rmax);
  for (const Data::Pair& p : data.pairs) {
    Index b = Index(p.dist / bin);
    if (b < Index(c.r.size())) {
      c.value[b] += 2.0 / double(data.positions.size());
      c.count[b] += 1.0;
    }
  }
  return c;
}

StochasticNetworkFit::Curve StochasticNetworkFit::MeanLog10J2(const Data& data,
                                                              Index carrier,
                                                              double bin,
                                                              double rmax) {
  Curve c = Bins(bin, rmax);
  for (const Data::Pair& p : data.pairs) {
    Index b = Index(p.dist / bin);
    if (b < Index(c.r.size()) && p.J2[carrier] > 0.0) {
      c.value[b] += Log10J2eV(p.J2[carrier]);
      c.count[b] += 1.0;
    }
  }
  for (std::size_t b = 0; b < c.r.size(); ++b) {
    if (c.count[b] > 0.0) {
      c.value[b] /= c.count[b];
    }
  }
  return c;
}

// ------------------------------------------------------------------- fit

StochasticNetworkFit::StochasticNetworkFit(const Data& data,
                                           const Settings& settings)
    : data_(data), settings_(settings) {
  if (settings_.carriers.empty()) {
    settings_.carriers = data_.carriers;
  }
  if (settings_.carriers != data_.carriers) {
    throw std::runtime_error(
        "networkfit: the data were read for other carriers");
  }
  if (!(settings_.bin > 0.0) || settings_.realizations <= 0) {
    throw std::runtime_error("networkfit: bin and realizations must be > 0");
  }
}

StochasticNetworkFit::PointFit StochasticNetworkFit::FitPoints() const {
  const Index n = Index(data_.positions.size());
  const Eigen::Vector3d& box = data_.box;
  const Index probes = std::max(Index(4000), n);
  const std::uint64_t probe_seed = 4711;
  PointFit fit;
  std::vector<double> nn = NearestDistances(data_.positions, box);
  std::vector<double> contact =
      ContactDistances(data_.positions, box, probes, probe_seed);
  const double step = 0.2 * settings_.bin;
  auto grid = [step](double rmax) {
    std::vector<double> r;
    for (double x = 0.0; x <= rmax; x += step) {
      r.push_back(x);
    }
    return r;
  };
  std::vector<double> rD = grid(1.5 * *std::max_element(nn.begin(), nn.end()));
  std::vector<double> rH =
      grid(1.5 * *std::max_element(contact.begin(), contact.end()));
  std::vector<double> sorted = nn;
  std::sort(sorted.begin(), sorted.end());
  const double q01 = sorted[std::size_t(0.01 * double(n - 1))];
  const double q50 = sorted[std::size_t(0.5 * double(n - 1))];
  const double rg = std::min(HalfBox(box), 4.0 * q50);
  fit.D = Cdf(nn, rD);
  fit.H = Cdf(contact, rH);
  fit.g = PairCorrelation(data_.positions, box, settings_.bin, rg);

  StochasticNetwork::PointProcess base;
  base.candidates = settings_.candidates;
  base.max_iterations = settings_.max_iterations;
  auto process = [&base](const std::array<double, 3>& x) {
    StochasticNetwork::PointProcess pp = base;
    pp.hard_core = std::abs(x[0]);
    pp.radius_mean = std::abs(x[1]);
    pp.radius_variance = pp.radius_mean > 0.0 ? x[2] * x[2] : 0.0;
    return pp;
  };
  auto msd = [](const Curve& a, const Curve& b) {
    double s = 0.0;
    for (std::size_t k = 0; k < a.value.size(); ++k) {
      s += (a.value[k] - b.value[k]) * (a.value[k] - b.value[k]);
    }
    return s / double(a.value.size());
  };
  auto evaluate = [&](const std::array<double, 3>& x, PointFit* store) {
    StochasticNetwork::PointProcess pp = process(x);
    Curve D = fit.D, H = fit.H, g = fit.g;
    for (Curve* c : {&D, &H, &g}) {
      std::fill(c->value.begin(), c->value.end(), 0.0);
    }
    const Index nr = settings_.realizations;
    std::vector<std::array<Curve, 3>> curves(nr);
    std::vector<char> saturated(nr, 0);
#pragma omp parallel for schedule(dynamic, 1)
    for (Index k = 0; k < nr; ++k) {
      std::vector<Eigen::Vector3d> pos;
      std::vector<double> radii;
      ModelRandom random(settings_.seed, 100 + std::uint64_t(k));
      try {
        StochasticNetwork::SamplePoints(pp, box, n, random, pos, radii);
      } catch (const std::runtime_error&) {
        saturated[k] = 1;  // radii too large
        continue;
      }
      curves[k] = {Cdf(NearestDistances(pos, box), rD),
                   Cdf(ContactDistances(pos, box, probes, probe_seed), rH),
                   PairCorrelation(pos, box, settings_.bin, rg)};
    }
    if (std::any_of(saturated.begin(), saturated.end(),
                    [](char s) { return s != 0; })) {
      return 1e3;
    }
    for (Index k = 0; k < nr; ++k) {
      const Curve& d = curves[k][0];
      const Curve& h = curves[k][1];
      const Curve& gg = curves[k][2];
      double f = 1.0 / double(nr);
      for (std::size_t q = 0; q < D.value.size(); ++q) {
        D.value[q] += f * d.value[q];
      }
      for (std::size_t q = 0; q < H.value.size(); ++q) {
        H.value[q] += f * h.value[q];
      }
      for (std::size_t q = 0; q < g.value.size(); ++q) {
        g.value[q] += f * gg.value[q];
      }
    }
    if (store != nullptr) {
      store->D_model = D;
      store->H_model = H;
      store->g_model = g;
    }
    return msd(D, fit.D) + msd(H, fit.H) + 0.25 * msd(g, fit.g);
  };

  // coarse grid, then Nelder-Mead
  std::array<double, 3> best{0.0, 0.0, 0.0};
  double fbest = std::numeric_limits<double>::infinity();
  Index evaluations = 0;
  for (double a : {0.6, 0.8, 1.0}) {
    for (double b : {0.25, 0.5, 0.75, 1.0}) {
      for (double c : {0.05, 0.15}) {
        double rh = a * q01;
        std::array<double, 3> x{rh, b * std::max(q50 - rh, 0.1 * q50), c * q50};
        double f = evaluate(x, nullptr);
        ++evaluations;
        if (f < fbest) {
          fbest = f;
          best = x;
        }
      }
    }
  }
  std::array<std::array<double, 3>, 4> s;
  std::array<double, 4> fs;
  s[0] = best;
  fs[0] = fbest;
  for (Index k = 0; k < 3; ++k) {
    s[k + 1] = best;
    s[k + 1][k] += 0.1 * q50;
    fs[k + 1] = evaluate(s[k + 1], nullptr);
    ++evaluations;
  }
  auto combine = [](const std::array<double, 3>& a,
                    const std::array<double, 3>& b, double t) {
    // a + t (b - a)
    return std::array<double, 3>{a[0] + t * (b[0] - a[0]),
                                 a[1] + t * (b[1] - a[1]),
                                 a[2] + t * (b[2] - a[2])};
  };
  while (evaluations < settings_.max_evaluations) {
    std::array<Index, 4> o{0, 1, 2, 3};
    std::sort(o.begin(), o.end(),
              [&fs](Index a, Index b) { return fs[a] < fs[b]; });
    double size = 0.0;
    for (Index k = 1; k < 4; ++k) {
      for (Index d = 0; d < 3; ++d) {
        size = std::max(size, std::abs(s[o[k]][d] - s[o[0]][d]));
      }
    }
    if (size < 1e-3 * q50) {
      break;
    }
    std::array<double, 3> centroid{0.0, 0.0, 0.0};
    for (Index k = 0; k < 3; ++k) {
      for (Index d = 0; d < 3; ++d) {
        centroid[d] += s[o[k]][d] / 3.0;
      }
    }
    const Index w = o[3];
    std::array<double, 3> xr = combine(centroid, s[w], -1.0);
    double fr = evaluate(xr, nullptr);
    ++evaluations;
    if (fr < fs[o[0]]) {
      std::array<double, 3> xe = combine(centroid, s[w], -2.0);
      double fe = evaluate(xe, nullptr);
      ++evaluations;
      if (fe < fr) {
        s[w] = xe;
        fs[w] = fe;
      } else {
        s[w] = xr;
        fs[w] = fr;
      }
    } else if (fr < fs[o[2]]) {
      s[w] = xr;
      fs[w] = fr;
    } else {
      bool outside = fr < fs[w];
      std::array<double, 3> xc = combine(centroid, outside ? xr : s[w], 0.5);
      double fc = evaluate(xc, nullptr);
      ++evaluations;
      if (fc < std::min(fr, fs[w])) {
        s[w] = xc;
        fs[w] = fc;
      } else {
        // shrink towards the best
        for (Index k = 1; k < 4; ++k) {
          s[o[k]] = combine(s[o[0]], s[o[k]], 0.5);
          fs[o[k]] = evaluate(s[o[k]], nullptr);
          ++evaluations;
        }
      }
    }
  }
  Index ib = Index(std::min_element(fs.begin(), fs.end()) - fs.begin());
  fit.process = process(s[ib]);
  fit.objective = evaluate(s[ib], &fit);
  fit.evaluations = evaluations + 1;
  return fit;
}

StochasticNetworkFit::Curve StochasticNetworkFit::FitConnection() const {
  const double bin = settings_.bin;
  double dmax = 0.0;
  for (const Data::Pair& p : data_.pairs) {
    dmax = std::max(dmax, p.dist);
  }
  const Index nb = Index(std::floor(dmax / bin)) + 1;
  const double rmax = double(nb) * bin;
  if (!(rmax < HalfBox(data_.box))) {
    throw std::runtime_error(
        "networkfit: the longest pair is longer than half the box; the "
        "connection probability cannot be sampled");
  }
  std::vector<double> pairs(nb, 0.0), all(nb, 0.0);
  PeriodicCells cells(data_.box, rmax, data_.positions);
  // distances of pairs and of all site pairs from the same positions
  for (const Data::Pair& p : data_.pairs) {
    double d =
        cells.Connection(data_.positions[p.i], data_.positions[p.j]).norm();
    pairs[std::min(Index(d / bin), nb - 1)] += 1.0;
  }
  for (Index i = 0; i < Index(data_.positions.size()); ++i) {
    cells.Within(data_.positions[i], rmax,
                 [&](Index j, const Eigen::Vector3d&, double dist) {
                   Index b = Index(dist / bin);
                   if (j > i && b < nb) {
                     all[b] += 1.0;
                   }
                 });
  }
  Curve c;
  Index first = 0;
  while (first < nb && all[first] == 0.0) {
    ++first;
  }
  for (Index b = 0; b < nb; ++b) {
    double f;
    if (b < first) {
      f = std::min(1.0, pairs[first] / all[first]);
    } else if (all[b] > 0.0) {
      f = std::min(1.0, pairs[b] / all[b]);
    } else {
      f = c.value.back();
    }
    c.r.push_back((double(b) + 0.5) * bin);
    c.value.push_back(f);
    c.count.push_back(all[b]);
  }
  c.r.push_back(rmax);
  c.value.push_back(0.0);
  c.count.push_back(0.0);
  return c;
}

StochasticNetworkFit::CouplingFit StochasticNetworkFit::FitCoupling(
    Index carrier, std::vector<double>& standardized) const {
  CouplingFit fit;
  const double nm = tools::conv::bohr2nm;
  std::vector<double> x, y;
  std::vector<Index> index;
  for (std::size_t k = 0; k < data_.pairs.size(); ++k) {
    const Data::Pair& p = data_.pairs[k];
    if (p.J2[carrier] > 0.0) {
      x.push_back(p.dist * nm);
      y.push_back(Log10J2eV(p.J2[carrier]));
      index.push_back(Index(k));
    } else {
      ++fit.excluded;
    }
  }
  fit.used = Index(x.size());
  standardized.assign(data_.pairs.size(),
                      std::numeric_limits<double>::quiet_NaN());
  if (fit.used < 3) {
    throw std::runtime_error("networkfit: carrier " +
                             data_.carriers[carrier].ToString() +
                             " has fewer than three nonzero couplings");
  }
  std::array<double, 2> m = Line(x, y);
  std::vector<double> e2(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) {
    double e = y[k] - m[0] - m[1] * x[k];
    e2[k] = e * e;
  }
  std::array<double, 2> v = Line(x, e2);
  double xmin = *std::min_element(x.begin(), x.end());
  double xmax = *std::max_element(x.begin(), x.end());
  if (!(v[0] + v[1] * xmin > 0.0) || !(v[0] + v[1] * xmax > 0.0)) {
    v = {Mean(e2), 0.0};  // constant variance
  }
  fit.model.mean_offset = m[0];
  fit.model.mean_slope = m[1];
  if (Mean(e2) < 1e-12) {
    // couplings are a function of the distance (model systems)
    fit.model.var_offset = 0.0;
    fit.model.var_slope = 0.0;
    return fit;
  }
  fit.model.var_offset = v[0];
  fit.model.var_slope = v[1];
  std::vector<double> u(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) {
    double var = std::max(v[0] + v[1] * x[k], 1e-300);
    u[k] = (y[k] - m[0] - m[1] * x[k]) / std::sqrt(var);
    standardized[index[k]] = u[k];
  }
  double mu = Mean(u);
  double var = Variance(u, mu);
  double m3 = 0.0, m4 = 0.0;
  for (double z : u) {
    m3 += std::pow(z - mu, 3);
    m4 += std::pow(z - mu, 4);
  }
  m3 /= double(u.size());
  m4 /= double(u.size());
  fit.skewness = var > 0.0 ? m3 / std::pow(var, 1.5) : 0.0;
  fit.excess_kurtosis = var > 0.0 ? m4 / (var * var) - 3.0 : 0.0;
  Index K = std::min(settings_.quantiles, fit.used / 10);
  if (K >= 2) {
    std::sort(u.begin(), u.end());
    for (Index k = 0; k < K; ++k) {
      double pos = (double(k) + 0.5) / double(K) * double(u.size()) - 0.5;
      pos = std::min(std::max(pos, 0.0), double(u.size() - 1));
      std::size_t lo = std::size_t(std::floor(pos));
      std::size_t hi = std::min(lo + 1, u.size() - 1);
      double t = pos - double(lo);
      fit.model.residual_quantiles.push_back(u[lo] + t * (u[hi] - u[lo]));
    }
  }
  return fit;
}

StochasticNetworkFit::Result StochasticNetworkFit::Fit() const {
  Result res;
  const Index n = Index(data_.positions.size());
  const std::size_t ncar = data_.carriers.size();
  StochasticNetwork::Options& opt = res.options;
  opt.box = data_.box;
  opt.sites = n;
  opt.seed = settings_.seed;

  if (settings_.fit_points) {
    res.points = FitPoints();
  } else {
    res.points.process = settings_.points;
  }
  opt.points = res.points.process;
  opt.points.candidates = settings_.candidates;
  opt.points.max_iterations = settings_.max_iterations;

  res.connection = FitConnection();
  opt.connection_r = res.connection.r;
  opt.connection_p = res.connection.value;

  // site energies
  double rc = settings_.correlation_cutoff > 0.0
                  ? std::min(settings_.correlation_cutoff, HalfBox(data_.box))
                  : HalfBox(data_.box);
  // trap sites are left out of the bulk statistics
  std::vector<char> trap_site(n, 0);
  for (const TrapInput& t : settings_.traps) {
    if (t.carrier < 0 || t.carrier >= Index(ncar) || t.sites.empty()) {
      throw std::runtime_error("networkfit: invalid trap input");
    }
    for (Index i : t.sites) {
      if (i < 0 || i >= n) {
        throw std::runtime_error("networkfit: trap site out of range");
      }
      trap_site[i] = 1;
    }
  }
  const Index nbulk =
      n - Index(std::count(trap_site.begin(), trap_site.end(), 1));
  if (nbulk < 2) {
    throw std::runtime_error("networkfit: no bulk sites left");
  }
  std::vector<std::vector<double>> bulk_energies(ncar);
  for (std::size_t c = 0; c < ncar; ++c) {
    for (Index i = 0; i < n; ++i) {
      if (!trap_site[i]) {
        bulk_energies[c].push_back(data_.energies[c][i]);
      }
    }
  }
  std::vector<std::vector<double>> z(ncar);
  bool any = false;
  for (std::size_t c = 0; c < ncar; ++c) {
    EnergyFit& e = res.energies.emplace_back();
    e.mean = Mean(bulk_energies[c]);
    e.sigma = std::sqrt(Variance(bulk_energies[c], e.mean));
    // below 1 micro-eV the energies are equal up to rounding
    e.fitted = e.sigma > 1e-6 * tools::conv::ev2hrt;
    if (!e.fitted) {
      e.sigma = 0.0;
    }
    if (e.fitted) {
      any = true;
      for (Index i = 0; i < n; ++i) {
        z[c].push_back(
            trap_site[i] ? 0.0 : (data_.energies[c][i] - e.mean) / e.sigma);
      }
    }
  }
  std::vector<Index> lb, lc, ls;
  for (Index l : settings_.l_b) {
    if (l >= 1 && l <= n) {
      lb.push_back(l);
    }
  }
  for (Index l : settings_.l_c) {
    if (l >= 1 && l <= n) {
      lc.push_back(l);
    }
  }
  if (any && (lb.empty() || lc.empty())) {
    throw std::runtime_error(
        "networkfit: the l_b and l_c grids have no values in [1, sites]");
  }
  std::vector<Curve> overlaps;
  if (any) {
    ls = lb;
    ls.insert(ls.end(), lc.begin(), lc.end());
    std::sort(ls.begin(), ls.end());
    ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
    overlaps = Overlaps(data_.positions, data_.box, ls, settings_.bin, rc);
  }
  auto overlap = [&](Index l) -> const Curve& {
    return overlaps[std::lower_bound(ls.begin(), ls.end(), l) - ls.begin()];
  };
  if (ncar > 1 && res.energies[0].fitted && res.energies[1].fitted) {
    double s = 0.0;
    for (Index i = 0; i < n; ++i) {
      s += z[0][i] * z[1][i];  // zero on trap sites
    }
    opt.energy_correlation = std::min(1.0, std::max(-1.0, s / double(nbulk)));
  }
  // The microscopic kappa is estimated with the sample mean and variance,
  // which biases it: for a moving average with mean correlation cbar over
  // all pairs (cbar = (w_a + w_b l_b + w_c l_c) / N), the estimate is about
  // (kappa - cbar) / (1 - cbar). The fit compares in this form, so that a
  // large model network has the unbiased kappa.
  std::vector<std::vector<double>> expected(ncar);  // unbiased model kappa
  double cbar0 = 0.0;
  for (std::size_t c = 0; c < ncar; ++c) {
    EnergyFit& e = res.energies[c];
    if (!e.fitted) {
      continue;
    }
    e.kappa = Correlation(data_.positions, data_.box, z[c], settings_.bin, rc,
                          trap_site);
    const Index nb = Index(e.kappa.value.size());
    const double rho = c == 1 ? opt.energy_correlation : 0.0;
    const double r2 = rho * rho;
    if (c == 1 && std::abs(rho) > 0.999) {
      const EnergyFit& f = res.energies[0];
      e.w_a = f.w_a;
      e.w_b = f.w_b;
      e.w_c = f.w_c;
      e.l_b = f.l_b;
      e.l_c = f.l_c;
      expected[c] = expected[0];
      e.kappa_model = f.kappa_model;
      e.cbar = f.cbar;
      e.rms = f.rms;
    } else {
      // the second carrier: kappa = rho^2 kappa_1 + (1 - rho^2) kappa_w
      std::vector<double> base(nb, 0.0);
      if (r2 > 0.0) {
        base = expected[0];
      }
      auto model = [&](Index b, Index cc, double wb, double wc,
                       std::vector<double>& m) {
        double wa = std::max(0.0, 1.0 - wb - wc);
        double cw = (wa + wb * double(b) + wc * double(cc)) / double(nbulk);
        m.resize(nb);
        for (Index q = 0; q < nb; ++q) {
          m[q] = r2 * base[q] + (1.0 - r2) * (wb * overlap(b).value[q] +
                                              wc * overlap(cc).value[q]);
        }
        return r2 * cbar0 + (1.0 - r2) * cw;  // cbar of the model
      };
      auto sample_residual = [&](const std::vector<double>& m, double cbar) {
        double s = 0.0, w = 0.0;
        for (Index q = 0; q < nb; ++q) {
          double d = (m[q] - cbar) / (1.0 - cbar) - e.kappa.value[q];
          s += e.kappa.count[q] * d * d;
          w += e.kappa.count[q];
        }
        return s / w;
      };
      double best = std::numeric_limits<double>::infinity();
      std::vector<double> m;
      Curve target = e.kappa;
      e.l_c_scan.r.assign(lc.begin(), lc.end());
      e.l_c_scan.value.assign(lc.size(),
                              std::numeric_limits<double>::infinity());
      for (Index b : lb) {
        for (std::size_t ic = 0; ic < lc.size(); ++ic) {
          const Index cc = lc[ic];
          if (cc <= b) {
            continue;
          }
          double wb = 0.0, wc = 0.0, cbar = 0.0;
          for (Index iter = 0; iter < 8; ++iter) {
            for (Index q = 0; q < nb; ++q) {
              double unbiased = cbar + (1.0 - cbar) * e.kappa.value[q];
              target.value[q] = (unbiased - r2 * base[q]) / (1.0 - r2);
            }
            FitWeights(target, overlap(b), overlap(cc), wb, wc);
            cbar = model(b, cc, wb, wc, m);
          }
          double f = sample_residual(m, cbar);
          e.l_c_scan.value[ic] = std::min(e.l_c_scan.value[ic], std::sqrt(f));
          if (f < best) {
            best = f;
            e.l_b = b;
            e.l_c = cc;
            e.w_b = wb;
            e.w_c = wc;
          }
        }
      }
      if (best == std::numeric_limits<double>::infinity()) {
        throw std::runtime_error(
            "networkfit: no l_c larger than l_b in the grids");
      }
      e.w_a = std::max(0.0, 1.0 - e.w_b - e.w_c);
      double cbar = model(e.l_b, e.l_c, e.w_b, e.w_c, expected[c]);
      if (c == 0) {
        cbar0 = cbar;
      }
      e.cbar = cbar;
      e.l_b_at_edge = lb.size() > 1 && e.l_b == lb.back();
      e.l_c_at_edge = lc.size() > 1 && e.l_c == lc.back();
      e.rms = std::sqrt(sample_residual(expected[c], cbar));
      e.kappa_model = e.kappa;
      for (Index q = 0; q < nb; ++q) {
        e.kappa_model.value[q] = (expected[c][q] - cbar) / (1.0 - cbar);
      }
    }
  }

  for (EnergyFit& e : res.energies) {
    e.sigma_model = e.fitted ? e.sigma / std::sqrt(1.0 - e.cbar) : 0.0;
  }

  // couplings
  std::vector<std::vector<double>> standardized(ncar);
  for (std::size_t c = 0; c < ncar; ++c) {
    res.couplings.push_back(FitCoupling(Index(c), standardized[c]));
  }
  if (ncar > 1) {
    std::vector<std::array<double, 2>> both;
    for (std::size_t k = 0; k < data_.pairs.size(); ++k) {
      if (std::isfinite(standardized[0][k]) &&
          std::isfinite(standardized[1][k])) {
        both.push_back({standardized[0][k], standardized[1][k]});
      }
    }
    if (both.size() > 2) {
      // normal scores of the ranks
      std::vector<std::vector<double>> score(2,
                                             std::vector<double>(both.size()));
      for (Index a = 0; a < 2; ++a) {
        std::vector<std::size_t> o(both.size());
        std::iota(o.begin(), o.end(), 0);
        std::sort(o.begin(), o.end(), [&](std::size_t p, std::size_t q) {
          return both[p][a] < both[q][a];
        });
        for (std::size_t r = 0; r < o.size(); ++r) {
          score[a][o[r]] =
              NormalQuantile((double(r) + 0.5) / double(both.size()));
        }
      }
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t k = 0; k < both.size(); ++k) {
        sxy += score[0][k] * score[1][k];
        sxx += score[0][k] * score[0][k];
        syy += score[1][k] * score[1][k];
      }
      opt.coupling_correlation =
          std::min(1.0, std::max(-1.0, sxy / std::sqrt(sxx * syy)));
    }
  }

  for (std::size_t c = 0; c < ncar; ++c) {
    StochasticNetwork::Carrier car;
    car.state = data_.carriers[c];
    const EnergyFit& e = res.energies[c];
    // the model energies are the site energies without the internal energy
    car.U_xX_nN = data_.U_xX_nN[c];
    car.U_nX_nN = data_.U_nX_nN[c];
    car.U_xN_xX = data_.U_xN_xX[c];
    car.energy.mean = e.mean - car.U_xX_nN;
    car.energy.sigma = e.sigma_model;
    car.energy.w_a = e.fitted ? e.w_a : 1.0;
    car.energy.w_b = e.fitted ? e.w_b : 0.0;
    car.energy.l_b = e.l_b;
    car.energy.l_c = e.l_c;
    car.coupling = res.couplings[c].model;
    opt.carriers.push_back(car);
  }
  // trap templates: energies relative to the bulk mean, all carriers
  for (const TrapInput& t : settings_.traps) {
    StochasticNetwork::Trap trap;
    trap.carrier = t.carrier;
    trap.density = 1.0 / double(n);
    trap.release_rate = t.release_rate;
    for (Index i : t.sites) {
      std::vector<double> o(ncar);
      for (std::size_t c = 0; c < ncar; ++c) {
        o[c] = data_.energies[c][i] - res.energies[c].mean;
      }
      trap.offsets.push_back(o);
    }
    opt.traps.push_back(trap);
  }
  return res;
}

}  // namespace xtp
}  // namespace votca
