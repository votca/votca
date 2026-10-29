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
#include <complex>
#include <numeric>
#include <stdexcept>

// Third party includes
#include <Eigen/Eigenvalues>

// Local VOTCA includes
#include "votca/xtp/trapfinder.h"

namespace votca {
namespace xtp {

namespace {
/// eigenvector scaled so that its entry of largest modulus is 1 (real for
/// a real eigenvalue), returned as real numbers
Eigen::VectorXd Normalized(const Eigen::VectorXcd& v) {
  Index imax = 0;
  v.cwiseAbs().maxCoeff(&imax);
  Eigen::VectorXcd w = v / v(imax);
  return w.real();
}
}  // namespace

TrapFinder::Result TrapFinder::Find(const ToFGraph& graph,
                                    const Settings& settings) {
  if (!(settings.threshold > 0.0 && settings.threshold < 1.0) ||
      !(settings.gap >= 1.0) || settings.modes < 3 ||
      !(settings.mass > 0.0 && settings.mass <= 1.0) || settings.max_size < 1) {
    throw std::runtime_error(
        "TrapFinder: needs 0 < threshold < 1, gap >= 1, modes >= 3, "
        "0 < mass <= 1 and max_size >= 1");
  }
  // sites with hops
  const Index n = graph.NumberOfSites();
  std::vector<Index> local(n, -1);
  std::vector<Index> sites;
  for (Index i = 0; i < n; ++i) {
    if (!graph.getSite(i).ignored && graph.OutEdges(i).size() > 0) {
      local[i] = Index(sites.size());
      sites.push_back(i);
    }
  }
  const Index m = Index(sites.size());
  if (m < 3) {
    throw std::runtime_error("TrapFinder: needs at least three sites");
  }
  if (m > settings.max_sites) {
    throw std::runtime_error(
        "TrapFinder: " + std::to_string(m) +
        " sites; the dense eigenvalue problem is limited to " +
        std::to_string(settings.max_sites));
  }
  // the network must be strongly connected (one closed class): all sites
  // are reached from the first one along hops and backwards
  {
    std::vector<std::vector<Index>> fwd(m), bwd(m);
    for (Index a = 0; a < m; ++a) {
      for (const ToFGraph::Edge& e : graph.OutEdges(sites[a])) {
        Index b = local[e.to];
        if (b >= 0 && e.rate > 0.0) {
          fwd[a].push_back(b);
          bwd[b].push_back(a);
        }
      }
    }
    for (const auto* g : {&fwd, &bwd}) {
      std::vector<char> seen(m, 0);
      std::vector<Index> stack{0};
      seen[0] = 1;
      Index count = 1;
      while (!stack.empty()) {
        Index a = stack.back();
        stack.pop_back();
        for (Index b : (*g)[a]) {
          if (!seen[b]) {
            seen[b] = 1;
            ++count;
            stack.push_back(b);
          }
        }
      }
      if (count != m) {
        throw std::runtime_error(
            "TrapFinder: the hopping network is not strongly connected (" +
            std::to_string(m - count) +
            " sites unreachable); traps are not defined");
      }
    }
  }
  // transposed generator D - W, scaled by the largest escape rate
  double scale = 0.0;
  for (Index i : sites) {
    scale = std::max(scale, graph.EscapeRate(i));
  }
  if (!(scale > 0.0)) {
    throw std::runtime_error("TrapFinder: no hops");
  }
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  for (Index a = 0; a < m; ++a) {
    for (const ToFGraph::Edge& e : graph.OutEdges(sites[a])) {
      Index b = local[e.to];
      if (b < 0 || e.rate <= 0.0) {
        continue;
      }
      At(b, a) -= e.rate / scale;
      At(a, a) += e.rate / scale;
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> es(At, true);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error("TrapFinder: eigenvalue problem failed");
  }
  const Eigen::VectorXcd& lambda = es.eigenvalues();
  std::vector<Index> order(m);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&lambda](Index a, Index b) {
    return lambda(a).real() < lambda(b).real();
  });
  Result res;
  const Index modes = std::min(settings.modes, m);
  for (Index k = 0; k < modes; ++k) {
    res.eigenvalues.push_back(lambda(order[k]).real() * scale);
  }
  // stationary distribution: density mode of the zero eigenvalue
  Eigen::VectorXd pi = Normalized(es.eigenvectors().col(order[0]));
  pi = pi.cwiseMax(0.0);
  pi /= pi.sum();
  res.occupation.assign(n, 0.0);
  for (Index a = 0; a < m; ++a) {
    res.occupation[sites[a]] = pi(a);
  }

  for (Index i = 1; i + 1 < modes; ++i) {
    double r = res.eigenvalues[i + 1] / res.eigenvalues[i];
    if (r > res.gap) {
      res.gap = r;
      res.slow_modes = i;
    }
  }

  // candidates: connected components of the localized parts of the slow
  // density modes
  // undirected adjacency of the hopping graph
  std::vector<std::vector<Index>> adj(m);
  for (Index a = 0; a < m; ++a) {
    for (const ToFGraph::Edge& e : graph.OutEdges(sites[a])) {
      Index b = local[e.to];
      if (b >= 0 && e.rate > 0.0) {
        adj[a].push_back(b);
        adj[b].push_back(a);
      }
    }
  }
  struct Candidate {
    std::vector<Index> sites;  // local indices, sorted
    Index mode;
  };
  std::vector<Candidate> candidates;
  std::vector<char> in(m, 0);
  std::vector<char> seen(m, 0);
  for (Index k = 1; k <= res.slow_modes; ++k) {
    Eigen::VectorXd v = Normalized(es.eigenvectors().col(order[k]));
    // A density mode sums to zero: it moves weight between the trap and the
    // rest. The trap is the side with fewer sites; when the trap holds most
    // of the equilibrium weight, the other side has entries of the same
    // size on the sites that hold the rest of the weight.
    Index positive = 0;
    Index negative = 0;
    const double tiny = 1e-12 * v.cwiseAbs().maxCoeff();
    for (Index a = 0; a < m; ++a) {
      positive += v(a) > tiny;
      negative += v(a) < -tiny;
    }
    const double sign = positive <= negative ? 1.0 : -1.0;
    Eigen::VectorXd side = (sign * v).cwiseMax(0.0);
    side /= side.maxCoeff();
    const double side_mass = side.sum();
    for (Index a = 0; a < m; ++a) {
      in[a] = side(a) > settings.threshold;
      seen[a] = 0;
    }
    for (Index a = 0; a < m; ++a) {
      if (!in[a] || seen[a]) {
        continue;
      }
      std::vector<Index> comp{a};
      seen[a] = 1;
      for (std::size_t q = 0; q < comp.size(); ++q) {
        for (Index b : adj[comp[q]]) {
          if (in[b] && !seen[b]) {
            seen[b] = 1;
            comp.push_back(b);
          }
        }
      }
      // a trap is localized: few sites carrying most of the side's weight
      double mass = 0.0;
      for (Index b : comp) {
        mass += side(b);
      }
      if (Index(comp.size()) > settings.max_size ||
          mass < settings.mass * side_mass) {
        continue;
      }
      std::sort(comp.begin(), comp.end());
      candidates.push_back({std::move(comp), k});
    }
  }
  // drop duplicates and components contained in others
  std::sort(
      candidates.begin(), candidates.end(),
      [](const Candidate& a, const Candidate& b) {
        return a.sites.size() > b.sites.size() ||
               (a.sites.size() == b.sites.size() &&
                (a.sites < b.sites || (a.sites == b.sites && a.mode < b.mode)));
      });
  std::vector<Candidate> traps;
  for (const Candidate& c : candidates) {
    bool contained =
        std::any_of(traps.begin(), traps.end(), [&c](const Candidate& t) {
          return std::includes(t.sites.begin(), t.sites.end(), c.sites.begin(),
                               c.sites.end());
        });
    if (!contained) {
      traps.push_back(c);
    }
  }
  for (const Candidate& t : traps) {
    Trap trap;
    std::vector<char> member(m, 0);
    for (Index a : t.sites) {
      member[a] = 1;
    }
    double out = 0.0;
    for (Index a : t.sites) {
      trap.weight += pi(a);
      for (const ToFGraph::Edge& e : graph.OutEdges(sites[a])) {
        Index b = local[e.to];
        if (b >= 0 && !member[b]) {
          out += pi(a) * e.rate;
        }
      }
      trap.sites.push_back(sites[a]);
    }
    trap.hop_rate = trap.weight > 0.0 ? out / trap.weight : 0.0;
    // two-state picture of the mode: relaxation rate = release + capture,
    // weight = capture / (release + capture)
    trap.relaxation_rate = res.eigenvalues[t.mode];
    trap.release_rate = trap.relaxation_rate * (1.0 - trap.weight);
    trap.capture_rate = trap.relaxation_rate * trap.weight;
    // a trap holds a carrier much longer than the network relaxes otherwise
    if (trap.release_rate * settings.gap >
        res.eigenvalues[res.slow_modes + 1]) {
      continue;
    }
    std::sort(trap.sites.begin(), trap.sites.end(), [&graph](Index a, Index b) {
      return graph.getSite(a).energy < graph.getSite(b).energy;
    });
    res.traps.push_back(std::move(trap));
  }
  std::sort(res.traps.begin(), res.traps.end(),
            [](const Trap& a, const Trap& b) {
              return a.release_rate < b.release_rate;
            });
  return res;
}

}  // namespace xtp
}  // namespace votca
