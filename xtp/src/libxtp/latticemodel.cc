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
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/segment.h"

namespace votca {
namespace xtp {

double LatticeModel::NormalGenerator::Uniform() {
  // 53 random bits -> (0,1]
  return (double(engine_() >> 11) + 1.0) * (1.0 / 9007199254740992.0);
}

double LatticeModel::NormalGenerator::operator()() {
  if (has_spare_) {
    has_spare_ = false;
    return spare_;
  }
  double u1 = Uniform();
  double u2 = Uniform();
  double r = std::sqrt(-2.0 * std::log(u1));
  double phi = 2.0 * tools::conv::Pi * u2;
  spare_ = r * std::sin(phi);
  has_spare_ = true;
  return r * std::cos(phi);
}

LatticeModel::Disorder LatticeModel::DisorderFromString(
    const std::string& name) {
  if (name == "none") {
    return Disorder::None;
  } else if (name == "uncorrelated") {
    return Disorder::Uncorrelated;
  } else if (name == "correlated") {
    return Disorder::Correlated;
  }
  throw std::runtime_error("LatticeModel: unknown disorder type '" + name +
                           "', use none, uncorrelated or correlated");
}

LatticeModel::LatticeModel(const Options& opt) : opt_(opt) {
  if (std::any_of(opt_.dims.begin(), opt_.dims.end(),
                  [](Index d) { return d < 1; })) {
    throw std::runtime_error("LatticeModel: dimensions must be positive");
  }
  if (opt_.spacing <= 0.0 || opt_.cutoff <= 0.0 || opt_.decay_length <= 0.0) {
    throw std::runtime_error(
        "LatticeModel: spacing, cutoff and decay length must be positive");
  }
  if (opt_.reorganization <= 0.0) {
    throw std::runtime_error(
        "LatticeModel: the reorganization energy must be positive");
  }
  box_ = Eigen::Matrix3d::Zero();
  for (Index d = 0; d < 3; ++d) {
    box_(d, d) = double(opt_.dims[d]) * opt_.spacing;
  }
  // a pair must not see two images of the same site
  double half_box = 0.5 * box_.diagonal().minCoeff();
  if (opt_.cutoff >= half_box) {
    throw std::runtime_error(
        "LatticeModel: the cutoff must be smaller than half the box");
  }

  positions_.reserve(opt_.dims[0] * opt_.dims[1] * opt_.dims[2]);
  for (Index ix = 0; ix < opt_.dims[0]; ++ix) {
    for (Index iy = 0; iy < opt_.dims[1]; ++iy) {
      for (Index iz = 0; iz < opt_.dims[2]; ++iz) {
        positions_.push_back(opt_.spacing * Eigen::Vector3d(double(ix) + 0.5,
                                                            double(iy) + 0.5,
                                                            double(iz) + 0.5));
      }
    }
  }
  GenerateEnergies();
}

Eigen::Vector3d LatticeModel::ShortestConnection(
    const Eigen::Vector3d& r1, const Eigen::Vector3d& r2) const {
  const Eigen::Array3d box = box_.diagonal();
  const Eigen::Array3d r12 = r2 - r1;
  return (r12 - box * (r12 / box).round()).matrix();
}

Index LatticeModel::SiteIndex(Index ix, Index iy, Index iz) const {
  auto wrap = [](Index i, Index n) { return ((i % n) + n) % n; };
  return (wrap(ix, opt_.dims[0]) * opt_.dims[1] + wrap(iy, opt_.dims[1])) *
             opt_.dims[2] +
         wrap(iz, opt_.dims[2]);
}

std::vector<std::pair<double, std::array<Index, 3>>> LatticeModel::Offsets(
    double radius) const {
  // distances from site 0 to all others; site j has the coordinates of the
  // offset (sites are numbered ix, iy, iz with iz fastest)
  std::vector<std::pair<double, std::array<Index, 3>>> result;
  for (Index ix = 0; ix < opt_.dims[0]; ++ix) {
    for (Index iy = 0; iy < opt_.dims[1]; ++iy) {
      for (Index iz = 0; iz < opt_.dims[2]; ++iz) {
        Index j = SiteIndex(ix, iy, iz);
        if (j == 0) {
          continue;
        }
        double d = ShortestConnection(positions_[0], positions_[j]).norm();
        if (d <= radius) {
          result.push_back({d, {ix, iy, iz}});
        }
      }
    }
  }
  return result;
}

std::vector<std::vector<Index>> LatticeModel::NearestNeighbours(Index l) const {
  Index n = NumberOfSites();
  if (l > n - 1) {
    throw std::runtime_error(
        "LatticeModel: more neighbours requested for the correlated energies "
        "(" +
        std::to_string(l) + ") than other sites (" + std::to_string(n - 1) +
        ")");
  }
  if (l <= 0) {
    return std::vector<std::vector<Index>>(n);
  }
  // the l-th smallest distance is the same for every site; all offsets up
  // to it (with ties) are the candidates, sorted per site exactly as a
  // search over all sites would: distances compared with a tolerance, ties
  // broken by site index, so the result does not depend on rounding
  double tol = 1e-9 * opt_.spacing;
  std::vector<std::pair<double, std::array<Index, 3>>> all =
      Offsets(std::numeric_limits<double>::max());
  std::vector<double> distances(all.size());
  std::transform(all.begin(), all.end(), distances.begin(),
                 [](const auto& o) { return o.first; });
  std::nth_element(distances.begin(), distances.begin() + (l - 1),
                   distances.end());
  std::vector<std::pair<double, std::array<Index, 3>>> offsets =
      Offsets(distances[l - 1] + tol);

  auto less = [tol](const std::pair<double, Index>& a,
                    const std::pair<double, Index>& b) {
    if (std::abs(a.first - b.first) > tol) {
      return a.first < b.first;
    }
    return a.second < b.second;
  };
  std::vector<std::vector<Index>> result(n);
  std::vector<std::pair<double, Index>> dist(offsets.size());
  for (Index ix = 0; ix < opt_.dims[0]; ++ix) {
    for (Index iy = 0; iy < opt_.dims[1]; ++iy) {
      for (Index iz = 0; iz < opt_.dims[2]; ++iz) {
        Index i = SiteIndex(ix, iy, iz);
        for (Index k = 0; k < Index(offsets.size()); ++k) {
          const std::array<Index, 3>& o = offsets[k].second;
          Index j = SiteIndex(ix + o[0], iy + o[1], iz + o[2]);
          dist[k] = {ShortestConnection(positions_[i], positions_[j]).norm(),
                     j};
        }
        std::partial_sort(dist.begin(), dist.begin() + l, dist.end(), less);
        result[i].reserve(l);
        for (Index m = 0; m < l; ++m) {
          result[i].push_back(dist[m].second);
        }
      }
    }
  }
  return result;
}

void LatticeModel::GenerateEnergies() {
  Index n = NumberOfSites();
  energies_.assign(n, opt_.mean_energy);
  if (opt_.disorder == Disorder::None || opt_.sigma == 0.0) {
    return;
  }
  NormalGenerator normal(opt_.seed);
  if (opt_.disorder == Disorder::Uncorrelated) {
    // std::generate evaluates the generator in order
    std::vector<double> noise(energies_.size());
    std::generate(noise.begin(), noise.end(), [&normal]() { return normal(); });
    std::transform(energies_.begin(), energies_.end(), noise.begin(),
                   energies_.begin(),
                   [this](double e, double z) { return e + opt_.sigma * z; });
    return;
  }

  // Correlated: eps_i = sqrt(eta_a) M^a_i + sqrt(eta_b/l_b) sum_{j in N_b(i)}
  // M^b_j + sqrt((1-eta_a-eta_b)/l_c) sum_{j in N_c(i)} M^c_j + mean, with
  // M ~ N(0, sigma^2) and N_l(i) the l sites nearest to i, i itself included
  // (include_self) or not, so that the variance is sigma^2.
  double eta_c = 1.0 - opt_.eta_a - opt_.eta_b;
  if (opt_.eta_a < 0.0 || opt_.eta_b < 0.0 || eta_c < -1e-12) {
    throw std::runtime_error(
        "LatticeModel: weights eta_a, eta_b must be non-negative with "
        "eta_a + eta_b <= 1");
  }
  eta_c = std::max(eta_c, 0.0);
  if (opt_.l_b < 1 || opt_.l_c < 1) {
    throw std::runtime_error(
        "LatticeModel: the neighbour numbers l_b and l_c must be positive");
  }
  std::vector<double> Ma(n), Mb(n), Mc(n);
  for (Index i = 0; i < n; ++i) {
    Ma[i] = opt_.sigma * normal();
  }
  for (Index i = 0; i < n; ++i) {
    Mb[i] = opt_.sigma * normal();
  }
  for (Index i = 0; i < n; ++i) {
    Mc[i] = opt_.sigma * normal();
  }
  // with the site itself, only l - 1 other sites are needed
  const Index self = opt_.include_self ? 1 : 0;
  std::vector<std::vector<Index>> neighbours =
      NearestNeighbours(std::max(opt_.l_b, opt_.l_c) - self);
  double wa = std::sqrt(opt_.eta_a);
  double wb = std::sqrt(opt_.eta_b / double(opt_.l_b));
  double wc = std::sqrt(eta_c / double(opt_.l_c));
  for (Index i = 0; i < n; ++i) {
    double sb = opt_.include_self ? Mb[i] : 0.0;
    for (Index m = 0; m < opt_.l_b - self; ++m) {
      sb += Mb[neighbours[i][m]];
    }
    double sc = opt_.include_self ? Mc[i] : 0.0;
    for (Index m = 0; m < opt_.l_c - self; ++m) {
      sc += Mc[neighbours[i][m]];
    }
    energies_[i] += wa * Ma[i] + wb * sb + wc * sc;
  }
}

void LatticeModel::BuildTopology(Topology& top) const {
  if (!top.Segments().empty() || top.NBList().size() > 0) {
    throw std::runtime_error("LatticeModel: topology must be empty");
  }
  top.setBox(box_);
  top.setStep(0);
  top.setTime(0.0);
  Index n = NumberOfSites();
  // half of the reorganization energy on each segment, so that
  // reorg12 = U_nX_nN(1) + U_xN_xX(2) = reorg21 = reorganization
  double half_lambda = 0.5 * opt_.reorganization;
  for (Index i = 0; i < n; ++i) {
    Segment& seg = top.AddSegment(opt_.segment_name);
    seg.push_back(Atom(0, "C", positions_[i]));
    for (QMStateType s : opt_.states) {
      seg.setEMpoles(s, energies_[i]);
      seg.setU_xX_nN(0.0, s);
      seg.setU_nX_nN(half_lambda, s);
      seg.setU_xN_xX(half_lambda, s);
    }
  }
  const std::vector<Segment>& segs = top.Segments();
  QMNBList& nblist = top.NBList();
  // candidates from the lattice offsets within the cutoff (O(n)); pairs are
  // added in the order of a search over all i < j, and the cutoff is tested
  // on the actual pair vector as before
  std::vector<std::pair<double, std::array<Index, 3>>> offsets =
      Offsets(opt_.cutoff + 1e-9 * opt_.spacing);
  std::vector<Index> partners;
  for (Index i = 0; i < n; ++i) {
    Index ix = i / (opt_.dims[1] * opt_.dims[2]);
    Index iy = (i / opt_.dims[2]) % opt_.dims[1];
    Index iz = i % opt_.dims[2];
    partners.clear();
    for (const auto& o : offsets) {
      Index j = SiteIndex(ix + o.second[0], iy + o.second[1], iz + o.second[2]);
      if (j > i) {
        partners.push_back(j);
      }
    }
    std::sort(partners.begin(), partners.end());
    partners.erase(std::unique(partners.begin(), partners.end()),
                   partners.end());
    for (Index j : partners) {
      Eigen::Vector3d r = ShortestConnection(positions_[i], positions_[j]);
      double dist = r.norm();
      if (dist >= opt_.cutoff) {
        continue;
      }
      QMPair& pair = nblist.Add(segs[i], segs[j], r);
      double J2 =
          opt_.J2 * std::exp(-(dist - opt_.spacing) / opt_.decay_length);
      for (QMStateType s : opt_.states) {
        pair.setJeff2(J2, s);
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca
