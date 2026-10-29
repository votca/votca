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
#include <queue>
#include <sstream>
#include <stdexcept>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/rate_engine.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/tofgraph.h"

namespace votca {
namespace xtp {

namespace {
double CarrierCharge(QMStateType type) {
  if (type == QMStateType::Electron) {
    return -1.0;
  } else if (type == QMStateType::Hole) {
    return 1.0;
  }
  return 0.0;
}

bool MatchesPattern(const std::string& pattern, const std::string& name) {
  return tools::wildcmp(pattern, name) != 0;
}

std::vector<std::string> SplitTypes(const std::string& list) {
  return tools::Tokenizer(list, ", \t\n").ToVector();
}
}  // namespace

ToFGraph::Options ToFGraph::BulkOptions(const std::string& injection_pattern,
                                        const std::string& ignore_types) {
  Options opt;
  opt.bulk = true;
  opt.source_thickness = 0.0;
  opt.sink_thickness = 0.0;
  opt.source_pattern = (injection_pattern == "*") ? "" : injection_pattern;
  opt.sink_pattern = "";
  opt.ignore_types = SplitTypes(ignore_types);
  return opt;
}

ToFGraph::Options ToFGraph::OptionsFromProperty(const tools::Property& tof) {
  Options opt;
  std::string axis = tof.get("axis").as<std::string>();
  if (axis == "a") {
    opt.axis = 0;
  } else if (axis == "b") {
    opt.axis = 1;
  } else if (axis == "c") {
    opt.axis = 2;
  } else {
    throw std::runtime_error("ToFGraph: axis must be a, b or c, not " + axis);
  }
  opt.source_thickness =
      tof.get("source.thickness").as<double>() * tools::conv::nm2bohr;
  opt.sink_thickness =
      tof.get("sink.thickness").as<double>() * tools::conv::nm2bohr;
  opt.source_pattern = tof.get("source.pattern").as<std::string>();
  opt.sink_pattern = tof.get("sink.pattern").as<std::string>();
  opt.ignore_types = SplitTypes(
      tof.ifExistsReturnElseReturnDefault<std::string>("ignoretypes", ""));
  return opt;
}

ToFGraph::ToFGraph(const Topology& top, const Options& opt, QMStateType carrier,
                   double temperature, const Eigen::Vector3d& field)
    : opt_(opt), carrier_(carrier), temperature_(temperature), field_(field) {
  if (!carrier_.isKMCState()) {
    throw std::runtime_error("ToFGraph: " + carrier_.ToLongString() +
                             " is not a hopping state");
  }
  if (temperature_ <= 0.0) {
    throw std::runtime_error("ToFGraph: the temperature must be positive");
  }
  if (opt_.axis < 0 || opt_.axis > 2) {
    throw std::runtime_error("ToFGraph: axis must be 0, 1 or 2");
  }
  bool source_criterion =
      opt_.source_thickness > 0.0 || !opt_.source_pattern.empty();
  bool sink_criterion = opt_.sink_thickness > 0.0 || !opt_.sink_pattern.empty();
  if (opt_.bulk) {
    if (opt_.source_thickness > 0.0 || sink_criterion) {
      throw std::runtime_error(
          "ToFGraph: a bulk graph has no slabs and no sinks");
    }
  } else if (!source_criterion || !sink_criterion) {
    throw std::runtime_error(
        "ToFGraph: source and sink each need a slab thickness and/or a "
        "segment pattern");
  }
  if (top.Segments().empty()) {
    throw std::runtime_error("ToFGraph: the topology contains no segments");
  }
  SetupSites(top);
  SetupEdges(top);
  if (opt_.bulk) {
    SetupBulkSources();
  } else {
    FindReachability();
  }
}

void ToFGraph::SetupSites(const Topology& top) {
  const Eigen::Matrix3d& box = top.getBox();  // columns: box vectors a, b, c
  double volume = std::abs(box.determinant());
  if (!(volume > 0.0)) {
    throw std::runtime_error(
        "ToFGraph: the box is degenerate, a periodic box is required");
  }
  if (!opt_.bulk) {
    const Eigen::Vector3d u = box.col((opt_.axis + 1) % 3);
    const Eigen::Vector3d v = box.col((opt_.axis + 2) % 3);
    normal_ = u.cross(v).normalized();
    if (normal_.dot(box.col(opt_.axis)) < 0.0) {
      normal_ = -normal_;
    }
    thickness_ = box.col(opt_.axis).dot(normal_);
  }
  if (!opt_.bulk && opt_.source_thickness + opt_.sink_thickness >= thickness_) {
    throw std::runtime_error(
        "ToFGraph: source and sink slabs together are thicker than the box");
  }

  const Eigen::Matrix3d inv = box.inverse();
  const double charge = CarrierCharge(carrier_);
  const std::vector<Segment>& segs = top.Segments();
  sites_.resize(segs.size());
  for (Index i = 0; i < Index(segs.size()); ++i) {
    const Segment& seg = segs[i];
    if (seg.getId() != i) {
      throw std::runtime_error(
          "ToFGraph: segment ids must be consecutive and start at 0");
    }
    Site& site = sites_[i];
    site.segment = i;
    Eigen::Vector3d frac = inv * seg.getPos();
    frac = frac.array() - frac.array().floor();
    // floor can return 1 for -1e-17
    for (Index d = 0; d < 3; ++d) {
      if (frac[d] >= 1.0) {
        frac[d] = 0.0;
      }
    }
    site.pos = box * frac;
    site.x = frac[opt_.axis] * thickness_;  // 0 for a bulk graph
    site.energy = seg.getSiteEnergy(carrier_);
    site.potential = site.energy - charge * field_.dot(site.pos);
    site.ignored =
        std::any_of(opt_.ignore_types.begin(), opt_.ignore_types.end(),
                    [&seg](const std::string& p) {
                      return MatchesPattern(p, seg.getType());
                    });
    if (site.ignored) {
      ++ignored_sites_;
      continue;  // role Bulk, no edges
    }

    bool source = true;
    bool sink = true;
    if (opt_.source_thickness > 0.0) {
      source = source && site.x < opt_.source_thickness;
    }
    if (!opt_.source_pattern.empty()) {
      source = source && MatchesPattern(opt_.source_pattern, seg.getType());
    }
    if (opt_.sink_thickness > 0.0) {
      sink = sink && site.x >= thickness_ - opt_.sink_thickness;
    }
    if (!opt_.sink_pattern.empty()) {
      sink = sink && MatchesPattern(opt_.sink_pattern, seg.getType());
    }
    sink = sink && !opt_.bulk;
    if (source && sink) {
      throw std::runtime_error("ToFGraph: segment " + std::to_string(i) +
                               " is selected as source and as sink");
    }
    if (source) {
      site.role = Role::Source;
    } else if (sink) {
      site.role = Role::Sink;
      sinks_.push_back(i);
    }
  }
  if (sinks_.empty() && !opt_.bulk) {
    throw std::runtime_error("ToFGraph: no sink sites selected");
  }
}

void ToFGraph::SetupEdges(const Topology& top) {
  const Eigen::Matrix3d& box = top.getBox();
  const Eigen::Matrix3d inv = box.inverse();
  Rate_Engine engine(temperature_, field_);
  const QMNBList& nblist = top.NBList();

  struct RawEdge {
    Index from;
    Edge edge;
  };
  std::vector<RawEdge> raw;
  raw.reserve(2 * nblist.size());
  db_deviation_ = 0.0;
  for (const QMPair* pair : nblist) {
    Index i = pair->Seg1()->getId();
    Index j = pair->Seg2()->getId();
    const Eigen::Vector3d& R = pair->R();
    // image shift between the pair vector and the wrapped positions
    Eigen::Vector3d shift = inv * (R - (sites_[j].pos - sites_[i].pos));
    Eigen::Vector3d rounded = shift.array().round();
    if ((shift - rounded).cwiseAbs().maxCoeff() > 1e-6) {
      throw std::runtime_error(
          "ToFGraph: the vector of pair " + std::to_string(pair->getId()) +
          " does not connect its segments or a periodic image of them");
    }
    if (sites_[i].ignored || sites_[j].ignored) {
      continue;
    }
    if (!opt_.bulk && rounded[opt_.axis] != 0.0) {
      ++cut_pairs_;
      continue;
    }
    Rate_Engine::PairRates rates = engine.Rate(*pair, carrier_);
    if (!std::isfinite(rates.rate12) || !std::isfinite(rates.rate21) ||
        rates.rate12 < 0.0 || rates.rate21 < 0.0) {
      throw std::runtime_error("ToFGraph: invalid rate for pair " +
                               std::to_string(pair->getId()));
    }
    if (rates.rate12 > 0.0 && rates.rate21 > 0.0) {
      double lhs = std::log(rates.rate12 / rates.rate21);
      double rhs = (sites_[i].potential - sites_[j].potential) / temperature_;
      db_deviation_ = std::max(db_deviation_, std::abs(lhs - rhs));
    } else if (rates.rate12 > 0.0 || rates.rate21 > 0.0) {
      db_deviation_ = std::numeric_limits<double>::infinity();
    }
    raw.push_back({i, Edge{j, rates.rate12, R, pair->getId()}});
    raw.push_back({j, Edge{i, rates.rate21, -R, pair->getId()}});
  }
  reversible_ = db_deviation_ < 1e-6;

  Index n = NumberOfSites();
  offsets_.assign(n + 1, 0);
  for (const RawEdge& e : raw) {
    ++offsets_[e.from + 1];
  }
  for (Index i = 0; i < n; ++i) {
    offsets_[i + 1] += offsets_[i];
  }
  edges_.resize(raw.size());
  std::vector<Index> fill(offsets_.begin(), offsets_.end() - 1);
  for (const RawEdge& e : raw) {
    edges_[fill[e.from]++] = e.edge;
  }
  escape_rates_.assign(n, 0.0);
  for (Index i = 0; i < n; ++i) {
    for (const Edge& e : OutEdges(i)) {
      escape_rates_[i] += e.rate;
    }
  }
}

void ToFGraph::FindReachability() {
  Index n = NumberOfSites();
  // incoming edges with positive rate
  std::vector<std::vector<Index>> incoming(n);
  for (Index i = 0; i < n; ++i) {
    for (const Edge& e : OutEdges(i)) {
      if (e.rate > 0.0) {
        incoming[e.to].push_back(i);
      }
    }
  }
  std::queue<Index> queue;
  for (Index s : sinks_) {
    sites_[s].reaches_sink = true;
    queue.push(s);
  }
  while (!queue.empty()) {
    Index j = queue.front();
    queue.pop();
    for (Index i : incoming[j]) {
      if (!sites_[i].reaches_sink) {
        sites_[i].reaches_sink = true;
        queue.push(i);
      }
    }
  }
  sites_without_sink_ = 0;
  double mean_x = 0.0;
  for (Index i = 0; i < n; ++i) {
    if (!sites_[i].reaches_sink) {
      ++sites_without_sink_;
    }
    if (sites_[i].role == Role::Source) {
      if (sites_[i].reaches_sink) {
        sources_.push_back(i);
        mean_x += sites_[i].x;
      } else {
        disconnected_sources_.push_back(i);
      }
    }
  }
  if (sources_.empty()) {
    throw std::runtime_error(
        "ToFGraph: no source site from which a sink can be reached (" +
        std::to_string(disconnected_sources_.size()) +
        " source sites without connection to a sink)");
  }
  mean_x /= double(sources_.size());
  double min_sink_x = std::accumulate(
      sinks_.begin(), sinks_.end(), std::numeric_limits<double>::max(),
      [this](double m, Index s) { return std::min(m, sites_[s].x); });
  length_ = min_sink_x - mean_x;
  if (length_ <= 0.0) {
    throw std::runtime_error(
        "ToFGraph: the sink sites are not further along the transport axis "
        "than the source sites");
  }
}

void ToFGraph::SetupBulkSources() {
  for (Index i = 0; i < NumberOfSites(); ++i) {
    if (sites_[i].role == Role::Source) {
      sources_.push_back(i);
    }
  }
  if (sources_.empty()) {
    throw std::runtime_error(
        "ToFGraph: no injection site matches the pattern '" +
        opt_.source_pattern + "'");
  }
  sites_without_sink_ = NumberOfSites();
}

std::string ToFGraph::Summary() const {
  std::stringstream s;
  double nm = tools::conv::bohr2nm;
  const char* axes[] = {"a", "b", "c"};
  if (opt_.bulk) {
    s << "Bulk hopping graph (periodic in a, b and c, no sinks) for "
      << carrier_.ToLongString() << "s\n";
    s << "  sites: " << NumberOfSites() << " (" << ignored_sites_
      << " of ignored types), hops: " << NumberOfEdges() << "\n";
    s << "  injection sites: " << sources_.size() << "\n";
    s << "  detailed balance: " << (reversible_ ? "yes" : "no")
      << " (max deviation of ln(w_ij/w_ji) " << db_deviation_ << ")";
    return s.str();
  }
  s << "ToF graph for " << carrier_.ToLongString() << "s\n";
  s << "  periodicity cut across box vector " << axes[opt_.axis] << ", normal ("
    << normal_.x() << " " << normal_.y() << " " << normal_.z()
    << "), thickness " << thickness_ * nm << " nm\n";
  s << "  sites: " << NumberOfSites() << " (" << ignored_sites_
    << " of ignored types), hops: " << NumberOfEdges()
    << ", pairs cut at the periodic boundary: " << cut_pairs_ << "\n";
  s << "  sources: " << sources_.size() << " connected, "
    << disconnected_sources_.size() << " without path to a sink (excluded)\n";
  s << "  sinks: " << sinks_.size()
    << ", sites without path to a sink: " << sites_without_sink_ << "\n";
  s << "  transport length (mean source to first sink): " << length_ * nm
    << " nm\n";
  s << "  detailed balance: " << (reversible_ ? "yes" : "no")
    << " (max deviation of ln(w_ij/w_ji) " << db_deviation_ << ")";
  return s.str();
}

}  // namespace xtp
}  // namespace votca
