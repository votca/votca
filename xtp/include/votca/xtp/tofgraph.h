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
#ifndef VOTCA_XTP_TOFGRAPH_H
#define VOTCA_XTP_TOFGRAPH_H

// Standard includes
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
 * \brief Hopping graph for time-of-flight (ToF) transport.
 *
 * Built from a periodic morphology (segments, box, neighbour list with
 * couplings). The periodicity is cut across one box vector (the transport
 * axis): pairs that connect two segments through a periodic image along that
 * vector are dropped, the other two directions stay periodic. Source and sink
 * sites are selected by slabs at the two faces normal to the transport axis
 * and/or by segment type patterns. Rates are Marcus rates from Rate_Engine at
 * the given temperature and field.
 *
 * The graph is the common input of the ToF kinetic Monte Carlo and of the
 * graph random walk (hitting time) calculations, so both see exactly the same
 * sites, edges and rates.
 *
 * With Options::bulk nothing is cut, the periodicity is kept in all three
 * directions and there are no sinks: the graph of the bulk KMC
 * (kmcmultiple). The ToF methods (GraphRandomWalk, ToFSimulation) refuse such
 * a graph; a ToF graph stays periodic in the two directions normal to the
 * transport axis. The
 * source sites are then the injection sites (source_pattern, all sites if
 * empty); normal, thickness, length and the coordinate x are zero. In this
 * mode a field makes U = E - q F.r multivalued (it changes by q F.(box
 * vector) around the torus), so the rates satisfy detailed balance only
 * without a field component along a box vector, and the potential is only
 * informative.
 *
 * Segments whose type matches one of Options::ignore_types take no part in
 * transport: they keep their index but have no edges and are neither
 * source nor sink.
 *
 * Coordinates along the transport axis are measured along the unit normal of
 * the cut face, from the face at fractional coordinate 0 (source side) to the
 * one at 1 (sink side). All quantities are in atomic units, rates in 1/s.
 */
class ToFGraph {
 public:
  struct Options {
    /// box vector across which the periodicity is cut (0: a, 1: b, 2: c)
    Index axis = 0;
    /// slab criteria: thickness of the source/sink slab along the normal
    /// (0: no slab criterion)
    double source_thickness = 0.0;
    double sink_thickness = 0.0;
    /// segment type wildcards (empty: no pattern criterion)
    /// if slab and pattern are both given, a site must satisfy both
    std::string source_pattern = "";
    std::string sink_pattern = "";
    /// bulk: no cut (periodic in all three directions), no sinks; slab
    /// thicknesses and the sink pattern must then be empty, axis is unused
    bool bulk = false;
    /// segment type wildcards of segments without edges
    std::vector<std::string> ignore_types;
  };

  /// bulk options: no cut, no sinks, injection on segments matching
  /// injection_pattern (empty or "*": all), ignore_types as for Options;
  /// the list is separated by commas or white space
  static Options BulkOptions(const std::string& injection_pattern,
                             const std::string& ignore_types);

  /// Options from the "tof" subpackage (subpackages/tof.xml), lengths in nm
  static Options OptionsFromProperty(const tools::Property& tof);

  enum class Role { Bulk, Source, Sink };

  struct Site {
    Index segment = -1;
    /// position with the fractional coordinates wrapped into [0,1)
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    /// coordinate along the transport normal, in [0, thickness)
    double x = 0.0;
    Role role = Role::Bulk;
    /// a sink can be reached from this site (always true for sinks)
    bool reaches_sink = false;
    /// site energy of the carrier
    double energy = 0.0;
    /// energy + field term, U = E - q F.r: hops release U_i - U_j
    double potential = 0.0;
    /// segment type in Options::ignore_types, no edges
    bool ignored = false;
  };

  struct Edge {
    Index to = -1;
    double rate = 0.0;
    Eigen::Vector3d dR = Eigen::Vector3d::Zero();  // pos(to) - pos(from)
    Index pair = -1;                               // id in the neighbour list
  };

  /// [begin, end) of the out-edges of one site
  class EdgeRange {
   public:
    EdgeRange(const Edge* b, const Edge* e) : begin_(b), end_(e) {}
    const Edge* begin() const { return begin_; }
    const Edge* end() const { return end_; }
    Index size() const { return Index(end_ - begin_); }

   private:
    const Edge* begin_;
    const Edge* end_;
  };

  /// temperature as kT in Hartree, field in Hartree/bohr (as Rate_Engine)
  ToFGraph(const Topology& top, const Options& opt, QMStateType carrier,
           double temperature, const Eigen::Vector3d& field);

  Index NumberOfSites() const { return Index(sites_.size()); }
  const std::vector<Site>& Sites() const { return sites_; }
  const Site& getSite(Index i) const { return sites_[i]; }

  EdgeRange OutEdges(Index i) const {
    return EdgeRange(edges_.data() + offsets_[i],
                     edges_.data() + offsets_[i + 1]);
  }
  Index NumberOfEdges() const { return Index(edges_.size()); }
  /// global index of the first out-edge of site i (edges are stored per
  /// site, OutEdges(i) covers [EdgeOffset(i), EdgeOffset(i+1)))
  Index EdgeOffset(Index i) const { return offsets_[i]; }
  const Edge& getEdge(Index e) const { return edges_[e]; }
  /// sum of the rates of all out-edges
  double EscapeRate(Index i) const { return escape_rates_[i]; }

  /// source sites from which a sink can be reached
  const std::vector<Index>& SourceSites() const { return sources_; }
  /// source sites from which no sink can be reached (excluded above)
  const std::vector<Index>& DisconnectedSources() const {
    return disconnected_sources_;
  }
  const std::vector<Index>& SinkSites() const { return sinks_; }

  Index NumberOfCutPairs() const { return cut_pairs_; }
  /// sites of ignored segment types
  Index NumberOfIgnoredSites() const { return ignored_sites_; }
  bool isBulk() const { return opt_.bulk; }

  /// unit normal of the cut face, pointing from source to sink side
  const Eigen::Vector3d& Normal() const { return normal_; }
  /// thickness of the box along the normal
  double Thickness() const { return thickness_; }
  /// transport length: distance along the normal from the mean position of
  /// the (connected) source sites to the closest sink site
  double Length() const { return length_; }

  /// true if the rates satisfy detailed balance with the potential U,
  /// pi_i w_ij = pi_j w_ji with pi_i ~ exp(-U_i/kT); requires equal
  /// reorganization energies for both directions of each pair and a field
  /// along the normal (a field component along a periodic direction is not a
  /// potential)
  bool isReversible() const { return reversible_; }
  /// largest |ln(w_ij/w_ji) - (U_i - U_j)/kT| over all pairs
  double DetailedBalanceDeviation() const { return db_deviation_; }

  QMStateType CarrierType() const { return carrier_; }
  double Temperature() const { return temperature_; }
  const Eigen::Vector3d& Field() const { return field_; }

  /// human readable summary
  std::string Summary() const;

 private:
  Options opt_;
  QMStateType carrier_;
  double temperature_;
  Eigen::Vector3d field_;

  std::vector<Site> sites_;
  std::vector<Edge> edges_;
  std::vector<Index> offsets_;
  std::vector<double> escape_rates_;
  std::vector<Index> sources_;
  std::vector<Index> disconnected_sources_;
  std::vector<Index> sinks_;
  Index cut_pairs_ = 0;
  Index ignored_sites_ = 0;
  Index sites_without_sink_ = 0;
  Eigen::Vector3d normal_ = Eigen::Vector3d::Zero();
  double thickness_ = 0.0;
  double length_ = 0.0;
  bool reversible_ = false;
  double db_deviation_ = 0.0;

  void SetupSites(const Topology& top);
  void SetupEdges(const Topology& top);
  void FindReachability();
  void SetupBulkSources();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_TOFGRAPH_H
