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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE tofgraph_test

// Standard includes
#include <cmath>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/rate_engine.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/tofgraph.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
// 300 K as kT in Hartree
const double kT = votca::tools::conv::kB * 300.0 * ev;
// V/m -> Hartree/bohr
const double Vm = ev / (1e9 * nm);

LatticeModel::Options PaperLattice(LatticeModel::Disorder disorder) {
  LatticeModel::Options opt;
  opt.dims = {8, 8, 8};
  opt.spacing = 1.0 * nm;
  opt.cutoff = 1.5 * nm;
  opt.J2 = 1.79e-4 * ev * ev;
  opt.decay_length = 1.0 * nm;
  opt.reorganization = 0.23 * ev;
  opt.mean_energy = -0.76 * ev;
  opt.sigma = 0.19 * ev;
  opt.disorder = disorder;
  return opt;
}

ToFGraph::Options Slabs(Index axis = 0) {
  ToFGraph::Options opt;
  opt.axis = axis;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 1.0 * nm;
  return opt;
}

void AddSite(Topology& top, const std::string& type, const Eigen::Vector3d& pos,
             double energy) {
  Segment& seg = top.AddSegment(type);
  seg.push_back(Atom(0, "C", pos));
  for (QMStateType s :
       {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
    seg.setEMpoles(s, energy);
    seg.setU_nX_nN(0.1 * ev, s);
    seg.setU_xN_xX(0.1 * ev, s);
  }
}

void AddPair(Topology& top, Index i, Index j, const Eigen::Vector3d& R) {
  QMPair& pair = top.NBList().Add(top.Segments()[i], top.Segments()[j], R);
  for (QMStateType s :
       {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
    pair.setJeff2(1e-4 * ev * ev, s);
  }
}
}  // namespace

BOOST_AUTO_TEST_SUITE(tofgraph_test)

BOOST_AUTO_TEST_CASE(paper_lattice_cut_and_regions) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None));
  Topology top;
  model.BuildTopology(top);

  for (Index axis = 0; axis < 3; ++axis) {
    ToFGraph graph(top, Slabs(axis), QMStateType::Hole, kT,
                   Eigen::Vector3d::Zero());
    // each of the 64 sites in the last layer has 1 + 4 neighbours in the
    // first layer through the periodic boundary
    BOOST_CHECK_EQUAL(graph.NumberOfCutPairs(), 64 * 5);
    BOOST_CHECK_EQUAL(graph.NumberOfEdges(), 2 * (512 * 9 - 320));
    BOOST_CHECK_EQUAL(graph.SourceSites().size(), 64);
    BOOST_CHECK_EQUAL(graph.SinkSites().size(), 64);
    BOOST_CHECK(graph.DisconnectedSources().empty());
    BOOST_CHECK_CLOSE(graph.Thickness(), 8.0 * nm, 1e-10);
    // sources at 0.5 a, first sinks at 7.5 a: L = (N-1) a as in the paper
    BOOST_CHECK_CLOSE(graph.Length(), 7.0 * nm, 1e-10);
    BOOST_CHECK_CLOSE(graph.Normal()[axis], 1.0, 1e-12);
    for (Index s : graph.SourceSites()) {
      BOOST_CHECK_CLOSE(top.Segments()[s].getPos()[axis], 0.5 * nm, 1e-10);
    }
    for (Index s : graph.SinkSites()) {
      BOOST_CHECK_CLOSE(top.Segments()[s].getPos()[axis], 7.5 * nm, 1e-10);
    }
    // no edge crosses the cut
    for (Index i = 0; i < graph.NumberOfSites(); ++i) {
      for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
        double dx = graph.getSite(e.to).pos[axis] - graph.getSite(i).pos[axis];
        BOOST_CHECK_SMALL(dx - e.dR[axis], 1e-8);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(rates_match_rate_engine) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated));
  Topology top;
  model.BuildTopology(top);
  Eigen::Vector3d field(3e7 * Vm, 0.0, 0.0);
  ToFGraph graph(top, Slabs(), QMStateType::Electron, kT, field);
  Rate_Engine engine(kT, field);
  for (Index i = 0; i < graph.NumberOfSites(); ++i) {
    double escape = 0.0;
    for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
      const QMPair* pair = top.NBList()[e.pair];
      Rate_Engine::PairRates r = engine.Rate(*pair, QMStateType::Electron);
      double expected = (pair->Seg1()->getId() == i) ? r.rate12 : r.rate21;
      BOOST_CHECK_CLOSE(e.rate, expected, 1e-12);
      escape += e.rate;
    }
    BOOST_CHECK_CLOSE(graph.EscapeRate(i), escape, 1e-10);
  }
}

BOOST_AUTO_TEST_CASE(detailed_balance_flag) {
  for (auto disorder : {LatticeModel::Disorder::Uncorrelated,
                        LatticeModel::Disorder::Correlated}) {
    LatticeModel model(PaperLattice(disorder));
    Topology top;
    model.BuildTopology(top);
    for (QMStateType carrier :
         {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
      // field along the transport normal: a potential, rates are reversible
      ToFGraph along(top, Slabs(0), carrier, kT,
                     Eigen::Vector3d(5e7 * Vm, 0, 0));
      BOOST_CHECK(along.isReversible());
      BOOST_CHECK_SMALL(along.DetailedBalanceDeviation(), 1e-8);
      // field component along a periodic direction: not a potential
      ToFGraph across(top, Slabs(0), carrier, kT,
                      Eigen::Vector3d(5e7 * Vm, 1e7 * Vm, 0));
      BOOST_CHECK(!across.isReversible());
    }
  }
  // different reorganization energies for the two hop directions
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated));
  Topology top;
  model.BuildTopology(top);
  top.getSegment(100).setU_nX_nN(0.2 * ev, QMStateType::Hole);
  ToFGraph graph(top, Slabs(0), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  BOOST_CHECK(!graph.isReversible());
}

BOOST_AUTO_TEST_CASE(potential_sign) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None));
  Topology top;
  model.BuildTopology(top);
  Eigen::Vector3d field(1e7 * Vm, 0, 0);
  ToFGraph holes(top, Slabs(), QMStateType::Hole, kT, field);
  ToFGraph electrons(top, Slabs(), QMStateType::Electron, kT, field);
  // holes drift along the field (lower potential at larger x), electrons
  // against it
  Index src = holes.SourceSites()[0];
  Index snk = holes.SinkSites()[0];
  BOOST_CHECK_GT(holes.getSite(src).potential, holes.getSite(snk).potential);
  BOOST_CHECK_LT(electrons.getSite(src).potential,
                 electrons.getSite(snk).potential);
}

BOOST_AUTO_TEST_CASE(patterns_and_reachability) {
  // chain along x in a 6 x 4 x 4 nm box: "in" - b - b - b - "out", plus an
  // isolated "in" segment without pairs
  Topology top;
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << 6 * nm, 4 * nm, 4 * nm;
  top.setBox(box);
  AddSite(top, "in", Eigen::Vector3d(0.5, 2, 2) * nm, 0.0);
  AddSite(top, "b", Eigen::Vector3d(1.5, 2, 2) * nm, 0.0);
  AddSite(top, "b", Eigen::Vector3d(2.5, 2, 2) * nm, 0.0);
  AddSite(top, "b", Eigen::Vector3d(3.5, 2, 2) * nm, 0.0);
  AddSite(top, "out", Eigen::Vector3d(4.5, 2, 2) * nm, 0.0);
  AddSite(top, "in", Eigen::Vector3d(0.5, 0.5, 0.5) * nm, 0.0);
  for (Index i = 0; i < 4; ++i) {
    AddPair(top, i, i + 1, Eigen::Vector3d(1, 0, 0) * nm);
  }
  ToFGraph::Options opt;
  opt.source_pattern = "i*";
  opt.sink_pattern = "out";
  ToFGraph graph(top, opt, QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(graph.SourceSites().size(), 1);
  BOOST_CHECK_EQUAL(graph.SourceSites()[0], 0);
  BOOST_REQUIRE_EQUAL(graph.DisconnectedSources().size(), 1);
  BOOST_CHECK_EQUAL(graph.DisconnectedSources()[0], 5);
  BOOST_CHECK_EQUAL(graph.SinkSites().size(), 1);
  BOOST_CHECK(!graph.getSite(5).reaches_sink);
  BOOST_CHECK_CLOSE(graph.Length(), 4.0 * nm, 1e-10);

  // pattern and slab together: the sink slab of 1 nm contains no "out"
  opt.sink_thickness = 1.0 * nm;
  BOOST_CHECK_THROW(
      ToFGraph(top, opt, QMStateType::Hole, kT, Eigen::Vector3d::Zero()),
      std::runtime_error);
  opt.sink_thickness = 2.0 * nm;
  ToFGraph both(top, opt, QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(both.SinkSites().size(), 1);

  // no criteria
  ToFGraph::Options none;
  BOOST_CHECK_THROW(
      ToFGraph(top, none, QMStateType::Hole, kT, Eigen::Vector3d::Zero()),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(triclinic_box) {
  // box vectors a = (4,0,0), b = (1.5,5,0), c = (0,0,4) nm; cut across b.
  // Two chains along b at fractional a = 0.2 and 0.8, five sites each;
  // pairs along each chain (one through the b boundary, to be cut) and
  // between the chains through the a boundary (to be kept).
  Topology top;
  Eigen::Matrix3d box;
  box.col(0) = Eigen::Vector3d(4, 0, 0) * nm;
  box.col(1) = Eigen::Vector3d(1.5, 5, 0) * nm;
  box.col(2) = Eigen::Vector3d(0, 0, 4) * nm;
  top.setBox(box);
  for (double fa : {0.2, 0.8}) {
    for (Index k = 0; k < 5; ++k) {
      Eigen::Vector3d frac(fa, (double(k) + 0.5) / 5.0, 0.5);
      AddSite(top, "s", box * frac, 0.0);
    }
  }
  Eigen::Vector3d step_b = box.col(1) / 5.0;
  Eigen::Vector3d across_a = box * Eigen::Vector3d(0.4, 0, 0);  // 0.8 -> 1.2
  for (Index chain = 0; chain < 2; ++chain) {
    for (Index k = 0; k < 5; ++k) {
      AddPair(top, chain * 5 + k, chain * 5 + (k + 1) % 5, step_b);
    }
  }
  for (Index k = 0; k < 5; ++k) {
    AddPair(top, 5 + k, k, across_a);
  }
  ToFGraph::Options opt;
  opt.axis = 1;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 1.0 * nm;
  ToFGraph graph(top, opt, QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(graph.NumberOfCutPairs(), 2);
  BOOST_CHECK_EQUAL(graph.NumberOfEdges(), 2 * (15 - 2));
  BOOST_CHECK_CLOSE(graph.Thickness(), 5.0 * nm, 1e-10);
  BOOST_CHECK_SMALL((graph.Normal() - Eigen::Vector3d(0, 1, 0)).norm(), 1e-12);
  // slabs by the coordinate along the normal, not by Cartesian x
  BOOST_CHECK_EQUAL(graph.SourceSites().size(), 2);
  BOOST_CHECK_EQUAL(graph.SinkSites().size(), 2);
  BOOST_CHECK_CLOSE(graph.Length(), 4.0 * nm, 1e-10);
}

BOOST_AUTO_TEST_CASE(bulk_graph) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Hole, kT,
                 Eigen::Vector3d::Zero());
  BOOST_CHECK(graph.isBulk());
  BOOST_CHECK_EQUAL(graph.NumberOfCutPairs(), 0);
  BOOST_CHECK_EQUAL(graph.NumberOfEdges(), 2 * 512 * 9);
  BOOST_CHECK_EQUAL(graph.SourceSites().size(), 512);
  BOOST_CHECK(graph.SinkSites().empty());
  BOOST_CHECK(graph.isReversible());
  // a field makes the potential multivalued around the torus
  ToFGraph driven(top, ToFGraph::BulkOptions("*", ""), QMStateType::Hole, kT,
                  Eigen::Vector3d(1e7, 0, 0) * Vm);
  BOOST_CHECK(!driven.isReversible());
  // slabs or sinks are not allowed
  ToFGraph::Options bad = ToFGraph::BulkOptions("*", "");
  bad.sink_pattern = "site";
  BOOST_CHECK_THROW(
      ToFGraph(top, bad, QMStateType::Hole, kT, Eigen::Vector3d::Zero()),
      std::runtime_error);
  BOOST_CHECK_THROW(ToFGraph(top, ToFGraph::BulkOptions("nomatch", ""),
                             QMStateType::Hole, kT, Eigen::Vector3d::Zero()),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ignored_types) {
  Topology top;
  Eigen::Matrix3d box = Eigen::Matrix3d::Identity() * 4.0 * nm;
  top.setBox(box);
  AddSite(top, "A", Eigen::Vector3d(0.5, 0.5, 0.5) * nm, 0.0);
  AddSite(top, "B1", Eigen::Vector3d(1.5, 0.5, 0.5) * nm, 0.0);
  AddSite(top, "A", Eigen::Vector3d(2.5, 0.5, 0.5) * nm, 0.0);
  AddSite(top, "C", Eigen::Vector3d(3.5, 0.5, 0.5) * nm, 0.0);
  AddPair(top, 0, 1, Eigen::Vector3d(1, 0, 0) * nm);
  AddPair(top, 1, 2, Eigen::Vector3d(1, 0, 0) * nm);
  AddPair(top, 2, 3, Eigen::Vector3d(1, 0, 0) * nm);
  AddPair(top, 3, 0, Eigen::Vector3d(1, 0, 0) * nm);
  ToFGraph graph(top, ToFGraph::BulkOptions("A", "B*, C"), QMStateType::Hole,
                 kT, Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(graph.NumberOfIgnoredSites(), 2);
  BOOST_CHECK_EQUAL(graph.NumberOfEdges(), 0);
  BOOST_CHECK(graph.getSite(1).ignored && graph.getSite(3).ignored);
  BOOST_CHECK_EQUAL(graph.SourceSites().size(), 2);
  ToFGraph one(top, ToFGraph::BulkOptions("*", "C"), QMStateType::Hole, kT,
               Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(one.NumberOfEdges(), 4);
  BOOST_CHECK_EQUAL(one.SourceSites().size(), 3);
}

BOOST_AUTO_TEST_CASE(options_from_property) {
  votca::tools::Property tof;
  tof.add("axis", "b");
  votca::tools::Property& source = tof.add("source", "");
  source.add("thickness", "1.5");
  source.add("pattern", "in*");
  votca::tools::Property& sink = tof.add("sink", "");
  sink.add("thickness", "0");
  sink.add("pattern", "out");
  ToFGraph::Options opt = ToFGraph::OptionsFromProperty(tof);
  BOOST_CHECK_EQUAL(opt.axis, 1);
  BOOST_CHECK_CLOSE(opt.source_thickness, 1.5 * nm, 1e-12);
  BOOST_CHECK_EQUAL(opt.sink_thickness, 0.0);
  BOOST_CHECK_EQUAL(opt.source_pattern, "in*");
  BOOST_CHECK_EQUAL(opt.sink_pattern, "out");
  BOOST_CHECK(opt.ignore_types.empty());
  tof.add("ignoretypes", "X*, Y");
  opt = ToFGraph::OptionsFromProperty(tof);
  BOOST_REQUIRE_EQUAL(opt.ignore_types.size(), 2);
  BOOST_CHECK_EQUAL(opt.ignore_types[1], "Y");
  tof.set("axis", "x");
  BOOST_CHECK_THROW(ToFGraph::OptionsFromProperty(tof), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(inconsistent_pair_vector) {
  Topology top;
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << 6 * nm, 4 * nm, 4 * nm;
  top.setBox(box);
  AddSite(top, "s", Eigen::Vector3d(0.5, 2, 2) * nm, 0.0);
  AddSite(top, "s", Eigen::Vector3d(5.5, 2, 2) * nm, 0.0);
  AddPair(top, 0, 1, Eigen::Vector3d(1.3, 0, 0) * nm);
  BOOST_CHECK_THROW(
      ToFGraph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero()),
      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
