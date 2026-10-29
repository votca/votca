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

#define BOOST_TEST_MODULE graphrandomwalk_test

// Standard includes
#include <cmath>
#include <map>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/segment.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT = votca::tools::conv::kB * 300.0 * ev;
const double Vm = ev / (1e9 * nm);

LatticeModel::Options PaperLattice(LatticeModel::Disorder disorder,
                                   Index n = 8) {
  LatticeModel::Options opt;
  opt.dims = {n, n, n};
  opt.spacing = 1.0 * nm;
  opt.cutoff = 1.5 * nm;
  opt.J2 = 1.79e-4 * ev * ev;
  opt.decay_length = 1.0 * nm;
  opt.reorganization = 0.23 * ev;
  opt.mean_energy = -0.76 * ev;
  opt.sigma = 0.19 * ev;
  opt.disorder = disorder;
  opt.seed = 5;
  return opt;
}

ToFGraph::Options Slabs() {
  ToFGraph::Options opt;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 1.0 * nm;
  return opt;
}

// chain of n sites along x with random energies, source = site 0, sink =
// site n-1
void BuildChain(Topology& top, Index n, double sigma) {
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << double(n + 2) * nm, 5 * nm, 5 * nm;
  top.setBox(box);
  LatticeModel::NormalGenerator normal(17);
  for (Index i = 0; i < n; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(
        Atom(0, "C", Eigen::Vector3d(0.5 + double(i), 2.5, 2.5) * nm));
    seg.setEMpoles(QMStateType::Hole, sigma * normal());
    seg.setU_nX_nN(0.1 * ev, QMStateType::Hole);
    seg.setU_xN_xX(0.1 * ev, QMStateType::Hole);
  }
  for (Index i = 0; i + 1 < n; ++i) {
    QMPair& pair = top.NBList().Add(top.Segments()[i], top.Segments()[i + 1],
                                    Eigen::Vector3d(1, 0, 0) * nm);
    pair.setJeff2(1e-4 * ev * ev, QMStateType::Hole);
  }
}

double Rate(const ToFGraph& g, Index from, Index to) {
  for (const ToFGraph::Edge& e : g.OutEdges(from)) {
    if (e.to == to) {
      return e.rate;
    }
  }
  return 0.0;
}

// dense reference: solve Q tau = 1 on the transient sites
Eigen::VectorXd DenseHittingTimes(const ToFGraph& g) {
  Index n = g.NumberOfSites();
  std::vector<Index> idx(n, -1);
  Index m = 0;
  for (Index i = 0; i < n; ++i) {
    if (g.getSite(i).role != ToFGraph::Role::Sink) {
      idx[i] = m++;
    }
  }
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(m, m);
  for (Index i = 0; i < n; ++i) {
    if (idx[i] < 0) {
      continue;
    }
    Q(idx[i], idx[i]) = g.EscapeRate(i);
    for (const ToFGraph::Edge& e : g.OutEdges(i)) {
      if (idx[e.to] >= 0) {
        Q(idx[i], idx[e.to]) -= e.rate;
      }
    }
  }
  // with iterative refinement: Q is badly conditioned with disorder
  Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
  Eigen::VectorXd one = Eigen::VectorXd::Ones(m);
  Eigen::VectorXd t = lu.solve(one);
  for (Index it = 0; it < 3; ++it) {
    t += lu.solve(one - Q * t);
  }
  Eigen::VectorXd tau = Eigen::VectorXd::Zero(n);
  for (Index i = 0; i < n; ++i) {
    if (idx[i] >= 0) {
      tau[i] = t[idx[i]];
    }
  }
  return tau;
}
// dense two-carrier reference, written independently of the pair ranking:
// states are sorted site pairs (and, for transit, single sites) of non-sink
// sites in a map; returns the mean time for every pair state
std::map<std::pair<Index, Index>, double> DenseTwoCarriers(const ToFGraph& g,
                                                           bool transit) {
  Index n = g.NumberOfSites();
  auto sink = [&g](Index i) {
    return g.getSite(i).role == ToFGraph::Role::Sink;
  };
  std::map<std::pair<Index, Index>, Index> pair_index;
  std::map<Index, Index> single_index;
  Index count = 0;
  for (Index a = 0; a < n; ++a) {
    for (Index b = a + 1; b < n; ++b) {
      if (!sink(a) && !sink(b)) {
        pair_index[{a, b}] = count++;
      }
    }
  }
  if (transit) {
    for (Index a = 0; a < n; ++a) {
      if (!sink(a)) {
        single_index[a] = count++;
      }
    }
  }
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(count, count);
  Eigen::VectorXd rhs = Eigen::VectorXd::Ones(count);
  for (const auto& [ab, s] : pair_index) {
    for (int c = 0; c < 2; ++c) {
      Index x = (c == 0) ? ab.first : ab.second;
      Index y = (c == 0) ? ab.second : ab.first;
      for (const ToFGraph::Edge& e : g.OutEdges(x)) {
        if (e.to == y) {
          continue;
        }
        Q(s, s) += e.rate;
        if (sink(e.to)) {
          if (transit) {
            Q(s, single_index[y]) -= e.rate;
          }
        } else {
          Q(s, pair_index[{std::min(e.to, y), std::max(e.to, y)}]) -= e.rate;
        }
      }
    }
  }
  for (const auto& [a, s] : single_index) {
    rhs[s] = 0.5;
    for (const ToFGraph::Edge& e : g.OutEdges(a)) {
      Q(s, s) += e.rate;
      if (!sink(e.to)) {
        Q(s, single_index[e.to]) -= e.rate;
      }
    }
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu(Q);
  Eigen::VectorXd t = lu.solve(rhs);
  for (Index it = 0; it < 3; ++it) {
    t += lu.solve(rhs - Q * t);
  }
  std::map<std::pair<Index, Index>, double> result;
  for (const auto& [ab, s] : pair_index) {
    result[ab] = t[s];
  }
  return result;
}

// the five-site graph of Fig. 2 in arXiv:2405.15836: edges 1-2, 1-3, 2-4,
// 3-4, 3-5, 4-5, sources 1 and 2, sink 5 (0-based here)
void BuildFiveSites(Topology& top) {
  Eigen::Matrix3d box = Eigen::Matrix3d::Zero();
  box.diagonal() << 3.5 * nm, 5 * nm, 5 * nm;
  top.setBox(box);
  std::vector<Eigen::Vector3d> pos = {{0.5, 2.5, 2.5},
                                      {0.5, 1.5, 2.5},
                                      {1.5, 2.5, 2.5},
                                      {1.5, 1.5, 2.5},
                                      {2.5, 2.0, 2.5}};
  std::vector<double> energy = {0.02, -0.05, 0.08, -0.11, 0.0};
  for (Index i = 0; i < 5; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(Atom(0, "C", pos[i] * nm));
    seg.setEMpoles(QMStateType::Hole, energy[i] * ev);
    seg.setU_nX_nN(0.1 * ev, QMStateType::Hole);
    seg.setU_xN_xX(0.1 * ev, QMStateType::Hole);
  }
  std::vector<std::array<Index, 2>> edges = {{0, 1}, {0, 2}, {1, 3},
                                             {2, 3}, {2, 4}, {3, 4}};
  double J2 = 1e-4;
  for (const auto& e : edges) {
    const Segment& a = top.Segments()[e[0]];
    const Segment& b = top.Segments()[e[1]];
    QMPair& pair = top.NBList().Add(a, b, b.getPos() - a.getPos());
    pair.setJeff2(J2 * ev * ev, QMStateType::Hole);
    J2 *= 0.7;
  }
}

ToFGraph::Options FiveSiteRegions() {
  ToFGraph::Options opt;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 1.2 * nm;
  return opt;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(graphrandomwalk_test)

BOOST_AUTO_TEST_CASE(chain_analytic) {
  // birth-death chain 0..n-1, reflecting at 0, absorbing at n-1:
  // tau_0 = sum_{k=0}^{n-2} sum_{j=0}^{k} (1/f_j) prod_{l=j+1}^{k} b_l/f_l
  Index n = 12;
  Topology top;
  BuildChain(top, n, 0.05 * ev);
  ToFGraph::Options opt;
  opt.source_thickness = 1.0 * nm;
  opt.sink_thickness = 2.6 * nm;  // box is n+2 nm: only the last site
  for (Eigen::Vector3d field :
       {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(2e7 * Vm, 0, 0),
        Eigen::Vector3d(-1e7 * Vm, 0, 0)}) {
    ToFGraph graph(top, opt, QMStateType::Hole, kT, field);
    BOOST_REQUIRE_EQUAL(graph.SourceSites().size(), 1);
    BOOST_REQUIRE_EQUAL(graph.SinkSites().size(), 1);
    BOOST_REQUIRE_EQUAL(graph.SinkSites()[0], n - 1);
    double expected = 0.0;
    for (Index k = 0; k <= n - 2; ++k) {
      for (Index j = 0; j <= k; ++j) {
        double term = 1.0 / Rate(graph, j, j + 1);
        for (Index l = j + 1; l <= k; ++l) {
          term *= Rate(graph, l, l - 1) / Rate(graph, l, l + 1);
        }
        expected += term;
      }
    }
    GraphRandomWalk grw(graph);
    for (auto solver :
         {GraphRandomWalk::Solver::Cholesky, GraphRandomWalk::Solver::LU}) {
      GraphRandomWalk::Result r = grw.SingleCarrier(solver);
      BOOST_CHECK_CLOSE(r.tau_mean, expected, 1e-9);
      BOOST_CHECK_CLOSE(r.tau_harmonic, expected, 1e-9);
      BOOST_CHECK_CLOSE(r.occupation_time.sum(), r.tau_mean, 1e-9);
      BOOST_CHECK_LT(r.backward_error, 1e-13);
    }
  }
}

BOOST_AUTO_TEST_CASE(paper_lattice_no_disorder) {
  // independent reference (scipy, same constants, single source (0,0,0); all
  // sources are equivalent here): 1.136718239319e-11 s, CV 0.8237544710
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Result r = grw.SingleCarrier();
  BOOST_CHECK(r.solver.find("cholesky") != std::string::npos);
  BOOST_CHECK_EQUAL(r.unknowns, 448);
  BOOST_CHECK_EQUAL(r.sources.size(), 64);
  BOOST_CHECK_CLOSE(r.tau_mean, 1.136718239319e-11, 1e-4);
  BOOST_CHECK_CLOSE(r.tau_harmonic, r.tau_mean, 1e-8);
  BOOST_CHECK_CLOSE(r.tau_std_injection / r.tau_mean, 0.8237544710, 1e-4);
  BOOST_CHECK_CLOSE(r.occupation_time.sum(), r.tau_mean, 1e-9);
  for (Index s : graph.SinkSites()) {
    BOOST_CHECK_EQUAL(r.occupation_time[s], 0.0);
  }
}

BOOST_AUTO_TEST_CASE(disorder_solvers_and_dense_reference) {
  for (auto disorder : {LatticeModel::Disorder::Uncorrelated,
                        LatticeModel::Disorder::Correlated}) {
    LatticeModel::Options lopt = PaperLattice(disorder, 8);
    LatticeModel model(lopt);
    Topology top;
    model.BuildTopology(top);
    for (QMStateType carrier :
         {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole)}) {
      // field along the normal: reversible; electrons drift against it
      ToFGraph graph(top, Slabs(), carrier, kT,
                     Eigen::Vector3d(4e7 * Vm, 0, 0));
      BOOST_REQUIRE(graph.isReversible());
      GraphRandomWalk grw(graph);
      GraphRandomWalk::Result chol =
          grw.SingleCarrier(GraphRandomWalk::Solver::Cholesky);
      GraphRandomWalk::Result lu =
          grw.SingleCarrier(GraphRandomWalk::Solver::LU);
      GraphRandomWalk::Result cg =
          grw.SingleCarrier(GraphRandomWalk::Solver::CG);
      GraphRandomWalk::Result bicg =
          grw.SingleCarrier(GraphRandomWalk::Solver::BiCGSTAB);
      BOOST_CHECK(cg.converged);
      BOOST_CHECK(bicg.converged);
      BOOST_CHECK_GT(cg.iterations, 0);
      Eigen::VectorXd dense = DenseHittingTimes(graph);
      // tolerances in percent; with sigma = 0.19 eV the rows of Q cancel
      // to about 1 part in D_i tau_i ~ 1e8, which limits every solver to a
      // relative accuracy of about 1e-9
      for (Index k = 0; k < Index(chol.sources.size()); ++k) {
        BOOST_CHECK_CLOSE(chol.tau[k], dense[chol.sources[k][0]], 1e-6);
        BOOST_CHECK_CLOSE(lu.tau[k], dense[lu.sources[k][0]], 1e-6);
        BOOST_CHECK_CLOSE(cg.tau[k], dense[cg.sources[k][0]], 1e-6);
        BOOST_CHECK_CLOSE(bicg.tau[k], dense[bicg.sources[k][0]], 1e-6);
        BOOST_CHECK_CLOSE(chol.tau_std[k], lu.tau_std[k], 1e-5);
        BOOST_CHECK_CLOSE(cg.tau_std[k], lu.tau_std[k], 1e-5);
      }
      BOOST_CHECK_CLOSE(chol.occupation_time.sum(), chol.tau_mean, 1e-6);
      BOOST_CHECK_CLOSE(lu.occupation_time.sum(), lu.tau_mean, 1e-6);
      BOOST_CHECK_LT((chol.occupation_time - lu.occupation_time).norm(),
                     1e-7 * lu.occupation_time.norm());
      BOOST_CHECK_LT((cg.occupation_time - lu.occupation_time).norm(),
                     1e-7 * lu.occupation_time.norm());
      BOOST_CHECK_LT((bicg.occupation_time - lu.occupation_time).norm(),
                     1e-7 * lu.occupation_time.norm());
      BOOST_CHECK_LE(chol.tau_harmonic, chol.tau_mean);
    }
  }
}

BOOST_AUTO_TEST_CASE(nonreversible_field) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated, 6));
  Topology top;
  model.BuildTopology(top);
  // field component along a periodic direction
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT,
                 Eigen::Vector3d(3e7 * Vm, 2e7 * Vm, 0));
  BOOST_REQUIRE(!graph.isReversible());
  GraphRandomWalk grw(graph);
  BOOST_CHECK_THROW(grw.SingleCarrier(GraphRandomWalk::Solver::Cholesky),
                    std::runtime_error);
  BOOST_CHECK_THROW(grw.SingleCarrier(GraphRandomWalk::Solver::CG),
                    std::runtime_error);
  GraphRandomWalk::Result bicg =
      grw.SingleCarrier(GraphRandomWalk::Solver::BiCGSTAB);
  GraphRandomWalk::Result r = grw.SingleCarrier();
  BOOST_CHECK_EQUAL(r.solver, "lu");
  Eigen::VectorXd dense = DenseHittingTimes(graph);
  for (Index k = 0; k < Index(r.sources.size()); ++k) {
    BOOST_CHECK_CLOSE(r.tau[k], dense[r.sources[k][0]], 1e-6);
    BOOST_CHECK_CLOSE(bicg.tau[k], dense[bicg.sources[k][0]], 1e-6);
  }
  BOOST_CHECK_CLOSE(r.occupation_time.sum(), r.tau_mean, 1e-6);
}

BOOST_AUTO_TEST_CASE(auto_selection_and_iteration_limit) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated, 6));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Settings settings;
  settings.direct_limit = 10;  // 180 unknowns: iterative
  grw.setSettings(settings);
  GraphRandomWalk::Result r = grw.SingleCarrier();
  BOOST_CHECK(r.solver.find("cg") == 0);
  BOOST_CHECK(r.converged);
  settings.max_iterations = 3;
  grw.setSettings(settings);
  GraphRandomWalk::Result truncated = grw.SingleCarrier();
  BOOST_CHECK(!truncated.converged);
}

BOOST_AUTO_TEST_CASE(unreachable_tolerance) {
  // a tolerance below rounding: CG stops at its floor, BiCGSTAB may report
  // NoConvergence depending on the Eigen version; both are accepted when the
  // refined result is accurate, not only when the iteration limit is hit
  LatticeModel model(PaperLattice(LatticeModel::Disorder::Uncorrelated, 6));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT,
                 Eigen::Vector3d(2e7 * Vm, 0, 0));
  Eigen::VectorXd dense = DenseHittingTimes(graph);
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Settings settings;
  settings.tolerance = 1e-19;
  settings.max_iterations = 20000;
  grw.setSettings(settings);
  for (auto solver :
       {GraphRandomWalk::Solver::CG, GraphRandomWalk::Solver::BiCGSTAB}) {
    GraphRandomWalk::Result r = grw.SingleCarrier(solver);
    BOOST_CHECK_LT(r.backward_error, 1e-10);
    for (Index k = 0; k < Index(r.sources.size()); ++k) {
      BOOST_CHECK_CLOSE(r.tau[k], dense[r.sources[k][0]], 1e-6);
    }
    if (solver == GraphRandomWalk::Solver::CG) {
      // stopped at the floor, far below the iteration limit
      BOOST_CHECK(r.converged);
      BOOST_CHECK_LT(r.iterations, 3 * 6 * 20000 / 10);
    }
  }
}

BOOST_AUTO_TEST_CASE(two_carriers_five_sites) {
  Topology top;
  BuildFiveSites(top);
  for (Eigen::Vector3d field :
       {Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(3e7 * Vm, 0, 0)}) {
    ToFGraph graph(top, FiveSiteRegions(), QMStateType::Hole, kT, field);
    BOOST_REQUIRE_EQUAL(graph.SourceSites().size(), 2);
    BOOST_REQUIRE_EQUAL(graph.SinkSites().size(), 1);
    GraphRandomWalk grw(graph);
    for (auto obs : {GraphRandomWalk::Observable::First,
                     GraphRandomWalk::Observable::Transit}) {
      bool transit = obs == GraphRandomWalk::Observable::Transit;
      auto dense = DenseTwoCarriers(graph, transit);
      for (auto solver :
           {GraphRandomWalk::Solver::Cholesky, GraphRandomWalk::Solver::LU,
            GraphRandomWalk::Solver::CG, GraphRandomWalk::Solver::BiCGSTAB}) {
        GraphRandomWalk::Result r = grw.TwoCarriers(obs, solver);
        BOOST_CHECK_EQUAL(r.unknowns, 6);  // pairs of the 4 non-sink sites
        BOOST_REQUIRE_EQUAL(r.sources.size(), 1);
        BOOST_CHECK_CLOSE(r.tau[0], dense.at({0, 1}), 1e-8);
        BOOST_CHECK_CLOSE(r.occupation_time.sum(), 2.0 * r.tau_mean, 1e-8);
        BOOST_CHECK(r.converged);
      }
    }
    // first arrival is earlier than the mean transit time
    BOOST_CHECK_LT(
        grw.TwoCarriers(GraphRandomWalk::Observable::First).tau_mean,
        grw.TwoCarriers(GraphRandomWalk::Observable::Transit).tau_mean);
  }
}

BOOST_AUTO_TEST_CASE(two_carriers_disordered_lattice) {
  // 4^3 lattice: 48 transient sites, 1128 pair states (dense reference)
  LatticeModel::Options lopt =
      PaperLattice(LatticeModel::Disorder::Uncorrelated, 4);
  lopt.cutoff = 1.5 * nm;
  LatticeModel model(lopt);
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Electron, kT,
                 Eigen::Vector3d(-2e7 * Vm, 0, 0));
  GraphRandomWalk grw(graph);
  for (auto obs : {GraphRandomWalk::Observable::First,
                   GraphRandomWalk::Observable::Transit}) {
    auto dense =
        DenseTwoCarriers(graph, obs == GraphRandomWalk::Observable::Transit);
    for (auto solver :
         {GraphRandomWalk::Solver::Auto, GraphRandomWalk::Solver::LU,
          GraphRandomWalk::Solver::CG, GraphRandomWalk::Solver::BiCGSTAB}) {
      GraphRandomWalk::Result r = grw.TwoCarriers(obs, solver);
      BOOST_CHECK_EQUAL(r.unknowns, 1128);
      BOOST_CHECK_EQUAL(r.sources.size(), 16 * 15 / 2);
      for (Index k = 0; k < Index(r.sources.size()); ++k) {
        Index a = std::min(r.sources[k][0], r.sources[k][1]);
        Index b = std::max(r.sources[k][0], r.sources[k][1]);
        BOOST_CHECK_CLOSE(r.tau[k], dense.at({a, b}), 1e-6);
      }
      BOOST_CHECK_CLOSE(r.occupation_time.sum(), 2.0 * r.tau_mean, 1e-6);
      BOOST_CHECK_LT(r.backward_error, 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(two_carriers_paper_lattice_no_disorder) {
  // independent reference (scipy, same constants): sources (0,0,0) and
  // (0,0,1) nm + 0.5 = segments 0 and 1, first arrival 6.651326554700e-12 s
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Result r =
      grw.TwoCarriers(GraphRandomWalk::Observable::First);
  BOOST_CHECK_EQUAL(r.unknowns, 448 * 447 / 2);
  BOOST_CHECK(r.solver.find("cg") == 0);
  bool found = false;
  for (Index k = 0; k < Index(r.sources.size()); ++k) {
    if (std::min(r.sources[k][0], r.sources[k][1]) == 0 &&
        std::max(r.sources[k][0], r.sources[k][1]) == 1) {
      BOOST_CHECK_CLOSE(r.tau[k], 6.651326554700e-12, 1e-3);
      found = true;
    }
  }
  BOOST_CHECK(found);
}

BOOST_AUTO_TEST_CASE(transit_symmetric_rates_identity) {
  // with symmetric rates (no site-energy disorder, no field) exclusion drops
  // out of the mean occupations (as in the symmetric simple exclusion
  // process): the mean transit time per carrier of a pair equals the
  // average of the two single-carrier hitting times
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None, 6));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Result one = grw.SingleCarrier();
  std::map<Index, double> single;
  for (Index k = 0; k < Index(one.sources.size()); ++k) {
    single[one.sources[k][0]] = one.tau[k];
  }
  GraphRandomWalk::Result two =
      grw.TwoCarriers(GraphRandomWalk::Observable::Transit);
  for (Index k = 0; k < Index(two.sources.size()); ++k) {
    double expected =
        0.5 * (single[two.sources[k][0]] + single[two.sources[k][1]]);
    BOOST_CHECK_CLOSE(two.tau[k], expected, 1e-7);
  }

  // with energetic disorder the identity fails: occupied traps speed up the
  // second carrier
  LatticeModel disordered(
      PaperLattice(LatticeModel::Disorder::Uncorrelated, 6));
  Topology top2;
  disordered.BuildTopology(top2);
  ToFGraph graph2(top2, Slabs(), QMStateType::Hole, kT,
                  Eigen::Vector3d::Zero());
  GraphRandomWalk grw2(graph2);
  double one_mean = grw2.SingleCarrier().tau_mean;
  double two_mean =
      grw2.TwoCarriers(GraphRandomWalk::Observable::Transit).tau_mean;
  BOOST_CHECK_LT(two_mean, 0.99 * one_mean);
}

BOOST_AUTO_TEST_CASE(two_carriers_state_limit) {
  LatticeModel model(PaperLattice(LatticeModel::Disorder::None, 6));
  Topology top;
  model.BuildTopology(top);
  ToFGraph graph(top, Slabs(), QMStateType::Hole, kT, Eigen::Vector3d::Zero());
  GraphRandomWalk grw(graph);
  GraphRandomWalk::Settings settings;
  settings.max_states = 1000;
  grw.setSettings(settings);
  BOOST_CHECK_THROW(grw.TwoCarriers(GraphRandomWalk::Observable::First),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(solver_names) {
  BOOST_CHECK(GraphRandomWalk::SolverFromString("lu") ==
              GraphRandomWalk::Solver::LU);
  BOOST_CHECK(GraphRandomWalk::SolverFromString("cg") ==
              GraphRandomWalk::Solver::CG);
  BOOST_CHECK(GraphRandomWalk::SolverFromString("bicgstab") ==
              GraphRandomWalk::Solver::BiCGSTAB);
  BOOST_CHECK_THROW(GraphRandomWalk::SolverFromString("qmr"),
                    std::runtime_error);
  BOOST_CHECK(GraphRandomWalk::ObservableFromString("transit") ==
              GraphRandomWalk::Observable::Transit);
  BOOST_CHECK_THROW(GraphRandomWalk::ObservableFromString("last"),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
