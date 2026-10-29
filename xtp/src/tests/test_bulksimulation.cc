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

#define BOOST_TEST_MODULE bulksimulation_test

// Standard includes
#include <cmath>
#include <cstdint>
#include <unordered_map>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/bulksimulation.h"
#include "votca/xtp/eigen.h"
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/latticemodel.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/tofsimulation.h"

using namespace votca::xtp;
using votca::Index;

namespace {
const double nm = votca::tools::conv::nm2bohr;
const double ev = votca::tools::conv::ev2hrt;
const double kT = votca::tools::conv::kB * 300.0 * ev;
const double Vm = ev / (1e9 * nm);

// 3x3x3 periodic cubic lattice, nearest neighbours, Gaussian disorder
void BuildLattice(Topology& top) {
  LatticeModel::Options opt;
  opt.dims = {3, 3, 3};
  opt.spacing = 1.0 * nm;
  opt.cutoff = 1.1 * nm;
  opt.J2 = 1.79e-4 * ev * ev;
  opt.decay_length = 1.0 * nm;
  opt.reorganization = 0.23 * ev;
  opt.mean_energy = -0.76 * ev;
  opt.sigma = 0.1 * ev;
  opt.disorder = LatticeModel::Disorder::Uncorrelated;
  opt.seed = 1;
  LatticeModel(opt).BuildTopology(top);
}

Eigen::MatrixXd Generator(const ToFGraph& g) {
  // (L f)_i = sum_j w_ij (f_j - f_i)
  Index n = g.NumberOfSites();
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
  for (Index i = 0; i < n; ++i) {
    for (const ToFGraph::Edge& e : g.OutEdges(i)) {
      L(i, e.to) += e.rate;
      L(i, i) -= e.rate;
    }
  }
  return L;
}

// stationary distribution of the chain with generator L (rows: from)
Eigen::VectorXd Stationary(const Eigen::MatrixXd& L) {
  Eigen::MatrixXd A = L.transpose();
  A.row(0).setOnes();
  Eigen::VectorXd b = Eigen::VectorXd::Zero(L.rows());
  b(0) = 1.0;
  return A.partialPivLu().solve(b);
}

// one carrier: exact drift velocity and diffusion tensor; with the
// corrector chi, L chi = v - b, r + chi(site) - v t is a martingale and
// D = 1/2 sum_i pi_i sum_j w_ij (R_ij + chi_j - chi_i)(...)^T
void SingleCarrierExact(const ToFGraph& g, Eigen::Vector3d& v,
                        Eigen::Matrix3d& D, Eigen::VectorXd& pi) {
  Index n = g.NumberOfSites();
  Eigen::MatrixXd L = Generator(g);
  pi = Stationary(L);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n, 3);
  for (Index i = 0; i < n; ++i) {
    for (const ToFGraph::Edge& e : g.OutEdges(i)) {
      b.row(i) += e.rate * e.dR.transpose();
    }
  }
  v = (pi.transpose() * b).transpose();
  Eigen::MatrixXd A(n + 1, n);
  A << L, pi.transpose();
  Eigen::MatrixXd rhs(n + 1, 3);
  rhs << (Eigen::VectorXd::Ones(n) * v.transpose() - b),
      Eigen::RowVector3d::Zero();
  Eigen::MatrixXd chi = A.colPivHouseholderQr().solve(rhs);
  D.setZero();
  for (Index i = 0; i < n; ++i) {
    for (const ToFGraph::Edge& e : g.OutEdges(i)) {
      Eigen::Vector3d d = e.dR + (chi.row(e.to) - chi.row(i)).transpose();
      D += 0.5 * pi(i) * e.rate * d * d.transpose();
    }
  }
}

// N carriers with exclusion: exact mean drift velocity per carrier from
// the stationary distribution over all configurations (n <= 32)
Eigen::Vector3d ExclusionExact(const ToFGraph& g, Index N) {
  const Index n = g.NumberOfSites();
  std::vector<std::uint32_t> configs;
  for (std::uint64_t c = 0; c < (std::uint64_t(1) << n); ++c) {
    if (Index(__builtin_popcountll(c)) == N) {
      configs.push_back(std::uint32_t(c));
    }
  }
  std::unordered_map<std::uint32_t, Index> index;
  for (Index k = 0; k < Index(configs.size()); ++k) {
    index[configs[k]] = k;
  }
  const Index M = Index(configs.size());
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(M, M);
  Eigen::MatrixXd drift = Eigen::MatrixXd::Zero(M, 3);
  for (Index k = 0; k < M; ++k) {
    std::uint32_t c = configs[k];
    for (Index s = 0; s < n; ++s) {
      if (!(c & (1u << s))) {
        continue;
      }
      for (const ToFGraph::Edge& e : g.OutEdges(s)) {
        if (c & (1u << e.to)) {
          continue;
        }
        std::uint32_t next = (c & ~(1u << s)) | (1u << e.to);
        L(k, index.at(next)) += e.rate;
        L(k, k) -= e.rate;
        drift.row(k) += e.rate * e.dR.transpose();
      }
    }
  }
  Eigen::VectorXd p = Stationary(L);
  return (p.transpose() * drift).transpose() / double(N);
}

BulkSimulation::Settings EventLimited(Index carriers, Index runs,
                                      Index events) {
  BulkSimulation::Settings s;
  s.seed = 11;
  s.runs = runs;
  s.carriers = carriers;
  s.limits.max_events = events;
  return s;
}

// |kmc - exact| within 4.5 standard errors and a meaningful error
void CheckValue(double kmc, double error, double exact, double max_rel) {
  BOOST_CHECK_GT(error, 0.0);
  BOOST_CHECK_LT(std::abs(kmc - exact), 4.5 * error);
  BOOST_CHECK_LT(error, max_rel * std::abs(exact));
}
}  // namespace

BOOST_AUTO_TEST_SUITE(bulksimulation_test)

BOOST_AUTO_TEST_CASE(single_carrier_drift_and_occupation) {
  Topology top;
  BuildLattice(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                 Eigen::Vector3d(5e7, 0, 0) * Vm);
  Eigen::Vector3d v;
  Eigen::Matrix3d D;
  Eigen::VectorXd pi;
  SingleCarrierExact(graph, v, D, pi);

  BulkSimulation sim(graph, EventLimited(1, 4, 2000000));
  BulkSimulation::Result r = sim.Simulate();
  // electrons drift against the field
  BOOST_CHECK_LT(v.x(), 0.0);
  CheckValue(r.velocity.x(), r.velocity_error.x(), v.x(), 0.02);
  BOOST_CHECK_CLOSE(r.occupation.sum(), 1.0, 1e-8);
  BOOST_CHECK_LT((r.occupation - pi).cwiseAbs().maxCoeff(), 0.02);
}

BOOST_AUTO_TEST_CASE(exclusion_drift) {
  // the previous kmcmultiple forced a hop after every waiting time drawn
  // from the full escape rates, which overestimated the velocity by 5 % (2
  // carriers) and 22 % (3 carriers) on this lattice
  Topology top;
  BuildLattice(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                 Eigen::Vector3d(5e7, 0, 0) * Vm);
  for (Index N : {2, 25}) {
    Eigen::Vector3d v = ExclusionExact(graph, N);
    BulkSimulation sim(graph, EventLimited(N, 4, 2000000));
    BulkSimulation::Result r = sim.Simulate();
    CheckValue(r.velocity.x(), r.velocity_error.x(), v.x(), 0.02);
    BOOST_CHECK_GT(r.null_events, 0);
    BOOST_CHECK_CLOSE(r.occupation.sum(), double(N), 1e-8);
  }
}

BOOST_AUTO_TEST_CASE(single_carrier_diffusion) {
  Topology top;
  BuildLattice(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                 Eigen::Vector3d::Zero());
  Eigen::Vector3d v;
  Eigen::Matrix3d D;
  Eigen::VectorXd pi;
  SingleCarrierExact(graph, v, D, pi);
  BOOST_CHECK_SMALL(v.norm(), 1e-6 * std::sqrt(D.trace()) / nm);

  BulkSimulation::Settings s = EventLimited(1, 8, 4000000);
  s.blocks = 100;
  BulkSimulation sim(graph, s);
  BulkSimulation::Result r = sim.Simulate();
  for (Index d = 0; d < 3; ++d) {
    CheckValue(r.diffusion(d, d), r.diffusion_error(d, d), D(d, d), 0.1);
    BOOST_CHECK_LT(std::abs(r.velocity(d)), 4.5 * r.velocity_error(d));
  }
  for (const BulkSimulation::RunData& run : r.runs) {
    BOOST_CHECK_GE(run.blocks, 100);
    BOOST_CHECK_LE(run.blocks, 200);
  }
}

BOOST_AUTO_TEST_CASE(exact_against_dense_reference) {
  Topology top;
  BuildLattice(top);
  for (double F : {0.0, 5e7}) {
    ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron,
                   kT, Eigen::Vector3d(F, 0, 0) * Vm);
    Eigen::Vector3d v;
    Eigen::Matrix3d D;
    Eigen::VectorXd pi;
    SingleCarrierExact(graph, v, D, pi);
    for (Index limit : {Index(100), Index(0)}) {  // direct and iterative
      BulkSimulation::ExactResult r = BulkSimulation::Exact(graph, limit);
      BOOST_REQUIRE(r.available);
      BOOST_CHECK(r.converged);
      BOOST_CHECK_EQUAL(r.sites, 27);
      // gross flux sum_i pi_i sum_j w_ij |R_ij|: the scale of v's rounding
      double vscale = 0.0;
      for (Index i = 0; i < graph.NumberOfSites(); ++i) {
        for (const ToFGraph::Edge& e : graph.OutEdges(i)) {
          vscale += pi(i) * e.rate * e.dR.norm();
        }
      }
      BOOST_CHECK_SMALL((r.velocity - v).norm() / vscale, 1e-8);
      BOOST_CHECK_SMALL((r.diffusion - D).norm() / D.norm(), 1e-8);
      BOOST_CHECK_SMALL((r.occupation - pi).cwiseAbs().maxCoeff(), 1e-8);
    }
  }
}

BOOST_AUTO_TEST_CASE(exact_needs_one_closed_class) {
  // two separate chains: two closed classes
  Topology top;
  Eigen::Matrix3d box = Eigen::Matrix3d::Identity() * 4.0 * nm;
  top.setBox(box);
  for (Index i = 0; i < 4; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(
        Atom(0, "C", Eigen::Vector3d(0.5 + double(i), 0.5, 0.5) * nm));
    seg.setEMpoles(QMStateType::Electron, 0.0);
    seg.setU_nX_nN(0.1 * ev, QMStateType::Electron);
    seg.setU_xN_xX(0.1 * ev, QMStateType::Electron);
  }
  for (Index i : {0, 2}) {
    QMPair& pair = top.NBList().Add(top.Segments()[i], top.Segments()[i + 1],
                                    Eigen::Vector3d(1, 0, 0) * nm);
    pair.setJeff2(1e-4 * ev * ev, QMStateType::Electron);
  }
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                 Eigen::Vector3d::Zero());
  BulkSimulation::ExactResult r = BulkSimulation::Exact(graph);
  BOOST_CHECK(!r.available);
  BOOST_CHECK(!r.message.empty());
}

BOOST_AUTO_TEST_CASE(time_limit_window_and_trajectory) {
  Topology top;
  BuildLattice(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Hole, kT,
                 Eigen::Vector3d(0, 2e7, 0) * Vm);
  BulkSimulation::Settings s;
  s.carriers = 3;
  s.runs = 2;
  s.limits.max_time = 1e-7;
  s.equilibration = 2e-8;
  s.sample_interval = 1e-8;
  BulkSimulation sim(graph, s);
  BulkSimulation::Result r = sim.Simulate();
  for (const BulkSimulation::RunData& run : r.runs) {
    BOOST_CHECK_EQUAL(run.start, 2e-8);
    BOOST_CHECK_EQUAL(run.end, 1e-7);
    BOOST_CHECK_CLOSE(double(run.blocks) * run.block_length, 8e-8, 5.0);
  }
  BOOST_REQUIRE_EQUAL(r.trajectory.size(), 10);
  BOOST_CHECK_EQUAL(r.trajectory.front().time, 0.0);
  BOOST_CHECK_CLOSE(r.trajectory.back().time, 9e-8, 1e-10);
  // the samples follow the carriers: unwrapped positions are continuous in
  // the sense that the displacement over the window matches the trajectory
  BOOST_CHECK_EQUAL(r.trajectory.front().positions.size(), 3);
  // holes drift along the field
  BOOST_CHECK_GT(r.velocity.y(), 0.0);
  BOOST_CHECK_CLOSE(r.occupation.sum(), 3.0, 1e-8);
}

BOOST_AUTO_TEST_CASE(reproducible) {
  Topology top;
  BuildLattice(top);
  ToFGraph graph(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                 Eigen::Vector3d(5e7, 0, 0) * Vm);
  BulkSimulation sim(graph, EventLimited(3, 3, 100000));
  BulkSimulation::Result a = sim.Simulate();
  BulkSimulation::Result b = sim.Simulate();
  BOOST_CHECK_EQUAL(a.hops, b.hops);
  BOOST_CHECK(a.velocity == b.velocity);
  BOOST_CHECK(a.diffusion == b.diffusion);
  BOOST_CHECK(a.runs[1].velocity != a.runs[2].velocity);
}

BOOST_AUTO_TEST_CASE(invalid_setups) {
  Topology top;
  BuildLattice(top);
  ToFGraph bulk(top, ToFGraph::BulkOptions("*", ""), QMStateType::Electron, kT,
                Eigen::Vector3d::Zero());
  ToFGraph::Options slabs;
  slabs.source_thickness = 1.0 * nm;
  slabs.sink_thickness = 1.0 * nm;
  ToFGraph tof(top, slabs, QMStateType::Electron, kT, Eigen::Vector3d::Zero());
  BOOST_CHECK_THROW(BulkSimulation(tof, EventLimited(1, 1, 10)),
                    std::runtime_error);
  BOOST_CHECK_THROW(BulkSimulation(bulk, EventLimited(28, 1, 10)),
                    std::runtime_error);
  BulkSimulation::Settings unlimited;
  BOOST_CHECK_THROW(BulkSimulation(bulk, unlimited), std::runtime_error);
  BOOST_CHECK_THROW(GraphRandomWalk{bulk}, std::runtime_error);
  ToFSimulation::Settings tofsettings;
  BOOST_CHECK_THROW(ToFSimulation(bulk, tofsettings), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
