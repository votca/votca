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
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

// Third party includes
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

// Local VOTCA includes
#include "votca/xtp/graphrandomwalk.h"
#include "votca/xtp/scaledsolve.h"

namespace votca {
namespace xtp {

namespace {
using SpMat = Eigen::SparseMatrix<double>;
using SpMatR = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using SolveFn = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;
using Solver = GraphRandomWalk::Solver;
using Settings = GraphRandomWalk::Settings;

// conjugate gradients with the Jacobi preconditioner; the hitting times are
// dominated by the slow modes (traps), so the residual has to be driven far
// down (default 1e-14) before they are resolved
Eigen::VectorXd PCG(const SpMatR& A, const Eigen::VectorXd& b,
                    const Settings& settings, Index& its, bool& converged) {
  Eigen::VectorXd dinv = A.diagonal().cwiseInverse();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(b.size());
  Eigen::VectorXd r = b;
  double bnorm = b.norm();
  its = 0;
  converged = true;
  if (bnorm == 0.0) {
    return x;
  }
  Eigen::VectorXd z = dinv.cwiseProduct(r);
  Eigen::VectorXd d = z;
  Eigen::VectorXd Ad(b.size());
  double rz = r.dot(z);
  const double threshold2 =
      settings.tolerance * settings.tolerance * bnorm * bnorm;
  // stop at a rounding floor just above the tolerance instead of running
  // into the iteration limit; the refinement decides. Only close to the
  // tolerance: further up, plateaus of the residual are normal while the
  // slow trap modes are resolved.
  const double floor2 = 1e6 * threshold2;  // 1000 x tolerance
  double best2 = r.squaredNorm();
  Index last_improvement = 0;
  while (r.squaredNorm() > threshold2) {
    if (its >= settings.max_iterations) {
      converged = false;
      break;
    }
    if (best2 < floor2 && its - last_improvement > 1000) {
      break;
    }
    Ad.noalias() = A * d;
    double alpha = rz / d.dot(Ad);
    x += alpha * d;
    r -= alpha * Ad;
    z = dinv.cwiseProduct(r);
    double rz_new = r.dot(z);
    d = z + (rz_new / rz) * d;
    rz = rz_new;
    ++its;
    double r2 = r.squaredNorm();
    if (r2 < 0.81 * best2) {  // residual down by 10 %
      best2 = r2;
      last_improvement = its;
    }
  }
  return x;
}

/**
 * Solves Q x = b and Q^T x = b for a generator Q (row major), with
 * symmetrization by exp(h) if h is given (detailed balance), the solver
 * chosen as requested or automatically, and refinement against Q.
 */
class GeneratorSolver {
 public:
  GeneratorSolver(const SpMatR& Q, const Eigen::VectorXd& h, Solver solver,
                  const Settings& settings)
      : Q_(Q), settings_(settings) {
    Qabs_ = Q_.cwiseAbs();
    Qt_ = Q_.transpose();
    Qtabs_ = Qt_.cwiseAbs();
    const Index m = Q_.rows();

    bool symmetric = (solver == Solver::Cholesky || solver == Solver::CG ||
                      (solver == Solver::Auto && h.size() == m));
    if ((solver == Solver::Cholesky || solver == Solver::CG) && h.size() != m) {
      throw std::runtime_error(
          "GraphRandomWalk: the cholesky and cg solvers need rates with "
          "detailed balance, use lu or bicgstab");
    }
    if (symmetric && h.cwiseAbs().maxCoeff() > 600.0) {
      if (solver != Solver::Auto) {
        throw std::runtime_error(
            "GraphRandomWalk: the potential range is too large for the "
            "symmetrized solvers, use lu or bicgstab");
      }
      symmetric = false;
    }
    if (solver == Solver::Auto) {
      bool direct = m <= settings_.direct_limit;
      if (symmetric) {
        solver = direct ? Solver::Cholesky : Solver::CG;
      } else {
        solver = direct ? Solver::LU : Solver::BiCGSTAB;
      }
    }
    iterative_ = (solver == Solver::CG || solver == Solver::BiCGSTAB);

    if (symmetric) {
      p_ = h.array().exp();
      pinv_ = p_.cwiseInverse();
      // S = P Q P^-1 is symmetric; Q = P^-1 S P, Q^T = P S P^-1
      SpMatR S = p_.asDiagonal() * Q_ * pinv_.asDiagonal();
      SpMatR St = S.transpose();
      S = 0.5 * (S + St);
      if (solver == Solver::Cholesky) {
        name = "cholesky (LDLT of the symmetrized generator)";
        ldlt_.compute(SpMat(S));
        if (ldlt_.info() != Eigen::Success) {
          throw std::runtime_error(
              "GraphRandomWalk: Cholesky factorization failed");
        }
        solveS_ = [this](const Eigen::VectorXd& b) {
          return Eigen::VectorXd(ldlt_.solve(b));
        };
      } else {
        name = "cg (Jacobi preconditioned, symmetrized generator)";
        S_ = std::move(S);
        solveS_ = [this](const Eigen::VectorXd& b) {
          Index its = 0;
          bool ok = true;
          Eigen::VectorXd y = PCG(S_, b, settings_, its, ok);
          iterations += its;
          converged = converged && ok;
          return y;
        };
      }
      solve_ = [this](const Eigen::VectorXd& b) {
        Eigen::VectorXd y = solveS_((p_.array() * b.array()).matrix());
        return Eigen::VectorXd(pinv_.array() * y.array());
      };
      solve_adjoint_ = [this](const Eigen::VectorXd& b) {
        Eigen::VectorXd y = solveS_((pinv_.array() * b.array()).matrix());
        return Eigen::VectorXd(p_.array() * y.array());
      };
    } else if (solver == Solver::LU) {
      name = "lu";
      Qc_ = Q_;
      lu_.analyzePattern(Qc_);
      lu_.factorize(Qc_);
      if (lu_.info() != Eigen::Success) {
        throw std::runtime_error("GraphRandomWalk: LU factorization failed: " +
                                 lu_.lastErrorMessage());
      }
      solve_ = [this](const Eigen::VectorXd& b) {
        return Eigen::VectorXd(lu_.solve(b));
      };
      solve_adjoint_ = [this](const Eigen::VectorXd& b) {
        return Eigen::VectorXd(lu_.transpose().solve(b));
      };
    } else {
      name = "bicgstab (incomplete LU preconditioned)";
      for (auto* b : {&bicg_, &bicg_adjoint_}) {
        b->setTolerance(settings_.tolerance);
        b->setMaxIterations(settings_.max_iterations);
        // a cheap incomplete factorization; accuracy comes from the
        // iterations and the refinement
        b->preconditioner().setDroptol(1e-4);
        b->preconditioner().setFillfactor(5);
      }
      // the iterative solvers keep a reference to the matrix
      Qc_ = Q_;
      Qtc_ = Qt_;
      bicg_.compute(Qc_);
      bicg_adjoint_.compute(Qtc_);
      auto make = [this](BiCG& b) {
        return [this, &b](const Eigen::VectorXd& rhs) {
          bool solved = false;
          Eigen::VectorXd x = ScaledSolve(b, rhs, &solved);
          if (!solved) {
            return x;  // zero right-hand side
          }
          iterations += b.iterations();
          // only the iteration limit counts here: Eigen versions differ in
          // reporting a residual floor just above the tolerance or a restart
          // as NoConvergence/NumericalIssue; the refinement against Q and the
          // backward error decide whether the result is usable
          converged = converged && b.iterations() < b.maxIterations();
          return x;
        };
      };
      solve_ = make(bicg_);
      solve_adjoint_ = make(bicg_adjoint_);
    }
  }

  Eigen::VectorXd Solve(const Eigen::VectorXd& b) {
    return Refine(Q_, Qabs_, solve_, b);
  }
  Eigen::VectorXd SolveAdjoint(const Eigen::VectorXd& b) {
    return Refine(Qt_, Qtabs_, solve_adjoint_, b);
  }

  std::string name;
  Index iterations = 0;
  bool converged = true;
  double backward_error = 0.0;
  Index refinement_steps = 0;

 private:
  using BiCG = Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>>;

  // Solve A x = b with an approximate solver and iterative refinement
  // against A; tracks the largest componentwise backward error
  // |b - A x|_i / (|A| |x| + |b|)_i.
  Eigen::VectorXd Refine(const SpMatR& A, const SpMatR& Aabs,
                         const SolveFn& solve, const Eigen::VectorXd& b) {
    // direct solvers: refine to rounding level; iterative solvers: 1e-12
    // suffices, each step is a full solve
    double target =
        iterative_ ? std::max(1e-12, 100.0 * settings_.tolerance) : 1e-14;
    auto error = [&A, &Aabs, &b](const Eigen::VectorXd& y) {
      Eigen::VectorXd r = b - A * y;
      Eigen::VectorXd scale = Aabs * y.cwiseAbs() + b.cwiseAbs();
      double e = 0.0;
      for (Index i = 0; i < r.size(); ++i) {
        if (scale[i] > 0.0) {
          e = std::max(e, std::abs(r[i]) / scale[i]);
        }
      }
      return std::make_pair(e, Eigen::VectorXd(r));
    };
    Eigen::VectorXd x = solve(b);
    Index steps = 0;
    auto [err, r] = error(x);
    // refine until the backward error is at the target or stops improving
    while (err > target && steps < 5) {
      Eigen::VectorXd x_new = x + solve(r);
      auto [err_new, r_new] = error(x_new);
      ++steps;
      if (err_new > 0.5 * err) {
        if (err_new < err) {
          x = x_new;
          err = err_new;
        }
        break;
      }
      x = x_new;
      err = err_new;
      r = r_new;
    }
    if (!x.allFinite()) {
      throw std::runtime_error("GraphRandomWalk: solution is not finite");
    }
    backward_error = std::max(backward_error, err);
    refinement_steps = std::max(refinement_steps, steps);
    return x;
  }

  const SpMatR& Q_;
  Settings settings_;
  SpMatR Qabs_;
  SpMatR Qt_;
  SpMatR Qtabs_;
  bool iterative_ = false;
  Eigen::VectorXd p_;
  Eigen::VectorXd pinv_;
  SpMatR S_;
  SpMat Qc_;  // column-major copies for LU and BiCGSTAB
  SpMat Qtc_;
  Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int>> ldlt_;
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu_;
  BiCG bicg_;
  BiCG bicg_adjoint_;
  SolveFn solveS_;
  SolveFn solve_;
  SolveFn solve_adjoint_;
};

// row-major sparse matrix from rows of (column, value), duplicates summed
class RowBuilder {
 public:
  explicit RowBuilder(Index rows) : rows_(rows) {
    outer_.reserve(rows + 1);
    outer_.push_back(0);
  }
  void Add(Index col, double value) { row_.emplace_back(col, value); }
  void FinishRow() {
    std::sort(row_.begin(), row_.end());
    for (std::size_t k = 0; k < row_.size(); ++k) {
      if (k > 0 && row_[k].first == row_[k - 1].first) {
        values_.back() += row_[k].second;
      } else {
        inner_.push_back(row_[k].first);
        values_.push_back(row_[k].second);
      }
    }
    row_.clear();
    outer_.push_back(Index(inner_.size()));
  }
  SpMatR Build(Index cols) const {
    if (Index(outer_.size()) != rows_ + 1) {
      throw std::runtime_error("RowBuilder: incomplete matrix");
    }
    SpMatR A(rows_, cols);
    A.reserve(Index(inner_.size()));
    for (Index i = 0; i < rows_; ++i) {
      A.startVec(i);
      for (Index k = outer_[i]; k < outer_[i + 1]; ++k) {
        A.insertBack(i, inner_[k]) = values_[k];
      }
    }
    A.finalize();
    return A;
  }

 private:
  Index rows_;
  std::vector<std::pair<Index, double>> row_;
  std::vector<Index> outer_;
  std::vector<Index> inner_;
  std::vector<double> values_;
};

inline Index PairRank(Index i, Index j) {
  if (i > j) {
    std::swap(i, j);
  }
  return j * (j - 1) / 2 + i;
}
}  // namespace

GraphRandomWalk::Solver GraphRandomWalk::SolverFromString(
    const std::string& name) {
  if (name == "auto") {
    return Solver::Auto;
  } else if (name == "cholesky") {
    return Solver::Cholesky;
  } else if (name == "lu") {
    return Solver::LU;
  } else if (name == "cg") {
    return Solver::CG;
  } else if (name == "bicgstab") {
    return Solver::BiCGSTAB;
  }
  throw std::runtime_error("GraphRandomWalk: unknown solver '" + name +
                           "', use auto, cholesky, lu, cg or bicgstab");
}

GraphRandomWalk::Observable GraphRandomWalk::ObservableFromString(
    const std::string& name) {
  if (name == "first") {
    return Observable::First;
  } else if (name == "transit") {
    return Observable::Transit;
  }
  throw std::runtime_error("GraphRandomWalk: unknown observable '" + name +
                           "', use first or transit");
}

GraphRandomWalk::Transient GraphRandomWalk::FindTransient() const {
  const Index n = graph_.NumberOfSites();
  const std::vector<ToFGraph::Site>& sites = graph_.Sites();
  Transient t;
  t.index_of.assign(n, -1);
  // transient sites: not a sink, sink reachable
  for (Index i = 0; i < n; ++i) {
    if (sites[i].role != ToFGraph::Role::Sink && sites[i].reaches_sink) {
      t.index_of[i] = Index(t.site_of.size());
      t.site_of.push_back(i);
    }
  }
  const Index m = Index(t.site_of.size());
  if (m == 0) {
    throw std::runtime_error("GraphRandomWalk: no transient sites");
  }
  for (Index i : t.site_of) {
    if (!(graph_.EscapeRate(i) > 0.0)) {
      throw std::runtime_error("GraphRandomWalk: site " + std::to_string(i) +
                               " has no hops");
    }
    // a carrier can enter a site from which no sink is reachable
    ToFGraph::EdgeRange out = graph_.OutEdges(i);
    const ToFGraph::Edge* trap =
        std::find_if(out.begin(), out.end(), [&](const ToFGraph::Edge& e) {
          return e.rate > 0.0 && sites[e.to].role != ToFGraph::Role::Sink &&
                 t.index_of[e.to] < 0;
        });
    if (trap != out.end()) {
      throw std::runtime_error(
          "GraphRandomWalk: site " + std::to_string(trap->to) +
          " traps carriers (no path to a sink), the time of flight is "
          "infinite");
    }
  }
  // symmetrization weights, only with detailed balance
  if (graph_.isReversible()) {
    Eigen::VectorXd U(m);
    for (Index k = 0; k < m; ++k) {
      U[k] = sites[t.site_of[k]].potential;
    }
    std::vector<double> sorted(U.data(), U.data() + m);
    std::nth_element(sorted.begin(), sorted.begin() + m / 2, sorted.end());
    double Uref = sorted[m / 2];
    t.half_log_weight = -(U.array() - Uref) / (2.0 * graph_.Temperature());
  }
  return t;
}

GraphRandomWalk::Result GraphRandomWalk::SingleCarrier(Solver solver) const {
  const Index n = graph_.NumberOfSites();
  const std::vector<ToFGraph::Site>& sites = graph_.Sites();
  Transient tr = FindTransient();
  const Index m = Index(tr.site_of.size());

  // Q = diag(D) - W on the transient sites
  RowBuilder rows(m);
  for (Index t = 0; t < m; ++t) {
    Index i = tr.site_of[t];
    rows.Add(t, graph_.EscapeRate(i));
    for (const ToFGraph::Edge& e : graph_.OutEdges(i)) {
      if (e.rate > 0.0 && sites[e.to].role != ToFGraph::Role::Sink) {
        rows.Add(tr.index_of[e.to], -e.rate);
      }
    }
    rows.FinishRow();
  }
  SpMatR Q = rows.Build(m);
  GeneratorSolver gs(Q, tr.half_log_weight, solver, settings_);

  Eigen::VectorXd tau = gs.Solve(Eigen::VectorXd::Ones(m));
  Eigen::VectorXd m2 = gs.Solve(2.0 * tau);
  const std::vector<Index>& sources = graph_.SourceSites();
  Eigen::VectorXd w = Eigen::VectorXd::Zero(m);
  for (Index s : sources) {
    w[tr.index_of[s]] = 1.0 / double(sources.size());
  }
  Eigen::VectorXd g = gs.SolveAdjoint(w);

  Result result;
  result.solver = gs.name;
  result.carriers = 1;
  result.unknowns = m;
  result.backward_error = gs.backward_error;
  result.refinement_steps = gs.refinement_steps;
  result.iterations = gs.iterations;
  result.converged = gs.converged && gs.backward_error < 1e-10;

  Index ns = Index(sources.size());
  result.tau.resize(ns);
  result.tau_std.resize(ns);
  double mean2 = 0.0;
  for (Index k = 0; k < ns; ++k) {
    Index t = tr.index_of[sources[k]];
    result.sources.push_back({sources[k]});
    result.tau[k] = tau[t];
    result.tau_std[k] = std::sqrt(std::max(m2[t] - tau[t] * tau[t], 0.0));
    mean2 += m2[t];
  }
  result.tau_mean = result.tau.mean();
  mean2 /= double(ns);
  result.tau_std_injection =
      std::sqrt(std::max(mean2 - result.tau_mean * result.tau_mean, 0.0));
  result.tau_harmonic = double(ns) / result.tau.cwiseInverse().sum();

  result.occupation_time = Eigen::VectorXd::Zero(n);
  for (Index t = 0; t < m; ++t) {
    result.occupation_time[tr.site_of[t]] = g[t];
  }
  return result;
}

GraphRandomWalk::Result GraphRandomWalk::TwoCarriers(Observable observable,
                                                     Solver solver) const {
  const Index n = graph_.NumberOfSites();
  const std::vector<ToFGraph::Site>& sites = graph_.Sites();
  Transient tr = FindTransient();
  const Index m = Index(tr.site_of.size());
  if (m < 2) {
    throw std::runtime_error(
        "GraphRandomWalk: two carriers need at least two transient sites");
  }
  const Index M = m * (m - 1) / 2;
  if (M > settings_.max_states || M > Index(std::numeric_limits<int>::max())) {
    throw std::runtime_error(
        "GraphRandomWalk: " + std::to_string(M) +
        " pair states exceed the limit of " +
        std::to_string(settings_.max_states) +
        " (maxstates); use kinetic Monte Carlo for this system");
  }
  const std::vector<Index>& sources = graph_.SourceSites();
  if (sources.size() < 2) {
    throw std::runtime_error(
        "GraphRandomWalk: two carriers need at least two source sites");
  }
  const bool transit = (observable == Observable::Transit);

  // pair block Q_PP; for Transit also W_PS (pair -> remaining single carrier)
  RowBuilder pp(M);
  RowBuilder ps(M);
  for (Index j = 1; j < m; ++j) {
    for (Index i = 0; i < j; ++i) {
      double D = 0.0;
      const std::array<std::array<Index, 2>, 2> movers = {
          {{i, j}, {j, i}}};  // {moving carrier, other carrier}
      for (const auto& mv : movers) {
        Index a = mv[0];
        Index b = mv[1];
        for (const ToFGraph::Edge& e : graph_.OutEdges(tr.site_of[a])) {
          if (e.rate <= 0.0 || e.to == tr.site_of[b]) {
            continue;  // blocked: site occupied by the other carrier
          }
          D += e.rate;
          if (sites[e.to].role == ToFGraph::Role::Sink) {
            if (transit) {
              ps.Add(b, e.rate);
            }
          } else {
            pp.Add(PairRank(tr.index_of[e.to], b), -e.rate);
          }
        }
      }
      pp.Add(PairRank(i, j), D);
      pp.FinishRow();
      ps.FinishRow();
    }
  }
  SpMatR Q = pp.Build(M);

  // symmetrization weights of a pair: product of the site weights
  Eigen::VectorXd h;
  if (tr.half_log_weight.size() == m) {
    h.resize(M);
    for (Index j = 1; j < m; ++j) {
      for (Index i = 0; i < j; ++i) {
        h[PairRank(i, j)] = tr.half_log_weight[i] + tr.half_log_weight[j];
      }
    }
  }
  GeneratorSolver gs(Q, h, solver, settings_);

  // injection: uniform over pairs of distinct source sites
  std::vector<Index> src(sources.size());
  std::transform(sources.begin(), sources.end(), src.begin(),
                 [&tr](Index s) { return tr.index_of[s]; });
  std::vector<std::array<Index, 2>> pairs;
  for (std::size_t a = 0; a < src.size(); ++a) {
    for (std::size_t b = a + 1; b < src.size(); ++b) {
      pairs.push_back({src[a], src[b]});
    }
  }
  Eigen::VectorXd w = Eigen::VectorXd::Zero(M);
  for (const auto& pr : pairs) {
    w[PairRank(pr[0], pr[1])] = 1.0 / double(pairs.size());
  }

  Result result;
  result.carriers = 2;
  result.observable = observable;
  result.unknowns = M;
  Eigen::VectorXd tau;
  Eigen::VectorXd m2;
  Eigen::VectorXd g;
  Eigen::VectorXd g_single;
  std::string single_solver;
  double single_error = 0.0;
  if (!transit) {
    tau = gs.Solve(Eigen::VectorXd::Ones(M));
    m2 = gs.Solve(2.0 * tau);
    g = gs.SolveAdjoint(w);
  } else {
    // remaining single carrier: Q1 tau_S = 1/2
    RowBuilder rows(m);
    for (Index t = 0; t < m; ++t) {
      Index i = tr.site_of[t];
      rows.Add(t, graph_.EscapeRate(i));
      for (const ToFGraph::Edge& e : graph_.OutEdges(i)) {
        if (e.rate > 0.0 && sites[e.to].role != ToFGraph::Role::Sink) {
          rows.Add(tr.index_of[e.to], -e.rate);
        }
      }
      rows.FinishRow();
    }
    SpMatR Q1 = rows.Build(m);
    SpMatR W = ps.Build(m);
    GeneratorSolver gs1(Q1, tr.half_log_weight, solver, settings_);
    Eigen::VectorXd tau_s = gs1.Solve(Eigen::VectorXd::Constant(m, 0.5));
    tau = gs.Solve(Eigen::VectorXd::Ones(M) + W * tau_s);
    g = gs.SolveAdjoint(w);
    // expected number of transitions into each single-carrier state
    Eigen::VectorXd inflow = W.transpose() * g;
    g_single = gs1.SolveAdjoint(inflow);
    single_solver = gs1.name;
    single_error = gs1.backward_error;
    result.iterations += gs1.iterations;
    result.converged = gs1.converged;
  }
  result.solver = gs.name;
  if (transit && single_solver != gs.name) {
    result.solver += " (remaining carrier: " + single_solver + ")";
  }
  result.backward_error = std::max(gs.backward_error, single_error);
  result.refinement_steps = gs.refinement_steps;
  result.iterations += gs.iterations;
  result.converged =
      result.converged && gs.converged && result.backward_error < 1e-10;

  Index np = Index(pairs.size());
  result.tau.resize(np);
  if (!transit) {
    result.tau_std.resize(np);
  }
  double mean2 = 0.0;
  for (Index k = 0; k < np; ++k) {
    Index r = PairRank(pairs[k][0], pairs[k][1]);
    result.sources.push_back(
        {tr.site_of[pairs[k][0]], tr.site_of[pairs[k][1]]});
    result.tau[k] = tau[r];
    if (!transit) {
      result.tau_std[k] = std::sqrt(std::max(m2[r] - tau[r] * tau[r], 0.0));
      mean2 += m2[r];
    }
  }
  result.tau_mean = result.tau.mean();
  if (!transit) {
    mean2 /= double(np);
    result.tau_std_injection =
        std::sqrt(std::max(mean2 - result.tau_mean * result.tau_mean, 0.0));
  }
  result.tau_harmonic = double(np) / result.tau.cwiseInverse().sum();

  // time each site is occupied: sum over the pair states containing it,
  // plus the single-carrier states for Transit
  result.occupation_time = Eigen::VectorXd::Zero(n);
  for (Index j = 1; j < m; ++j) {
    for (Index i = 0; i < j; ++i) {
      double t = g[PairRank(i, j)];
      result.occupation_time[tr.site_of[i]] += t;
      result.occupation_time[tr.site_of[j]] += t;
    }
  }
  if (transit) {
    for (Index t = 0; t < m; ++t) {
      result.occupation_time[tr.site_of[t]] += g_single[t];
    }
  }
  return result;
}

}  // namespace xtp
}  // namespace votca
