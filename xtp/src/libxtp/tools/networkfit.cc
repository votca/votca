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
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/bulksimulation.h"
#include "votca/xtp/statesaver.h"
#include "votca/xtp/tofgraph.h"
#include "votca/xtp/trapfinder.h"

// Local private VOTCA includes
#include "networkfit.h"

namespace votca {
namespace xtp {

namespace {
std::vector<Index> Grid(const std::vector<Index>& g, const std::string& name) {
  if (g.size() != 3 || g[0] < 1 || g[1] < g[0] || g[2] < 1) {
    throw std::runtime_error("networkfit: " + name +
                             " needs first, last, step with 1 <= first <= "
                             "last and step >= 1");
  }
  std::vector<Index> r;
  for (Index l = g[0]; l <= g[1]; l += g[2]) {
    r.push_back(l);
  }
  return r;
}

double MeanOf(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0) / double(x.size());
}

/// value of a binned curve at distance r (bohr), NaN if the bin is empty
double At(const StochasticNetworkFit::Curve& c, double r) {
  for (std::size_t k = 0; k < c.r.size(); ++k) {
    double half = k + 1 < c.r.size() ? 0.5 * (c.r[k + 1] - c.r[k])
                                     : 0.5 * (c.r[k] - c.r[k - 1]);
    if (std::abs(r - c.r[k]) <= half) {
      return (c.count.empty() || c.count[k] > 0.0)
                 ? c.value[k]
                 : std::numeric_limits<double>::quiet_NaN();
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

struct Log10Stats {
  double mean = 0.0;
  double sd = 0.0;
};
Log10Stats LogJ(const StochasticNetworkFit::Data& d, Index c) {
  std::vector<double> y;
  for (const auto& p : d.pairs) {
    if (p.J2[c] > 0.0) {
      y.push_back(
          std::log10(p.J2[c] * tools::conv::hrt2ev * tools::conv::hrt2ev));
    }
  }
  Log10Stats s;
  if (y.empty()) {
    return s;
  }
  s.mean = MeanOf(y);
  double v = 0.0;
  for (double x : y) {
    v += (x - s.mean) * (x - s.mean);
  }
  s.sd = std::sqrt(v / double(y.size()));
  return s;
}

double Sigma(const std::vector<double>& e) {
  double m = MeanOf(e);
  double v = 0.0;
  for (double x : e) {
    v += (x - m) * (x - m);
  }
  return std::sqrt(v / double(e.size()));
}
}  // namespace

void NetworkFit::ParseOptions(const tools::Property& options) {
  const double nm = tools::conv::nm2bohr;
  statefile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".statefile", job_name_ + ".hdf5");
  frame_ = options.get(".frame").as<Index>();
  parameterfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".parameterfile", job_name_ + "_network.xml");
  energyfile_ =
      options.ifExistsReturnElseReturnDefault<std::string>(".energyfile", "");
  statsfile_ =
      options.ifExistsReturnElseReturnDefault<std::string>(".statsfile", "");
  settings_.carriers.clear();
  for (const std::string& s :
       options.get(".carriers").as<std::vector<std::string>>()) {
    QMStateType type(s);
    if (!type.isKMCState()) {
      throw std::runtime_error("networkfit: " + s +
                               " is not a hopping state (e, h, s, t)");
    }
    settings_.carriers.push_back(type);
  }
  settings_.bin = options.get(".bin").as<double>() * nm;
  settings_.correlation_cutoff =
      options.get(".correlationcutoff").as<double>() * nm;
  settings_.l_b = Grid(options.get(".l_b").as<std::vector<Index>>(), "l_b");
  settings_.l_c = Grid(options.get(".l_c").as<std::vector<Index>>(), "l_c");
  settings_.quantiles = options.get(".quantiles").as<Index>();
  settings_.fit_points = options.get(".points.fit").as<bool>();
  settings_.realizations = options.get(".points.realizations").as<Index>();
  settings_.max_evaluations = options.get(".points.evaluations").as<Index>();
  settings_.candidates = options.get(".points.candidates").as<double>();
  settings_.max_iterations = options.get(".points.maxiterations").as<Index>();
  settings_.points.hard_core =
      options.get(".points.hardcore").as<double>() * nm;
  settings_.points.radius_mean =
      options.get(".points.radiusmean").as<double>() * nm;
  settings_.points.radius_variance =
      options.get(".points.radiusvariance").as<double>() * nm * nm;
  settings_.seed = options.get(".seed").as<std::uint64_t>();
  validate_ = options.get(".validate").as<Index>();
  find_traps_ = options.get(".traps.find").as<bool>();
  traps_.threshold = options.get(".traps.threshold").as<double>();
  traps_.gap = options.get(".traps.gap").as<double>();
  traps_.max_sites = options.get(".traps.maxsites").as<Index>();
  traps_.mass = options.get(".traps.mass").as<double>();
  traps_.max_size = options.get(".traps.maxsize").as<Index>();
  temperature_ = options.get(".temperature").as<double>();
}

bool NetworkFit::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");
  const double nm = tools::conv::bohr2nm;
  const double ev = tools::conv::hrt2ev;
  const double kT = tools::conv::kB * temperature_;

  StateSaver reader(statefile_);
  Index frame = frame_;
  if (frame < 0) {
    std::vector<Index> frames = reader.getFrames();
    if (frames.empty()) {
      throw std::runtime_error("networkfit: " + statefile_ + " has no frames");
    }
    frame = frames.front();
  }
  Topology top = reader.ReadFrame(frame);
  StochasticNetworkFit::Data data =
      StochasticNetworkFit::FromTopology(top, settings_.carriers);
  const Index n = Index(data.positions.size());
  XTP_LOG(Log::error, log_)
      << "Read " << statefile_ << " frame " << frame << ": " << n << " sites, "
      << data.pairs.size() << " pairs, box " << data.box.x() * nm << " x "
      << data.box.y() * nm << " x " << data.box.z() * nm << " nm^3"
      << std::flush;

  // traps of the microscopic system (zero field)
  auto find_traps = [this](const Topology& t, Index c) {
    ToFGraph graph(t, ToFGraph::BulkOptions("", ""), settings_.carriers[c],
                   temperature_ * tools::conv::kB * tools::conv::ev2hrt,
                   Eigen::Vector3d::Zero());
    return TrapFinder::Find(graph, traps_);
  };
  settings_.traps.clear();
  if (find_traps_) {
    for (std::size_t c = 0; c < settings_.carriers.size(); ++c) {
      TrapFinder::Result tr = find_traps(top, Index(c));
      std::ostringstream s;
      s << "Traps for " << settings_.carriers[c].ToString()
        << " (slow modes of the master equation at " << temperature_
        << " K): largest eigenvalue ratio " << std::setprecision(3) << tr.gap
        << " after " << tr.slow_modes << " slow mode(s)";
      if (tr.traps.empty()) {
        s << ": no traps";
      }
      for (const TrapFinder::Trap& t : tr.traps) {
        s << "\n    trap of " << t.sites.size() << " site(s):";
        for (Index i : t.sites) {
          double e = data.energies[c][i];
          double m = 0.0;
          for (double x : data.energies[c]) {
            m += x / double(n);
          }
          s << " " << i << " (" << (e - m) * tools::conv::hrt2ev << " eV)";
        }
        s << "\n      equilibrium weight " << t.weight << "; release rate "
          << t.release_rate << " 1/s, capture rate " << t.capture_rate
          << " 1/s (slow mode " << t.relaxation_rate << " 1/s); hops out "
          << t.hop_rate << " 1/s";
        settings_.traps.push_back({Index(c), t.sites, t.release_rate});
      }
      XTP_LOG(Log::error, log_) << s.str() << std::flush;
    }
  }

  StochasticNetworkFit fitter(data, settings_);
  StochasticNetworkFit::Result res = fitter.Fit();
  const StochasticNetwork::Options& opt = res.options;

  const StochasticNetwork::PointProcess& pp = opt.points;
  if (settings_.fit_points) {
    XTP_LOG(Log::error, log_)
        << "Point process: r_h = " << pp.hard_core * nm
        << " nm, radius mean = " << pp.radius_mean * nm
        << " nm, variance = " << pp.radius_variance * nm * nm
        << " nm^2 (objective " << res.points.objective << ", "
        << res.points.evaluations << " evaluations)" << std::flush;
  }
  XTP_LOG(Log::error, log_)
      << "Connection probability tabulated up to "
      << opt.connection_r.back() * nm << " nm" << std::flush;
  for (std::size_t c = 0; c < opt.carriers.size(); ++c) {
    const auto& e = res.energies[c];
    const auto& car = opt.carriers[c];
    const auto& j = res.couplings[c];
    std::ostringstream s;
    s << "Carrier " << car.state.ToString() << ": site energies mean "
      << e.mean * ev << " eV, sigma " << e.sigma * ev << " eV (sample), "
      << e.sigma_model * ev << " eV = " << e.sigma_model * ev / kT
      << " kT (model, corrected for the finite system)";
    if (e.fitted) {
      s << "\n    moving average w_a = " << e.w_a << ", w_b = " << e.w_b
        << " (l_b = " << e.l_b << "), w_c = " << e.w_c << " (l_c = " << e.l_c
        << "), rms of kappa " << e.rms << ", mean correlation over all pairs "
        << e.cbar;
      if (e.l_b_at_edge) {
        s << "\n    WARNING: l_b is at the end of its grid, extend the grid";
      }
      if (e.l_c_at_edge) {
        s << "\n    WARNING: l_c is at the end of its grid: the range of "
             "the correlation is not determined by this system size (see "
             "the l_c scan in the stats file); extend the grid or use a "
             "larger microscopic system";
      }
    } else {
      s << "; no energy disorder (all site energies equal)";
    }
    s << "\n    log10(J^2/eV^2) = " << j.model.mean_offset << " + "
      << j.model.mean_slope << " r/nm, variance " << j.model.var_offset << " + "
      << j.model.var_slope << " r/nm; residual skewness " << j.skewness
      << ", excess kurtosis " << j.excess_kurtosis << "; " << j.used
      << " couplings";
    if (j.excluded > 0) {
      s << ", " << j.excluded << " zero couplings ignored";
    }
    s << "\n    lambda = " << (car.U_nX_nN + car.U_xN_xX) * ev
      << " eV (U_nX_nN " << car.U_nX_nN * ev << ", U_xN_xX " << car.U_xN_xX * ev
      << "), internal energy U_xX_nN " << car.U_xX_nN * ev << " eV";
    double s2 = (e.sigma_model * ev / kT) * (e.sigma_model * ev / kT);
    s << "\n    the density of occupied states at " << temperature_
      << " K needs about " << std::setprecision(2) << std::scientific
      << std::exp((s2 + 5.7) / 1.05) << " sites" << std::defaultfloat
      << std::setprecision(6);
    XTP_LOG(Log::error, log_) << s.str() << std::flush;
  }
  if (opt.carriers.size() > 1) {
    XTP_LOG(Log::error, log_)
        << "Correlation of the site energies of "
        << opt.carriers[0].state.ToString() << " and "
        << opt.carriers[1].state.ToString() << ": " << opt.energy_correlation
        << ", of the coupling residuals: " << opt.coupling_correlation
        << std::flush;
  }

  {
    std::ofstream out(parameterfile_);
    out << StochasticNetwork::ToXML(opt);
  }
  XTP_LOG(Log::error, log_) << "Wrote " << parameterfile_ << std::flush;

  if (!energyfile_.empty()) {
    const std::size_t ncar = settings_.carriers.size();
    // per carrier: 1 if the site belongs to a trap of that carrier
    std::vector<std::vector<int>> trap_of(ncar, std::vector<int>(n, 0));
    for (const auto& t : settings_.traps) {
      for (Index i : t.sites) {
        trap_of[t.carrier][i] = 1;
      }
    }
    std::ofstream out(energyfile_);
    out << "# segment id, x, y, z [nm] (wrapped), site energies [eV] of";
    for (const auto& c : settings_.carriers) {
      out << " " << c.ToString();
    }
    out << ", trap site of";
    for (const auto& c : settings_.carriers) {
      out << " " << c.ToString();
    }
    out << " (1: yes)\n" << std::setprecision(12);
    for (Index i = 0; i < n; ++i) {
      Eigen::Vector3d p = data.positions[i] * nm;
      out << i << " " << p.x() << " " << p.y() << " " << p.z();
      for (std::size_t c = 0; c < ncar; ++c) {
        out << " " << data.energies[c][i] * ev;
      }
      for (std::size_t c = 0; c < ncar; ++c) {
        out << " " << trap_of[c][i];
      }
      out << "\n";
    }
    XTP_LOG(Log::error, log_) << "Wrote " << energyfile_ << std::flush;
  }

  // validation: model networks of the same size, statistics averaged
  StochasticNetworkFit::Data model;  // the first, for the stats file
  const bool have_model = validate_ > 0;
  if (have_model) {
    const double half = 0.4999 * data.box.minCoeff();
    const double rk = std::min(half, 2.5 * tools::conv::nm2bohr);
    const std::size_t ncar = settings_.carriers.size();
    auto statistics = [&](const StochasticNetworkFit::Data& d) {
      std::vector<std::pair<std::string, double>> rows;
      rows.emplace_back("pairs per site",
                        2.0 * double(d.pairs.size()) / double(n));
      rows.emplace_back(
          "mean nearest neighbour distance (nm)",
          MeanOf(StochasticNetworkFit::NearestDistances(d.positions, d.box)) *
              nm);
      for (std::size_t c = 0; c < ncar; ++c) {
        const std::string name = settings_.carriers[c].ToString();
        double sd = Sigma(d.energies[c]);
        rows.emplace_back("sigma " + name + " (eV)", sd * ev);
        Log10Stats l = LogJ(d, Index(c));
        rows.emplace_back("mean log10(J^2/eV^2) " + name, l.mean);
        rows.emplace_back("sd log10(J^2/eV^2) " + name, l.sd);
        if (res.energies[c].fitted) {
          double md = MeanOf(d.energies[c]);
          std::vector<double> z;
          for (double x : d.energies[c]) {
            z.push_back((x - md) / sd);
          }
          auto k = StochasticNetworkFit::Correlation(d.positions, d.box, z,
                                                     settings_.bin, rk);
          for (double r : {0.75, 1.0, 1.5, 2.0}) {
            std::ostringstream q;
            q << "kappa " << name << " at " << std::fixed
              << std::setprecision(2) << r << " nm";
            rows.emplace_back(q.str(), At(k, r * tools::conv::nm2bohr));
          }
        }
      }
      if (ncar > 1 && res.energies[0].fitted && res.energies[1].fitted) {
        const std::vector<double>& a = d.energies[0];
        const std::vector<double>& b = d.energies[1];
        double ma = MeanOf(a), mb = MeanOf(b);
        double sab = 0.0, saa = 0.0, sbb = 0.0;
        for (Index i = 0; i < n; ++i) {
          sab += (a[i] - ma) * (b[i] - mb);
          saa += (a[i] - ma) * (a[i] - ma);
          sbb += (b[i] - mb) * (b[i] - mb);
        }
        rows.emplace_back("site energy correlation",
                          sab / std::sqrt(saa * sbb));
      }
      return rows;
    };
    // transport and traps of a network: exact zero-field diffusion and the
    // number of traps found, per carrier
    auto trap_rows = [&](const Topology& t,
                         std::vector<std::pair<std::string, double>>& rows) {
      for (std::size_t c = 0; c < ncar; ++c) {
        const std::string name = settings_.carriers[c].ToString();
        // transport: exact zero-field diffusion of one carrier
        ToFGraph graph(t, ToFGraph::BulkOptions("", ""), settings_.carriers[c],
                       temperature_ * tools::conv::kB * tools::conv::ev2hrt,
                       Eigen::Vector3d::Zero());
        BulkSimulation::ExactResult ex = BulkSimulation::Exact(graph);
        if (ex.available) {
          double d = ex.diffusion.trace() / 3.0 * tools::conv::bohr2nm *
                     tools::conv::bohr2nm * 1e-18;  // m^2/s
          rows.emplace_back("log10 D zero field " + name + " (m^2/s)",
                            std::log10(d));
        }
        if (!find_traps_) {
          continue;
        }
        TrapFinder::Result tr = TrapFinder::Find(graph, traps_);
        rows.emplace_back("traps found " + name, double(tr.traps.size()));
      }
    };
    auto micro = statistics(data);
    trap_rows(top, micro);
    std::vector<double> mean(micro.size(), 0.0), var(micro.size(), 0.0);
    Index inserted = 0;
    for (Index k = 0; k < validate_; ++k) {
      StochasticNetwork::Options o = opt;
      o.seed = opt.seed + std::uint64_t(k);
      StochasticNetwork net(o);
      inserted += Index(net.Traps().size());
      StochasticNetworkFit::Data d =
          StochasticNetworkFit::FromNetwork(net, settings_.carriers, data.box);
      auto rows = statistics(d);
      {
        Topology t;
        net.BuildTopology(t);
        trap_rows(t, rows);
      }
      for (std::size_t q = 0; q < rows.size(); ++q) {
        mean[q] += rows[q].second / double(validate_);
        var[q] += rows[q].second * rows[q].second / double(validate_);
      }
      if (k == 0) {
        model = std::move(d);
      }
    }
    std::ostringstream s;
    s << std::setprecision(4);
    s << "Validation: " << validate_ << " model network(s) of " << n
      << " sites (seeds " << opt.seed << " to "
      << opt.seed + std::uint64_t(validate_) - 1
      << "), mean and standard deviation over the networks; " << inserted
      << " trap(s) inserted in total:\n";
    s << "    " << std::left << std::setw(38) << "quantity" << std::setw(14)
      << "microscopic" << std::setw(14) << "model" << "sd\n";
    for (std::size_t q = 0; q < micro.size(); ++q) {
      double sd = std::sqrt(std::max(0.0, var[q] - mean[q] * mean[q]));
      s << "    " << std::left << std::setw(38) << micro[q].first
        << std::setw(14) << micro[q].second << std::setw(14) << mean[q] << sd
        << "\n";
    }
    XTP_LOG(Log::error, log_) << s.str() << std::flush;
  }

  if (!statsfile_.empty()) {
    std::ofstream out(statsfile_);
    out << std::setprecision(8);
    using Curve = StochasticNetworkFit::Curve;
    auto block = [&out](const std::string& title, const std::string& columns,
                        const Curve& a, const Curve* b) {
      out << "# " << title << "\n# " << columns << "\n";
      for (std::size_t k = 0; k < a.r.size(); ++k) {
        bool empty = !a.count.empty() && a.count[k] <= 0.0;
        out << a.r[k] * tools::conv::bohr2nm << " "
            << (empty ? std::numeric_limits<double>::quiet_NaN() : a.value[k]);
        if (b != nullptr) {
          bool bempty = !b->count.empty() && b->count[k] <= 0.0;
          out << " "
              << (bempty ? std::numeric_limits<double>::quiet_NaN()
                         : b->value[k]);
        }
        out << "\n";
      }
      out << "\n\n";
    };
    if (settings_.fit_points) {
      block("nearest neighbour distance distribution D(r)",
            "r [nm], microscopic, model", res.points.D, &res.points.D_model);
      block("spherical contact distribution H(r)", "r [nm], microscopic, model",
            res.points.H, &res.points.H_model);
      block("pair correlation function g(r)", "r [nm], microscopic, model",
            res.points.g, &res.points.g_model);
    }
    block("connection probability f(r)", "r [nm], f", res.connection, nullptr);
    double rmax = res.connection.r.back();
    if (have_model) {
      Curve pd = StochasticNetworkFit::PairsPerSite(data, settings_.bin, rmax);
      Curve pm = StochasticNetworkFit::PairsPerSite(model, settings_.bin, rmax);
      pd.count.clear();
      pm.count.clear();
      block("pairs per site in distance bins", "r [nm], microscopic, model", pd,
            &pm);
    }
    for (std::size_t c = 0; c < settings_.carriers.size(); ++c) {
      const std::string name = settings_.carriers[c].ToString();
      if (res.energies[c].fitted) {
        block("energy correlation kappa(r), carrier " + name,
              "r [nm], microscopic, fitted moving average",
              res.energies[c].kappa, &res.energies[c].kappa_model);
        out << "# rms of the best kappa fit for each l_c, carrier " << name
            << "\n# l_c, rms\n";
        const Curve& scan = res.energies[c].l_c_scan;
        for (std::size_t k = 0; k < scan.r.size(); ++k) {
          out << scan.r[k] << " " << scan.value[k] << "\n";
        }
        out << "\n\n";
      }
      Curve jd = StochasticNetworkFit::MeanLog10J2(data, Index(c),
                                                   settings_.bin, rmax);
      if (have_model) {
        Curve jm = StochasticNetworkFit::MeanLog10J2(model, Index(c),
                                                     settings_.bin, rmax);
        block("mean log10(J^2/eV^2), carrier " + name,
              "r [nm], microscopic, model", jd, &jm);
      } else {
        block("mean log10(J^2/eV^2), carrier " + name, "r [nm], microscopic",
              jd, nullptr);
      }
    }
    XTP_LOG(Log::error, log_) << "Wrote " << statsfile_ << std::flush;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
