#include "sigma_ri_reduced.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

Sigma_RI_Reduced::Sigma_RI_Reduced(TCMatrix_gwbse& Mmn, RPA& rpa)
    : Sigma_base(Mmn, rpa), rpa_red_(Mmn) {}

void Sigma_RI_Reduced::PrepareScreening() {
  rpa_red_.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa_red_.configure_qp_window(opt_.qpmin, opt_.qpmax);
  rpa_red_.setRPAInputEnergies(rpa_.getRPAInputEnergies());

  RPA_RI_Reduced::options ropt;
  ropt.imag_omega_max = opt_red_.imag_omega_max;
  ropt.imag_omega_points = opt_red_.imag_omega_points;
  ropt.basis_threshold = opt_red_.basis_threshold;
  ropt.max_rank = opt_red_.max_rank;
  ropt.sigma_aware_basis = opt_red_.sigma_aware_basis;
  ropt.sigma_mix = opt_red_.sigma_mix;
  ropt.normalize_metric_components = opt_red_.normalize_metric_components;
  rpa_red_.configure_reduced(ropt);

  rpa_red_.BuildReducedBasis();

  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  reduced_couplings_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    reduced_couplings_[gw_level] =
        Mmn_[gw_level + qpoffset] * rpa_red_.BasisU();
  }

  BuildReducedTransitionModel();
  BuildPoleExpansion();

  if (opt_red_.run_pole_diagnostics) {
    const auto diags = RunPoleDiagnostics();

    double worst_rel2 = 0.0;
    double worst_rel4 = 0.0;
    for (const auto& d : diags) {
      worst_rel2 = std::max(worst_rel2, d.rel_diff_pole2);
      worst_rel4 = std::max(worst_rel4, d.rel_diff_pole4);
    }

    std::ostringstream oss;
    oss << std::setprecision(10);
    oss << "[RI-REDUCED POLE DIAGNOSTIC] rank = " << rpa_red_.rank()
        << "  npoles = " << rpa_omegas_.size() << "\n";
    for (const auto& d : diags) {
      oss << "  omega = " << d.omega
          << "  ||Wdir||_F = " << d.norm_direct
          << "  rel(2Ω) = " << d.rel_diff_pole2
          << "  rel(4Ω) = " << d.rel_diff_pole4 << "\n";
    }
    std::cout << oss.str() << std::flush;

    const double best_rel = std::min(worst_rel2, worst_rel4);
    if (best_rel > opt_red_.pole_reconstruction_tol) {
      std::ostringstream err;
      err << "Sigma_RI_Reduced: reduced pole reconstruction failed. "
          << "worst rel(2Ω) = " << worst_rel2
          << ", worst rel(4Ω) = " << worst_rel4
          << ". The extracted poles do not yet reproduce BuildReducedWcImag.";
      throw std::runtime_error(err.str());
    }
  }

  const Index npoles = rpa_omegas_.size();

  // Raw Wc amplitudes
  wc_amplitudes_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    wc_amplitudes_[gw_level] = reduced_couplings_[gw_level] * pole_vectors_;
    if (wc_amplitudes_[gw_level].cols() != npoles) {
      throw std::runtime_error("Sigma_RI_Reduced: wc_amplitudes shape mismatch.");
    }
  }

  // Sigma-style residues
  residues_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    residues_[gw_level] = wc_amplitudes_[gw_level];

    for (Index s = 0; s < npoles; ++s) {
      const double omega_s = rpa_omegas_(s);
      if (!(omega_s > 0.0)) {
        throw std::runtime_error(
            "Sigma_RI_Reduced: encountered non-positive screening pole.");
      }
      residues_[gw_level].col(s) /= std::sqrt(2.0 * omega_s);
    }
  }

  if (opt_red_.run_contracted_pole_diagnostics) {
    RunContractedWcDiagnostics();
  }
  if (opt_red_.run_full_vs_reduced_contracted_diagnostics) {
    RunFullVsReducedContractedDiagnostics();
  }
  if (opt_red_.run_sigma_diagnostics) {
    RunSigmaDiagnostics();
  }
    if (opt_red_.run_pole_weight_diagnostics) {
    RunPoleWeightDiagnostics();
  }
}

void Sigma_RI_Reduced::RunPoleWeightDiagnostics() const {
  if (rpa_omegas_.size() == 0 || residues_.empty()) {
    std::cout << "[RI-REDUCED POLE-WEIGHT DIAGNOSTIC] no poles or residues available\n";
    return;
  }

  std::vector<Index> gw_levels_to_check;
  gw_levels_to_check.push_back(std::max<Index>(0, opt_.homo - opt_.qpmin));
  if (opt_.homo + 1 >= opt_.qpmin && opt_.homo + 1 <= opt_.qpmax) {
    gw_levels_to_check.push_back(opt_.homo + 1 - opt_.qpmin);
  }

  const Index topn = std::max<Index>(1, opt_red_.pole_weight_topn);

  for (Index gw_level : gw_levels_to_check) {
    const Eigen::MatrixXd& R = residues_[gw_level];

    if (R.cols() != rpa_omegas_.size()) {
      throw std::runtime_error(
          "Sigma_RI_Reduced::RunPoleWeightDiagnostics: residue/pole dimension mismatch.");
    }

    struct PoleInfo {
      Index s;
      double omega;
      double weight;
      double max_abs_residue;
      Index max_m;
    };

    std::vector<PoleInfo> poles;
    poles.reserve(static_cast<std::size_t>(rpa_omegas_.size()));

    double total_weight = 0.0;

    for (Index s = 0; s < rpa_omegas_.size(); ++s) {
      const Eigen::VectorXd col = R.col(s);
      const double weight = col.squaredNorm();
      total_weight += weight;

      Index max_m = 0;
      double max_abs_res = 0.0;
      for (Index m = 0; m < col.size(); ++m) {
        const double val = std::abs(col(m));
        if (val > max_abs_res) {
          max_abs_res = val;
          max_m = m;
        }
      }

      poles.push_back({s, rpa_omegas_(s), weight, max_abs_res, max_m});
    }

    std::sort(poles.begin(), poles.end(),
              [](const PoleInfo& a, const PoleInfo& b) {
                return a.weight > b.weight;
              });

    const Index nprint =
        std::min<Index>(topn, static_cast<Index>(poles.size()));

    std::cout << "[RI-REDUCED POLE-WEIGHT DIAGNOSTIC] gw_level = "
              << gw_level + opt_.qpmin
              << "  npoles = " << rpa_omegas_.size()
              << "  total_weight = " << total_weight
              << "\n";

    double cumulative = 0.0;
    for (Index k = 0; k < nprint; ++k) {
      const PoleInfo& p = poles[static_cast<std::size_t>(k)];
      cumulative += p.weight;
      const double frac =
          (total_weight > 1e-16) ? (p.weight / total_weight) : 0.0;
      const double cumfrac =
          (total_weight > 1e-16) ? (cumulative / total_weight) : 0.0;

      std::cout << "  rank = " << (k + 1)
                << "  pole = " << p.s
                << "  Omega = " << p.omega
                << "  weight = " << p.weight
                << "  frac = " << frac
                << "  cum_frac = " << cumfrac
                << "  max|R_im| = " << p.max_abs_residue
                << "  at m = " << p.max_m
                << "\n";
    }

    // Also print the poles in energy order for the dominant subset
    std::vector<PoleInfo> dominant_by_energy(
        poles.begin(), poles.begin() + static_cast<std::size_t>(nprint));
    std::sort(dominant_by_energy.begin(), dominant_by_energy.end(),
              [](const PoleInfo& a, const PoleInfo& b) {
                return a.omega < b.omega;
              });

    std::cout << "  dominant poles in ascending Omega:\n";
    for (const PoleInfo& p : dominant_by_energy) {
      const double frac =
          (total_weight > 1e-16) ? (p.weight / total_weight) : 0.0;
      std::cout << "    pole = " << p.s
                << "  Omega = " << p.omega
                << "  weight = " << p.weight
                << "  frac = " << frac
                << "  max|R_im| = " << p.max_abs_residue
                << "  at m = " << p.max_m
                << "\n";
    }
  }
}

void Sigma_RI_Reduced::BuildReducedTransitionModel() {
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index n_trans = n_occ * n_unocc;
  const Eigen::VectorXd& eps = rpa_.getRPAInputEnergies();
  const Eigen::MatrixXd& U = rpa_red_.BasisU();
  const Index rank = U.cols();

  transition_energies_.resize(n_trans);
  transition_couplings_.resize(n_trans, rank);

  Index t = 0;
  for (Index v = 0; v < n_occ; ++v) {
    const double eps_v = eps(v);

    const Eigen::MatrixXd Mvc_red =
        Mmn_[v].middleRows(n_occ, n_unocc) * U;

    for (Index c = 0; c < n_unocc; ++c) {
      const double delta = eps(n_occ + c) - eps_v;
      transition_energies_(t) = delta;
      transition_couplings_.row(t) = Mvc_red.row(c);
      ++t;
    }
  }

  std::vector<Index> order(static_cast<std::size_t>(n_trans));
  for (Index i = 0; i < n_trans; ++i) {
    order[static_cast<std::size_t>(i)] = i;
  }

  std::sort(order.begin(), order.end(),
            [&](Index a, Index b) { return transition_energies_(a) < transition_energies_(b); });

  Eigen::VectorXd delta_sorted(n_trans);
  Eigen::MatrixXd G_sorted(n_trans, rank);

  for (Index i = 0; i < n_trans; ++i) {
    const Index old = order[static_cast<std::size_t>(i)];
    delta_sorted(i) = transition_energies_(old);
    G_sorted.row(i) = transition_couplings_.row(old);
  }

  transition_energies_ = std::move(delta_sorted);
  transition_couplings_ = std::move(G_sorted);
}

Eigen::MatrixXd Sigma_RI_Reduced::BuildAReal(double omega) const {
  const Index rank = transition_couplings_.cols();
  Eigen::MatrixXd A = Eigen::MatrixXd::Identity(rank, rank);

  const double omega2 = omega * omega;

  for (Index t = 0; t < transition_energies_.size(); ++t) {
    const double delta = transition_energies_(t);
    const double denom = omega2 - delta * delta;
    const double coeff = 4.0 * delta / denom;

    const Eigen::VectorXd g = transition_couplings_.row(t).transpose();
    A.noalias() -= coeff * (g * g.transpose());
  }

  return A;
}

Eigen::MatrixXd Sigma_RI_Reduced::BuildARealDerivative(double omega) const {
  const Index rank = transition_couplings_.cols();
  Eigen::MatrixXd dA = Eigen::MatrixXd::Zero(rank, rank);

  const double omega2 = omega * omega;

  for (Index t = 0; t < transition_energies_.size(); ++t) {
    const double delta = transition_energies_(t);
    const double denom = omega2 - delta * delta;
    const double coeff = 8.0 * omega * delta / (denom * denom);

    const Eigen::VectorXd g = transition_couplings_.row(t).transpose();
    dA.noalias() += coeff * (g * g.transpose());
  }

  return dA;
}

Index Sigma_RI_Reduced::CountNegativeEigenvalues(const Eigen::MatrixXd& A) const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  if (es.info() != Eigen::Success) {
    throw std::runtime_error("Sigma_RI_Reduced: failed to diagonalize reduced dielectric matrix.");
  }

  const double scale = std::max(1.0, A.norm());
  const double tol = 1e-10 * scale;

  Index nneg = 0;
  for (Index i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) < -tol) {
      ++nneg;
    }
  }
  return nneg;
}

void Sigma_RI_Reduced::IsolateRootsInInterval(
    double left, double right, Index nneg_left, Index nneg_right,
    std::vector<std::pair<double, double>>& brackets) const {

  struct Node {
    double left;
    double right;
    Index nneg_left;
    Index nneg_right;
  };

  std::vector<Node> stack;
  stack.push_back({left, right, nneg_left, nneg_right});

  const double width_tol = std::max(opt_red_.pole_tol, 1e-12);

  while (!stack.empty()) {
    const Node node = stack.back();
    stack.pop_back();

    if (node.nneg_left == node.nneg_right) {
      continue;
    }

    if (node.nneg_right - node.nneg_left == 1) {
      brackets.emplace_back(node.left, node.right);
      continue;
    }

    const double width = node.right - node.left;
    if (!(width > width_tol)) {
      brackets.emplace_back(node.left, node.right);
      continue;
    }

    const double mid = 0.5 * (node.left + node.right);
    if (!(mid > node.left && mid < node.right)) {
      brackets.emplace_back(node.left, node.right);
      continue;
    }

    const Index nneg_mid = CountNegativeEigenvalues(BuildAReal(mid));
    const bool no_progress =
        (nneg_mid == node.nneg_left && nneg_mid == node.nneg_right);

    if (no_progress) {
      brackets.emplace_back(node.left, node.right);
      continue;
    }

    if (nneg_mid != node.nneg_right) {
      stack.push_back({mid, node.right, nneg_mid, node.nneg_right});
    }
    if (nneg_mid != node.nneg_left) {
      stack.push_back({node.left, mid, node.nneg_left, nneg_mid});
    }
  }
}

double Sigma_RI_Reduced::RefineRoot(double left, double right, Index nneg_left) const {
  double a = left;
  double b = right;

  for (Index iter = 0; iter < opt_red_.max_bisection; ++iter) {
    const double mid = 0.5 * (a + b);

    if (!(mid > a && mid < b) || (b - a) < opt_red_.pole_tol) {
      break;
    }

    const Index nneg_mid = CountNegativeEigenvalues(BuildAReal(mid));

    if (nneg_mid > nneg_left) {
      b = mid;
    } else {
      a = mid;
    }
  }

  auto minabs_eval = [&](double x) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(BuildAReal(x));
    if (es.info() != Eigen::Success) {
      return std::numeric_limits<double>::infinity();
    }
    double v = std::abs(es.eigenvalues()(0));
    for (Index i = 1; i < es.eigenvalues().size(); ++i) {
      v = std::min(v, std::abs(es.eigenvalues()(i)));
    }
    return v;
  };

  const double x1 = a;
  const double x2 = 0.5 * (a + b);
  const double x3 = b;

  double best_x = x2;
  double best_v = minabs_eval(x2);

  const double v1 = minabs_eval(x1);
  if (v1 < best_v) {
    best_v = v1;
    best_x = x1;
  }

  const double v3 = minabs_eval(x3);
  if (v3 < best_v) {
    best_x = x3;
  }

  return best_x;
}

std::vector<double> Sigma_RI_Reduced::UniqueSortedTransitionEnergies() const {
  std::vector<double> bare(static_cast<std::size_t>(transition_energies_.size()));
  for (Index i = 0; i < transition_energies_.size(); ++i) {
    bare[static_cast<std::size_t>(i)] = transition_energies_(i);
  }

  std::sort(bare.begin(), bare.end());

  std::vector<double> uniq;
  uniq.reserve(bare.size());
  for (double x : bare) {
    if (uniq.empty() || std::abs(x - uniq.back()) > opt_red_.pole_shift) {
      uniq.push_back(x);
    }
  }
  return uniq;
}

void Sigma_RI_Reduced::BuildPoleExpansion() {
  const std::vector<double> bare = UniqueSortedTransitionEnergies();
  if (bare.empty()) {
    throw std::runtime_error("Sigma_RI_Reduced: no transition energies available.");
  }

  std::vector<std::pair<double, double>> brackets;

  {
    const double left = std::max(1e-12, bare.front() * 0.5);
    const double right = bare.front() - opt_red_.pole_shift;
    if (right > left) {
      const Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
      const Index nneg_right = CountNegativeEigenvalues(BuildAReal(right));
      IsolateRootsInInterval(left, right, nneg_left, nneg_right, brackets);
    }
  }

  for (std::size_t i = 0; i + 1 < bare.size(); ++i) {
    const double left = bare[i] + opt_red_.pole_shift;
    const double right = bare[i + 1] - opt_red_.pole_shift;
    if (right <= left) {
      continue;
    }

    const Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
    const Index nneg_right = CountNegativeEigenvalues(BuildAReal(right));
    IsolateRootsInInterval(left, right, nneg_left, nneg_right, brackets);
  }

  {
    const double left = bare.back() + opt_red_.pole_shift;
    double right = std::max(2.0 * bare.back(), bare.back() + 1.0);

    const Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
    Index nneg_right = CountNegativeEigenvalues(BuildAReal(right));

    while (nneg_right > 0) {
      right *= 2.0;
      nneg_right = CountNegativeEigenvalues(BuildAReal(right));
      if (right > 1e6) {
        throw std::runtime_error("Sigma_RI_Reduced: failed to bracket high-energy reduced screening poles.");
      }
    }

    IsolateRootsInInterval(left, right, nneg_left, nneg_right, brackets);
  }

  std::vector<double> roots;
  std::vector<Eigen::VectorXd> zvecs;

  roots.reserve(brackets.size());
  zvecs.reserve(brackets.size());

  for (const auto& br : brackets) {
    const double left_i = br.first;
    const double right_i = br.second;
    const Index nneg_left_i = CountNegativeEigenvalues(BuildAReal(left_i));
    const double omega = RefineRoot(left_i, right_i, nneg_left_i);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(BuildAReal(omega));
    if (es.info() != Eigen::Success) {
      throw std::runtime_error("Sigma_RI_Reduced: failed to diagonalize A(omega) at pole.");
    }

    Index imin = 0;
    double minabs = std::abs(es.eigenvalues()(0));
    for (Index i = 1; i < es.eigenvalues().size(); ++i) {
      const double val = std::abs(es.eigenvalues()(i));
      if (val < minabs) {
        minabs = val;
        imin = i;
      }
    }

    Eigen::VectorXd p = es.eigenvectors().col(imin);
    p.normalize();

    const Eigen::MatrixXd dA = BuildARealDerivative(omega);
    const double alpha = p.dot(dA * p);

    if (!(alpha > 0.0)) {
      continue;
    }

    roots.push_back(omega);
    zvecs.push_back(p / std::sqrt(alpha));
  }

  if (roots.empty()) {
    throw std::runtime_error("Sigma_RI_Reduced: no valid reduced screening poles found.");
  }

  std::vector<Index> order(roots.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = static_cast<Index>(i);
  }
  std::sort(order.begin(), order.end(),
            [&](Index a, Index b) { return roots[a] < roots[b]; });

  std::vector<double> roots_unique;
  std::vector<Eigen::VectorXd> zvecs_unique;

  const double root_merge_tol = std::max(10.0 * opt_red_.pole_tol, 1e-8);

  for (Index idx : order) {
    const double omega = roots[idx];
    const Eigen::VectorXd& z = zvecs[idx];

    if (roots_unique.empty() ||
        std::abs(omega - roots_unique.back()) >
            root_merge_tol * std::max(1.0, std::abs(omega))) {
      roots_unique.push_back(omega);
      zvecs_unique.push_back(z);
    } else {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_old(BuildAReal(roots_unique.back()));
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_new(BuildAReal(omega));

      double old_min = std::abs(es_old.eigenvalues()(0));
      for (Index i = 1; i < es_old.eigenvalues().size(); ++i) {
        old_min = std::min(old_min, std::abs(es_old.eigenvalues()(i)));
      }

      double new_min = std::abs(es_new.eigenvalues()(0));
      for (Index i = 1; i < es_new.eigenvalues().size(); ++i) {
        new_min = std::min(new_min, std::abs(es_new.eigenvalues()(i)));
      }

      if (new_min < old_min) {
        roots_unique.back() = omega;
        zvecs_unique.back() = z;
      }
    }
  }

  const Index npoles = static_cast<Index>(roots_unique.size());
  rpa_omegas_.resize(npoles);
  pole_vectors_.resize(rpa_red_.rank(), npoles);

  for (Index s = 0; s < npoles; ++s) {
    rpa_omegas_(s) = roots_unique[static_cast<std::size_t>(s)];
    pole_vectors_.col(s) = zvecs_unique[static_cast<std::size_t>(s)];
  }

  std::cout << "[RI-REDUCED POLES] raw_brackets = " << brackets.size()
            << "  unique_poles = " << npoles
            << "  rank = " << rpa_red_.rank() << "\n";
}

Eigen::MatrixXd Sigma_RI_Reduced::BuildReducedWcImagFromPoles(
    double omega, double prefactor) const {
  const Index rank = rpa_red_.rank();
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(rank, rank);

  const double omega2 = omega * omega;
  for (Index s = 0; s < rpa_omegas_.size(); ++s) {
    const double Omega = rpa_omegas_(s);
    const double coeff = -prefactor * Omega / (omega2 + Omega * Omega);
    const Eigen::VectorXd z = pole_vectors_.col(s);
    W.noalias() += coeff * (z * z.transpose());
  }

  return W;
}

std::vector<Sigma_RI_Reduced::pole_diagnostic>
Sigma_RI_Reduced::RunPoleDiagnostics() const {
  std::vector<double> grid;

  const Index npts = std::max<Index>(3, opt_red_.imag_omega_points);
  grid.reserve(static_cast<std::size_t>(npts));

  if (npts == 1) {
    grid.push_back(0.0);
  } else {
    const double denom = static_cast<double>(npts - 1);
    for (Index i = 0; i < npts; ++i) {
      grid.push_back(opt_red_.imag_omega_max * static_cast<double>(i) / denom);
    }
  }

  std::vector<pole_diagnostic> out;
  out.reserve(grid.size());

  for (double omega : grid) {
    const Eigen::MatrixXd W_direct = rpa_red_.BuildReducedWcImag(omega);
    const Eigen::MatrixXd W_pole2 = BuildReducedWcImagFromPoles(omega, 2.0);
    const Eigen::MatrixXd W_pole4 = BuildReducedWcImagFromPoles(omega, 4.0);

    pole_diagnostic d;
    d.omega = omega;

    d.norm_direct = W_direct.norm();

    d.norm_pole2 = W_pole2.norm();
    d.abs_diff_pole2 = (W_direct - W_pole2).norm();
    d.rel_diff_pole2 =
        (d.norm_direct > 1e-14) ? (d.abs_diff_pole2 / d.norm_direct)
                                : d.abs_diff_pole2;

    d.norm_pole4 = W_pole4.norm();
    d.abs_diff_pole4 = (W_direct - W_pole4).norm();
    d.rel_diff_pole4 =
        (d.norm_direct > 1e-14) ? (d.abs_diff_pole4 / d.norm_direct)
                                : d.abs_diff_pole4;

    out.push_back(d);
  }

  return out;
}

double Sigma_RI_Reduced::CalcCorrelationDiagElement(Index gw_level,
                                                    double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;

  double sigma = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); ++s) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();

    Eigen::ArrayXd temp = -rpa_.getRPAInputEnergies().array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;

    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    sigma += (res_12 * temp / denom).sum();
  }

  return 2.0 * sigma;
}

double Sigma_RI_Reduced::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;

  double dsigma_domega = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); ++s) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();

    Eigen::ArrayXd temp = -rpa_.getRPAInputEnergies().array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;

    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    dsigma_domega += ((eta2 - temp.abs2()) * res_12 / denom.abs2()).sum();
  }

  return 2.0 * dsigma_domega;
}

double Sigma_RI_Reduced::CalcCorrelationOffDiagElement(Index gw_level1,
                                                       Index gw_level2,
                                                       double frequency1,
                                                       double frequency2) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;

  double sigma_c = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); ++s) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::VectorXd& res1 = residues_[gw_level1].col(s);
    const Eigen::VectorXd& res2 = residues_[gw_level2].col(s);
    const Eigen::VectorXd res_12 = res1.cwiseProduct(res2);

    Eigen::ArrayXd temp1 = -rpa_.getRPAInputEnergies().array();
    temp1.segment(0, n_occ) += eigenvalue;
    temp1.segment(n_occ, n_unocc) -= eigenvalue;

    const Eigen::ArrayXd temp2 = temp1 + frequency2;
    temp1 += frequency1;

    const Eigen::ArrayXd numer1 = res_12.array() * temp1;
    const Eigen::ArrayXd numer2 = res_12.array() * temp2;
    const Eigen::ArrayXd denom1 = temp1.abs2() + eta2;
    const Eigen::ArrayXd denom2 = temp2.abs2() + eta2;

    sigma_c += 0.5 * ((numer1 / denom1) + (numer2 / denom2)).sum();
  }

  return 2.0 * sigma_c;
}

double Sigma_RI_Reduced::ContractedFullWcImag(Index gw_level, Index m,
                                              double omega) const {
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd Wfull = rpa_red_.BuildWcImag(omega);
  const Eigen::VectorXd cfull = Mmn_[gw_level + qpoffset].row(m).transpose();
  return cfull.dot(Wfull * cfull);
}

double Sigma_RI_Reduced::ContractedProjectedReducedWcImag(Index gw_level, Index m,
                                                          double omega) const {
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd Wproj = rpa_red_.BuildProjectedReducedWcImag(omega);
  const Eigen::VectorXd cfull = Mmn_[gw_level + qpoffset].row(m).transpose();
  return cfull.dot(Wproj * cfull);
}

double Sigma_RI_Reduced::ContractedDirectWcImag(Index gw_level, Index m,
                                                double omega) const {
  const Eigen::MatrixXd Wred = rpa_red_.BuildReducedWcImag(omega);
  const Eigen::VectorXd c = reduced_couplings_[gw_level].row(m).transpose();
  return c.dot(Wred * c);
}

double Sigma_RI_Reduced::ContractedPoleWcImag(Index gw_level, Index m,
                                              double omega) const {
  double val = 0.0;
  const double omega2 = omega * omega;

  for (Index s = 0; s < rpa_omegas_.size(); ++s) {
    const double Omega = rpa_omegas_(s);
    const double Aim = wc_amplitudes_[gw_level](m, s);
    val += -2.0 * Omega * Aim * Aim / (omega2 + Omega * Omega);
  }

  return val;
}

void Sigma_RI_Reduced::RunContractedWcDiagnostics() const {
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpatotal = n_occ + n_unocc;

  std::vector<Index> gw_levels_to_check;
  gw_levels_to_check.push_back(std::max<Index>(0, opt_.homo - opt_.qpmin));
  if (opt_.homo + 1 >= opt_.qpmin && opt_.homo + 1 <= opt_.qpmax) {
    gw_levels_to_check.push_back(opt_.homo + 1 - opt_.qpmin);
  }

  std::vector<Index> m_to_check = {0, n_occ - 1, n_occ, rpatotal - 1};
  std::vector<double> omegas = {0.0, 0.5, 1.0};

  for (Index gw_level : gw_levels_to_check) {
    std::cout << "[RI-REDUCED CONTRACTED-W DIAGNOSTIC] gw_level = "
              << gw_level + opt_.qpmin << "\n";

    for (Index m : m_to_check) {
      if (m < 0 || m >= rpatotal) {
        continue;
      }

      for (double omega : omegas) {
        const double direct = ContractedDirectWcImag(gw_level, m, omega);
        const double pole = ContractedPoleWcImag(gw_level, m, omega);
        const double abs_diff = std::abs(direct - pole);
        const double rel_diff =
            (std::abs(direct) > 1e-14) ? abs_diff / std::abs(direct) : abs_diff;

        std::cout << "  m = " << m
                  << "  omega = " << omega
                  << "  direct = " << direct
                  << "  pole = " << pole
                  << "  abs_diff = " << abs_diff
                  << "  rel_diff = " << rel_diff << "\n";
      }
    }
  }
}

void Sigma_RI_Reduced::RunFullVsReducedContractedDiagnostics() const {
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpatotal = n_occ + n_unocc;

  std::vector<Index> gw_levels_to_check;
  gw_levels_to_check.push_back(std::max<Index>(0, opt_.homo - opt_.qpmin));
  if (opt_.homo + 1 >= opt_.qpmin && opt_.homo + 1 <= opt_.qpmax) {
    gw_levels_to_check.push_back(opt_.homo + 1 - opt_.qpmin);
  }

  std::vector<Index> m_to_check = {0, n_occ - 1, n_occ, rpatotal - 1};
  std::vector<double> omegas = {0.0, 0.5, 1.0};

  for (Index gw_level : gw_levels_to_check) {
    std::cout << "[RI-REDUCED FULL-vs-PROJECTED CONTRACTED-W DIAGNOSTIC] gw_level = "
              << gw_level + opt_.qpmin << "\n";

    for (Index m : m_to_check) {
      if (m < 0 || m >= rpatotal) {
        continue;
      }

      for (double omega : omegas) {
        const double full = ContractedFullWcImag(gw_level, m, omega);
        const double proj = ContractedProjectedReducedWcImag(gw_level, m, omega);
        const double red = ContractedDirectWcImag(gw_level, m, omega);

        const double abs_diff_proj = std::abs(full - proj);
        const double rel_diff_proj =
            (std::abs(full) > 1e-14) ? abs_diff_proj / std::abs(full) : abs_diff_proj;

        const double abs_diff_red = std::abs(full - red);
        const double rel_diff_red =
            (std::abs(full) > 1e-14) ? abs_diff_red / std::abs(full) : abs_diff_red;

        std::cout << "  m = " << m
                  << "  omega = " << omega
                  << "  full = " << full
                  << "  proj = " << proj
                  << "  red = " << red
                  << "  abs(full-proj) = " << abs_diff_proj
                  << "  rel(full-proj) = " << rel_diff_proj
                  << "  abs(full-red) = " << abs_diff_red
                  << "  rel(full-red) = " << rel_diff_red
                  << "\n";
      }
    }
  }
}

double Sigma_RI_Reduced::CalcCorrelationDiagElementDirectReduced(
    Index gw_level, double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpatotal = n_occ + n_unocc;

  double sigma = 0.0;

  // simple imaginary-axis quadrature for diagnostics only
  const Index nw = std::max<Index>(16, opt_red_.imag_omega_points);
  const double wmax = std::max(4.0, opt_red_.imag_omega_max);

  for (Index iw = 0; iw < nw; ++iw) {
    const double x = static_cast<double>(iw) / static_cast<double>(nw - 1);
    const double omega = wmax * x;
    const double weight =
        (iw == 0 || iw == nw - 1) ? 0.5 * wmax / static_cast<double>(nw - 1)
                                  :       wmax / static_cast<double>(nw - 1);

    const Eigen::MatrixXd Wred = rpa_red_.BuildReducedWcImag(omega);

    for (Index m = 0; m < rpatotal; ++m) {
      const Eigen::VectorXd c = reduced_couplings_[gw_level].row(m).transpose();
      const double Wim = c.dot(Wred * c);

      double denom = frequency - rpa_.getRPAInputEnergies()(m);
      if (m < n_occ) {
        denom += omega;
      } else {
        denom -= omega;
      }

      sigma += weight * Wim * denom / (denom * denom + eta2);
    }
  }

  return sigma / M_PI;
}

double Sigma_RI_Reduced::CalcCorrelationOffDiagElementDirectReduced(
    Index gw_level1, Index gw_level2, double frequency1, double frequency2) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpatotal = n_occ + n_unocc;

  double sigma = 0.0;

  const Index nw = std::max<Index>(16, opt_red_.imag_omega_points);
  const double wmax = std::max(4.0, opt_red_.imag_omega_max);

  for (Index iw = 0; iw < nw; ++iw) {
    const double x = static_cast<double>(iw) / static_cast<double>(nw - 1);
    const double omega = wmax * x;
    const double weight =
        (iw == 0 || iw == nw - 1) ? 0.5 * wmax / static_cast<double>(nw - 1)
                                  :       wmax / static_cast<double>(nw - 1);

    const Eigen::MatrixXd Wred = rpa_red_.BuildReducedWcImag(omega);

    for (Index m = 0; m < rpatotal; ++m) {
      const Eigen::VectorXd c1 = reduced_couplings_[gw_level1].row(m).transpose();
      const Eigen::VectorXd c2 = reduced_couplings_[gw_level2].row(m).transpose();
      const double Wim = c1.dot(Wred * c2);

      double d1 = frequency1 - rpa_.getRPAInputEnergies()(m);
      double d2 = frequency2 - rpa_.getRPAInputEnergies()(m);

      if (m < n_occ) {
        d1 += omega;
        d2 += omega;
      } else {
        d1 -= omega;
        d2 -= omega;
      }

      sigma += weight * 0.5 *
               (Wim * d1 / (d1 * d1 + eta2) + Wim * d2 / (d2 * d2 + eta2));
    }
  }

  return sigma / M_PI;
}

void Sigma_RI_Reduced::RunSigmaDiagnostics() const {
  std::vector<Index> gw_levels_to_check;
  gw_levels_to_check.push_back(std::max<Index>(0, opt_.homo - opt_.qpmin));
  if (opt_.homo + 1 >= opt_.qpmin && opt_.homo + 1 <= opt_.qpmax) {
    gw_levels_to_check.push_back(opt_.homo + 1 - opt_.qpmin);
  }

  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  const Eigen::VectorXd& eps = rpa_.getRPAInputEnergies();

  for (Index gw_level : gw_levels_to_check) {
    const double freq = eps(gw_level + qpoffset);

    const double sigma_pole = CalcCorrelationDiagElement(gw_level, freq);
    const double sigma_direct = CalcCorrelationDiagElementDirectReduced(gw_level, freq);

    const double abs_diff = std::abs(sigma_pole - sigma_direct);
    const double rel_diff =
        (std::abs(sigma_direct) > 1e-14) ? abs_diff / std::abs(sigma_direct) : abs_diff;

    std::cout << "[RI-REDUCED SIGMA DIAGNOSTIC] gw_level = "
              << gw_level + opt_.qpmin
              << "  freq = " << freq
              << "  sigma_pole = " << sigma_pole
              << "  sigma_direct = " << sigma_direct
              << "  abs_diff = " << abs_diff
              << "  rel_diff = " << rel_diff
              << "\n";
  }
}

}  // namespace xtp
}  // namespace votca