#include "sigma_ri_reduced.h"

#include <algorithm>
#include <cmath>
#include <limits>
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
  // Build reduced RI basis U from imaginary-axis Pi(0)
  rpa_red_.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa_red_.setRPAInputEnergies(rpa_.getRPAInputEnergies());

  RPA_RI_Reduced::options ropt;
  ropt.imag_omega_max = opt_red_.imag_omega_max;
  ropt.imag_omega_points = opt_red_.imag_omega_points;
  ropt.basis_threshold = opt_red_.basis_threshold;
  ropt.max_rank = opt_red_.max_rank;
  rpa_red_.configure_reduced(ropt);

  rpa_red_.BuildReducedBasis();

  // Build reduced couplings C_i(m,a) for all QP states i
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  reduced_couplings_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    reduced_couplings_[gw_level] =
        Mmn_[gw_level + qpoffset] * rpa_red_.BasisU();
  }

  // Build reduced transition model G_ta and Delta_t
  BuildReducedTransitionModel();

  // Extract exact pole expansion within reduced model
  BuildPoleExpansion();

  // Convert reduced pole vectors z_s into Sigma_Exact-style residues
  const Index npoles = rpa_omegas_.size();
  residues_.resize(qptotal_);

#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; ++gw_level) {
    residues_[gw_level] = reduced_couplings_[gw_level] * pole_vectors_;
    if (residues_[gw_level].cols() != npoles) {
      throw std::runtime_error("Sigma_RI_Reduced: residue shape mismatch.");
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

    // rows c=0..n_unocc-1 correspond to virtual states
    const Eigen::MatrixXd Mvc_red =
        Mmn_[v].middleRows(n_occ, n_unocc) * U;

    for (Index c = 0; c < n_unocc; ++c) {
      const double delta = eps(n_occ + c) - eps_v;
      transition_energies_(t) = delta;
      transition_couplings_.row(t) = Mvc_red.row(c);
      ++t;
    }
  }

  // Sort transitions by energy
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

  Index nneg = 0;
  for (Index i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) < 0.0) {
      ++nneg;
    }
  }
  return nneg;
}

void Sigma_RI_Reduced::IsolateRootsInInterval(
    double left, double right, Index nneg_left, Index nneg_right,
    std::vector<std::pair<double, double>>& brackets) const {

  if (nneg_left <= nneg_right) {
    return;  // no root in this interval
  }

  const Index nroots = nneg_left - nneg_right;
  if (nroots == 1 && std::abs(right - left) < 100.0 * opt_red_.pole_tol) {
    brackets.emplace_back(left, right);
    return;
  }

  const double mid = 0.5 * (left + right);
  const Index nneg_mid = CountNegativeEigenvalues(BuildAReal(mid));

  if (mid == left || mid == right) {
    brackets.emplace_back(left, right);
    return;
  }

  IsolateRootsInInterval(left, mid, nneg_left, nneg_mid, brackets);
  IsolateRootsInInterval(mid, right, nneg_mid, nneg_right, brackets);
}

double Sigma_RI_Reduced::RefineRoot(double left, double right, Index nneg_left) const {
  double a = left;
  double b = right;

  for (Index iter = 0; iter < opt_red_.max_bisection; ++iter) {
    const double mid = 0.5 * (a + b);
    const Index nneg_mid = CountNegativeEigenvalues(BuildAReal(mid));

    if (nneg_mid < nneg_left) {
      b = mid;
    } else {
      a = mid;
    }

    if (std::abs(b - a) < opt_red_.pole_tol * (1.0 + std::abs(mid))) {
      return 0.5 * (a + b);
    }
  }

  return 0.5 * (a + b);
}

std::vector<double> Sigma_RI_Reduced::UniqueSortedTransitionEnergies() const {
  std::vector<double> vals;
  vals.reserve(static_cast<std::size_t>(transition_energies_.size()));

  for (Index i = 0; i < transition_energies_.size(); ++i) {
    vals.push_back(transition_energies_(i));
  }

  std::sort(vals.begin(), vals.end());

  const double tol = 1e-12;
  std::vector<double> uniq;
  for (double x : vals) {
    if (uniq.empty() || std::abs(x - uniq.back()) > tol) {
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

  // First interval: (0, Delta_1)
  {
    const double left = 0.0;
    const double right = bare.front() - opt_red_.pole_shift;
    if (right > left) {
      const Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
      const Index nneg_right = CountNegativeEigenvalues(BuildAReal(right));
      IsolateRootsInInterval(left, right, nneg_left, nneg_right, brackets);
    }
  }

  // Intervals between bare transitions
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

  // Final interval: (Delta_max, infinity). Increase upper bound until A > 0.
  {
    const double left = bare.back() + opt_red_.pole_shift;
    double right = std::max(2.0 * bare.back(), bare.back() + 1.0);

    Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
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

  if (brackets.empty()) {
    throw std::runtime_error("Sigma_RI_Reduced: no reduced screening poles found.");
  }

  const Index npoles = static_cast<Index>(brackets.size());
  rpa_omegas_.resize(npoles);
  pole_vectors_.resize(rpa_red_.rank(), npoles);

  for (Index s = 0; s < npoles; ++s) {
    const double left = brackets[static_cast<std::size_t>(s)].first;
    const double right = brackets[static_cast<std::size_t>(s)].second;
    const Index nneg_left = CountNegativeEigenvalues(BuildAReal(left));
    const double omega = RefineRoot(left, right, nneg_left);

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
      throw std::runtime_error("Sigma_RI_Reduced: non-positive residue normalization encountered.");
    }

    rpa_omegas_(s) = omega;
    pole_vectors_.col(s) = p / std::sqrt(alpha);
  }
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

  // restricted spin factor, same as Sigma_Exact
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

}  // namespace xtp
}  // namespace votca