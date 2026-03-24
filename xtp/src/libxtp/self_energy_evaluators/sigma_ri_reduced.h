#pragma once
#ifndef VOTCA_XTP_SIGMA_RI_REDUCED_H
#define VOTCA_XTP_SIGMA_RI_REDUCED_H

#include <utility>
#include <vector>

#include "votca/xtp/rpa_ri_reduced.h"
#include "votca/xtp/sigma_base.h"

namespace votca {
namespace xtp {

class Sigma_RI_Reduced : public Sigma_base {
 public:
  Sigma_RI_Reduced(TCMatrix_gwbse& Mmn, RPA& rpa);

  struct reduced_options {
    double imag_omega_max = 2.0;
    Index imag_omega_points = 12;
    double basis_threshold = 1e-8;
    Index max_rank = -1;

    // pole search
    double pole_shift = 1e-8;
    double pole_tol = 1e-10;
    Index max_bisection = 200;

    // diagnostics
    bool run_pole_diagnostics = true;
    double pole_reconstruction_tol = 1e-8;
  };

  struct pole_diagnostic {
    double omega = 0.0;

    double norm_direct = 0.0;

    double norm_pole2 = 0.0;
    double abs_diff_pole2 = 0.0;
    double rel_diff_pole2 = 0.0;

    double norm_pole4 = 0.0;
    double abs_diff_pole4 = 0.0;
    double rel_diff_pole4 = 0.0;
  };

  void configure_reduced(reduced_options opt) { opt_red_ = opt; }

  void PrepareScreening() final;

  double CalcCorrelationDiagElement(Index gw_level,
                                    double frequency) const final;

  double CalcCorrelationDiagElementDerivative(Index gw_level,
                                              double frequency) const final;

  double CalcCorrelationOffDiagElement(Index gw_level1, Index gw_level2,
                                       double frequency1,
                                       double frequency2) const final;

 private:
  reduced_options opt_red_;
  RPA_RI_Reduced rpa_red_;

  // Reduced transition model
  Eigen::VectorXd transition_energies_;
  Eigen::MatrixXd transition_couplings_;

  // Pole expansion of reduced Wc
  Eigen::VectorXd rpa_omegas_;
  Eigen::MatrixXd pole_vectors_;

  // Sigma_Exact-style residues
  std::vector<Eigen::MatrixXd> residues_;

  // Per-QP reduced couplings C_i(m,a) = sum_mu M_im^mu U_mu,a
  std::vector<Eigen::MatrixXd> reduced_couplings_;

  void BuildReducedTransitionModel();
  void BuildPoleExpansion();

  Eigen::MatrixXd BuildAReal(double omega) const;
  Eigen::MatrixXd BuildARealDerivative(double omega) const;

  Index CountNegativeEigenvalues(const Eigen::MatrixXd& A) const;
  void IsolateRootsInInterval(double left, double right, Index nneg_left,
                              Index nneg_right,
                              std::vector<std::pair<double, double>>& brackets) const;
  double RefineRoot(double left, double right, Index nneg_left) const;

  std::vector<double> UniqueSortedTransitionEnergies() const;

  // diagnostics
  Eigen::MatrixXd BuildReducedWcImagFromPoles(double omega,
                                              double prefactor) const;
  std::vector<pole_diagnostic> RunPoleDiagnostics() const;

struct contracted_w_diagnostic {
  Index gw_level = 0;
  Index m = 0;
  double omega = 0.0;
  double direct = 0.0;
  double pole = 0.0;
  double abs_diff = 0.0;
  double rel_diff = 0.0;
};

double ContractedDirectWcImag(Index gw_level, Index m, double omega) const;
double ContractedPoleWcImag(Index gw_level, Index m, double omega) const;
void RunContractedWcDiagnostics() const;

};

}  // namespace xtp
}  // namespace votca

#endif