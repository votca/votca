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
    double pole_shift = 1e-8;      // shift away from bare transition poles
    double pole_tol = 1e-10;       // bisection tolerance
    Index max_bisection = 200;
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
  Eigen::VectorXd transition_energies_;   // Delta_t
  Eigen::MatrixXd transition_couplings_;  // rows: t, cols: reduced index a

  // Pole expansion of reduced Wc
  Eigen::VectorXd rpa_omegas_;            // Omega_s
  Eigen::MatrixXd pole_vectors_;          // cols: z_s in reduced RI basis

  // Sigma_Exact-style residues: one matrix per GW state, size (rpatotal x npoles)
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
};

}  // namespace xtp
}  // namespace votca

#endif