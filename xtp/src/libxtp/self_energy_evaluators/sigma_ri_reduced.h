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

    bool sigma_aware_basis = false;
    double sigma_mix = 0.25;
    bool normalize_metric_components = true;

    // New targeted sigma metric options
    bool sigma_targeted_basis = false;
    double sigma_target_delta = 0.10;  // Hartree
    std::vector<Index> sigma_target_levels;  // absolute MO indices

    double pole_shift = 1e-8;
    double pole_tol = 1e-10;
    Index max_bisection = 200;

    bool run_pole_diagnostics = true;
    double pole_reconstruction_tol = 1e-8;

    bool run_contracted_pole_diagnostics = true;
    bool run_full_vs_reduced_contracted_diagnostics = true;

    bool run_sigma_diagnostics = true;

    bool run_pole_weight_diagnostics = true;
    Index pole_weight_topn = 12;

    // New m-resolved residue-weight diagnostic
    bool run_m_weight_diagnostics = true;
    Index m_weight_topn = 12;
    // Automatically choose dominant m channels from the residue-weight
    // diagnostic when running contracted-W comparisons.
    bool contracted_use_top_m_weights = true;
    Index contracted_top_m = 8;

    // Frequencies (Hartree) used in contracted-W diagnostics.
    std::vector<double> contracted_diag_omegas = {0.0, 0.5, 1.0};

        bool run_sigma_term_diagnostics = true;
    Index sigma_term_topn = 20;

     // Frequencies for Σ-term diagnostics.
    // If empty, use the input orbital energy of the target level.
    std::vector<double> sigma_term_omegas;

        bool run_sigma_partial_sum_diagnostics = true;
    std::vector<Index> sigma_partial_sum_ns = {1, 2, 5, 10, 20, 50};
    
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

  // Raw Wc amplitudes: wc_amplitudes_(m,s) = (C_i z_s)_m
  std::vector<Eigen::MatrixXd> wc_amplitudes_;

  // Sigma-style residues used in CalcCorrelation...
  // residues_(m,s) = wc_amplitudes_(m,s) / sqrt(2 Omega_s)
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

  Eigen::MatrixXd BuildReducedWcImagFromPoles(double omega,
                                              double prefactor) const;
  std::vector<pole_diagnostic> RunPoleDiagnostics() const;

  double ContractedFullWcImag(Index gw_level, Index m, double omega) const;
  double ContractedProjectedReducedWcImag(Index gw_level, Index m, double omega) const;
  double ContractedDirectWcImag(Index gw_level, Index m, double omega) const;
  double ContractedPoleWcImag(Index gw_level, Index m, double omega) const;

  void RunContractedWcDiagnostics() const;
  void RunFullVsReducedContractedDiagnostics() const;

  double CalcCorrelationDiagElementDirectReduced(Index gw_level,
                                                 double frequency) const;
  double CalcCorrelationOffDiagElementDirectReduced(Index gw_level1,
                                                    Index gw_level2,
                                                    double frequency1,
                                                    double frequency2) const;
  void RunSigmaDiagnostics() const;

  void RunPoleWeightDiagnostics() const;
  void RunMResolvedResidueWeightDiagnostics() const;

    std::vector<Index> GetTopMChannelsForLevel(Index gw_level, Index topn) const;
  std::vector<double> GetContractedDiagnosticOmegas() const;

    struct SigmaTermEntry {
    Index m = -1;
    Index pole = -1;
    double eps_m = 0.0;
    double omega_pole = 0.0;
    double residue2 = 0.0;
    double denom = 0.0;
    double term = 0.0;
    bool occupied = false;
  };

  std::vector<double> GetSigmaTermDiagnosticOmegas(Index gw_level) const;
  std::vector<SigmaTermEntry> BuildSigmaTermEntries(Index gw_level,
                                                    double frequency) const;
  void RunSigmaTermDiagnostics() const;

    void RunSigmaPartialSumDiagnostics() const;
};

}  // namespace xtp
}  // namespace votca

#endif