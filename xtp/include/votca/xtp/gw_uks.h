/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 */
#pragma once
#ifndef VOTCA_XTP_GW_UKS_H
#define VOTCA_XTP_GW_UKS_H

#include <memory>
#include <unordered_set>

#include "logger.h"
#include "orbitals.h"
#include "qp_solver_utils.h"
#include "rpa_uks.h"
#include "sigma_base_uks.h"
#include "threecenter.h"
#include "votca/xtp/ppm.h"

namespace votca {
namespace xtp {

class GW_UKS {

  using EvalStage = qp_solver::EvalStage;
  using QPStats = qp_solver::Stats;
  using QPRootCandidate = qp_solver::RootCandidate;
  using QPWindowDiagnostics = qp_solver::WindowDiagnostics;

 public:
  struct options {
    Index homo_alpha;
    Index homo_beta;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double eta;
    double g_sc_limit;
    Index g_sc_max_iterations;
    double gw_sc_limit;
    Index gw_sc_max_iterations;
    double shift = 0;
    double ScaHFX = 0.0;
    std::string sigma_integration;
    Index reset_3c;
    std::string qp_solver;
    double qp_solver_alpha = 0.75;
    Index qp_grid_steps;
    double qp_grid_spacing;
    Index gw_mixing_order;
    double gw_mixing_alpha;
    std::string quadrature_scheme;
    Index order;
    double alpha;
    bool qp_restrict_search = true;
    double qp_zero_margin = 1e-6;
    double qp_virtual_min_energy = -0.1;
    std::string qp_root_finder = "bisection";
    std::string qp_grid_search_mode = "adaptive_with_dense_fallback";
  };

  GW_UKS(Logger& log, TCMatrix_gwbse_spin& Mmn,
         const Eigen::MatrixXd& vxc_alpha, const Eigen::MatrixXd& vxc_beta,
         const Eigen::VectorXd& dft_energies_alpha,
         const Eigen::VectorXd& dft_energies_beta);

  void configure(const options& opt);
  void CalculateGWPerturbation();
  void CalculateHQP();

  Eigen::VectorXd getGWAResultsAlpha() const;
  Eigen::VectorXd getGWAResultsBeta() const;
  const Eigen::VectorXd& RPAInputEnergiesAlpha() const;
  const Eigen::VectorXd& RPAInputEnergiesBeta() const;
  Eigen::MatrixXd getHQPAlpha() const;
  Eigen::MatrixXd getHQPBeta() const;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagonalizeQPHamiltonianAlpha()
      const;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagonalizeQPHamiltonianBeta()
      const;

 private:
  enum class Spin { Alpha, Beta };

  class QPFunc {
   public:
    QPFunc(Index gw_level, const Sigma_base_UKS& sigma, double offset)
        : gw_level_(gw_level), offset_(offset), sigma_c_func_(sigma) {}

    std::pair<double, double> operator()(double frequency) const {
      std::pair<double, double> result;
      result.first = value(frequency, EvalStage::Other);
      result.second = deriv(frequency);
      return result;
    }

    double sigma(double frequency, EvalStage stage = EvalStage::Other) const {
      const std::uint64_t key = FrequencyKey(frequency);

      auto insert_result = seen_frequencies_.insert(key);
      if (!insert_result.second) {
        ++stats_.sigma_repeat_calls;
      } else {
        ++stats_.sigma_unique_frequencies;
      }

      CountSigmaStage(stage);
      return sigma_c_func_.CalcCorrelationDiagElement(gw_level_, frequency);
    }

    double value(double frequency, EvalStage stage = EvalStage::Other) const {
      return sigma(frequency, stage) + offset_ - frequency;
    }

    double deriv(double frequency) const {
      ++stats_.deriv_calls;
      return sigma_c_func_.CalcCorrelationDiagElementDerivative(gw_level_,
                                                                frequency) -
             1.0;
    }

    const QPStats& GetStats() const { return stats_; }

   private:
    static std::uint64_t FrequencyKey(double x) {
      std::uint64_t key = 0;
      static_assert(sizeof(double) == sizeof(std::uint64_t),
                    "Unexpected double size");
      std::memcpy(&key, &x, sizeof(double));
      return key;
    }

    void CountSigmaStage(EvalStage stage) const {
      switch (stage) {
        case EvalStage::Scan:
          ++stats_.sigma_scan_calls;
          break;
        case EvalStage::Refine:
          ++stats_.sigma_refine_calls;
          break;
        case EvalStage::Derivative:
          ++stats_.sigma_derivative_calls;
          break;
        case EvalStage::Other:
        default:
          ++stats_.sigma_other_calls;
          break;
      }
    }

    Index gw_level_;
    double offset_;
    const Sigma_base_UKS& sigma_c_func_;

    mutable std::unordered_set<std::uint64_t> seen_frequencies_;
    mutable QPStats stats_;
  };

  PPM ppm_;
  const Eigen::VectorXd& DftEnergies(Spin spin) const;
  const Eigen::MatrixXd& Vxc(Spin spin) const;
  Eigen::MatrixXd& SigmaX(Spin spin);
  Eigen::MatrixXd& SigmaC(Spin spin);
  const Eigen::MatrixXd& SigmaX(Spin spin) const;
  const Eigen::MatrixXd& SigmaC(Spin spin) const;
  Sigma_base_UKS& SigmaEvaluator(Spin spin);
  const Sigma_base_UKS& SigmaEvaluator(Spin spin) const;
  Index Homo(Spin spin) const;
  const char* SpinName(Spin spin) const;

  Eigen::VectorXd ScissorShift_DFTlevel(const Eigen::VectorXd& dft_energies,
                                        Index homo) const;
  double CalcSpinHomoLumoShift(const Eigen::VectorXd& frequencies,
                               Spin spin) const;
  void PrintGWA_Energies(Spin spin) const;
  void PrintQP_Energies(Spin spin,
                        const Eigen::VectorXd& qp_diag_energies) const;
  Eigen::VectorXd SolveQP(Spin spin, const Eigen::VectorXd& frequencies) const;

  boost::optional<double> SolveQP_Grid(Spin spin, double intercept0,
                                       double frequency0, Index gw_level,
                                       QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed(
      Spin spin, double intercept0, double frequency0, Index gw_level,
      double left_limit, double right_limit, QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed_Adaptive(
      Spin spin, double intercept0, double frequency0, Index gw_level,
      double left_limit, double right_limit, QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Grid_Windowed_Dense(
      Spin spin, double intercept0, double frequency0, Index gw_level,
      double left_limit, double right_limit, QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_FixedPoint(Spin spin, double intercept0,
                                             double frequency0, Index gw_level,
                                             QPStats* stats = nullptr) const;

  boost::optional<double> SolveQP_Linearisation(Spin spin, double intercept0,
                                                double frequency0,
                                                Index gw_level,
                                                QPStats* stats = nullptr) const;

  bool Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                 double epsilon) const;

  boost::optional<QPRootCandidate> RefineQPInterval(
      double lowerbound, double f_lowerbound, double upperbound,
      double f_upperbound, const QPFunc& f, double reference) const;

  std::string LevelLabel(Spin spin, Index level) const;
  const char* OccupationTag(Spin spin, Index level) const;

  Index qptotal_ = 0;
  options opt_;
  Logger& log_;
  Index gw_sc_iteration_ = 0;
  TCMatrix_gwbse_spin& Mmn_;
  const Eigen::MatrixXd& vxc_alpha_;
  const Eigen::MatrixXd& vxc_beta_;
  const Eigen::VectorXd& dft_energies_alpha_;
  const Eigen::VectorXd& dft_energies_beta_;
  RPA_UKS rpa_;
  std::unique_ptr<Sigma_base_UKS> sigma_alpha_;
  std::unique_ptr<Sigma_base_UKS> sigma_beta_;
  Eigen::MatrixXd Sigma_x_alpha_;
  Eigen::MatrixXd Sigma_x_beta_;
  Eigen::MatrixXd Sigma_c_alpha_;
  Eigen::MatrixXd Sigma_c_beta_;
};

}  // namespace xtp
}  // namespace votca

#endif
