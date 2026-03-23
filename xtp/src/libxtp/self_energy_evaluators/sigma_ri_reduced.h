#pragma once
#ifndef VOTCA_XTP_SIGMA_RI_REDUCED_H
#define VOTCA_XTP_SIGMA_RI_REDUCED_H

#include "votca/xtp/sigma_base.h"
#include "votca/xtp/rpa_ri_reduced.h"

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
  };

  void configure_reduced(reduced_options opt) { opt_red_ = opt; }

void PrepareScreening() final;

double CalcCorrelationDiagElement(Index gw_level,
                                  double frequency) const final;

double CalcCorrelationDiagElementDerivative(Index gw_level,
                                            double frequency) const final;

double CalcCorrelationOffDiagElement(Index, Index, double, double) const final {
  throw std::runtime_error("Off-diagonal Sigma not implemented yet");
}

 private:
  std::vector<double> ImagFrequencyGrid() const;
  std::vector<double> ImagFrequencyWeights() const;

  RPA_RI_Reduced rpa_red_;
  reduced_options opt_red_;

  std::vector<Eigen::MatrixXd> reduced_couplings_;  // per QP state
};

}  // namespace xtp
}  // namespace votca

#endif