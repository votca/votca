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

#pragma once
#ifndef VOTCA_XTP_CHARGERESPONSESOLVER_H
#define VOTCA_XTP_CHARGERESPONSESOLVER_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "classicalsegment.h"
#include "eeinteractor.h"
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Induced charges and dipoles of a polar region with charge-response
 * segments (Thole segments may be mixed in).
 *
 * Unknowns x = (dq of the charge-response sites, mu of all sites), sources
 * c = (-phi, E) from the potentials and fields accumulated in the sites
 * (Phi + Phi_noE, -(V + V_noE)), M = blockdiag(Pi of each charge-response
 * segment, alpha of each site). With T the couplings of the induced
 * charges and dipoles ((T x)_q = induced potential, (T x)_mu = minus the
 * induced field; nothing inside a charge-response segment, the Thole
 * couplings inside a Thole segment):
 *
 *   x = M (c - T x)  <=>  (1 + M^1/2 T M^1/2) y = M^1/2 c,  x = M^1/2 y,
 *
 * solved by conjugate gradients on the symmetric operator. The null space
 * of Pi (charge conservation) is handled by M^1/2. The kernels are those of
 * eeInteractor (FillDipoleInteraction, GaussianKernel).
 */
class ChargeResponseSolver {
 public:
  ChargeResponseSolver(std::vector<PolarSegment>& segments,
                       const eeInteractor& interactor);

  struct Info {
    Index iterations = 0;
    double error = 0.0;  // relative residual
    bool converged = false;
  };

  /// Solves and writes dq and mu to the sites. With use_previous the
  /// induced values in the sites are the initial guess, else the response
  /// without couplings.
  Info Solve(double tolerance, Index max_iter, bool use_previous);

  /// 1/2 dq^T Pi^+ dq summed over the charge-response segments, the
  /// internal energy of the induced charges (that of the dipoles is in
  /// PolarSite::InternalEnergy).
  double InternalChargeEnergy() const;

  Index size() const { return dim_; }

  // y + M^1/2 T M^1/2 y (public for tests)
  Eigen::VectorXd Apply(const Eigen::VectorXd& y) const;

  /// The operator 1 + M^1/2 T M^1/2 as a dense matrix (for many right-hand
  /// sides at once, EnvironmentScreening).
  Eigen::MatrixXd DenseOperator() const;
  /// M^1/2 X for a matrix of unknown-space columns
  Eigen::MatrixXd ApplyMhalf(const Eigen::MatrixXd& X) const;
  /// Index of the induced charge of site k (the sites of all segments in
  /// order) among the unknowns, -1 for a Thole site; dipole of site k at
  /// NumCharges() + 3 k.
  Index ChargeIndex(Index k) const { return sites_[std::size_t(k)].q; }
  Index NumCharges() const { return nq_; }

 private:
  struct SiteRef {
    Index segment;
    Index site;
    Index q;  // index of its induced charge in x, -1 for Thole sites
  };
  // the couplings of sites a < b (FillDipoleInteraction / Thole inside a
  // Thole segment, Gaussian kernel for the charges); none inside a CR
  // segment (skip)
  struct Coupling {
    bool skip = false;
    Eigen::Matrix3d T;  // grad phi at b of the dipole at a
    bool charges = false;
    Eigen::Vector3d r;  // pos b - pos a
    double B0 = 0.0;
    double B1 = 0.0;
  };
  Coupling PairCoupling(Index a, Index b) const;
  Eigen::VectorXd ApplyT(const Eigen::VectorXd& x) const;
  Eigen::VectorXd ApplyMhalf(const Eigen::VectorXd& x, bool inverse) const;
  Eigen::VectorXd Sources() const;
  Eigen::VectorXd ReadFromSites() const;
  void WriteToSites(const Eigen::VectorXd& x);

  std::vector<PolarSegment>& segments_;
  const eeInteractor& interactor_;
  std::vector<SiteRef> sites_;
  Index nq_ = 0;
  Index dim_ = 0;
  // per charge-response segment: its first charge index, Pi^1/2 and the
  // pseudo-inverse square root Pi^+1/2
  struct SegmentBlock {
    Index segment;
    Index q0;
    Eigen::MatrixXd sqrt;
    Eigen::MatrixXd inv_sqrt;
  };
  std::vector<SegmentBlock> blocks_;
  // per site alpha^1/2 and alpha^-1/2
  std::vector<Eigen::Matrix3d> alpha_sqrt_;
  std::vector<Eigen::Matrix3d> alpha_inv_sqrt_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CHARGERESPONSESOLVER_H
