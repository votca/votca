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

// Local VOTCA includes
#include "votca/xtp/chargeresponsesolver.h"
#include "votca/xtp/gaussiankernel.h"

namespace votca {
namespace xtp {

namespace {
// symmetric square root and pseudo-inverse square root of a positive
// semidefinite matrix; eigenvalues below rel_cutoff * largest are null
template <class Matrix>
void SqrtAndInvSqrt(const Matrix& A, Matrix& sqrt, Matrix& inv_sqrt,
                    double rel_cutoff = 1e-12) {
  Eigen::SelfAdjointEigenSolver<Matrix> es(0.5 * (A + A.transpose()));
  const auto& w = es.eigenvalues();
  const double cutoff = rel_cutoff * std::max(w.maxCoeff(), 0.0);
  auto s = w.array().max(0.0).sqrt().matrix().eval();
  auto si = s.eval();
  for (Index i = 0; i < Index(w.size()); ++i) {
    si(i) = (w(i) > cutoff && w(i) > 0.0) ? 1.0 / s(i) : 0.0;
  }
  sqrt = es.eigenvectors() * s.asDiagonal() * es.eigenvectors().transpose();
  inv_sqrt =
      es.eigenvectors() * si.asDiagonal() * es.eigenvectors().transpose();
}
}  // namespace

ChargeResponseSolver::ChargeResponseSolver(std::vector<PolarSegment>& segments,
                                           const eeInteractor& interactor)
    : segments_(segments), interactor_(interactor) {
  for (Index s = 0; s < Index(segments_.size()); ++s) {
    const PolarSegment& seg = segments_[s];
    if (seg.hasChargeResponse()) {
      SegmentBlock block;
      block.segment = s;
      block.q0 = nq_;
      SqrtAndInvSqrt(seg.getChargeResponse().Pi(), block.sqrt, block.inv_sqrt);
      blocks_.push_back(std::move(block));
    }
    for (Index i = 0; i < seg.size(); ++i) {
      if (seg.hasChargeResponse() != seg[i].isChargeResponse()) {
        throw std::runtime_error(
            "ChargeResponseSolver: the sites of segment " + seg.getType() +
            " are not set up for its charge response (PolarRegion)");
      }
      const Index q = seg.hasChargeResponse() ? nq_++ : -1;
      sites_.push_back({s, i, q});
      Eigen::Matrix3d sq, isq;
      SqrtAndInvSqrt(Eigen::Matrix3d(seg[i].getpolarization()), sq, isq);
      alpha_sqrt_.push_back(sq);
      alpha_inv_sqrt_.push_back(isq);
    }
  }
  dim_ = nq_ + 3 * Index(sites_.size());
}

ChargeResponseSolver::Coupling ChargeResponseSolver::PairCoupling(
    Index a, Index b) const {
  const SiteRef& ra = sites_[std::size_t(a)];
  const SiteRef& rb = sites_[std::size_t(b)];
  const PolarSite& sa = segments_[ra.segment][ra.site];
  const PolarSite& sb = segments_[rb.segment][rb.site];
  Coupling c;
  if (ra.segment == rb.segment) {
    if (ra.q >= 0) {
      c.skip = true;  // nothing inside a charge-response segment
      return c;
    }
    // inside a Thole segment: the Thole couplings of its dipoles
    c.T = interactor_.FillTholeInteraction(sa, sb);
    return c;
  }
  c.T = interactor_.FillDipoleInteraction(sa, sb);
  if (ra.q >= 0 || rb.q >= 0) {
    c.charges = true;
    c.r = sb.getPos() - sa.getPos();
    const GaussianKernel k(c.r.norm(), eeInteractor::PairWidth(sa, sb));
    c.B0 = k.B0;
    c.B1 = k.B1;
  }
  return c;
}

Eigen::VectorXd ChargeResponseSolver::ApplyT(const Eigen::VectorXd& x) const {
  const Index ns = Index(sites_.size());
  Eigen::VectorXd result = Eigen::VectorXd::Zero(dim_);
#pragma omp parallel for schedule(dynamic) reduction(+ : result)
  for (Index a = 0; a < ns; ++a) {
    const Index qa_idx = sites_[std::size_t(a)].q;
    const Index ma = nq_ + 3 * a;
    for (Index b = a + 1; b < ns; ++b) {
      const Coupling c = PairCoupling(a, b);
      if (c.skip) {
        continue;
      }
      const Index qb_idx = sites_[std::size_t(b)].q;
      const Index mb = nq_ + 3 * b;
      // grad phi at b of the dipole at a, and vice versa
      result.segment<3>(mb) += c.T * x.segment<3>(ma);
      result.segment<3>(ma) += c.T.transpose() * x.segment<3>(mb);
      if (!c.charges) {
        continue;
      }
      if (qa_idx >= 0) {
        const double qa = x(qa_idx);
        result.segment<3>(mb) -= qa * c.B1 * c.r;  // grad phi at b
        result(qa_idx) -= c.B1 * x.segment<3>(mb).dot(c.r);
        if (qb_idx >= 0) {
          result(qb_idx) += qa * c.B0;
        }
      }
      if (qb_idx >= 0) {
        const double qb = x(qb_idx);
        result.segment<3>(ma) += qb * c.B1 * c.r;  // grad phi at a
        result(qb_idx) += c.B1 * x.segment<3>(ma).dot(c.r);
        if (qa_idx >= 0) {
          result(qa_idx) += qb * c.B0;
        }
      }
    }
  }
  return result;
}

Eigen::MatrixXd ChargeResponseSolver::DenseOperator() const {
  // T as a dense matrix, entry by entry as ApplyT applies it
  const Index ns = Index(sites_.size());
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(dim_, dim_);
#pragma omp parallel for schedule(dynamic)
  for (Index a = 0; a < ns; ++a) {
    const Index qa = sites_[std::size_t(a)].q;
    const Index ma = nq_ + 3 * a;
    for (Index b = a + 1; b < ns; ++b) {
      const Coupling c = PairCoupling(a, b);
      if (c.skip) {
        continue;
      }
      const Index qb = sites_[std::size_t(b)].q;
      const Index mb = nq_ + 3 * b;
      // rows of a and of b are written by this thread only for this pair;
      // block (ma, mb) and (mb, ma) belong to the pair a < b
      T.block<3, 3>(mb, ma) = c.T;
      T.block<3, 3>(ma, mb) = c.T.transpose();
      if (!c.charges) {
        continue;
      }
      if (qa >= 0) {
        T.block<3, 1>(mb, qa) = -c.B1 * c.r;
        T.block<1, 3>(qa, mb) = -c.B1 * c.r.transpose();
      }
      if (qb >= 0) {
        T.block<3, 1>(ma, qb) = c.B1 * c.r;
        T.block<1, 3>(qb, ma) = c.B1 * c.r.transpose();
      }
      if (qa >= 0 && qb >= 0) {
        T(qa, qb) = c.B0;
        T(qb, qa) = c.B0;
      }
    }
  }
  Eigen::MatrixXd A = ApplyMhalf(ApplyMhalf(T).transpose());
  A = 0.5 * (A + A.transpose()).eval();
  A.diagonal().array() += 1.0;
  return A;
}

Eigen::MatrixXd ChargeResponseSolver::ApplyMhalf(
    const Eigen::MatrixXd& X) const {
  Eigen::MatrixXd result(X.rows(), X.cols());
  for (const SegmentBlock& block : blocks_) {
    const Index n = Index(block.sqrt.rows());
    result.middleRows(block.q0, n) = block.sqrt * X.middleRows(block.q0, n);
  }
  for (Index a = 0; a < Index(sites_.size()); ++a) {
    const Index m = nq_ + 3 * a;
    result.middleRows<3>(m) = alpha_sqrt_[std::size_t(a)] * X.middleRows<3>(m);
  }
  return result;
}

Eigen::VectorXd ChargeResponseSolver::ApplyMhalf(const Eigen::VectorXd& x,
                                                 bool inverse) const {
  Eigen::VectorXd result(dim_);
  for (const SegmentBlock& block : blocks_) {
    const Index n = Index(block.sqrt.rows());
    result.segment(block.q0, n) =
        (inverse ? block.inv_sqrt : block.sqrt) * x.segment(block.q0, n);
  }
  for (Index a = 0; a < Index(sites_.size()); ++a) {
    const Index m = nq_ + 3 * a;
    result.segment<3>(m) =
        (inverse ? alpha_inv_sqrt_[a] : alpha_sqrt_[a]) * x.segment<3>(m);
  }
  return result;
}

Eigen::VectorXd ChargeResponseSolver::Apply(const Eigen::VectorXd& y) const {
  return y + ApplyMhalf(ApplyT(ApplyMhalf(y, false)), false);
}

Eigen::VectorXd ChargeResponseSolver::Sources() const {
  Eigen::VectorXd c(dim_);
  for (Index a = 0; a < Index(sites_.size()); ++a) {
    const PolarSite& site = segments_[sites_[a].segment][sites_[a].site];
    if (sites_[a].q >= 0) {
      c(sites_[a].q) = -(site.Phi() + site.Phi_noE());
    }
    c.segment<3>(nq_ + 3 * a) = -(site.V() + site.V_noE());
  }
  return c;
}

Eigen::VectorXd ChargeResponseSolver::ReadFromSites() const {
  Eigen::VectorXd x(dim_);
  for (Index a = 0; a < Index(sites_.size()); ++a) {
    const PolarSite& site = segments_[sites_[a].segment][sites_[a].site];
    if (sites_[a].q >= 0) {
      x(sites_[a].q) = site.Induced_Charge();
    }
    x.segment<3>(nq_ + 3 * a) = site.Induced_Dipole();
  }
  return x;
}

void ChargeResponseSolver::WriteToSites(const Eigen::VectorXd& x) {
  for (Index a = 0; a < Index(sites_.size()); ++a) {
    PolarSite& site = segments_[sites_[a].segment][sites_[a].site];
    if (sites_[a].q >= 0) {
      site.setInduced_Charge(x(sites_[a].q));
    }
    site.setInduced_Dipole(x.segment<3>(nq_ + 3 * a));
  }
}

ChargeResponseSolver::Info ChargeResponseSolver::Solve(double tolerance,
                                                       Index max_iter,
                                                       bool use_previous) {
  const Eigen::VectorXd rhs = ApplyMhalf(Sources(), false);
  Info info;
  const double rhs_norm = rhs.norm();
  if (rhs_norm == 0.0) {
    WriteToSites(Eigen::VectorXd::Zero(dim_));
    info.converged = true;
    return info;
  }
  Eigen::VectorXd y =
      use_previous ? ApplyMhalf(ReadFromSites(), true) : Eigen::VectorXd(rhs);
  // conjugate gradients on the symmetric positive definite operator
  Eigen::VectorXd r = rhs - Apply(y);
  Eigen::VectorXd p = r;
  double rr = r.squaredNorm();
  for (info.iterations = 0; info.iterations < max_iter; ++info.iterations) {
    info.error = std::sqrt(rr) / rhs_norm;
    if (info.error < tolerance) {
      info.converged = true;
      break;
    }
    const Eigen::VectorXd Ap = Apply(p);
    const double pAp = p.dot(Ap);
    if (!(pAp > 0.0)) {
      throw std::runtime_error(
          "ChargeResponseSolver: the operator 1 + M^1/2 T M^1/2 is not "
          "positive definite (polarization catastrophe); check the widths "
          "and polarizabilities");
    }
    const double alpha = rr / pAp;
    y += alpha * p;
    r -= alpha * Ap;
    const double rr_new = r.squaredNorm();
    p = r + (rr_new / rr) * p;
    rr = rr_new;
  }
  if (!info.converged) {
    info.error = std::sqrt(rr) / rhs_norm;
    info.converged = info.error < tolerance;
  }
  WriteToSites(ApplyMhalf(y, false));
  return info;
}

double ChargeResponseSolver::InternalChargeEnergy() const {
  double e = 0.0;
  for (const SegmentBlock& block : blocks_) {
    const PolarSegment& seg = segments_[block.segment];
    Eigen::VectorXd dq(seg.size());
    for (Index i = 0; i < seg.size(); ++i) {
      dq(i) = seg[i].Induced_Charge();
    }
    e += 0.5 * (block.inv_sqrt * dq).squaredNorm();
  }
  return e;
}

}  // namespace xtp
}  // namespace votca
