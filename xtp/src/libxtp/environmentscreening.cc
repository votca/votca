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
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/threecenter.h"
#include <votca/tools/constants.h>

// libint2 last, as in libint2_calls.cc, otherwise it overrides Eigen. And,
// as in libint2_derivative_calls.cc, WITHOUT <libint2/statics_definition.h>:
// that header defines storage for libint2's static tables and may appear in
// exactly one translation unit of the library, which is libint2_calls.cc.
#include "votca/xtp/make_libint_work.h"
#define LIBINT2_CONSTEXPR_STATICS 0
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#pragma GCC diagnostic ignored "-Wcpp"
#endif
#include <libint2.hpp>
#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace votca {
namespace xtp {

Eigen::MatrixXd EnvironmentScreening::AuxFieldAtPoints(
    const AOBasis& auxbasis, const std::vector<Eigen::Vector3d>& points,
    const std::vector<double>& widths) {
  if (!widths.empty() && widths.size() != points.size()) {
    throw std::runtime_error("EnvironmentScreening::AuxFieldAtPoints: " +
                             std::to_string(widths.size()) + " widths for " +
                             std::to_string(points.size()) + " points.");
  }

  // Same shells, same function offsets, same normalization as
  // TCMatrix_gwbse::Fill3cMO -- which is what makes the rows of F line up
  // with the auxiliary index of the stored three-centre integrals.
  const std::vector<libint2::Shell> shells = auxbasis.GenerateLibintBasis();
  const std::vector<Index> shell2bf = auxbasis.getMapToBasisFunctions();
  const Index n_points = Index(points.size());

  Eigen::MatrixXd F =
      Eigen::MatrixXd::Zero(auxbasis.AOBasisSize(), 3 * n_points);
  if (n_points == 0 || shells.empty()) {
    return F;
  }

  const Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] =
      libint2::Engine(libint2::Operator::nuclear, int(auxbasis.getMaxNprim()),
                      int(auxbasis.getMaxL()), 0);
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }
  // Smeared sites: (chi_Q | g) with g a single s primitive, as the
  // auxiliary metric is computed (xs_xs, unit shells on the other side).
  std::vector<libint2::Engine> coulomb(nthreads);
  coulomb[0] =
      libint2::Engine(libint2::Operator::coulomb, int(auxbasis.getMaxNprim()),
                      int(auxbasis.getMaxL()), 0);
  coulomb[0].set(libint2::BraKet::xs_xs);
  for (Index i = 1; i < nthreads; ++i) {
    coulomb[i] = coulomb[0];
  }

  // Four-point central difference for d/dC, error O(h^4):
  //   f'(C) = [-f(C+2h) + 8 f(C+h) - 8 f(C-h) + f(C-2h)] / (12 h)
  const double h = 1e-3;
  const std::array<double, 4> offset = {2.0 * h, h, -h, -2.0 * h};
  const std::array<double, 4> weight = {-1.0 / (12.0 * h), 8.0 / (12.0 * h),
                                        -8.0 / (12.0 * h), 1.0 / (12.0 * h)};

  // Sign. libint2's nuclear operator with a unit charge at C is
  // -1/|r - C|, so the integral against Shell::unit() is -phi_P(C), minus
  // the potential of the density chi_P at C (verified against the closed
  // form for an s-Gaussian). The field is E = -grad_C phi = +grad_C of
  // that integral, so the stencil is applied to the raw integral as is.
#pragma omp parallel for schedule(dynamic)
  for (Index p = 0; p < n_points; ++p) {
    const double width = widths.empty() ? 0.0 : widths[std::size_t(p)];
    if (width > 0.0) {
      // Field of chi_Q averaged over g: E = -grad_C (chi_Q | g_C), with
      // (chi_Q | g_C) the potential of chi_Q averaged over g (positive
      // for a positive chi_Q), so the stencil is applied with a minus.
      // libint2 normalizes the primitive to unit L2 norm, N = (2b/pi)^3/4;
      // dividing by its charge N (pi/b)^3/2 makes it a unit charge.
      libint2::Engine& engine = coulomb[OPENMP::getThreadId()];
      const libint2::Engine::target_ptr_vec& buf = engine.results();
      const double beta = 1.0 / (width * width);
      const double charge =
          std::pow(2.0 * beta / M_PI, 0.75) * std::pow(M_PI / beta, 1.5);
      for (Index k = 0; k < 3; ++k) {
        for (std::size_t s = 0; s < offset.size(); ++s) {
          Eigen::Vector3d C = points[std::size_t(p)];
          C(k) += offset[s];
          const libint2::Shell g{
              {beta}, {{0, false, {1.0}}}, {{C.x(), C.y(), C.z()}}};
          for (std::size_t sh = 0; sh < shells.size(); ++sh) {
            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs,
                            0>(shells[sh], libint2::Shell::unit(), g,
                               libint2::Shell::unit());
            if (buf[0] == nullptr) {
              continue;
            }
            const Index start = shell2bf[sh];
            for (std::size_t f = 0; f < shells[sh].size(); ++f) {
              F(start + Index(f), 3 * p + k) -= weight[s] * buf[0][f] / charge;
            }
          }
        }
      }
      continue;
    }
    libint2::Engine& engine = engines[OPENMP::getThreadId()];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    for (Index k = 0; k < 3; ++k) {
      for (std::size_t s = 0; s < offset.size(); ++s) {
        Eigen::Vector3d C = points[std::size_t(p)];
        C(k) += offset[s];
        std::vector<std::pair<double, std::array<double, 3>>> charge{
            {1.0, {C.x(), C.y(), C.z()}}};
        engine.set_params(charge);
        for (std::size_t sh = 0; sh < shells.size(); ++sh) {
          engine.compute(shells[sh], libint2::Shell::unit());
          if (buf[0] == nullptr) {
            continue;  // screened out by libint2: exactly zero
          }
          const Index start = shell2bf[sh];
          for (std::size_t f = 0; f < shells[sh].size(); ++f) {
            F(start + Index(f), 3 * p + k) += weight[s] * buf[0][f];
          }
        }
      }
    }
  }
  return F;
}

namespace {
// The polarizable sites in the order ReactionFieldKernel and ShellKernel
// index them, and the order AuxFieldAtPoints must be called with.
std::vector<const PolarSite*> SitesOf(
    const std::vector<PolarSegment>& segments) {
  std::vector<const PolarSite*> sites;
  for (const PolarSegment& seg : segments) {
    for (const PolarSite& site : seg) {
      sites.push_back(&site);
    }
  }
  return sites;
}

void CheckShape(const Eigen::MatrixXd& F, Index n_sites,
                const std::string& who) {
  if (F.cols() != 3 * n_sites) {
    throw std::runtime_error(
        "EnvironmentScreening::" + who + ": F has " + std::to_string(F.cols()) +
        " columns for " + std::to_string(n_sites) +
        " polarizable sites; expected " + std::to_string(3 * n_sites) +
        ". F must be AuxFieldAtPoints at exactly these sites, in order.");
  }
}
// A = alpha^-1 + T_Thole, block for block as
// DipoleDipoleInteraction::multiply applies it: PInv on the diagonal,
// FillTholeInteraction(i, j) above and its transpose below. Every pair,
// intramolecular ones included, exactly as PolarRegion's CG sees it.
Eigen::MatrixXd AssembleThole(const std::vector<const PolarSite*>& sites,
                              double exp_damp) {
  const Index n = Index(sites.size());
  const eeInteractor interactor(exp_damp);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3 * n, 3 * n);
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < n; ++i) {
    A.block<3, 3>(3 * i, 3 * i) = sites[std::size_t(i)]->getPInv();
    for (Index j = i + 1; j < n; ++j) {
      const Eigen::Matrix3d block = interactor.FillTholeInteraction(
          *sites[std::size_t(i)], *sites[std::size_t(j)]);
      A.block<3, 3>(3 * i, 3 * j) = block;
      A.block<3, 3>(3 * j, 3 * i) = block.transpose();
    }
  }
  return A;
}
}  // namespace

Eigen::MatrixXd EnvironmentScreening::ReactionFieldKernel(
    const Eigen::MatrixXd& F, const std::vector<PolarSegment>& segments,
    double exp_damp) {
  const std::vector<const PolarSite*> sites = SitesOf(segments);
  const Index n = Index(sites.size());
  CheckShape(F, n, "ReactionFieldKernel");
  if (n == 0) {
    return Eigen::MatrixXd::Zero(F.rows(), F.rows());
  }

  const Eigen::MatrixXd A = AssembleThole(sites, exp_damp);
  const Eigen::LLT<Eigen::MatrixXd> llt(A);
  if (llt.info() != Eigen::Success) {
    throw std::runtime_error(
        "EnvironmentScreening::ReactionFieldKernel: the Thole interaction "
        "matrix of the polar region (" +
        std::to_string(n) +
        " sites) is not positive definite. That is a polarization "
        "catastrophe -- sites too close together for their damping -- and "
        "there is no stable induced-dipole response for it to screen with.");
  }

  // B = -F A^-1 F^T = -(L^-1 F^T)^T (L^-1 F^T): symmetric and negative
  // semidefinite by construction.
  const Eigen::MatrixXd Y = llt.matrixL().solve(Eigen::MatrixXd(F.transpose()));
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(F.rows(), F.rows());
  B.selfadjointView<Eigen::Lower>().rankUpdate(Y.transpose(), -1.0);
  return B.selfadjointView<Eigen::Lower>();
}

Eigen::MatrixXd EnvironmentScreening::ShellKernel(
    const Eigen::MatrixXd& F, const std::vector<PolarSegment>& segments,
    double epsilon) {
  const std::vector<const PolarSite*> sites = SitesOf(segments);
  const Index n = Index(sites.size());
  CheckShape(F, n, "ShellKernel");
  if (!(epsilon >= 1.0)) {
    throw std::runtime_error(
        "EnvironmentScreening::ShellKernel: epsilon must be at least 1 "
        "(it screens the QM field; below 1 it would amplify it). Got " +
        std::to_string(epsilon) + ".");
  }

  // B_shell = - sum_j F_j (alpha_j / epsilon) F_j^T. Written as -Z Z^T with
  // Z_j = F_j (alpha_j / epsilon)^(1/2), so it is symmetric and negative
  // semidefinite by construction, like the explicit-region kernel.
  Eigen::MatrixXd Z(F.rows(), 3 * n);
  for (Index j = 0; j < n; ++j) {
    const Eigen::Matrix3d alpha = sites[std::size_t(j)]->getPInv().inverse();
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(alpha / epsilon);
    const Eigen::Matrix3d sqrt_alpha = es.operatorSqrt();
    Z.middleCols<3>(3 * j) = F.middleCols<3>(3 * j) * sqrt_alpha;
  }
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(F.rows(), F.rows());
  B.selfadjointView<Eigen::Lower>().rankUpdate(Z, -1.0);
  return B.selfadjointView<Eigen::Lower>();
}

Eigen::MatrixXd EnvironmentScreening::SymmetrizedReactionField(
    const Eigen::MatrixXd& B, const Eigen::MatrixXd& T) {
  if (B.rows() != T.rows() || B.cols() != T.rows()) {
    throw std::runtime_error(
        "EnvironmentScreening::SymmetrizedReactionField: B is " +
        std::to_string(B.rows()) + "x" + std::to_string(B.cols()) +
        " but the metric T has " + std::to_string(T.rows()) +
        " rows. Both must be over the same auxiliary basis.");
  }
  Eigen::MatrixXd R = T.transpose() * B * T;
  // Exactly symmetric, so the eigensolvers downstream see what they expect.
  R = 0.5 * (R + R.transpose()).eval();

  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      R, Eigen::EigenvaluesOnly);
  const double lowest = es.eigenvalues().minCoeff();
  if (!(lowest > -1.0)) {
    throw std::runtime_error(
        "EnvironmentScreening::SymmetrizedReactionField: the lowest "
        "eigenvalue of R is " +
        std::to_string(lowest) +
        ", so 1 + R, the effective interaction v + v_reac in this metric, "
        "is not positive definite. The environment would screen a charge "
        "fluctuation by more than the fluctuation itself.");
  }
  return R;
}

std::vector<double> EnvironmentScreening::SiteWidths(
    const std::vector<PolarSegment>& segments, double site_width) {
  if (!(site_width >= 0.0)) {
    throw std::runtime_error(
        "EnvironmentScreening::SiteWidths: site_width must be >= 0, got " +
        std::to_string(site_width) + ".");
  }
  std::vector<double> widths;
  for (const PolarSite* site : SitesOf(segments)) {
    const double alpha_iso = site->getPInv().inverse().trace() / 3.0;
    widths.push_back(site_width * std::cbrt(alpha_iso));
  }
  return widths;
}

Eigen::MatrixXd EnvironmentScreening::Kernel(const AOBasis& auxbasis,
                                             const ScreeningEnvironment& env) {
  auto positions = [](const std::vector<PolarSegment>& segs) {
    std::vector<Eigen::Vector3d> pos;
    for (const PolarSite* site : SitesOf(segs)) {
      pos.push_back(site->getPos());
    }
    return pos;
  };
  const Index naux = auxbasis.AOBasisSize();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(naux, naux);
  if (!env.explicit_segments.empty()) {
    const Eigen::MatrixXd F =
        AuxFieldAtPoints(auxbasis, positions(env.explicit_segments),
                         SiteWidths(env.explicit_segments, env.site_width));
    B += ReactionFieldKernel(F, env.explicit_segments, env.exp_damp);
  }
  if (!env.shell_segments.empty()) {
    const Eigen::MatrixXd F =
        AuxFieldAtPoints(auxbasis, positions(env.shell_segments),
                         SiteWidths(env.shell_segments, env.site_width));
    B += ShellKernel(F, env.shell_segments, env.shell_dielectric);
  }
  return B;
}

Eigen::MatrixXd EnvironmentScreening::DressingMatrix(const Eigen::MatrixXd& R) {
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(R);
  const double lowest = es.eigenvalues().minCoeff();
  if (!(lowest > -1.0)) {
    throw std::runtime_error(
        "EnvironmentScreening::DressingMatrix: 1 + R is not positive "
        "definite (lowest eigenvalue of R " +
        std::to_string(lowest) + ").");
  }
  return es.eigenvectors() *
         (1.0 + es.eigenvalues().array()).sqrt().matrix().asDiagonal() *
         es.eigenvectors().transpose();
}

Eigen::MatrixXd EnvironmentScreening::Metric(const AOBasis& auxbasis) {
  AOCoulomb coulomb;
  coulomb.Fill(auxbasis);
  AOOverlap overlap;
  overlap.Fill(auxbasis);
  return coulomb.Pseudo_InvSqrt_GWBSE(overlap,
                                      TCMatrix_gwbse::metric_tolerance);
}

namespace {
// Net charge of every auxiliary function, integral chi_Q d^3r: an overlap
// integral against the unit shell, in the order and normalization of the
// three-centre integrals (same shells, same offsets).
Eigen::VectorXd AuxCharges(const AOBasis& auxbasis) {
  const std::vector<libint2::Shell> shells = auxbasis.GenerateLibintBasis();
  const std::vector<Index> shell2bf = auxbasis.getMapToBasisFunctions();
  Eigen::VectorXd q = Eigen::VectorXd::Zero(auxbasis.AOBasisSize());
  libint2::Engine engine(libint2::Operator::overlap,
                         int(auxbasis.getMaxNprim()), int(auxbasis.getMaxL()),
                         0);
  const libint2::Engine::target_ptr_vec& buf = engine.results();
  for (std::size_t sh = 0; sh < shells.size(); ++sh) {
    engine.compute(shells[sh], libint2::Shell::unit());
    if (buf[0] == nullptr) {
      continue;
    }
    for (std::size_t f = 0; f < shells[sh].size(); ++f) {
      q(shell2bf[sh] + Index(f)) = buf[0][f];
    }
  }
  return q;
}

struct SiteInfo {
  const PolarSite* site;
  Index segment;
  bool shell;
  double alpha;   // isotropic, bohr^3
  double d_qm;    // to the nearest QM atom, bohr
  Index qm_atom;  // that atom's position in the QM molecule
};
}  // namespace

ScreeningCheck EnvironmentScreening::Check(const AOBasis& auxbasis,
                                           const QMMolecule& atoms,
                                           const ScreeningEnvironment& env,
                                           const Eigen::MatrixXd& T,
                                           Index n_modes) {
  const double b2a = tools::conv::bohr2ang;
  ScreeningCheck out;
  std::ostringstream rep;
  rep << std::fixed;

  // Every site, explicit first then shell, with its distance to the QM
  // molecule -- the quantity the failure mode depends on.
  std::vector<SiteInfo> info;
  auto collect = [&](const std::vector<PolarSegment>& segs, bool shell) {
    for (const PolarSegment& seg : segs) {
      for (const PolarSite& site : seg) {
        SiteInfo si{&site,
                    seg.getId(),
                    shell,
                    3.0 / site.getPInv().trace(),
                    std::numeric_limits<double>::max(),
                    -1};
        for (Index a = 0; a < atoms.size(); ++a) {
          const double d = (atoms[a].getPos() - site.getPos()).norm();
          if (d < si.d_qm) {
            si.d_qm = d;
            si.qm_atom = a;  // position, as AOShell::getAtomIndex counts
          }
        }
        info.push_back(si);
      }
    }
  };
  collect(env.explicit_segments, false);
  collect(env.shell_segments, true);
  const Index n_explicit = Index(SitesOf(env.explicit_segments).size());

  auto describe_site = [&](const SiteInfo& si) {
    std::ostringstream o;
    o << std::fixed << std::setprecision(2) << "segment " << si.segment << " "
      << si.site->getElement() << " (alpha " << si.alpha << " bohr^3"
      << (si.shell ? ", shell" : "") << ") at " << si.d_qm * b2a
      << " A from QM atom " << si.qm_atom << " "
      << atoms[si.qm_atom].getElement();
    return o.str();
  };

  // Geometry: how close does the environment come?
  {
    std::vector<std::size_t> order(info.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
      order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      return info[a].d_qm < info[b].d_qm;
    });
    rep << "  " << info.size() << " polar sites around " << atoms.size()
        << " QM atoms; sites within";
    for (double r : {1.5, 2.0, 2.5, 3.0}) {
      const Index n = std::count_if(info.begin(), info.end(), [&](auto& si) {
        return si.d_qm * b2a < r;
      });
      rep << std::setprecision(1) << " " << r << " A: " << n << ",";
    }
    rep << "\n  closest sites:\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(5, order.size()); ++k) {
      rep << "    " << describe_site(info[order[k]]) << "\n";
    }
  }

  // B and R exactly as Kernel and SymmetrizedReactionField build them,
  // keeping the pieces needed to take the lowest modes apart.
  const Index naux = auxbasis.AOBasisSize();
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(naux, naux);
  Eigen::MatrixXd F_exp, F_sh;
  Eigen::LLT<Eigen::MatrixXd> llt;
  {
    std::vector<Eigen::Vector3d> pos;
    for (const SiteInfo& si : info) {
      if (!si.shell) {
        pos.push_back(si.site->getPos());
      }
    }
    if (!pos.empty()) {
      F_exp = AuxFieldAtPoints(
          auxbasis, pos, SiteWidths(env.explicit_segments, env.site_width));
      llt.compute(AssembleThole(SitesOf(env.explicit_segments), env.exp_damp));
      if (llt.info() != Eigen::Success) {
        out.lowest = -std::numeric_limits<double>::infinity();
        out.report = rep.str() +
                     "  The Thole matrix of the explicit polar sites is not "
                     "positive definite: a polarization catastrophe inside "
                     "the environment itself, independent of the QM region.\n";
        return out;
      }
      const Eigen::MatrixXd Y =
          llt.matrixL().solve(Eigen::MatrixXd(F_exp.transpose()));
      B.selfadjointView<Eigen::Lower>().rankUpdate(Y.transpose(), -1.0);
      B = B.selfadjointView<Eigen::Lower>();
    }
  }
  if (!env.shell_segments.empty()) {
    std::vector<Eigen::Vector3d> pos;
    for (const SiteInfo& si : info) {
      if (si.shell) {
        pos.push_back(si.site->getPos());
      }
    }
    F_sh = AuxFieldAtPoints(auxbasis, pos,
                            SiteWidths(env.shell_segments, env.site_width));
    B += ShellKernel(F_sh, env.shell_segments, env.shell_dielectric);
  }

  Eigen::MatrixXd R = T.transpose() * B * T;
  R = 0.5 * (R + R.transpose()).eval();
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(R);
  const Eigen::VectorXd& lam = es.eigenvalues();
  out.lowest = lam(0);
  out.n_unstable = (lam.array() <= -1.0).count();
  rep << std::setprecision(2) << "  site_width " << env.site_width
      << (env.site_width > 0.0 ? " alpha^(1/3)" : " (point sites)") << "\n";
  rep << std::setprecision(4) << "  lowest eigenvalues of R (must be > -1):";
  for (Index m = 0; m < std::min<Index>(6, lam.size()); ++m) {
    rep << " " << lam(m);
  }
  rep << "\n  " << out.n_unstable << " at or below -1, "
      << (lam.array() < -0.5).count() << " below -0.5, of " << lam.size()
      << "\n";

  // The lowest modes. Column w of the eigenvectors is in the metric of
  // the three-centre integrals; c = T w are its coefficients over the
  // auxiliary functions as charge densities, rho = sum_Q c_Q chi_Q, with
  // c^T V c = 1.
  const Eigen::VectorXd q_aux = AuxCharges(auxbasis);
  for (Index m = 0; m < std::min<Index>(n_modes, lam.size()); ++m) {
    const Eigen::VectorXd c = T * es.eigenvectors().col(m);
    rep << std::setprecision(4) << "  mode " << m << ": lambda " << lam(m)
        << ", net charge " << c.dot(q_aux)
        << " (in units where its Coulomb self energy is 1/2)\n";

    // Where the charge sits: squared coefficients by QM atom, and the
    // three auxiliary shells carrying most of them.
    const double norm = c.squaredNorm();
    std::vector<double> by_atom(atoms.size(), 0.0);
    std::vector<std::pair<double, std::string>> by_shell;
    for (const AOShell& shell : auxbasis) {
      const double w =
          c.segment(shell.getStartIndex(), shell.getNumFunc()).squaredNorm() /
          norm;
      by_atom[std::size_t(shell.getAtomIndex())] += w;
      std::ostringstream o;
      o << std::setprecision(3) << "atom " << shell.getAtomIndex() << " "
        << atoms[shell.getAtomIndex()].getElement() << " "
        << EnumToString(shell.getL()) << " exp " << shell.getMinDecay();
      by_shell.push_back({w, o.str()});
    }
    std::sort(by_shell.rbegin(), by_shell.rend());
    rep << "    aux shells:";
    for (std::size_t k = 0; k < std::min<std::size_t>(3, by_shell.size());
         ++k) {
      rep << std::setprecision(2) << " " << by_shell[k].second << " ("
          << 100 * by_shell[k].first << "%);";
    }
    rep << "\n";

    // Which sites carry its reaction: g = F^T c is the mode's field at the
    // sites, mu its induced dipoles, and -c^T B c = sum_j g_j . mu_j.
    std::vector<double> e_site(info.size(), 0.0);
    double e_total = 0.0;
    if (F_exp.size() > 0) {
      const Eigen::VectorXd g = F_exp.transpose() * c;
      const Eigen::VectorXd mu = llt.solve(g);
      for (Index j = 0; j < n_explicit; ++j) {
        e_site[std::size_t(j)] = g.segment<3>(3 * j).dot(mu.segment<3>(3 * j));
      }
    }
    if (F_sh.size() > 0) {
      const Eigen::VectorXd g = F_sh.transpose() * c;
      for (Index j = 0; j < Index(info.size()) - n_explicit; ++j) {
        const SiteInfo& si = info[std::size_t(n_explicit + j)];
        const Eigen::Matrix3d alpha = si.site->getPInv().inverse();
        e_site[std::size_t(n_explicit + j)] =
            g.segment<3>(3 * j).dot(alpha * g.segment<3>(3 * j)) /
            env.shell_dielectric;
      }
    }
    double e_near = 0.0;
    for (std::size_t j = 0; j < info.size(); ++j) {
      e_total += e_site[j];
      if (info[j].d_qm * b2a < 3.0) {
        e_near += e_site[j];
      }
    }
    std::vector<std::size_t> order(info.size());
    for (std::size_t j = 0; j < order.size(); ++j) {
      order[j] = j;
    }
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      return std::abs(e_site[a]) > std::abs(e_site[b]);
    });
    rep << std::setprecision(1) << "    " << 100 * e_near / e_total
        << "% of its reaction from sites within 3 A of the QM atoms; "
           "largest:\n";
    for (std::size_t k = 0; k < std::min<std::size_t>(5, order.size()); ++k) {
      rep << "      " << std::setprecision(1)
          << 100 * e_site[order[k]] / e_total << "%  "
          << describe_site(info[order[k]]) << "\n";
    }
  }
  out.report = rep.str();
  return out;
}

}  // namespace xtp
}  // namespace votca
