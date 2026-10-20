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
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/threecenter.h"
#include <votca/tools/constants.h>

// Local VOTCA includes (charge response)
#include "votca/xtp/chargeresponsesolver.h"

namespace votca {
namespace xtp {

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
bool HasChargeResponse(const std::vector<PolarSegment>& segments) {
  return std::any_of(
      segments.begin(), segments.end(),
      [](const PolarSegment& seg) { return seg.hasChargeResponse(); });
}

// Positions and widths of the charge-response sites, in order
void ChargeResponsePoints(const std::vector<PolarSegment>& segments,
                          std::vector<Eigen::Vector3d>& points,
                          std::vector<double>& widths) {
  for (const PolarSite* site : SitesOf(segments)) {
    if (site->isChargeResponse()) {
      points.push_back(site->getPos());
      widths.push_back(site->CRWidth());
    }
  }
}

// The pieces of ChargeResponseKernel: M^1/2 C and the Cholesky factor of
// A = 1 + M^1/2 T M^1/2 (ok() false if A is not positive definite)
struct ChargeResponseSystem {
  Eigen::MatrixXd MhC;
  CholeskyFactor llt;
};

ChargeResponseSystem BuildChargeResponseSystem(
    const Eigen::MatrixXd& F, const Eigen::MatrixXd& Phi,
    const std::vector<PolarSegment>& segments, double exp_damp) {
  std::vector<PolarSegment> segs = segments;  // the solver needs its own
  const eeInteractor interactor(exp_damp);
  ChargeResponseSolver solver(segs, interactor);
  const Index nq = solver.NumCharges();
  const Index ns = Index(SitesOf(segments).size());
  if (Phi.cols() != nq || Phi.rows() != F.rows()) {
    throw std::runtime_error(
        "EnvironmentScreening::ChargeResponseKernel: Phi is " +
        std::to_string(Phi.rows()) + "x" + std::to_string(Phi.cols()) +
        ", expected " + std::to_string(F.rows()) + "x" + std::to_string(nq) +
        " (AuxPotentialAtPoints at the charge-response sites, in order).");
  }
  // sources c = (-phi, E) per auxiliary function, unknowns in the solver's
  // order: charges first (site order), then 3 dipole components per site
  Eigen::MatrixXd C(solver.size(), F.rows());
  C.topRows(nq) = -Phi.transpose();
  C.bottomRows(3 * ns) = F.transpose();
  ChargeResponseSystem sys;
  sys.MhC = solver.ApplyMhalf(C);
  sys.llt.compute(solver.DenseOperator());
  return sys;
}
}  // namespace

Eigen::MatrixXd EnvironmentScreening::ChargeResponseKernel(
    const Eigen::MatrixXd& F, const Eigen::MatrixXd& Phi,
    const std::vector<PolarSegment>& segments, double exp_damp) {
  const Index n = Index(SitesOf(segments).size());
  CheckShape(F, n, "ChargeResponseKernel");
  if (n == 0) {
    return Eigen::MatrixXd::Zero(F.rows(), F.rows());
  }
  const ChargeResponseSystem sys =
      BuildChargeResponseSystem(F, Phi, segments, exp_damp);
  if (!sys.llt.ok()) {
    throw std::runtime_error(
        "EnvironmentScreening::ChargeResponseKernel: 1 + M^1/2 T M^1/2 of "
        "the polar region (" +
        std::to_string(n) +
        " sites, charge response) is not positive definite: a polarization "
        "catastrophe, no stable response to screen with.");
  }
  // B = -(L^-1 M^1/2 C)^T (L^-1 M^1/2 C)
  const Eigen::MatrixXd Y = sys.llt.SolveL(sys.MhC);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(F.rows(), F.rows());
  B.selfadjointView<Eigen::Lower>().rankUpdate(Y.transpose(), -1.0);
  return B.selfadjointView<Eigen::Lower>();
}

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
  const CholeskyFactor llt(A);
  if (!llt.ok()) {
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
  const Eigen::MatrixXd Y = llt.SolveL(F.transpose());
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

  const SymmetricEigenSystem es = SymmetricEigen(R, false);
  const double lowest = es.values.minCoeff();
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
    if (site->isChargeResponse()) {
      widths.push_back(site->CRWidth());
      continue;
    }
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
    if (HasChargeResponse(env.explicit_segments)) {
      std::vector<Eigen::Vector3d> cr_points;
      std::vector<double> cr_widths;
      ChargeResponsePoints(env.explicit_segments, cr_points, cr_widths);
      const Eigen::MatrixXd Phi =
          AuxPotentialAtPoints(auxbasis, cr_points, cr_widths);
      B += ChargeResponseKernel(F, Phi, env.explicit_segments, env.exp_damp);
    } else {
      B += ReactionFieldKernel(F, env.explicit_segments, env.exp_damp);
    }
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
  const SymmetricEigenSystem es = SymmetricEigen(R);
  const double lowest = es.values.minCoeff();
  if (!(lowest > -1.0)) {
    throw std::runtime_error(
        "EnvironmentScreening::DressingMatrix: 1 + R is not positive "
        "definite (lowest eigenvalue of R " +
        std::to_string(lowest) + ").");
  }
  return es.vectors * (1.0 + es.values.array()).sqrt().matrix().asDiagonal() *
         es.vectors.transpose();
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
  CholeskyFactor llt;
  ChargeResponseSystem cr_sys;
  bool have_cr = false;
  {
    std::vector<Eigen::Vector3d> pos;
    for (const SiteInfo& si : info) {
      if (!si.shell) {
        pos.push_back(si.site->getPos());
      }
    }
    if (!pos.empty() && HasChargeResponse(env.explicit_segments)) {
      F_exp = AuxFieldAtPoints(
          auxbasis, pos, SiteWidths(env.explicit_segments, env.site_width));
      std::vector<Eigen::Vector3d> cr_points;
      std::vector<double> cr_widths;
      ChargeResponsePoints(env.explicit_segments, cr_points, cr_widths);
      cr_sys = BuildChargeResponseSystem(
          F_exp, AuxPotentialAtPoints(auxbasis, cr_points, cr_widths),
          env.explicit_segments, env.exp_damp);
      if (!cr_sys.llt.ok()) {
        out.lowest = -std::numeric_limits<double>::infinity();
        out.report = rep.str() +
                     "  1 + M^1/2 T M^1/2 of the explicit polar sites (charge "
                     "response) is not positive definite: a polarization "
                     "catastrophe inside the environment itself.\n";
        return out;
      }
      const Eigen::MatrixXd Y = cr_sys.llt.SolveL(cr_sys.MhC);
      B.selfadjointView<Eigen::Lower>().rankUpdate(Y.transpose(), -1.0);
      B = B.selfadjointView<Eigen::Lower>();
      have_cr = true;
    } else if (!pos.empty()) {
      F_exp = AuxFieldAtPoints(
          auxbasis, pos, SiteWidths(env.explicit_segments, env.site_width));
      llt.compute(AssembleThole(SitesOf(env.explicit_segments), env.exp_damp));
      if (!llt.ok()) {
        out.lowest = -std::numeric_limits<double>::infinity();
        out.report = rep.str() +
                     "  The Thole matrix of the explicit polar sites is not "
                     "positive definite: a polarization catastrophe inside "
                     "the environment itself, independent of the QM region.\n";
        return out;
      }
      const Eigen::MatrixXd Y = llt.SolveL(F_exp.transpose());
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
  const SymmetricEigenSystem es = SymmetricEigen(R);
  const Eigen::VectorXd& lam = es.values;
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
    const Eigen::VectorXd c = T * es.vectors.col(m);
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
    if (have_cr) {
      // -c^T B c = |L^-1 M^1/2 C c|^2 = sum_i (M^1/2 C c)_i (A^-1 M^1/2 C c)_i
      // over the unknowns; summed per site (its charge and dipole)
      const Eigen::VectorXd g = cr_sys.MhC * c;
      const Eigen::VectorXd y = cr_sys.llt.Solve(g);
      const Eigen::VectorXd gy = g.cwiseProduct(y);
      const Index nq = gy.size() - 3 * n_explicit;
      Index q = 0;
      for (Index j = 0; j < n_explicit; ++j) {
        e_site[std::size_t(j)] = gy.segment<3>(nq + 3 * j).sum();
        if (info[std::size_t(j)].site->isChargeResponse()) {
          e_site[std::size_t(j)] += gy(q++);
        }
      }
    } else if (F_exp.size() > 0) {
      const Eigen::VectorXd g = F_exp.transpose() * c;
      const Eigen::VectorXd mu = llt.Solve(g);
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
