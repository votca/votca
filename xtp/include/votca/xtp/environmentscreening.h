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
#ifndef VOTCA_XTP_ENVIRONMENTSCREENING_H
#define VOTCA_XTP_ENVIRONMENTSCREENING_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "aobasis.h"
#include "classicalsegment.h"
#include "eigen.h"
#include "qmmolecule.h"

/**
 * \brief Screening of a GW-BSE calculation by a classical polarizable
 * environment, folded into the screened interaction rather than iterated
 * against it.
 *
 * Follows Li, D'Avino, Duchemin, Beljonne and Blase, Phys. Rev. B (2018)
 * (arXiv:1801.01755), with a Thole induced-dipole environment in place of
 * their charge-response model. The bare Coulomb interaction of the QM
 * subsystem is replaced by
 *
 *   v -> v + v_reac,   v_reac = v_12 chi^(2) v_21,
 *
 * the static response of the environment to a change in the QM charge
 * density, felt back by the QM subsystem. For inducible point dipoles,
 * mu = A^-1 E with A = alpha^-1 + T_Thole -- the same operator PolarRegion
 * solves -- so
 *
 *   v_reac = -F A^-1 F^T,
 *
 * with F mapping a charge density to the field it produces at the sites.
 * It is negative semidefinite: it screens.
 *
 * This class provides the pieces, each separately testable:
 *
 *   AuxFieldAtPoints        F, per auxiliary function
 *   ReactionFieldKernel     B = -F A^-1 F^T for the explicit Thole region
 *   ShellKernel             the cheap tail beyond it (legacy radial_dielectric)
 *   SymmetrizedReactionField  R = T^T B T, in the metric of the stored
 *                           three-centre integrals
 *
 * B lives between auxiliary functions taken as charge densities. In the
 * RI representation used by TCMatrix_gwbse, M = (mn|Q) T with
 * T = TCMatrix_gwbse::InvSqrt(), a density rho_mn has fit coefficients
 * c = V^-1 (Q|mn), and
 *
 *   (mn|v_reac|kl) = c_mn^T B c_kl = M_mn (T^T B T) M_kl^T,
 *
 * using V^-1 = T T^T. That is what R = T^T B T is for: it acts on M exactly
 * as the bare interaction does, which is the identity in that metric.
 */

namespace votca {
namespace xtp {

/**
 * \brief The polarizable environment a GW-BSE calculation is screened by:
 * what a QM/MM job hands to GWBSE (GWBSE::setScreeningEnvironment).
 *
 * explicit_segments respond as one coupled Thole system; shell_segments,
 * outside them, with alpha/epsilon and no coupling (legacy
 * radial_dielectric). Positions and polarizabilities are all that is
 * used; the permanent moments belong to the ground state, where the
 * QM/MM loop has already put them.
 */
struct ScreeningEnvironment {
  std::vector<PolarSegment> explicit_segments;
  double exp_damp = 0.39;
  std::vector<PolarSegment> shell_segments;
  double shell_dielectric = 4.0;
  bool include_kreac = true;
  // Width of the polar sites as seen by the QM charge density, in units of
  // alpha^(1/3): each site responds to the QM field averaged over a
  // unit-charge Gaussian exp(-r^2/R^2)/(pi^(3/2) R^3), R = site_width *
  // alpha_iso^(1/3). 0 gives point sites. See SiteWidths.
  double site_width = 0.5;

  bool empty() const {
    return explicit_segments.empty() && shell_segments.empty();
  }
};

/**
 * \brief Outcome of EnvironmentScreening::Check: whether 1 + R is positive
 * definite, and a human-readable account of what limits it.
 */
struct ScreeningCheck {
  // Lowest eigenvalue of R = T^T B T. Must be > -1.
  double lowest = 0.0;
  // Number of eigenvalues at or below -1.
  Index n_unstable = 0;
  bool ok() const { return lowest > -1.0; }
  // Multi-line report: spectrum, the geometry of the closest contacts, and
  // for the lowest modes where their charge sits on the QM side and which
  // polar sites carry their reaction. Always filled.
  std::string report;
};

class EnvironmentScreening {
 public:
  /**
   * \brief F: the electric field at each point produced by each auxiliary
   * basis function, taken as a charge density.
   *
   * Returns an N_aux x 3 N_points matrix; column 3p+k is field component k
   * at point p. Rows follow the auxiliary basis in exactly the order and
   * normalization TCMatrix_gwbse::Fill uses, because both go through
   * AOBasis::GenerateLibintBasis and getMapToBasisFunctions.
   *
   * The potential of a single auxiliary function at a point is a libint2
   * nuclear-attraction integral against Shell::unit(), which needs no
   * derivative support from libint2. The field is its gradient with
   * respect to the point, taken by a four-point central difference with
   * h = 1e-3 bohr. Measured against the closed-form field of a diffuse
   * s-Gaussian (exponent 0.12) at 2 to 15 bohr, that is accurate to about
   * 1e-12 relative; the two-point stencil at the same h is only 1e-8.
   *
   * The points are assumed to lie outside the auxiliary functions'
   * significant extent, as the method itself assumes -- the environment
   * does not overlap the QM subsystem. Inside a function the field is
   * still correct, but the reaction field would no longer mean anything.
   */
  static Eigen::MatrixXd AuxFieldAtPoints(
      const AOBasis& auxbasis, const std::vector<Eigen::Vector3d>& points,
      const std::vector<double>& widths = {});
  // widths (bohr, one per point, or empty for all zero): a point with
  // width R > 0 is a unit-charge Gaussian exp(-r^2/R^2)/(pi^(3/2) R^3)
  // rather than a point, and F is minus the gradient of the Coulomb
  // interaction (chi_Q | g_R) with respect to its centre -- a libint2
  // two-centre Coulomb integral, same stencil. Outside the extent of the
  // auxiliary functions a spherical charge acts as a point (Newton), so
  // this changes F only where a site sits inside their tails -- which is
  // exactly where point sites over-respond (see Check). R = 0 is the
  // nuclear-attraction path above, unchanged.

  /**
   * \brief R per site, in bohr: site_width * alpha_iso^(1/3), with
   * alpha_iso = tr(alpha)/3.
   *
   * Why smear at all: for a charge density rho and inducible points
   * outside it, <rho|v_reac|rho> is bounded by the dielectric limit,
   * above -<rho|v|rho>. A point dipole inside the tail of a diffuse
   * auxiliary function is not bounded, and 1 + R can lose positive
   * definiteness without any unphysical geometry: measured on a
   * production QM/MM job (55 QM atoms, def2-tzvp/aux-def2-tzvp, 3900 Thole
   * sites, closest contact 2.62 A), a pi-stacked neighbour's carbons with
   * alpha = 18.8 bohr^3 at 3.05-3.3 A carried 86% of a mode at lambda =
   * -1.056. Smearing each site over a Gaussian of width proportional to
   * alpha^(1/3) -- the length scale Thole damping itself uses -- removes
   * that: at site_width 0.5 the same job has lambda_min = -0.755, while
   * 1/2 <rho|v_reac|rho> for the HOMO, LUMO and HOMO-LUMO densities moves
   * by 1e-4 relative (0.8: 2e-3, 1.0: 7e-3).
   */
  static std::vector<double> SiteWidths(
      const std::vector<PolarSegment>& segments, double site_width);

  /**
   * \brief B = -F A^-1 F^T for an explicit Thole region.
   *
   * F must be AuxFieldAtPoints at the positions of the sites of `segments`,
   * in iteration order. A is assembled densely, block for block as
   * DipoleDipoleInteraction::multiply builds it, from the same eeInteractor
   * PolarRegion uses, so the environment responds exactly as the polar
   * region does in the iterative scheme.
   *
   * Dense on purpose. B needs A^-1 applied to N_aux right-hand sides at
   * once. A matrix-free CG re-evaluates every Thole block on every
   * iteration of every solve; assembling A once and factorizing it costs
   * one pass over the pairs plus O((3N)^3), which is far cheaper for any
   * polar region that fits in memory. Formed as -Y^T Y with Y = L^-1 F^T
   * from the Cholesky factor, so B is symmetric and negative semidefinite
   * by construction rather than up to round-off.
   *
   * Throws if A is not positive definite -- a Thole polarization
   * catastrophe, i.e. sites too close together for their damping -- since
   * then there is no stable induced-dipole response to speak of.
   */
  static Eigen::MatrixXd ReactionFieldKernel(
      const Eigen::MatrixXd& F, const std::vector<PolarSegment>& segments,
      double exp_damp);

  /**
   * \brief The tail beyond the explicit region: legacy's radial_dielectric.
   *
   * Legacy (Ewald3DnD::EvaluateRadialCorrection) took a shell of segments
   * outside the polar cutoff, applied the QM field screened by 1/epsilon,
   * induced directly -- no mutual induction -- and took the unscreened
   * interaction back with the QM region. That is exactly an effective
   * polarizability alpha_j / epsilon per site with no T coupling, so its
   * susceptibility is block-diagonal and
   *
   *   B_shell = - sum_j F_j (alpha_j / epsilon) F_j^T,
   *
   * which needs no linear solve. The 1/epsilon stands in for the mutual
   * induction it leaves out. Because the sum runs over real sites, the
   * geometry is whatever the morphology is -- a slab simply has no sites
   * in the vacuum -- which is why this, rather than an analytic Born term,
   * is the tail treatment.
   */
  static Eigen::MatrixXd ShellKernel(const Eigen::MatrixXd& F,
                                     const std::vector<PolarSegment>& segments,
                                     double epsilon);

  /**
   * \brief R = T^T B T, the reaction field in the metric of the stored
   * three-centre integrals, where the bare interaction is the identity.
   *
   * Throws unless 1 + R is positive definite. R itself is negative
   * semidefinite, and 1 + R is the effective interaction u = v + v_reac in
   * this metric; an eigenvalue of R at or below -1 would mean the
   * environment screens a charge fluctuation by more than the fluctuation
   * itself, which no stable environment does.
   */
  static Eigen::MatrixXd SymmetrizedReactionField(const Eigen::MatrixXd& B,
                                                  const Eigen::MatrixXd& T);

  /**
   * \brief S = (1 + R)^(1/2), the dressing of the auxiliary index that
   * turns the bare interaction into u = v + v_reac in the RPA and the
   * correlation self-energy (TCMatrix_gwbse::DressAuxIndex).
   *
   * Symmetric positive definite. Throws unless 1 + R is.
   */
  static Eigen::MatrixXd DressingMatrix(const Eigen::MatrixXd& R);

  /**
   * \brief B for a whole environment: ReactionFieldKernel of the explicit
   * segments plus ShellKernel of the shell, over the auxiliary basis.
   * The two do not couple to each other, as in legacy.
   */
  static Eigen::MatrixXd Kernel(const AOBasis& auxbasis,
                                const ScreeningEnvironment& env);

  /**
   * \brief T for an auxiliary basis on its own: the metric
   * TCMatrix_gwbse::Fill folds into the three-centre integrals, built the
   * same way (Pseudo_InvSqrt_GWBSE with the same tolerance), so R from it
   * is the R the GW-BSE run will see.
   */
  static Eigen::MatrixXd Metric(const AOBasis& auxbasis);

  /**
   * \brief Builds R for env and reports whether 1 + R is positive
   * definite, and why not.
   *
   * Everything R depends on -- the QM atoms, the auxiliary basis on them,
   * the polar sites and their polarizabilities -- is fixed before the
   * QM/MM loop starts, so this is what a job runs up front rather than
   * discovering the problem after the ground state has converged.
   *
   * The eigenvalues of R are those of B c = lambda V c: for the charge
   * distribution rho = sum_Q c_Q chi_Q, lambda is its reaction energy
   * relative to its own Coulomb energy. For inducible points outside rho
   * that ratio is bounded by the dielectric limit, above -1. What can
   * push it lower is a site inside the tail of rho, where the undamped
   * point-dipole response is no longer bounded -- so the report lists,
   * for the n_modes lowest modes, their charge, the QM atoms and auxiliary
   * shells carrying them, and the polar sites carrying their reaction,
   * with the distance of each to the nearest QM atom.
   */
  static ScreeningCheck Check(const AOBasis& auxbasis, const QMMolecule& atoms,
                              const ScreeningEnvironment& env,
                              const Eigen::MatrixXd& T, Index n_modes = 3);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ENVIRONMENTSCREENING_H
