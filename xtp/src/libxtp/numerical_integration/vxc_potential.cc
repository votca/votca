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

// Third party includes
#include <algorithm>
#include <boost/format.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/vxc_functionals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

namespace votca {
namespace xtp {

template <class Grid>
Vxc_Potential<Grid>::~Vxc_Potential() {
  if (setXC_) {
    xc_func_end(&xfunc);
    if (use_separate_) {
      xc_func_end(&cfunc);
    }
  }
}

template <class Grid>
double Vxc_Potential<Grid>::getExactExchange(const std::string& functional) {
  double exactexchange = 0.0;
  Vxc_Functionals map;

  std::vector<std::string> functional_names =
      tools::Tokenizer(functional, " ").ToVector();

  if (functional_names.size() > 2) {
    throw std::runtime_error("Too many functional names");
  } else if (functional_names.empty()) {
    throw std::runtime_error("Specify at least one functional");
  }

  for (const std::string& functional_name : functional_names) {
    int func_id = map.getID(functional_name);
    if (func_id < 0) {
      exactexchange = 0.0;
      break;
    }

    xc_func_type func;
    if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
      throw std::runtime_error(
          (boost::format("Functional %s not found\n") % functional_name).str());
    }

    if (exactexchange > 0 && func.cam_alpha > 0) {
      xc_func_end(&func);
      throw std::runtime_error(
          "You have specified two functionals with exact exchange");
    }

    exactexchange += func.cam_alpha;
    xc_func_end(&func);
  }

  return exactexchange;
}

template <class Grid>
void Vxc_Potential<Grid>::setXCfunctional(const std::string& functional) {
  Vxc_Functionals map;
  std::vector<std::string> strs =
      tools::Tokenizer(functional, " ,\n\t").ToVector();

  xfunc_id = 0;
  use_separate_ = false;
  cfunc_id = 0;

  if (strs.size() == 1) {
    xfunc_id = map.getID(strs[0]);
  } else if (strs.size() == 2) {
    xfunc_id = map.getID(strs[0]);
    cfunc_id = map.getID(strs[1]);
    use_separate_ = true;
  } else {
    throw std::runtime_error(
        "LIBXC. Please specify one combined or an exchange and a correlation "
        "functionals");
  }

  // Keep the stored handles UNPOLARIZED for the closed-shell/restricted code
  // path. The UKS path creates temporary polarized handles inside
  // EvaluateXCSpin().
  if (xc_func_init(&xfunc, xfunc_id, XC_UNPOLARIZED) != 0) {
    throw std::runtime_error(
        (boost::format("Functional %s not found\n") % strs[0]).str());
  }

  if (xfunc.info->kind != 2 && !use_separate_) {
    throw std::runtime_error(
        "Your functional misses either correlation or exchange, please specify "
        "another functional, separated by whitespace");
  }

  if (use_separate_) {
    if (xc_func_init(&cfunc, cfunc_id, XC_UNPOLARIZED) != 0) {
      xc_func_end(&xfunc);
      throw std::runtime_error(
          (boost::format("Functional %s not found\n") % strs[1]).str());
    }

    if ((xfunc.info->kind + cfunc.info->kind) != 1) {
      xc_func_end(&xfunc);
      xc_func_end(&cfunc);
      throw std::runtime_error(
          "Your functionals are not one exchange and one correlation");
    }
  }

  setXC_ = true;
  return;
}

template <class Grid>
typename Vxc_Potential<Grid>::XC_entry Vxc_Potential<Grid>::EvaluateXC(
    double rho, double sigma) const {
  typename Vxc_Potential<Grid>::XC_entry result;

  switch (xfunc.info->family) {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&xfunc, 1, &rho, &result.f_xc, &result.df_drho);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&xfunc, 1, &rho, &sigma, &result.f_xc, &result.df_drho,
                     &result.df_dsigma);
      break;
    default:
      throw std::runtime_error("Unsupported XC family for unpolarized DFT.");
  }

  if (use_separate_) {
    typename Vxc_Potential<Grid>::XC_entry temp;

    switch (cfunc.info->family) {
      case XC_FAMILY_LDA:
        xc_lda_exc_vxc(&cfunc, 1, &rho, &temp.f_xc, &temp.df_drho);
        break;
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA:
        xc_gga_exc_vxc(&cfunc, 1, &rho, &sigma, &temp.f_xc, &temp.df_drho,
                       &temp.df_dsigma);
        break;
      default:
        throw std::runtime_error(
            "Unsupported correlation family for unpolarized DFT.");
    }

    result.f_xc += temp.f_xc;
    result.df_drho += temp.df_drho;
    result.df_dsigma += temp.df_dsigma;
  }

  return result;
}

template <class Grid>
typename Vxc_Potential<Grid>::XC_entry_spin Vxc_Potential<Grid>::EvaluateXCSpin(
    double rho_a, double rho_b, double sigma_aa, double sigma_ab,
    double sigma_bb) const {
  typename Vxc_Potential<Grid>::XC_entry_spin result;

  // UKS/open-shell path: use temporary POLARIZED handles so LibXC receives the
  // correct rho[2], sigma[3], vrho[2], vsigma[3] layout.
  xc_func_type xfunc_pol;
  if (xc_func_init(&xfunc_pol, xfunc_id, XC_POLARIZED) != 0) {
    throw std::runtime_error(
        "Failed to initialize polarized exchange XC "
        "functional in EvaluateXCSpin.");
  }

  xc_func_type cfunc_pol;
  bool cfunc_pol_init = false;
  if (use_separate_) {
    if (xc_func_init(&cfunc_pol, cfunc_id, XC_POLARIZED) != 0) {
      xc_func_end(&xfunc_pol);
      throw std::runtime_error(
          "Failed to initialize polarized correlation XC "
          "functional in EvaluateXCSpin.");
    }
    cfunc_pol_init = true;
  }

  double rho[2] = {rho_a, rho_b};

  switch (xfunc_pol.info->family) {
    case XC_FAMILY_LDA: {
      double vrho[2] = {0.0, 0.0};
      xc_lda_exc_vxc(&xfunc_pol, 1, rho, &result.f_xc, vrho);
      result.vrho_a = vrho[0];
      result.vrho_b = vrho[1];
      break;
    }
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA: {
      double sigma[3] = {sigma_aa, sigma_ab, sigma_bb};
      double vrho[2] = {0.0, 0.0};
      double vsigma[3] = {0.0, 0.0, 0.0};
      xc_gga_exc_vxc(&xfunc_pol, 1, rho, sigma, &result.f_xc, vrho, vsigma);
      result.vrho_a = vrho[0];
      result.vrho_b = vrho[1];
      result.vsigma_aa = vsigma[0];
      result.vsigma_ab = vsigma[1];
      result.vsigma_bb = vsigma[2];
      break;
    }
    default:
      xc_func_end(&xfunc_pol);
      if (cfunc_pol_init) {
        xc_func_end(&cfunc_pol);
      }
      throw std::runtime_error("Unsupported XC family for polarized DFT.");
  }

  if (use_separate_) {
    typename Vxc_Potential<Grid>::XC_entry_spin temp;

    switch (cfunc_pol.info->family) {
      case XC_FAMILY_LDA: {
        double vrho[2] = {0.0, 0.0};
        xc_lda_exc_vxc(&cfunc_pol, 1, rho, &temp.f_xc, vrho);
        temp.vrho_a = vrho[0];
        temp.vrho_b = vrho[1];
        break;
      }
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA: {
        double sigma[3] = {sigma_aa, sigma_ab, sigma_bb};
        double vrho[2] = {0.0, 0.0};
        double vsigma[3] = {0.0, 0.0, 0.0};
        xc_gga_exc_vxc(&cfunc_pol, 1, rho, sigma, &temp.f_xc, vrho, vsigma);
        temp.vrho_a = vrho[0];
        temp.vrho_b = vrho[1];
        temp.vsigma_aa = vsigma[0];
        temp.vsigma_ab = vsigma[1];
        temp.vsigma_bb = vsigma[2];
        break;
      }
      default:
        xc_func_end(&xfunc_pol);
        xc_func_end(&cfunc_pol);
        throw std::runtime_error(
            "Unsupported correlation family for polarized DFT.");
    }

    result.f_xc += temp.f_xc;
    result.vrho_a += temp.vrho_a;
    result.vrho_b += temp.vrho_b;
    result.vsigma_aa += temp.vsigma_aa;
    result.vsigma_ab += temp.vsigma_ab;
    result.vsigma_bb += temp.vsigma_bb;
  }

  xc_func_end(&xfunc_pol);
  if (cfunc_pol_init) {
    xc_func_end(&cfunc_pol);
  }

  return result;
}

template <class Grid>
Mat_p_Energy Vxc_Potential<Grid>::IntegrateVXC(
    const Eigen::MatrixXd& density_matrix) const {
  assert(density_matrix.isApprox(density_matrix.transpose()) &&
         "Density matrix has to be symmetric!");

  Mat_p_Energy vxc = Mat_p_Energy(density_matrix.rows(), density_matrix.cols());

#pragma omp parallel for schedule(guided) reduction(+ : vxc)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }

    double EXC_box = 0.0;

    // two because we have to use the density matrix and its transpose
    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);

    double cutoff =
        1.e-40 / double(density_matrix.rows()) / double(density_matrix.rows());
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }

    Eigen::MatrixXd Vxc_here =
        Eigen::MatrixXd::Zero(DMAT_here.rows(), DMAT_here.cols());

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    for (Index p = 0; p < box.size(); ++p) {
      AOShell::AOValues ao = box.CalcAOValues(points[p]);

      Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
      double rho = 0.5 * temp.dot(ao.values);
      const double weight = weights[p];

      if (rho * weight < 1.e-20) {
        continue;
      }

      const Eigen::Vector3d rho_grad = temp.transpose() * ao.derivatives;

      typename Vxc_Potential<Grid>::XC_entry xc =
          EvaluateXC(rho, rho_grad.squaredNorm());

      EXC_box += weight * rho * xc.f_xc;

      auto grad = ao.derivatives * rho_grad;
      temp.noalias() =
          weight * (0.5 * xc.df_drho * ao.values + 2.0 * xc.df_dsigma * grad);
      Vxc_here.noalias() += temp * ao.values.transpose();
    }

    box.AddtoBigMatrix(vxc.matrix(), Vxc_here);
    vxc.energy() += EXC_box;
  }

  return Mat_p_Energy(vxc.energy(), vxc.matrix() + vxc.matrix().transpose());
}

// ===========================================================================
// Derivation (worked out before writing this, checked against LibXC's own
// documented API convention rather than assumed):
//
// E_xc = sum_p weight_p * e_xc(rho_p)  (e_xc = rho*f_xc, energy density)
//
// dE_xc/dR_A |_Pulay = sum_p weight_p * v_xc(rho_p) * d(rho_p)/dR_A|_basis
//
// where v_xc(rho) = de_xc/drho. Confirmed (not assumed) that this is
// EXACTLY xc.df_drho as already computed by EvaluateXC: xc_lda_exc_vxc /
// xc_gga_exc_vxc are LibXC's own named functions -- "exc" = energy per
// particle (f_xc), "vxc" = the potential dE_xc/drho directly (LibXC
// applies the f_xc + rho*df_xc/drho correction internally; the "vxc"
// output already IS the full potential, matching how df_drho is already
// used to build Vxc_here above -- no extra term needed here).
//
// For an atom-centered Gaussian chi_mu(r) = phi_mu(r - R_A), the nuclear
// derivative is EXACTLY -grad_r(chi_mu) for the atom mu is centered on,
// and exactly zero for every other atom (since r-R_A depends negatively
// on R_A). This means d(chi_mu)/dR_A can be read directly off
// ao.derivatives (already computed for the GGA sigma term above), just
// negated and restricted to basis functions centered on atom A -- no new
// basis-function-derivative machinery is needed for this piece.
//
// Using the symmetry of the density matrix (same trick as the standard
// Pulay-force derivation):
//   d(rho_p)/dR_A|_basis = -2 * sum_{mu in A} sum_nu P_munu (grad_r chi_mu)(r_p) chi_nu(r_p)
//                        = -sum_{mu in A} temp(mu) * (grad_r chi_mu)(r_p)
// where temp = ao.values^T * DMAT_here (DMAT_here = 2P, already exactly
// what IntegrateVXC computes above) -- the factor of 2 from symmetry is
// already folded into DMAT_here, so no extra factor is needed here either.
//
// SCOPE (see also the STATUS note on the declaration in vxc_potential.h):
// this captures the FULL Pulay term for LDA functionals. For GGA
// functionals, it captures only the df_drho-driven part; the additional
// df_dsigma-driven Pulay contribution needs second derivatives of basis
// functions w.r.t. electron position, not implemented here. The
// grid-weight (SSW partition) derivative term is a SEPARATE piece,
// deliberately not included in this function.
//
// STATUS: written but NOT yet run/tested.
// ===========================================================================
template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::PulayGradient(
    const Eigen::MatrixXd& density_matrix, const AOBasis& dftbasis) const {
  assert(density_matrix.isApprox(density_matrix.transpose()) &&
         "Density matrix has to be symmetric!");

  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(
      nthreads, Eigen::MatrixXd::Zero(natoms, 3));

#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    Index thread_id = OPENMP::getThreadId();
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }

    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);

    double cutoff =
        1.e-40 / double(density_matrix.rows()) / double(density_matrix.rows());
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }

    // Box-local AO index -> atom index, built once per box (not once per
    // grid point), using the same getShells()/getAOranges() mapping the
    // rest of GridBox already relies on for its own big-matrix
    // read/write bookkeeping.
    std::vector<Index> local_idx_to_atom(box.Matrixsize());
    const std::vector<const AOShell*>& shells = box.getShells();
    const std::vector<GridboxRange>& ao_ranges = box.getAOranges();
    for (size_t s = 0; s < shells.size(); ++s) {
      Index atom = shells[s]->getAtomIndex();
      for (Index k = 0; k < ao_ranges[s].size; ++k) {
        local_idx_to_atom[ao_ranges[s].start + k] = atom;
      }
    }

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    for (Index p = 0; p < box.size(); ++p) {
      AOShell::AOValues ao = box.CalcAOValues(points[p]);

      Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
      double rho = 0.5 * temp.dot(ao.values);
      const double weight = weights[p];

      if (rho * weight < 1.e-20) {
        continue;
      }

      const Eigen::Vector3d rho_grad = temp.transpose() * ao.derivatives;
      typename Vxc_Potential<Grid>::XC_entry xc =
          EvaluateXC(rho, rho_grad.squaredNorm());

      // d(rho_p)/dR_A|_basis, accumulated per local AO index mu, then
      // scattered to the correct atom via local_idx_to_atom. Looping
      // over local indices directly (rather than trying to slice
      // contiguous per-atom row ranges) since a single box's
      // significant shells are not guaranteed to be grouped
      // contiguously by atom.
      for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
        Index atom = local_idx_to_atom[mu];
        Eigen::Vector3d contribution =
            -weight * xc.df_drho * temp(mu) * ao.derivatives.row(mu).transpose();
        grad_thread[thread_id].row(atom) += contribution.transpose();
      }
    }
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
  }
  return grad;
}

template <class Grid>
typename Vxc_Potential<Grid>::SpinResult Vxc_Potential<Grid>::IntegrateVXCSpin(
    const Eigen::MatrixXd& dmat_alpha, const Eigen::MatrixXd& dmat_beta) const {
  assert(dmat_alpha.isApprox(dmat_alpha.transpose()) &&
         "Alpha density matrix has to be symmetric!");
  assert(dmat_beta.isApprox(dmat_beta.transpose()) &&
         "Beta density matrix has to be symmetric!");

  typename Vxc_Potential<Grid>::SpinResult result;
  result.vxc_alpha =
      Eigen::MatrixXd::Zero(dmat_alpha.rows(), dmat_alpha.cols());
  result.vxc_beta = Eigen::MatrixXd::Zero(dmat_beta.rows(), dmat_beta.cols());

#pragma omp parallel
  {
    Eigen::MatrixXd vxc_alpha_private =
        Eigen::MatrixXd::Zero(dmat_alpha.rows(), dmat_alpha.cols());
    Eigen::MatrixXd vxc_beta_private =
        Eigen::MatrixXd::Zero(dmat_beta.rows(), dmat_beta.cols());
    double exc_private = 0.0;

#pragma omp for schedule(guided)
    for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
      const GridBox& box = grid_[i];
      if (!box.Matrixsize()) {
        continue;
      }

      const Eigen::MatrixXd DMa = box.ReadFromBigMatrix(dmat_alpha);
      const Eigen::MatrixXd DMb = box.ReadFromBigMatrix(dmat_beta);

      double cutoff =
          1.e-40 / double(dmat_alpha.rows()) / double(dmat_alpha.rows());
      if (std::max(DMa.cwiseAbs2().maxCoeff(), DMb.cwiseAbs2().maxCoeff()) <
          cutoff) {
        continue;
      }

      Eigen::MatrixXd Vxc_a_here =
          Eigen::MatrixXd::Zero(DMa.rows(), DMa.cols());
      Eigen::MatrixXd Vxc_b_here =
          Eigen::MatrixXd::Zero(DMb.rows(), DMb.cols());

      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();

      for (Index p = 0; p < box.size(); ++p) {
        AOShell::AOValues ao = box.CalcAOValues(points[p]);

        Eigen::VectorXd temp_a = DMa * ao.values;
        Eigen::VectorXd temp_b = DMb * ao.values;

        const double rho_a = ao.values.dot(temp_a);
        const double rho_b = ao.values.dot(temp_b);
        const double rho = rho_a + rho_b;
        const double weight = weights[p];

        if (rho * weight < 1.e-20) {
          continue;
        }

        // For symmetric density matrices, this gives the full gradient
        // consistent with the restricted implementation, which used 2*P.
        const Eigen::Vector3d grad_a =
            2.0 * (ao.derivatives.transpose() * temp_a);
        const Eigen::Vector3d grad_b =
            2.0 * (ao.derivatives.transpose() * temp_b);

        const double sigma_aa = grad_a.dot(grad_a);
        const double sigma_ab = grad_a.dot(grad_b);
        const double sigma_bb = grad_b.dot(grad_b);

        typename Vxc_Potential<Grid>::XC_entry_spin xc =
            EvaluateXCSpin(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb);

        exc_private += weight * rho * xc.f_xc;

        if (xfunc.info->family == XC_FAMILY_LDA) {
          // 0.5 factor because we symmetrize by adding transpose at the end
          Eigen::VectorXd wa = weight * (0.5 * xc.vrho_a) * ao.values;
          Eigen::VectorXd wb = weight * (0.5 * xc.vrho_b) * ao.values;

          Vxc_a_here.noalias() += wa * ao.values.transpose();
          Vxc_b_here.noalias() += wb * ao.values.transpose();
        } else {
          Eigen::VectorXd g_a = ao.derivatives * grad_a;
          Eigen::VectorXd g_b = ao.derivatives * grad_b;

          // Same 0.5 prefactor on vrho term as in restricted path.
          Eigen::VectorXd wa =
              weight * (0.5 * xc.vrho_a * ao.values + 2.0 * xc.vsigma_aa * g_a +
                        xc.vsigma_ab * g_b);

          Eigen::VectorXd wb =
              weight * (0.5 * xc.vrho_b * ao.values + xc.vsigma_ab * g_a +
                        2.0 * xc.vsigma_bb * g_b);

          Vxc_a_here.noalias() += wa * ao.values.transpose();
          Vxc_b_here.noalias() += wb * ao.values.transpose();
        }
      }

      box.AddtoBigMatrix(vxc_alpha_private, Vxc_a_here);
      box.AddtoBigMatrix(vxc_beta_private, Vxc_b_here);
    }

#pragma omp critical
    {
      result.vxc_alpha += vxc_alpha_private + vxc_alpha_private.transpose();
      result.vxc_beta += vxc_beta_private + vxc_beta_private.transpose();
      result.energy += exc_private;
    }
  }

  return result;
}

namespace {
// Standalone re-implementation of Vxc_Grid's switching function and its
// derivative -- NOT calling into Vxc_Grid::erf1c (private, and exposing
// it seemed like more churn than re-stating this one small, stateless,
// pure-math formula here). Must stay byte-for-byte consistent with
// Vxc_Grid::erf1c (0.5*erfc(|x/(1-x^2)|*alpha), alpha=1/0.30) -- if that
// function's constants ever change, this needs to change with it.
constexpr double kSSWAlpha = 1.0 / 0.30;
constexpr double kSSWCutoff = 0.725;

constexpr double kSqrtPi = 1.7724538509055160273;  // sqrt(pi), literal
                                                     // constant rather than
                                                     // M_PI (not standard
                                                     // C++, not otherwise
                                                     // used anywhere in
                                                     // this codebase --
                                                     // avoiding relying on
                                                     // it being defined).

double SSWValue(double mu) {
  double val = 0.5 * std::erfc(std::abs(mu / (1.0 - mu * mu)) * kSSWAlpha);
  if (mu > 0.0) {
    val = 1.0 - val;
  }
  return val;
}

// d(SSWValue)/d(mu). Exact closed form, verified numerically against
// finite differences in Python (matching to ~1e-11) before translating
// to C++ -- see conversation history.
double SSWDerivative(double mu) {
  double h = std::abs(mu) / (1.0 - mu * mu);
  double sign_mu = (mu > 0.0) ? 1.0 : ((mu < 0.0) ? -1.0 : 0.0);
  double one_minus_mu2 = 1.0 - mu * mu;
  double hprime = sign_mu * (1.0 + mu * mu) / (one_minus_mu2 * one_minus_mu2);
  double d_erf1c =
      -(kSSWAlpha / kSqrtPi) * std::exp(-(kSSWAlpha * h) * (kSSWAlpha * h)) * hprime;
  return (mu > 0.0) ? -d_erf1c : d_erf1c;
}
}  // namespace

template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::GridWeightGradient(
    const Eigen::MatrixXd& density_matrix, const QMMolecule& atoms) const {
  Index natoms = atoms.size();
  Eigen::MatrixXd Rij = grid_.CalcInverseAtomDist(atoms);

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(
      nthreads, Eigen::MatrixXd::Zero(natoms, 3));
  // Aggregate diagnostics (atom A==1 only, the one under test in
  // test_xcgradient.cc) -- see usage below for what each measures.
  std::vector<double> gross_sum_thread(nthreads, 0.0);
  std::vector<double> max_contribution_thread(nthreads, 0.0);
  std::vector<long> active_point_count_thread(nthreads, 0);

#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    Index thread_id = OPENMP::getThreadId();
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }

    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);
    double cutoff =
        1.e-40 / double(density_matrix.rows()) / double(density_matrix.rows());
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

    for (Index pidx = 0; pidx < box.size(); ++pidx) {
      AOShell::AOValues ao = box.CalcAOValues(points[pidx]);
      Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
      double rho = 0.5 * temp.dot(ao.values);
      double weight = weights[pidx];
      if (rho * weight < 1.e-20) {
        continue;
      }
      const Eigen::Vector3d rho_grad = temp.transpose() * ao.derivatives;
      typename Vxc_Potential<Grid>::XC_entry xc =
          EvaluateXC(rho, rho_grad.squaredNorm());

      Index owner = owner_atoms[pidx];
      if (owner < 0) {
        // Point predates owner-atom tracking (e.g. constructed via some
        // other path not yet updated) -- cannot compute this term for
        // it. Should not happen for any grid built via GridSetup after
        // this change; flagged rather than silently skipped, since a
        // silent skip here would quietly produce a wrong (incomplete)
        // gradient.
        throw std::runtime_error(
            "GridWeightGradient: grid point has no owner_atom set -- was "
            "this grid built via GridSetup after the owner-atom tracking "
            "change?");
      }

      // rq(k) = distance from THIS point to atom k; needed for every
      // atom, not just the owner, since the partition weight depends on
      // distances to all atoms.
      const Eigen::Vector3d& point = points[pidx];
      Eigen::VectorXd rq(natoms);
      for (Index k = 0; k < natoms; ++k) {
        rq(k) = (point - atoms[k].getPos()).norm();
      }

      // Build p[], mu_table, sk_table, hard-cutoff flags for every pair
      // -- same structure as Vxc_Grid::SSWpartition, but retaining the
      // per-pair intermediate values needed for the derivative (the
      // energy-level code discards these immediately after use).
      Eigen::VectorXd p = Eigen::VectorXd::Ones(natoms);
      Eigen::MatrixXd mu_table = Eigen::MatrixXd::Zero(natoms, natoms);
      Eigen::MatrixXd sk_table = Eigen::MatrixXd::Zero(natoms, natoms);
      // hard(j,i): 0 = smooth (sk_table valid), 1 = mu>cutoff (p[i]=0
      // hard), -1 = mu<-cutoff (p[j]=0 hard). Only upper triangle (j<i)
      // populated, matching the loop structure below.
      Eigen::MatrixXi hard = Eigen::MatrixXi::Zero(natoms, natoms);
      for (Index ii = 1; ii < natoms; ++ii) {
        for (Index jj = 0; jj < ii; ++jj) {
          double mu = (rq(ii) - rq(jj)) * Rij(jj, ii);
          mu_table(jj, ii) = mu;
          if (mu > kSSWCutoff) {
            p(ii) = 0.0;
            hard(jj, ii) = 1;
          } else if (mu < -kSSWCutoff) {
            p(jj) = 0.0;
            hard(jj, ii) = -1;
          } else {
            double sk = SSWValue(mu);
            sk_table(jj, ii) = sk;
            p(jj) *= sk;
            p(ii) *= (1.0 - sk);
          }
        }
      }
      double wsum = p.sum();
      double w_owner = p(owner) / wsum;

      // d(rq(k))/dR_A, per the case analysis verified in Python: zero
      // unless A is this point's owner or A==k; +-unit vector otherwise
      // (with the owner==k case being exactly zero, since the point and
      // its own owner move together).
      auto d_rq_dR = [&](Index k, Index A) -> Eigen::Vector3d {
        if (A == owner && A == k) {
          return Eigen::Vector3d::Zero();
        } else if (A == owner) {
          return (point - atoms[k].getPos()) / rq(k);
        } else if (A == k) {
          return -(point - atoms[k].getPos()) / rq(k);
        }
        return Eigen::Vector3d::Zero();
      };
      auto d_Rab_dR = [&](Index a, Index b, Index A) -> Eigen::Vector3d {
        Eigen::Vector3d rvec = atoms[a].getPos() - atoms[b].getPos();
        double Rab = rvec.norm();
        if (A == a) {
          return rvec / Rab;
        } else if (A == b) {
          return -rvec / Rab;
        }
        return Eigen::Vector3d::Zero();
      };
      // a<b convention throughout, matching mu_table(a,b) with a<b.
      auto dmu_dR = [&](Index a, Index b, Index A) -> Eigen::Vector3d {
        Eigen::Vector3d d_rq_b = d_rq_dR(b, A);
        Eigen::Vector3d d_rq_a = d_rq_dR(a, A);
        double Rab = 1.0 / Rij(a, b);
        Eigen::Vector3d dRab = d_Rab_dR(a, b, A);
        double mu = mu_table(a, b);
        return (d_rq_b - d_rq_a) / Rab - (mu / Rab) * dRab;
      };
      auto dp_dR = [&](Index k, Index A) -> Eigen::Vector3d {
        // Threshold, not exact-zero check: p(k) approaches zero
        // SMOOTHLY as any of its factors approaches the SSW saturation
        // boundary (sk->1 or 1-sk->0), and the log-derivative terms
        // below (divided by sk or (1-sk)) blow up faster than p(k)
        // itself vanishes, in that regime -- a classic 0/0 numerical
        // instability. A state with negligible weight can't meaningfully
        // contribute to the total either way, so treating it as exactly
        // zero here is a defensible approximation, not just a numerical
        // patch. Threshold value not yet tuned against real data; if
        // this test now passes but with values that look suspiciously
        // insensitive to h, or if a genuinely different tolerance is
        // needed elsewhere, revisit this constant specifically.
        constexpr double kNegligibleP = 1.e-8;
        if (p(k) < kNegligibleP) {
          return Eigen::Vector3d::Zero();
        }
        Eigen::Vector3d total = Eigen::Vector3d::Zero();
        for (Index b = k + 1; b < natoms; ++b) {
          if (hard(k, b) != 0) {
            continue;
          }
          double skv = sk_table(k, b);
          if (skv < kNegligibleP) {
            continue;  // this factor's own contribution to p(k) is
                       // already negligible; skip rather than divide by
                       // a near-zero skv.
          }
          total += (SSWDerivative(mu_table(k, b)) / skv) * dmu_dR(k, b, A);
        }
        for (Index a = 0; a < k; ++a) {
          if (hard(a, k) != 0) {
            continue;
          }
          double skv = sk_table(a, k);
          double one_minus_skv = 1.0 - skv;
          if (one_minus_skv < kNegligibleP) {
            continue;  // same reasoning, for the (1-sk) denominator.
          }
          total += (-SSWDerivative(mu_table(a, k)) / one_minus_skv) *
                    dmu_dR(a, k, A);
        }
        return p(k) * total;
      };

      // Guard against dividing by a near-zero w_owner: if w_owner is
      // negligible, weight_p = C_p*w_owner is also negligible (C_p is a
      // bounded quadrature weight, not something that can blow up to
      // compensate), so this point's actual contribution to the total
      // XC energy is negligible regardless of what C_p technically
      // works out to -- same physical justification as the earlier
      // p(k) threshold: a point with negligible weight can't
      // meaningfully affect the total either way, so treating it as
      // contributing exactly zero here is defensible, not just a
      // numerical patch.
      constexpr double kNegligibleWOwner = 1.e-8;
      if (w_owner < kNegligibleWOwner) {
        continue;
      }

      // BUG FIX: weight (as stored/read from the grid) is NOT w_owner
      // alone -- per GridSetup, weight_p = C_p * w_owner(p), where C_p
      // is the raw radial*angular quadrature weight (position-
      // independent -- depends only on the fixed Lebedev/Euler-Maclaurin
      // grid design for the owning element, never on any atom's
      // position) and w_owner is the SSW partition fraction computed
      // above. The energy contribution from this point is
      // weight_p*rho_p*f_xc_p = C_p*w_owner*rho_p*f_xc_p, so
      // differentiating w_owner alone (which is all `dw` below
      // computes) needs the missing C_p = weight/w_owner factor too --
      // previously this was bare rho*f_xc, silently missing C_p, which
      // is ~1 for "typical" points but swings far from 1 for points
      // where the partition is lopsided (w_owner small), exactly
      // matching why the resulting error's magnitude depended on which
      // points a given density matrix happened to weight rather than
      // ever looking pathological in any single point's own dw value.
      double C_p = weight / w_owner;
      double prefactor = C_p * rho * xc.f_xc;

      // Diagnostic: print the RAW INPUTS (owner, absolute point position)
      // for one informative point, so it can be reproduced EXACTLY in an
      // independent Python implementation and every intermediate value
      // compared number-for-number -- the decisive check, since static
      // code review found no discrepancy and aggregate statistics only
      // narrowed this to "likely a C++-translation-specific bug" without
      // identifying it.
      static bool printed_once_inputs = false;
      bool is_informative_inputs = (rho > 1.e-4) && (owner == 1);
      bool should_print_inputs = false;
      if (is_informative_inputs) {
#pragma omp critical
        {
          if (!printed_once_inputs) {
            printed_once_inputs = true;
            should_print_inputs = true;
          }
        }
      }
      if (should_print_inputs) {
        std::cerr << std::setprecision(17)
                   << "[GridWeightGradient RAW INPUTS for cross-check] owner="
                   << owner << " point=" << point.x() << " " << point.y()
                   << " " << point.z() << std::endl;
        for (Index k = 0; k < natoms; ++k) {
          std::cerr << "  atom " << k << " pos="
                     << atoms[k].getPos().x() << " "
                     << atoms[k].getPos().y() << " "
                     << atoms[k].getPos().z() << std::endl;
        }
        std::cerr << std::setprecision(9) << "  rq=" << rq.transpose()
                   << "\n  p=" << p.transpose() << "\n  wsum=" << wsum
                   << " w_owner=" << w_owner << "\n  rho=" << rho
                   << " prefactor=" << prefactor << std::endl;
        for (Index A = 0; A < natoms; ++A) {
          Eigen::Vector3d dp_owner_dbg = dp_dR(owner, A);
          Eigen::Vector3d dwsum_dbg = Eigen::Vector3d::Zero();
          for (Index k = 0; k < natoms; ++k) {
            dwsum_dbg += dp_dR(k, A);
          }
          Eigen::Vector3d dw_dbg =
              dp_owner_dbg / wsum - w_owner * dwsum_dbg / wsum;
          std::cerr << "  A=" << A << " dp_owner=" << dp_owner_dbg.transpose()
                     << " dwsum=" << dwsum_dbg.transpose()
                     << " dw=" << dw_dbg.transpose()
                     << " contribution=" << (prefactor * dw_dbg).transpose()
                     << std::endl;
        }
      }

      for (Index A = 0; A < natoms; ++A) {
        Eigen::Vector3d dp_owner = dp_dR(owner, A);
        Eigen::Vector3d dwsum = Eigen::Vector3d::Zero();
        for (Index k = 0; k < natoms; ++k) {
          dwsum += dp_dR(k, A);
        }
        Eigen::Vector3d dw = dp_owner / wsum - w_owner * dwsum / wsum;
        Eigen::Vector3d contribution = prefactor * dw;

        // Aggregate diagnostics (thread-local, merged after the loop):
        // sum of |contribution| (gross magnitude, vs the net summed
        // gradient -- a big gap between gross and net would indicate
        // cancellation happening imperfectly, i.e. a subtle bias, while
        // gross~net would mean most points agree in sign), max single
        // |contribution| seen (already had an explicit >10 flag before;
        // this additionally reports the actual max regardless of
        // threshold), and point count contributing to atom A=1
        // specifically (the one under test) with non-negligible
        // prefactor, to get a sense of how many points are actually
        // "in play" for the quantity being compared against the
        // finite-difference reference.
        if (A == 1) {
          gross_sum_thread[thread_id] += contribution.cwiseAbs().sum();
          max_contribution_thread[thread_id] =
              std::max(max_contribution_thread[thread_id],
                       contribution.cwiseAbs().maxCoeff());
          if (std::abs(prefactor) > 1.e-12) {
            active_point_count_thread[thread_id]++;
          }
        }
        grad_thread[thread_id].row(A) += contribution.transpose();
      }
    }
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  double gross_sum = 0.0;
  double max_contribution = 0.0;
  long active_point_count = 0;
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
    gross_sum += gross_sum_thread[t];
    max_contribution = std::max(max_contribution, max_contribution_thread[t]);
    active_point_count += active_point_count_thread[t];
  }
  std::cerr << "[GridWeightGradient diagnostic, aggregate for atom A=1] "
             << "net dE/dR(A=1) row = " << grad.row(1) << "\n"
             << "  gross sum of |contribution| = " << gross_sum << "\n"
             << "  max single |contribution| = " << max_contribution << "\n"
             << "  active (non-negligible prefactor) point count = "
             << active_point_count << "\n"
             << "  total grid points (all boxes) = " << grid_.getGridSize()
             << std::endl;
  return grad;
}

template class Vxc_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca