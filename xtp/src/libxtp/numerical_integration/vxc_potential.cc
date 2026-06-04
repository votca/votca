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
#include <stdexcept>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
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

template class Vxc_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca