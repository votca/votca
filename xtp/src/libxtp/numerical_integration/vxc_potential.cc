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
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/gridboxaccumulator.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/vxc_functionals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

namespace votca {
namespace xtp {

namespace {
// Range-separated hybrids (HSE, CAM-B3LYP, ...) need erf-attenuated exact
// exchange, which neither the DFT nor GW parts implement: the code would use
// cam_alpha as a global exchange fraction and run a different functional
// without notice. Refuse them.
void RejectRangeSeparated(const xc_func_type& func, const std::string& name) {
  if (func.cam_beta != 0.0 || func.cam_omega != 0.0) {
    throw std::runtime_error(
        "Functional " + name +
        " is a range-separated hybrid; range-separated exact exchange is not "
        "implemented in xtp (DFT and GW). Use a global hybrid or a semilocal "
        "functional.");
  }
}
}  // namespace

template <class Grid>
Vxc_Potential<Grid>::~Vxc_Potential() {
  if (setXC_) {
    xc_func_end(&xfunc);
    xc_func_end(&xfunc_pol_);
    if (use_separate_) {
      xc_func_end(&cfunc);
      xc_func_end(&cfunc_pol_);
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

    try {
      RejectRangeSeparated(func, functional_name);
    } catch (...) {
      xc_func_end(&func);
      throw;
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

  try {
    RejectRangeSeparated(xfunc, strs[0]);
  } catch (...) {
    xc_func_end(&xfunc);
    throw;
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

    try {
      RejectRangeSeparated(cfunc, strs[1]);
    } catch (...) {
      xc_func_end(&xfunc);
      xc_func_end(&cfunc);
      throw;
    }

    if ((xfunc.info->kind + cfunc.info->kind) != 1) {
      xc_func_end(&xfunc);
      xc_func_end(&cfunc);
      throw std::runtime_error(
          "Your functionals are not one exchange and one correlation");
    }
  }

  // Polarized twins for the UKS path.
  if (xc_func_init(&xfunc_pol_, xfunc_id, XC_POLARIZED) != 0) {
    throw std::runtime_error(
        "Failed to initialize polarized exchange XC "
        "functional.");
  }
  if (use_separate_ && xc_func_init(&cfunc_pol_, cfunc_id, XC_POLARIZED) != 0) {
    xc_func_end(&xfunc_pol_);
    throw std::runtime_error(
        "Failed to initialize polarized correlation XC "
        "functional.");
  }

  auto is_gga = [](const xc_func_type& f) {
    return f.info->family == XC_FAMILY_GGA ||
           f.info->family == XC_FAMILY_HYB_GGA;
  };
  is_gga_ = is_gga(xfunc) || (use_separate_ && is_gga(cfunc));

  setXC_ = true;
  return;
}

namespace {
// One libxc call for n points with any supported family.
void CallLibxc(const xc_func_type& f, Index n, const double* rho,
               const double* sigma, double* exc, double* vrho, double* vsigma) {
  switch (f.info->family) {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&f, std::size_t(n), rho, exc, vrho);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&f, std::size_t(n), rho, sigma, exc, vrho, vsigma);
      break;
    default:
      throw std::runtime_error("Unsupported XC functional family.");
  }
}
}  // namespace

template <class Grid>
void Vxc_Potential<Grid>::EvaluateXCBatch(Index n, const double* rho,
                                          const double* sigma, double* exc,
                                          double* vrho, double* vsigma) const {
  std::fill(vsigma, vsigma + n, 0.0);
  CallLibxc(xfunc, n, rho, sigma, exc, vrho, vsigma);
  if (use_separate_) {
    std::vector<double> e(n), vr(n), vs(n, 0.0);
    CallLibxc(cfunc, n, rho, sigma, e.data(), vr.data(), vs.data());
    for (Index i = 0; i < n; ++i) {
      exc[i] += e[i];
      vrho[i] += vr[i];
      vsigma[i] += vs[i];
    }
  }
}

template <class Grid>
void Vxc_Potential<Grid>::EvaluateXCSpinBatch(Index n, const double* rho,
                                              const double* sigma, double* exc,
                                              double* vrho,
                                              double* vsigma) const {
  std::fill(vsigma, vsigma + 3 * n, 0.0);
  CallLibxc(xfunc_pol_, n, rho, sigma, exc, vrho, vsigma);
  if (use_separate_) {
    std::vector<double> e(n), vr(2 * n), vs(3 * n, 0.0);
    CallLibxc(cfunc_pol_, n, rho, sigma, e.data(), vr.data(), vs.data());
    for (Index i = 0; i < n; ++i) {
      exc[i] += e[i];
    }
    for (Index i = 0; i < 2 * n; ++i) {
      vrho[i] += vr[i];
    }
    for (Index i = 0; i < 3 * n; ++i) {
      vsigma[i] += vs[i];
    }
  }
}

template <class Grid>
typename Vxc_Potential<Grid>::XC_entry Vxc_Potential<Grid>::EvaluateXC(
    double rho, double sigma) const {
  typename Vxc_Potential<Grid>::XC_entry result;
  EvaluateXCBatch(1, &rho, &sigma, &result.f_xc, &result.df_drho,
                  &result.df_dsigma);
  return result;
}

template <class Grid>
typename Vxc_Potential<Grid>::XC_entry_spin Vxc_Potential<Grid>::EvaluateXCSpin(
    double rho_a, double rho_b, double sigma_aa, double sigma_ab,
    double sigma_bb) const {
  typename Vxc_Potential<Grid>::XC_entry_spin result;
  const double rho[2] = {rho_a, rho_b};
  const double sigma[3] = {sigma_aa, sigma_ab, sigma_bb};
  double vrho[2];
  double vsigma[3];
  EvaluateXCSpinBatch(1, rho, sigma, &result.f_xc, vrho, vsigma);
  result.vrho_a = vrho[0];
  result.vrho_b = vrho[1];
  result.vsigma_aa = vsigma[0];
  result.vsigma_ab = vsigma[1];
  result.vsigma_bb = vsigma[2];
  return result;
}

// Points are processed in blocks inside each box. With Phi the box's AO
// values (functions x points) and dPhi their gradients, the work that scales
// with the square of the box's basis size is two matrix products per block:
//   Y = D Phi                    (density and its gradient at every point)
//   Vxc += A Phi^T,  A_p = w_p (vrho_p/2 phi_p + 2 vsigma_p dPhi_p grad rho_p)
// and libxc is called once per block. Per point, only O(functions) work
// remains. Same formulas as the former point-by-point loop.
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

    // two because we have to use the density matrix and its transpose
    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);

    double cutoff =
        1.e-40 / double(density_matrix.rows()) / double(density_matrix.rows());
    if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
      continue;
    }

    const Index nf = box.Matrixsize();
    const Index block = std::min(kPointBlock, box.size());
    Eigen::MatrixXd Vxc_here = Eigen::MatrixXd::Zero(nf, nf);
    Eigen::MatrixXd phi(nf, block);
    Eigen::MatrixXd dphi(nf, 3 * block);
    Eigen::MatrixXd y(nf, block);
    Eigen::MatrixXd a(nf, block);
    std::vector<Index> kept;
    std::vector<double> rho, sigma, exc, vrho, vsigma;
    std::vector<Eigen::Vector3d> grad;
    double EXC_box = 0.0;

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    for (Index start = 0; start < box.size(); start += block) {
      const Index nb = std::min(block, box.size() - start);
      for (Index p = 0; p < nb; ++p) {
        box.CalcAOValues(points[start + p], phi, dphi, p);
      }
      y.leftCols(nb).noalias() = DMAT_here * phi.leftCols(nb);

      kept.clear();
      rho.clear();
      sigma.clear();
      grad.clear();
      for (Index p = 0; p < nb; ++p) {
        const double r = 0.5 * phi.col(p).dot(y.col(p));
        if (r * weights[start + p] < 1.e-20) {
          continue;
        }
        const Eigen::Vector3d g =
            dphi.middleCols(3 * p, 3).transpose() * y.col(p);
        kept.push_back(p);
        rho.push_back(r);
        sigma.push_back(g.squaredNorm());
        grad.push_back(g);
      }
      const Index nk = Index(kept.size());
      if (nk == 0) {
        continue;
      }
      exc.resize(nk);
      vrho.resize(nk);
      vsigma.resize(nk);
      EvaluateXCBatch(nk, rho.data(), sigma.data(), exc.data(), vrho.data(),
                      vsigma.data());

      // Kept points are compacted to the front: column k receives point
      // kept[k] >= k, which has been read by then.
      for (Index k = 0; k < nk; ++k) {
        const Index p = kept[k];
        const double w = weights[start + p];
        EXC_box += w * rho[k] * exc[k];
        a.col(k) = (w * 0.5 * vrho[k]) * phi.col(p);
        if (is_gga_) {
          a.col(k).noalias() +=
              (w * 2.0 * vsigma[k]) * (dphi.middleCols(3 * p, 3) * grad[k]);
        }
        if (k != p) {
          phi.col(k) = phi.col(p);
        }
      }
      Vxc_here.noalias() += a.leftCols(nk) * phi.leftCols(nk).transpose();
    }

    box.AddtoBigMatrix(vxc.matrix(), Vxc_here);
    vxc.energy() += EXC_box;
  }

  return Mat_p_Energy(vxc.energy(), vxc.matrix() + vxc.matrix().transpose());
}

namespace {
// first and second derivatives of one functional (unpolarized)
void CallLibxcFXC(const xc_func_type& f, Index n, const double* rho,
                  const double* sigma, double* vrho, double* vsigma,
                  double* v2rho2, double* v2rhosigma, double* v2sigma2) {
  if (!(f.info->flags & XC_FLAGS_HAVE_FXC)) {
    throw std::runtime_error(std::string("Functional ") + f.info->name +
                             " has no second derivatives in libxc");
  }
  switch (f.info->family) {
    case XC_FAMILY_LDA:
      xc_lda_vxc_fxc(&f, std::size_t(n), rho, vrho, v2rho2);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      xc_gga_vxc_fxc(&f, std::size_t(n), rho, sigma, vrho, vsigma, v2rho2,
                     v2rhosigma, v2sigma2);
      break;
    default:
      throw std::runtime_error("Unsupported XC functional family.");
  }
}
}  // namespace

template <class Grid>
void Vxc_Potential<Grid>::EvaluateFXCBatch(Index n, const double* rho,
                                           const double* sigma, double* vsigma,
                                           double* v2rho2, double* v2rhosigma,
                                           double* v2sigma2) const {
  std::vector<double> vrho(n);
  EvaluateVXCFXCBatch(n, rho, sigma, vrho.data(), vsigma, v2rho2, v2rhosigma,
                      v2sigma2);
}

template <class Grid>
void Vxc_Potential<Grid>::EvaluateVXCFXCBatch(Index n, const double* rho,
                                              const double* sigma, double* vrho,
                                              double* vsigma, double* v2rho2,
                                              double* v2rhosigma,
                                              double* v2sigma2) const {
  std::fill(vsigma, vsigma + n, 0.0);
  std::fill(v2rhosigma, v2rhosigma + n, 0.0);
  std::fill(v2sigma2, v2sigma2 + n, 0.0);
  CallLibxcFXC(xfunc, n, rho, sigma, vrho, vsigma, v2rho2, v2rhosigma,
               v2sigma2);
  if (use_separate_) {
    std::vector<double> vr(n), vs(n, 0.0), rr(n), rs(n, 0.0), ss(n, 0.0);
    CallLibxcFXC(cfunc, n, rho, sigma, vr.data(), vs.data(), rr.data(),
                 rs.data(), ss.data());
    for (Index i = 0; i < n; ++i) {
      vrho[i] += vr[i];
      vsigma[i] += vs[i];
      v2rho2[i] += rr[i];
      v2rhosigma[i] += rs[i];
      v2sigma2[i] += ss[i];
    }
  }
}

template <class Grid>
void Vxc_Potential<Grid>::EvaluateFXCSpinBatch(Index n, const double* rho,
                                               const double* sigma,
                                               double* vsigma, double* v2rho2,
                                               double* v2rhosigma,
                                               double* v2sigma2) const {
  std::vector<double> vrho(2 * n);
  std::fill(vsigma, vsigma + 3 * n, 0.0);
  std::fill(v2rhosigma, v2rhosigma + 6 * n, 0.0);
  std::fill(v2sigma2, v2sigma2 + 6 * n, 0.0);
  CallLibxcFXC(xfunc_pol_, n, rho, sigma, vrho.data(), vsigma, v2rho2,
               v2rhosigma, v2sigma2);
  if (use_separate_) {
    std::vector<double> vr(2 * n), vs(3 * n, 0.0), rr(3 * n), rs(6 * n, 0.0),
        ss(6 * n, 0.0);
    CallLibxcFXC(cfunc_pol_, n, rho, sigma, vr.data(), vs.data(), rr.data(),
                 rs.data(), ss.data());
    for (Index i = 0; i < 3 * n; ++i) {
      vsigma[i] += vs[i];
      v2rho2[i] += rr[i];
    }
    for (Index i = 0; i < 6 * n; ++i) {
      v2rhosigma[i] += rs[i];
      v2sigma2[i] += ss[i];
    }
  }
}

// First-order change of IntegrateVXCSpin. Per point, with g_s = grad rho_s,
// sigma = (g_a.g_a, g_a.g_b, g_b.g_b) and the first-order changes drho_s,
// dg_s = grad drho_s:
//   dsigma = (2 g_a.dg_a, g_a.dg_b + dg_a.g_b, 2 g_b.dg_b)
//   dv_rho_s   = sum_t v2rho2_st drho_t + sum_c v2rhosigma_sc dsigma_c
//   dv_sigma_c = sum_t v2rhosigma_tc drho_t + sum_d v2sigma2_cd dsigma_d
//   dV_a = int dv_rho_a chi chi + (2 dv_aa g_a + 2 v_aa dg_a + dv_ab g_b
//          + v_ab dg_b).grad(chi chi)
// and dV_b likewise (a <-> b, aa <-> bb); assembled as in IntegrateVXCSpin.
template <class Grid>
std::vector<std::array<Eigen::MatrixXd, 2>>
    Vxc_Potential<Grid>::IntegrateFXCSpin(
        const Eigen::MatrixXd& dmat_alpha, const Eigen::MatrixXd& dmat_beta,
        const std::vector<std::array<Eigen::MatrixXd, 2>>& response) const {
  const Index n_ao = dmat_alpha.rows();
  const Index n_resp = Index(response.size());
  std::vector<std::array<Eigen::MatrixXd, 2>> result(
      static_cast<std::size_t>(n_resp),
      {Eigen::MatrixXd::Zero(n_ao, n_ao), Eigen::MatrixXd::Zero(n_ao, n_ao)});
  if (n_resp == 0) {
    return result;
  }
  auto wrong_size = [n_ao](const std::array<Eigen::MatrixXd, 2>& r) {
    return r[0].rows() != n_ao || r[0].cols() != n_ao || r[1].rows() != n_ao ||
           r[1].cols() != n_ao;
  };
  if (std::any_of(response.begin(), response.end(), wrong_size)) {
    throw std::runtime_error(
        "IntegrateFXCSpin: response matrix does not match the density "
        "matrix");
  }
  // fail before the grid loop for functionals without fxc
  {
    const double rho[2] = {0.1, 0.05};
    const double sigma[3] = {0.01, 0.005, 0.004};
    double vs[3], rr[3], rs[6], ss[6];
    EvaluateFXCSpinBatch(1, rho, sigma, vs, rr, rs, ss);
  }

  // targets in the order alpha_0, beta_0, alpha_1, beta_1, ...
  std::vector<Eigen::MatrixXd*> targets;
  for (std::array<Eigen::MatrixXd, 2>& m : result) {
    targets.push_back(&m[0]);
    targets.push_back(&m[1]);
  }
  GridBoxAccumulator acc(targets);
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }
    const Eigen::MatrixXd DMa = box.ReadFromBigMatrix(dmat_alpha);
    const Eigen::MatrixXd DMb = box.ReadFromBigMatrix(dmat_beta);
    std::vector<Eigen::MatrixXd> dDa(static_cast<std::size_t>(n_resp));
    std::vector<Eigen::MatrixXd> dDb(static_cast<std::size_t>(n_resp));
    for (Index k = 0; k < n_resp; ++k) {
      dDa[std::size_t(k)] = box.ReadFromBigMatrix(response[std::size_t(k)][0]);
      dDb[std::size_t(k)] = box.ReadFromBigMatrix(response[std::size_t(k)][1]);
    }
    const Index nf = box.Matrixsize();
    const Index block = std::min(kPointBlock, box.size());
    std::vector<Eigen::MatrixXd> fxc_a(static_cast<std::size_t>(n_resp),
                                       Eigen::MatrixXd::Zero(nf, nf));
    std::vector<Eigen::MatrixXd> fxc_b(static_cast<std::size_t>(n_resp),
                                       Eigen::MatrixXd::Zero(nf, nf));
    Eigen::MatrixXd phi(nf, block);
    Eigen::MatrixXd dphi(nf, 3 * block);
    Eigen::MatrixXd ya(nf, block);
    Eigen::MatrixXd yb(nf, block);
    Eigen::MatrixXd yka(nf, block);
    Eigen::MatrixXd ykb(nf, block);
    Eigen::MatrixXd aa(nf, block);
    Eigen::MatrixXd ab(nf, block);
    std::vector<Index> kept;
    std::vector<double> rho, sigma, vsigma, v2rho2, v2rhosigma, v2sigma2;
    std::vector<Eigen::Vector3d> grad_a, grad_b;
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    for (Index start = 0; start < box.size(); start += block) {
      const Index nb = std::min(block, box.size() - start);
      for (Index p = 0; p < nb; ++p) {
        box.CalcAOValues(points[start + p], phi, dphi, p);
      }
      ya.leftCols(nb).noalias() = DMa * phi.leftCols(nb);
      yb.leftCols(nb).noalias() = DMb * phi.leftCols(nb);
      kept.clear();
      rho.clear();
      sigma.clear();
      grad_a.clear();
      grad_b.clear();
      for (Index p = 0; p < nb; ++p) {
        const double rho_a = phi.col(p).dot(ya.col(p));
        const double rho_b = phi.col(p).dot(yb.col(p));
        if ((rho_a + rho_b) * weights[start + p] < 1.e-20) {
          continue;
        }
        const Eigen::Vector3d ga =
            2.0 * (dphi.middleCols(3 * p, 3).transpose() * ya.col(p));
        const Eigen::Vector3d gb =
            2.0 * (dphi.middleCols(3 * p, 3).transpose() * yb.col(p));
        kept.push_back(p);
        rho.push_back(rho_a);
        rho.push_back(rho_b);
        sigma.push_back(ga.dot(ga));
        sigma.push_back(ga.dot(gb));
        sigma.push_back(gb.dot(gb));
        grad_a.push_back(ga);
        grad_b.push_back(gb);
      }
      const Index nk = Index(kept.size());
      if (nk == 0) {
        continue;
      }
      vsigma.resize(3 * nk);
      v2rho2.resize(3 * nk);
      v2rhosigma.resize(6 * nk);
      v2sigma2.resize(6 * nk);
      EvaluateFXCSpinBatch(nk, rho.data(), sigma.data(), vsigma.data(),
                           v2rho2.data(), v2rhosigma.data(), v2sigma2.data());
      Eigen::MatrixXd phik(nf, nk);
      for (Index q = 0; q < nk; ++q) {
        phik.col(q) = phi.col(kept[q]);
      }

      for (Index k = 0; k < n_resp; ++k) {
        yka.leftCols(nb).noalias() = dDa[std::size_t(k)] * phi.leftCols(nb);
        ykb.leftCols(nb).noalias() = dDb[std::size_t(k)] * phi.leftCols(nb);
        for (Index q = 0; q < nk; ++q) {
          const Index p = kept[q];
          const double w = weights[start + p];
          const double dra = phi.col(p).dot(yka.col(p));
          const double drb = phi.col(p).dot(ykb.col(p));
          const double* r2 = &v2rho2[3 * q];
          if (!is_gga_) {
            const double dva = r2[0] * dra + r2[1] * drb;
            const double dvb = r2[1] * dra + r2[2] * drb;
            aa.col(q) = (w * 0.5 * dva) * phi.col(p);
            ab.col(q) = (w * 0.5 * dvb) * phi.col(p);
            continue;
          }
          const auto d = dphi.middleCols(3 * p, 3);
          const Eigen::Vector3d dga = 2.0 * (d.transpose() * yka.col(p));
          const Eigen::Vector3d dgb = 2.0 * (d.transpose() * ykb.col(p));
          const Eigen::Vector3d& ga = grad_a[q];
          const Eigen::Vector3d& gb = grad_b[q];
          const double ds[3] = {2.0 * ga.dot(dga), ga.dot(dgb) + dga.dot(gb),
                                2.0 * gb.dot(dgb)};
          const double* rs = &v2rhosigma[6 * q];
          const double* s2 = &v2sigma2[6 * q];
          const double* vs = &vsigma[3 * q];
          // v2sigma2 as a symmetric 3x3 matrix
          const double S[3][3] = {{s2[0], s2[1], s2[2]},
                                  {s2[1], s2[3], s2[4]},
                                  {s2[2], s2[4], s2[5]}};
          const double dva = r2[0] * dra + r2[1] * drb + rs[0] * ds[0] +
                             rs[1] * ds[1] + rs[2] * ds[2];
          const double dvb = r2[1] * dra + r2[2] * drb + rs[3] * ds[0] +
                             rs[4] * ds[1] + rs[5] * ds[2];
          double dvs[3];
          for (int c = 0; c < 3; ++c) {
            dvs[c] = rs[c] * dra + rs[3 + c] * drb + S[c][0] * ds[0] +
                     S[c][1] * ds[1] + S[c][2] * ds[2];
          }
          const Eigen::Vector3d vec_a =
              2.0 * dvs[0] * ga + 2.0 * vs[0] * dga + dvs[1] * gb + vs[1] * dgb;
          const Eigen::Vector3d vec_b =
              2.0 * dvs[2] * gb + 2.0 * vs[2] * dgb + dvs[1] * ga + vs[1] * dga;
          aa.col(q) = (w * 0.5 * dva) * phi.col(p) + w * (d * vec_a);
          ab.col(q) = (w * 0.5 * dvb) * phi.col(p) + w * (d * vec_b);
        }
        fxc_a[std::size_t(k)].noalias() += aa.leftCols(nk) * phik.transpose();
        fxc_b[std::size_t(k)].noalias() += ab.leftCols(nk) * phik.transpose();
      }
    }
    std::vector<const Eigen::MatrixXd*> blocks;
    for (Index k = 0; k < n_resp; ++k) {
      blocks.push_back(&fxc_a[std::size_t(k)]);
      blocks.push_back(&fxc_b[std::size_t(k)]);
    }
    acc.Add(box, blocks);
  }
  acc.Finish();
  for (std::array<Eigen::MatrixXd, 2>& m : result) {
    m[0] = (m[0] + m[0].transpose()).eval();
    m[1] = (m[1] + m[1].transpose()).eval();
  }
  return result;
}

// First-order change of IntegrateVXC. With sigma = |grad rho|^2, at a point
//   drho, dg = grad drho, dsigma = 2 g.dg
//   dv_rho   = v2rho2 drho + v2rhosigma dsigma
//   dv_sigma = v2rhosigma drho + v2sigma2 dsigma
//   dV = int dv_rho chi chi + 2 (dv_sigma g + vsigma dg).grad(chi chi)
// assembled like IntegrateVXC: a = w (dv_rho/2 phi + 2 dphi.(...)),
// dV = a phi^T + transpose. The points IntegrateVXC skips are skipped.
template <class Grid>
typename Vxc_Potential<Grid>::FXCKernel Vxc_Potential<Grid>::PrepareFXC(
    const Eigen::MatrixXd& density_matrix) const {
  // fail before the grid loop for functionals without fxc
  {
    const double rho = 0.1;
    const double sigma = 0.01;
    double vs, rr, rs, ss;
    EvaluateFXCBatch(1, &rho, &sigma, &vs, &rr, &rs, &ss);
  }
  FXCKernel kernel;
  kernel.n_ao = density_matrix.rows();
  kernel.boxes.resize(std::size_t(grid_.getBoxesSize()));
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }
    typename FXCKernel::Box& kb = kernel.boxes[std::size_t(i)];
    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);
    const Index nf = box.Matrixsize();
    const Index block = std::min(kPointBlock, box.size());
    Eigen::MatrixXd phi(nf, block);
    Eigen::MatrixXd dphi(nf, 3 * block);
    Eigen::MatrixXd y(nf, block);
    std::vector<double> rho, sigma, w;
    std::vector<Eigen::Vector3d> grad;
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    for (Index start = 0; start < box.size(); start += block) {
      const Index nb = std::min(block, box.size() - start);
      for (Index p = 0; p < nb; ++p) {
        box.CalcAOValues(points[start + p], phi, dphi, p);
      }
      y.leftCols(nb).noalias() = DMAT_here * phi.leftCols(nb);
      for (Index p = 0; p < nb; ++p) {
        const double r = 0.5 * phi.col(p).dot(y.col(p));
        if (r * weights[start + p] < 1.e-20) {
          continue;
        }
        const Eigen::Vector3d g =
            dphi.middleCols(3 * p, 3).transpose() * y.col(p);
        kb.points.push_back(start + p);
        rho.push_back(r);
        sigma.push_back(g.squaredNorm());
        grad.push_back(g);
        w.push_back(weights[start + p]);
      }
    }
    const Index nk = Index(kb.points.size());
    kb.w = Eigen::Map<const Eigen::VectorXd>(w.data(), nk);
    kb.vsigma.resize(nk);
    kb.v2rho2.resize(nk);
    kb.v2rhosigma.resize(nk);
    kb.v2sigma2.resize(nk);
    if (nk > 0) {
      EvaluateFXCBatch(nk, rho.data(), sigma.data(), kb.vsigma.data(),
                       kb.v2rho2.data(), kb.v2rhosigma.data(),
                       kb.v2sigma2.data());
    }
    kb.grad.resize(3, nk);
    for (Index q = 0; q < nk; ++q) {
      kb.grad.col(q) = grad[std::size_t(q)];
    }
  }
  return kernel;
}

template <class Grid>
std::vector<Eigen::MatrixXd> Vxc_Potential<Grid>::ApplyFXC(
    const FXCKernel& kernel,
    const std::vector<Eigen::MatrixXd>& response_matrices) const {
  const Index n_ao = kernel.n_ao;
  const Index n_resp = Index(response_matrices.size());
  std::vector<Eigen::MatrixXd> result(static_cast<std::size_t>(n_resp),
                                      Eigen::MatrixXd::Zero(n_ao, n_ao));
  if (n_resp == 0) {
    return result;
  }
  if (std::any_of(response_matrices.begin(), response_matrices.end(),
                  [n_ao](const Eigen::MatrixXd& m) {
                    return m.rows() != n_ao || m.cols() != n_ao;
                  })) {
    throw std::runtime_error(
        "IntegrateFXC: response matrix does not match the density matrix");
  }
  std::vector<Eigen::MatrixXd*> targets;
  for (Eigen::MatrixXd& m : result) {
    targets.push_back(&m);
  }
  GridBoxAccumulator acc(targets);
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    const typename FXCKernel::Box& kb = kernel.boxes[std::size_t(i)];
    const Index nk_box = Index(kb.points.size());
    if (!box.Matrixsize() || nk_box == 0) {
      continue;
    }
    const Index nf = box.Matrixsize();
    std::vector<Eigen::MatrixXd> dmat_k(static_cast<std::size_t>(n_resp));
    for (Index k = 0; k < n_resp; ++k) {
      dmat_k[std::size_t(k)] =
          2 * box.ReadFromBigMatrix(response_matrices[std::size_t(k)]);
    }
    std::vector<Eigen::MatrixXd> fxc_here(static_cast<std::size_t>(n_resp),
                                          Eigen::MatrixXd::Zero(nf, nf));
    const Index block = std::min(kPointBlock, nk_box);
    Eigen::MatrixXd phi_buf(nf, block);
    Eigen::MatrixXd dphi_buf(nf, 3 * block);
    Eigen::MatrixXd yk(nf, block);
    Eigen::MatrixXd a(nf, block);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    for (Index start = 0; start < nk_box; start += block) {
      const Index nb = std::min(block, nk_box - start);
      for (Index q = 0; q < nb; ++q) {
        box.CalcAOValues(points[kb.points[std::size_t(start + q)]], phi_buf,
                         dphi_buf, q);
      }
      const auto phi = phi_buf.leftCols(nb);
      const auto dphi = dphi_buf.leftCols(3 * nb);
      for (Index k = 0; k < n_resp; ++k) {
        yk.leftCols(nb).noalias() = dmat_k[std::size_t(k)] * phi;
        for (Index q = 0; q < nb; ++q) {
          const Index kq = start + q;
          const double w = kb.w(kq);
          const double drho = 0.5 * phi.col(q).dot(yk.col(q));
          if (is_gga_) {
            const Eigen::Vector3d dg =
                dphi.middleCols(3 * q, 3).transpose() * yk.col(q);
            const Eigen::Vector3d g = kb.grad.col(kq);
            const double dsigma = 2.0 * g.dot(dg);
            const double dv_rho =
                kb.v2rho2(kq) * drho + kb.v2rhosigma(kq) * dsigma;
            const double dv_sigma =
                kb.v2rhosigma(kq) * drho + kb.v2sigma2(kq) * dsigma;
            const Eigen::Vector3d vec =
                2.0 * (dv_sigma * g + kb.vsigma(kq) * dg);
            a.col(q) = (w * 0.5 * dv_rho) * phi.col(q) +
                       w * (dphi.middleCols(3 * q, 3) * vec);
          } else {
            a.col(q) = (w * 0.5 * kb.v2rho2(kq) * drho) * phi.col(q);
          }
        }
        fxc_here[std::size_t(k)].noalias() += a.leftCols(nb) * phi.transpose();
      }
    }
    std::vector<const Eigen::MatrixXd*> blocks;
    for (const Eigen::MatrixXd& m : fxc_here) {
      blocks.push_back(&m);
    }
    acc.Add(box, blocks);
  }
  acc.Finish();
  for (Eigen::MatrixXd& m : result) {
    m = (m + m.transpose()).eval();
  }
  return result;
}

template <class Grid>
std::vector<Eigen::MatrixXd> Vxc_Potential<Grid>::IntegrateFXC(
    const Eigen::MatrixXd& density_matrix,
    const std::vector<Eigen::MatrixXd>& response_matrices) const {
  if (std::any_of(response_matrices.begin(), response_matrices.end(),
                  [&density_matrix](const Eigen::MatrixXd& m) {
                    return m.rows() != density_matrix.rows() ||
                           m.cols() != density_matrix.cols();
                  })) {
    throw std::runtime_error(
        "IntegrateFXC: response matrix does not match the density matrix");
  }
  if (response_matrices.empty()) {
    return {};
  }
  return ApplyFXC(PrepareFXC(density_matrix), response_matrices);
}

// ===========================================================================
// XC gradient at fixed grid weights: basis-function (Pulay) and grid-point
// translation terms, LDA and GGA. The grid-weight (SSW partition) term is
// in GridWeightGradient.
//
// E_xc = sum_p weight_p * e_xc(rho_p, sigma_p)  (e_xc = rho*f_xc)
//
// dE_xc/dR_A = sum_p weight_p [v_rho d(rho_p)/dR_A + v_sigma d(sigma_p)/dR_A]
//
// with v_rho = xc.df_drho, v_sigma = xc.df_dsigma (LibXC vrho/vsigma, the
// full derivatives of e_xc).
//
// Basis term: for chi_mu(r) = phi_mu(r - R_A), d(chi_mu)/dR_A =
// -grad_r chi_mu for mu on A and zero otherwise, so with
// temp = ao.values^T * DMAT_here (DMAT_here = 2P)
//   d(rho_p)/dR_A|_basis = -sum_{mu in A} temp(mu) (grad_r chi_mu)(r_p).
// ===========================================================================
template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::PulayGradient(
    const Eigen::MatrixXd& density_matrix, const AOBasis& dftbasis) const {
  assert(density_matrix.isApprox(density_matrix.transpose()) &&
         "Density matrix has to be symmetric!");

  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(nthreads,
                                           Eigen::MatrixXd::Zero(natoms, 3));

  std::exception_ptr eptr_pulay = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    try {
      Index thread_id = OPENMP::getThreadId();
      const GridBox& box = grid_[i];
      if (!box.Matrixsize()) {
        continue;
      }

      const Eigen::MatrixXd DMAT_here =
          2 * box.ReadFromBigMatrix(density_matrix);

      double cutoff = 1.e-40 / double(density_matrix.rows()) /
                      double(density_matrix.rows());
      if (DMAT_here.cwiseAbs2().maxCoeff() < cutoff) {
        continue;
      }

      // box-local AO index -> atom index
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
      // atom each grid point moves with, for the translation terms
      const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

      for (Index p = 0; p < box.size(); ++p) {
        AOShell::AOValuesHessian ao = box.CalcAOValuesHessian(points[p]);

        Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
        double rho = 0.5 * temp.dot(ao.values);
        const double weight = weights[p];

        if (rho * weight < 1.e-20) {
          continue;
        }

        const Eigen::Vector3d rho_grad = temp.transpose() * ao.derivatives;
        typename Vxc_Potential<Grid>::XC_entry xc =
            EvaluateXC(rho, rho_grad.squaredNorm());

        // ===========================================================
        // GGA sigma terms (need the AO Hessians, EvalAOspaceHessian).
        // They vanish for LDA, where xc.df_dsigma stays 0.
        //
        // sigma_p = |grad_r(rho_p)|^2 = rho_grad . rho_grad
        //
        // Gmat(mu,l) = sum_nu DMAT_mu,nu grad_nu,l  (DMAT_here*ao.derivatives)
        // s_vec(mu) = sum_l Gmat(mu,l) rho_grad_l  (Gmat*rho_grad)
        //
        // Basis term (mu in A):
        //   dsigma_p/dR_A|_basis = -2 sum_{mu in A} [grad_mu s_vec(mu)
        //                          + temp_mu (Hessian_mu . rho_grad)]
        //
        // Translation term (A == owner(p) only), with the density Hessian:
        //   Hessian_rho = sum_mu temp_mu*Hessian_mu
        //                 + symmetrized(ao.derivatives^T * Gmat)
        //   dsigma_p/dR_owner|_translation = 2 * (Hessian_rho * rho_grad)
        //
        // Both are contracted with weight_p * v_sigma.
        Eigen::MatrixX3d Gmat = DMAT_here * ao.derivatives;  // (Matrixsize x 3)
        Eigen::VectorXd s_vec = Gmat * rho_grad;             // (Matrixsize)

        Eigen::Matrix3d Hessian_rho = Eigen::Matrix3d::Zero();
        for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
          Hessian_rho += temp(mu) * ao.hessians[mu];
        }
        Eigen::Matrix3d M = ao.derivatives.transpose() * Gmat;
        Hessian_rho += 0.5 * (M + M.transpose());

        Index owner_of_point_sigma = owner_atoms[p];
        Eigen::Vector3d dsigma_translation = 2.0 * (Hessian_rho * rho_grad);
        grad_thread[thread_id].row(owner_of_point_sigma) +=
            (weight * xc.df_dsigma * dsigma_translation).transpose();
        // ===========================================================

        // ===========================================================
        // Grid-point translation term. Points move rigidly with their
        // owner atom, r_p = R_owner + const, so
        //   d(rho_p)/dR_A|_translation = grad_r(rho_p) if A == owner(p),
        //                                0 otherwise,
        //   dE_xc/dR_A|_translation = sum_{p: owner(p)==A}
        //                             weight_p v_rho(rho_p) grad_r(rho_p).
        Index owner_of_point = owner_atoms[p];
        grad_thread[thread_id].row(owner_of_point) +=
            (weight * xc.df_drho * rho_grad).transpose();
        // ===========================================================

        // Basis terms, per local AO index mu, added to the atom of mu (the
        // shells of a box are not necessarily grouped by atom).
        for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
          Index atom = local_idx_to_atom[mu];
          Eigen::Vector3d contribution = -weight * xc.df_drho * temp(mu) *
                                         ao.derivatives.row(mu).transpose();
          // GGA sigma basis term, see above
          Eigen::Vector3d dsigma_basis_mu =
              -2.0 * (ao.derivatives.row(mu).transpose() * s_vec(mu) +
                      temp(mu) * (ao.hessians[mu] * rho_grad));
          contribution += weight * xc.df_dsigma * dsigma_basis_mu;
          grad_thread[thread_id].row(atom) += contribution.transpose();
        }
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr_pulay) {
          eptr_pulay = std::current_exception();
        }
      }
    }
  }
  if (eptr_pulay) {
    std::rethrow_exception(eptr_pulay);
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
  }
  return grad;
}

// Blocked like IntegrateVXC, per spin channel. UKS density matrices are not
// pre-doubled: rho_s = phi^T D_s phi, grad rho_s = 2 dPhi^T D_s phi.
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

    std::vector<Index> kept;
    std::vector<double> rho, sigma, exc, vrho, vsigma;
    std::vector<Eigen::Vector3d> grad_a, grad_b;

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

      const Index nf = box.Matrixsize();
      const Index block = std::min(kPointBlock, box.size());
      Eigen::MatrixXd Vxc_a_here = Eigen::MatrixXd::Zero(nf, nf);
      Eigen::MatrixXd Vxc_b_here = Eigen::MatrixXd::Zero(nf, nf);
      Eigen::MatrixXd phi(nf, block);
      Eigen::MatrixXd dphi(nf, 3 * block);
      Eigen::MatrixXd ya(nf, block);
      Eigen::MatrixXd yb(nf, block);
      Eigen::MatrixXd aa(nf, block);
      Eigen::MatrixXd ab(nf, block);

      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();

      for (Index start = 0; start < box.size(); start += block) {
        const Index nb = std::min(block, box.size() - start);
        for (Index p = 0; p < nb; ++p) {
          box.CalcAOValues(points[start + p], phi, dphi, p);
        }
        ya.leftCols(nb).noalias() = DMa * phi.leftCols(nb);
        yb.leftCols(nb).noalias() = DMb * phi.leftCols(nb);

        kept.clear();
        rho.clear();
        sigma.clear();
        grad_a.clear();
        grad_b.clear();
        for (Index p = 0; p < nb; ++p) {
          const double rho_a = phi.col(p).dot(ya.col(p));
          const double rho_b = phi.col(p).dot(yb.col(p));
          if ((rho_a + rho_b) * weights[start + p] < 1.e-20) {
            continue;
          }
          const Eigen::Vector3d ga =
              2.0 * (dphi.middleCols(3 * p, 3).transpose() * ya.col(p));
          const Eigen::Vector3d gb =
              2.0 * (dphi.middleCols(3 * p, 3).transpose() * yb.col(p));
          kept.push_back(p);
          rho.push_back(rho_a);
          rho.push_back(rho_b);
          sigma.push_back(ga.dot(ga));
          sigma.push_back(ga.dot(gb));
          sigma.push_back(gb.dot(gb));
          grad_a.push_back(ga);
          grad_b.push_back(gb);
        }
        const Index nk = Index(kept.size());
        if (nk == 0) {
          continue;
        }
        exc.resize(nk);
        vrho.resize(2 * nk);
        vsigma.resize(3 * nk);
        EvaluateXCSpinBatch(nk, rho.data(), sigma.data(), exc.data(),
                            vrho.data(), vsigma.data());

        for (Index k = 0; k < nk; ++k) {
          const Index p = kept[k];
          const double w = weights[start + p];
          exc_private += w * (rho[2 * k] + rho[2 * k + 1]) * exc[k];
          // 0.5 because the result is symmetrized by adding the transpose
          aa.col(k) = (w * 0.5 * vrho[2 * k]) * phi.col(p);
          ab.col(k) = (w * 0.5 * vrho[2 * k + 1]) * phi.col(p);
          if (is_gga_) {
            const auto d = dphi.middleCols(3 * p, 3);
            const Eigen::VectorXd g_a = d * grad_a[k];
            const Eigen::VectorXd g_b = d * grad_b[k];
            const double vs_aa = vsigma[3 * k];
            const double vs_ab = vsigma[3 * k + 1];
            const double vs_bb = vsigma[3 * k + 2];
            aa.col(k).noalias() += w * (2.0 * vs_aa * g_a + vs_ab * g_b);
            ab.col(k).noalias() += w * (vs_ab * g_a + 2.0 * vs_bb * g_b);
          }
          if (k != p) {
            phi.col(k) = phi.col(p);
          }
        }
        Vxc_a_here.noalias() += aa.leftCols(nk) * phi.leftCols(nk).transpose();
        Vxc_b_here.noalias() += ab.leftCols(nk) * phi.leftCols(nk).transpose();
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
// SSW switching function and its derivative, restated from the private
// Vxc_Grid::erf1c (0.5*erfc(|x/(1-x^2)|*alpha), alpha = 1/0.30); must be kept
// consistent with it.
constexpr double kSSWAlpha = 1.0 / 0.30;
constexpr double kSSWCutoff = 0.725;

constexpr double kSqrtPi = 1.7724538509055160273;  // M_PI is not standard

double SSWValue(double mu) {
  double val = 0.5 * std::erfc(std::abs(mu / (1.0 - mu * mu)) * kSSWAlpha);
  if (mu > 0.0) {
    val = 1.0 - val;
  }
  return val;
}

// d(SSWValue)/d(mu), closed form.
double SSWDerivative(double mu) {
  double h = std::abs(mu) / (1.0 - mu * mu);
  double sign_mu = (mu > 0.0) ? 1.0 : ((mu < 0.0) ? -1.0 : 0.0);
  double one_minus_mu2 = 1.0 - mu * mu;
  double hprime = sign_mu * (1.0 + mu * mu) / (one_minus_mu2 * one_minus_mu2);
  double d_erf1c = -(kSSWAlpha / kSqrtPi) *
                   std::exp(-(kSSWAlpha * h) * (kSSWAlpha * h)) * hprime;
  return (mu > 0.0) ? -d_erf1c : d_erf1c;
}

// Adds e * C_p * d(w_owner)/dR_A for all atoms A to grad (natoms x 3), with
// the stored weight = C_p * w_owner, C_p the geometry-independent radial and
// angular quadrature weight and w_owner the SSW partition weight of the
// point's owner atom; e is the energy density at the point. The point moves
// with its owner. O(Natoms^2) per point: d(mu_ab)/dR_A is non-zero only for
// A in {owner, a, b}.
void AddSSWWeightDerivative(const Eigen::Vector3d& point, Index owner,
                            double weight, double e, const QMMolecule& atoms,
                            const Eigen::MatrixXd& Rij, Eigen::MatrixXd& grad) {
  const Index natoms = atoms.size();
  Eigen::VectorXd rq(natoms);
  for (Index k = 0; k < natoms; ++k) {
    rq(k) = (point - atoms[k].getPos()).norm();
  }
  // p[], mu_table, sk_table and hard-cutoff flags as in
  // Vxc_Grid::SSWpartition, keeping the per-pair values for the derivative.
  // hard(j,i): 0 = smooth, 1 = mu > cutoff (p[i] = 0), -1 = mu < -cutoff
  // (p[j] = 0); upper triangle (j < i) only.
  Eigen::VectorXd p = Eigen::VectorXd::Ones(natoms);
  Eigen::MatrixXd mu_table = Eigen::MatrixXd::Zero(natoms, natoms);
  Eigen::MatrixXd sk_table = Eigen::MatrixXd::Zero(natoms, natoms);
  Eigen::MatrixXi hard = Eigen::MatrixXi::Zero(natoms, natoms);
  for (Index ii = 1; ii < natoms; ++ii) {
    for (Index jj = 0; jj < ii; ++jj) {
      const double mu = (rq(ii) - rq(jj)) * Rij(jj, ii);
      mu_table(jj, ii) = mu;
      if (mu > kSSWCutoff) {
        p(ii) = 0.0;
        hard(jj, ii) = 1;
      } else if (mu < -kSSWCutoff) {
        p(jj) = 0.0;
        hard(jj, ii) = -1;
      } else {
        const double sk = SSWValue(mu);
        sk_table(jj, ii) = sk;
        p(jj) *= sk;
        p(ii) *= (1.0 - sk);
      }
    }
  }
  const double wsum = p.sum();
  const double w_owner = p(owner) / wsum;
  // points with negligible w_owner have negligible weight (C_p is bounded)
  constexpr double kNegligibleWOwner = 1.e-8;
  if (w_owner < kNegligibleWOwner) {
    return;
  }
  // d(rq(k))/dR_A: +-unit vector for A == owner or A == k, zero if both
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
  // coef * d(mu_ab)/dR_A (a < b) added to out for A in {owner, a, b}
  auto add_dmu = [&](Index a, Index b, double coef, Eigen::MatrixXd& out) {
    const double Rab = 1.0 / Rij(a, b);
    const Eigen::Vector3d rvec = atoms[a].getPos() - atoms[b].getPos();
    const double mu = mu_table(a, b);
    const Index cand[3] = {owner, a, b};
    for (Index c = 0; c < 3; ++c) {
      const Index A = cand[c];
      if ((c > 0 && A == cand[0]) || (c > 1 && A == cand[1])) {
        continue;
      }
      Eigen::Vector3d dRab = Eigen::Vector3d::Zero();
      if (A == a) {
        dRab = rvec / Rab;
      } else if (A == b) {
        dRab = -rvec / Rab;
      }
      const Eigen::Vector3d dmu =
          (d_rq_dR(b, A) - d_rq_dR(a, A)) / Rab - (mu / Rab) * dRab;
      out.row(A) += coef * dmu.transpose();
    }
  };
  // d p_k/dR = p_k sum over the smooth pairs of k of the log-derivatives of
  // s or (1-s); these divide by s or 1-s, so pairs (and p_k) near the SSW
  // saturation are treated as zero (threshold not tuned)
  constexpr double kNegligibleP = 1.e-8;
  Eigen::MatrixXd dp_owner = Eigen::MatrixXd::Zero(natoms, 3);
  Eigen::MatrixXd dp_sum = Eigen::MatrixXd::Zero(natoms, 3);
  Eigen::MatrixXd dp_k(natoms, 3);
  for (Index k = 0; k < natoms; ++k) {
    if (p(k) < kNegligibleP) {
      continue;
    }
    dp_k.setZero();
    for (Index b = k + 1; b < natoms; ++b) {
      if (hard(k, b) != 0) {
        continue;
      }
      const double skv = sk_table(k, b);
      if (skv < kNegligibleP) {
        continue;
      }
      add_dmu(k, b, SSWDerivative(mu_table(k, b)) / skv, dp_k);
    }
    for (Index a = 0; a < k; ++a) {
      if (hard(a, k) != 0) {
        continue;
      }
      const double one_minus_skv = 1.0 - sk_table(a, k);
      if (one_minus_skv < kNegligibleP) {
        continue;
      }
      add_dmu(a, k, -SSWDerivative(mu_table(a, k)) / one_minus_skv, dp_k);
    }
    dp_sum += p(k) * dp_k;
    if (k == owner) {
      dp_owner = p(k) * dp_k;
    }
  }
  const double C_p = weight / w_owner;
  grad += (C_p * e / wsum) * (dp_owner - w_owner * dp_sum);
}
}  // namespace

template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::GridWeightGradient(
    const Eigen::MatrixXd& density_matrix, const QMMolecule& atoms) const {
  Index natoms = atoms.size();
  Eigen::MatrixXd Rij = grid_.CalcInverseAtomDist(atoms);

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(nthreads,
                                           Eigen::MatrixXd::Zero(natoms, 3));

  std::exception_ptr eptr_weight = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    try {
      Index thread_id = OPENMP::getThreadId();
      const GridBox& box = grid_[i];
      if (!box.Matrixsize()) {
        continue;
      }

      const Eigen::MatrixXd DMAT_here =
          2 * box.ReadFromBigMatrix(density_matrix);
      double cutoff = 1.e-40 / double(density_matrix.rows()) /
                      double(density_matrix.rows());
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
          // Owner atoms are set by GridSetup; skipping the point would give
          // an incomplete gradient.
          throw std::runtime_error(
              "GridWeightGradient: grid point has no owner atom set; the "
              "grid must be built with Vxc_Grid::GridSetup.");
        }

        AddSSWWeightDerivative(points[pidx], owner, weight, rho * xc.f_xc,
                               atoms, Rij, grad_thread[thread_id]);
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr_weight) {
          eptr_weight = std::current_exception();
        }
      }
    }
  }
  if (eptr_weight) {
    std::rethrow_exception(eptr_weight);
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
  }
  return grad;
}

// ===========================================================================
// d/dt [PulayGradient + GridWeightGradient](D + t P) at t = 0.
//
// With T = 2D and T' = 2P on the box, per grid point (phi, dphi, H_mu the AO
// values, gradients and Hessians):
//   temp = phi^T T, g = dphi^T temp, Gm = T dphi, s = Gm g,
//   Hrho = sum_mu temp_mu H_mu + sym(dphi^T Gm),
// and the primed objects from T'. With drho = temp'.phi / 2,
// dsigma = 2 g.g', dvr = v2rho2 drho + v2rhosigma dsigma,
// dvs = v2rhosigma drho + v2sigma2 dsigma, s' = Gm' g + Gm g':
//   translation (owner):
//     w [dvr g + vr g' + 2 dvs Hrho g + 2 vs (Hrho' g + Hrho g')]
//   basis (mu on A):
//     -w (dvr temp_mu + vr temp'_mu) dphi_mu
//     -2w [dvs (dphi_mu s_mu + temp_mu H_mu g)
//          + vs (dphi_mu s'_mu + temp'_mu H_mu g + temp_mu H_mu g')]
//   grid weights: the SSW derivative with energy density
//     d(rho f_xc) = vr drho + vs dsigma.
// ===========================================================================
template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::GradientResponse(
    const Eigen::MatrixXd& density_matrix, const Eigen::MatrixXd& direction,
    const QMMolecule& atoms) const {
  const Index natoms = atoms.size();
  const Eigen::MatrixXd Rij = grid_.CalcInverseAtomDist(atoms);
  // fail before the grid loop for functionals without fxc
  {
    const double rho = 0.1;
    const double sigma = 0.01;
    double vr, vs, rr, rs, ss;
    EvaluateVXCFXCBatch(1, &rho, &sigma, &vr, &vs, &rr, &rs, &ss);
  }
  const Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(nthreads,
                                           Eigen::MatrixXd::Zero(natoms, 3));
  std::exception_ptr eptr = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    try {
      const Index thread_id = OPENMP::getThreadId();
      Eigen::MatrixXd& gt = grad_thread[thread_id];
      const GridBox& box = grid_[i];
      if (!box.Matrixsize()) {
        continue;
      }
      const Eigen::MatrixXd T = 2 * box.ReadFromBigMatrix(density_matrix);
      const Eigen::MatrixXd Tp = 2 * box.ReadFromBigMatrix(direction);
      const double cutoff = 1.e-40 / double(density_matrix.rows()) /
                            double(density_matrix.rows());
      if (T.cwiseAbs2().maxCoeff() < cutoff) {
        continue;
      }
      const Index nf = box.Matrixsize();
      std::vector<Index> local_idx_to_atom(nf);
      const std::vector<const AOShell*>& shells = box.getShells();
      const std::vector<GridboxRange>& ao_ranges = box.getAOranges();
      for (std::size_t sh = 0; sh < shells.size(); ++sh) {
        const Index atom = shells[sh]->getAtomIndex();
        for (Index k = 0; k < ao_ranges[sh].size; ++k) {
          local_idx_to_atom[ao_ranges[sh].start + k] = atom;
        }
      }
      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();
      const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

      for (Index pt = 0; pt < box.size(); ++pt) {
        const double w = weights[pt];
        const Index owner = owner_atoms[pt];
        if (owner < 0) {
          throw std::runtime_error(
              "GradientResponse: grid point has no owner atom set; the grid "
              "must be built with Vxc_Grid::GridSetup.");
        }
        AOShell::AOValuesHessian ao = is_gga_
                                          ? box.CalcAOValuesHessian(points[pt])
                                          : AOShell::AOValuesHessian(0);
        AOShell::AOValues aov(0);
        if (!is_gga_) {
          aov = box.CalcAOValues(points[pt]);
        }
        const Eigen::VectorXd& phi = is_gga_ ? ao.values : aov.values;
        const Eigen::MatrixX3d& dphi =
            is_gga_ ? ao.derivatives : aov.derivatives;

        const Eigen::VectorXd temp = T * phi;
        const double rho = 0.5 * temp.dot(phi);
        if (rho * w < 1.e-20) {
          continue;
        }
        const Eigen::VectorXd tempP = Tp * phi;
        const double drho = 0.5 * tempP.dot(phi);
        const Eigen::Vector3d g = dphi.transpose() * temp;
        const Eigen::Vector3d gp = dphi.transpose() * tempP;
        const double sigma = g.squaredNorm();
        const double dsigma = 2.0 * g.dot(gp);
        double vr, vs, rr, rs, ss;
        EvaluateVXCFXCBatch(1, &rho, &sigma, &vr, &vs, &rr, &rs, &ss);
        const double dvr = rr * drho + rs * dsigma;
        const double dvs = rs * drho + ss * dsigma;

        // LDA part: translation and basis terms
        gt.row(owner) += (w * (dvr * g + vr * gp)).transpose();
        for (Index mu = 0; mu < nf; ++mu) {
          gt.row(local_idx_to_atom[mu]) -=
              w * (dvr * temp(mu) + vr * tempP(mu)) * dphi.row(mu);
        }
        if (is_gga_) {
          const Eigen::MatrixX3d Gm = T * dphi;
          const Eigen::MatrixX3d Gp = Tp * dphi;
          const Eigen::VectorXd sv = Gm * g;
          const Eigen::VectorXd svp = Gp * g + Gm * gp;
          Eigen::Matrix3d Hrho = Eigen::Matrix3d::Zero();
          Eigen::Matrix3d Hp = Eigen::Matrix3d::Zero();
          for (Index mu = 0; mu < nf; ++mu) {
            Hrho += temp(mu) * ao.hessians[mu];
            Hp += tempP(mu) * ao.hessians[mu];
          }
          const Eigen::Matrix3d M = dphi.transpose() * Gm;
          const Eigen::Matrix3d Mp = dphi.transpose() * Gp;
          Hrho += 0.5 * (M + M.transpose());
          Hp += 0.5 * (Mp + Mp.transpose());
          gt.row(owner) +=
              (2.0 * w * (dvs * (Hrho * g) + vs * (Hp * g + Hrho * gp)))
                  .transpose();
          for (Index mu = 0; mu < nf; ++mu) {
            const Eigen::Matrix3d& H = ao.hessians[mu];
            const Eigen::Vector3d d = dphi.row(mu).transpose();
            const Eigen::Vector3d c =
                dvs * (d * sv(mu) + temp(mu) * (H * g)) +
                vs * (d * svp(mu) + tempP(mu) * (H * g) + temp(mu) * (H * gp));
            gt.row(local_idx_to_atom[mu]) -= (2.0 * w * c).transpose();
          }
        }
        AddSSWWeightDerivative(points[pt], owner, w, vr * drho + vs * dsigma,
                               atoms, Rij, gt);
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr) {
          eptr = std::current_exception();
        }
      }
    }
  }
  if (eptr) {
    std::rethrow_exception(eptr);
  }
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (const Eigen::MatrixXd& gtt : grad_thread) {
    grad += gtt;
  }
  return grad;
}

template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::PulayGradientUKS(
    const Eigen::MatrixXd& dmat_alpha, const Eigen::MatrixXd& dmat_beta,
    const AOBasis& dftbasis) const {
  // UKS version of PulayGradient (LDA and GGA). The three sigma terms
  // (sigma_aa, sigma_ab, sigma_bb) reduce to the restricted per-spin
  // expressions (Gmat_s, Hessian_rho_s), applied once per spin with the
  // effective gradient vectors
  //   V_a = 2*vsigma_aa*rho_a_grad + vsigma_ab*rho_b_grad
  //   V_b = 2*vsigma_bb*rho_b_grad + vsigma_ab*rho_a_grad.
  Index natoms = static_cast<Index>(dftbasis.getFuncPerAtom().size());
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(nthreads,
                                           Eigen::MatrixXd::Zero(natoms, 3));

  std::exception_ptr eptr_pulay_uks = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    try {
      Index thread_id = OPENMP::getThreadId();
      const GridBox& box = grid_[i];
      if (!box.Matrixsize()) {
        continue;
      }

      const Eigen::MatrixXd DMa = box.ReadFromBigMatrix(dmat_alpha);
      const Eigen::MatrixXd DMb = box.ReadFromBigMatrix(dmat_beta);

      // Box-local AO index -> atom index -- MUST use getAOranges() (start/
      // size per shell), not a naive sequential push_back per shell: a
      // single box's significant shells are not guaranteed to be grouped
      // contiguously by atom. Matches the exact construction already
      // validated in the restricted PulayGradient.
      std::vector<Index> local_idx_to_atom(box.Matrixsize());
      const std::vector<const AOShell*>& shells = box.getShells();
      const std::vector<GridboxRange>& ao_ranges = box.getAOranges();
      for (size_t s = 0; s < shells.size(); ++s) {
        Index atom = shells[s]->getAtomIndex();
        for (Index k = 0; k < ao_ranges[s].size; ++k) {
          local_idx_to_atom[ao_ranges[s].start + k] = atom;
        }
      }
      const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();

      for (Index p = 0; p < box.size(); ++p) {
        AOShell::AOValuesHessian ao = box.CalcAOValuesHessian(points[p]);

        Eigen::VectorXd temp_a = DMa * ao.values;
        Eigen::VectorXd temp_b = DMb * ao.values;
        const double rho_a = ao.values.dot(temp_a);
        const double rho_b = ao.values.dot(temp_b);
        const double rho = rho_a + rho_b;
        const double weight = weights[p];

        if (rho * weight < 1.e-20) {
          continue;
        }

        const Eigen::Vector3d rho_a_grad =
            2.0 * (ao.derivatives.transpose() * temp_a);
        const Eigen::Vector3d rho_b_grad =
            2.0 * (ao.derivatives.transpose() * temp_b);
        const double sigma_aa = rho_a_grad.dot(rho_a_grad);
        const double sigma_ab = rho_a_grad.dot(rho_b_grad);
        const double sigma_bb = rho_b_grad.dot(rho_b_grad);

        typename Vxc_Potential<Grid>::XC_entry_spin xc =
            EvaluateXCSpin(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb);

        // translation term (A == owner(p) only), rho part
        Index owner_of_point = owner_atoms[p];
        grad_thread[thread_id].row(owner_of_point) +=
            (weight * (xc.vrho_a * rho_a_grad + xc.vrho_b * rho_b_grad))
                .transpose();

        // GGA sigma terms: basis and translation terms as in the restricted
        // case, with Gmat_s/Hessian_rho_s per spin (from the spin density
        // matrices instead of DMAT_here), contracted with V_a/V_b.
        Eigen::MatrixX3d Gmat_a = 2.0 * (DMa * ao.derivatives);
        Eigen::MatrixX3d Gmat_b = 2.0 * (DMb * ao.derivatives);
        Eigen::Vector3d V_a =
            2.0 * xc.vsigma_aa * rho_a_grad + xc.vsigma_ab * rho_b_grad;
        Eigen::Vector3d V_b =
            2.0 * xc.vsigma_bb * rho_b_grad + xc.vsigma_ab * rho_a_grad;

        // Translation-type sigma contribution: Hessian_rho_s * V_s,
        // summed over both spins.
        Eigen::Matrix3d Hessian_rho_a = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d Hessian_rho_b = Eigen::Matrix3d::Zero();
        for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
          Hessian_rho_a += 2.0 * temp_a(mu) * ao.hessians[mu];
          Hessian_rho_b += 2.0 * temp_b(mu) * ao.hessians[mu];
        }
        Eigen::Matrix3d Ma = ao.derivatives.transpose() * Gmat_a;
        Eigen::Matrix3d Mb = ao.derivatives.transpose() * Gmat_b;
        Hessian_rho_a += 0.5 * (Ma + Ma.transpose());
        Hessian_rho_b += 0.5 * (Mb + Mb.transpose());

        grad_thread[thread_id].row(owner_of_point) +=
            (weight * (Hessian_rho_a * V_a + Hessian_rho_b * V_b)).transpose();

        // Basis-type term (LDA + sigma), accumulated per local AO index
        // mu, scattered to the correct atom via local_idx_to_atom.
        Eigen::VectorXd Gmat_a_dot_Va = Gmat_a * V_a;
        Eigen::VectorXd Gmat_b_dot_Vb = Gmat_b * V_b;
        for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
          Index atom = local_idx_to_atom[mu];
          Eigen::Vector3d lda_contribution =
              -2.0 * weight *
              (xc.vrho_a * temp_a(mu) + xc.vrho_b * temp_b(mu)) *
              ao.derivatives.row(mu).transpose();
          Eigen::Vector3d sigma_contribution =
              -weight *
              (ao.derivatives.row(mu).transpose() * Gmat_a_dot_Va(mu) +
               2.0 * temp_a(mu) * (ao.hessians[mu] * V_a) +
               ao.derivatives.row(mu).transpose() * Gmat_b_dot_Vb(mu) +
               2.0 * temp_b(mu) * (ao.hessians[mu] * V_b));
          grad_thread[thread_id].row(atom) +=
              (lda_contribution + sigma_contribution).transpose();
        }
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr_pulay_uks) {
          eptr_pulay_uks = std::current_exception();
        }
      }
    }
  }
  if (eptr_pulay_uks) {
    std::rethrow_exception(eptr_pulay_uks);
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
  }
  return grad;
}

template <class Grid>
Eigen::MatrixXd Vxc_Potential<Grid>::GridWeightGradientUKS(
    const Eigen::MatrixXd& dmat_alpha, const Eigen::MatrixXd& dmat_beta,
    const QMMolecule& atoms) const {
  // UKS version of GridWeightGradient (LDA and GGA). The SSW weight
  // derivative (AddSSWWeightDerivative) is purely geometric and identical; only
  // the prefactor differs: rho_a, rho_b from the spin density matrices (not
  // pre-doubled, as in IntegrateVXCSpin) and f_xc from EvaluateXCSpin.
  Index natoms = atoms.size();
  Eigen::MatrixXd Rij = grid_.CalcInverseAtomDist(atoms);

  Index nthreads = OPENMP::getMaxThreads();
  std::vector<Eigen::MatrixXd> grad_thread(nthreads,
                                           Eigen::MatrixXd::Zero(natoms, 3));

  std::exception_ptr eptr_weight_uks = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    try {
      Index thread_id = OPENMP::getThreadId();
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

      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();
      const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

      for (Index pidx = 0; pidx < box.size(); ++pidx) {
        AOShell::AOValues ao = box.CalcAOValues(points[pidx]);
        Eigen::VectorXd temp_a = DMa * ao.values;
        Eigen::VectorXd temp_b = DMb * ao.values;
        const double rho_a = ao.values.dot(temp_a);
        const double rho_b = ao.values.dot(temp_b);
        const double rho = rho_a + rho_b;
        double weight = weights[pidx];
        if (rho * weight < 1.e-20) {
          continue;
        }
        const Eigen::Vector3d rho_a_grad =
            2.0 * (ao.derivatives.transpose() * temp_a);
        const Eigen::Vector3d rho_b_grad =
            2.0 * (ao.derivatives.transpose() * temp_b);
        const double sigma_aa = rho_a_grad.dot(rho_a_grad);
        const double sigma_ab = rho_a_grad.dot(rho_b_grad);
        const double sigma_bb = rho_b_grad.dot(rho_b_grad);
        typename Vxc_Potential<Grid>::XC_entry_spin xc =
            EvaluateXCSpin(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb);

        Index owner = owner_atoms[pidx];
        if (owner < 0) {
          throw std::runtime_error(
              "GridWeightGradientUKS: grid point has no owner atom set; the "
              "grid must be built with Vxc_Grid::GridSetup.");
        }

        AddSSWWeightDerivative(points[pidx], owner, weight, rho * xc.f_xc,
                               atoms, Rij, grad_thread[thread_id]);
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr_weight_uks) {
          eptr_weight_uks = std::current_exception();
        }
      }
    }
  }
  if (eptr_weight_uks) {
    std::rethrow_exception(eptr_weight_uks);
  }

  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);
  for (Index t = 0; t < nthreads; ++t) {
    grad += grad_thread[t];
  }
  return grad;
}

template class Vxc_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca