/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include "votca/xtp/vc2index.h"
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_exact.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

void Sigma_Exact::PrepareScreening() {
  _rpa_solution = _rpa.Diagonalize_H2p();
  _residues = std::vector<Eigen::MatrixXd>(_qptotal);
#pragma omp parallel for
  for (Index gw_level = 0; gw_level < _qptotal; gw_level++) {
    _residues[gw_level] = CalcResidues(gw_level);
  }
  return;
}

double Sigma_Exact::CalcCorrelationDiagElement(Index gw_level,
                                               double frequency) const {
  const double eta = _opt.eta;
  const Index lumo = _opt.homo + 1;
  const Index n_occ = lumo - _opt.rpamin;
  const Index n_unocc = _opt.rpamax - _opt.homo;
  const Index rpasize = _rpa_solution.omega.size();
  double sigma_c = 0.0;
  for (Index s = 0; s < rpasize; s++) {
    const double eigenvalue = _rpa_solution.omega(s);
    const Eigen::VectorXd res_12 = _residues[gw_level].col(s).cwiseAbs2();
    Eigen::ArrayXd temp = -_rpa.getRPAInputEnergies().array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;
    const Eigen::ArrayXd numer = res_12.array() * temp;
    const Eigen::ArrayXd denom = temp.abs2() + eta * eta;
    sigma_c += (numer / denom).sum();
  }
  // Multiply with factor 2.0 to sum over both (identical) spin states
  return 2.0 * sigma_c;
}

double Sigma_Exact::CalcCorrelationOffDiagElement(Index gw_level1,
                                                  Index gw_level2,
                                                  double frequency1,
                                                  double frequency2) const {
  const double eta = _opt.eta;
  const Index lumo = _opt.homo + 1;
  const Index n_occ = lumo - _opt.rpamin;
  const Index n_unocc = _opt.rpamax - _opt.homo;
  const Index rpasize = _rpa_solution.omega.size();
  double sigma_c = 0.0;
  for (Index s = 0; s < rpasize; s++) {
    const double eigenvalue = _rpa_solution.omega(s);
    const Eigen::VectorXd& res1 = _residues[gw_level1].col(s);
    const Eigen::VectorXd& res2 = _residues[gw_level2].col(s);
    const Eigen::VectorXd res_12 = res1.cwiseProduct(res2);
    Eigen::ArrayXd temp1 = -_rpa.getRPAInputEnergies().array();
    temp1.segment(0, n_occ) += eigenvalue;
    temp1.segment(n_occ, n_unocc) -= eigenvalue;
    const Eigen::ArrayXd temp2 = temp1 + frequency2;
    temp1 += frequency1;
    const Eigen::ArrayXd numer1 = res_12.array() * temp1;
    const Eigen::ArrayXd numer2 = res_12.array() * temp2;
    const Eigen::ArrayXd denom1 = temp1.abs2() + eta * eta;
    const Eigen::ArrayXd denom2 = temp2.abs2() + eta * eta;
    sigma_c += 0.5 * ((numer1 / denom1) + (numer2 / denom2)).sum();
  }
  // Multiply with factor 2.0 to sum over both (identical) spin states
  return 2.0 * sigma_c;
}

Eigen::MatrixXd Sigma_Exact::CalcResidues(Index gw_level) const {
  const Index lumo = _opt.homo + 1;
  const Index n_occ = lumo - _opt.rpamin;
  const Index n_unocc = _opt.rpamax - _opt.homo;
  const Index rpasize = n_occ * n_unocc;
  const Index qpoffset = _opt.qpmin - _opt.rpamin;
  const Index auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);
  const Eigen::MatrixXd Mmn_i_T = _Mmn[gw_level + qpoffset].transpose();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_rpatotal, rpasize);
  for (Index v = 0; v < n_occ; v++) {  // Sum over v
    const Eigen::MatrixXd& Mmn_v = _Mmn[v].block(n_occ, 0, n_unocc, auxsize);
    const Eigen::MatrixXd fc = Mmn_v * Mmn_i_T;  // Sum over chi
    const Eigen::MatrixXd& XpY_v =
        _rpa_solution.XpY.block(vc.I(v, 0), 0, n_unocc, rpasize);
    res += fc.transpose() * XpY_v;  // Sum over c
  }
  return res;
}

}  // namespace xtp
}  // namespace votca
