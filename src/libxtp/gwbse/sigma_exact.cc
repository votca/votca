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
  _residues = CalcResidues();
  return;
}

Eigen::VectorXd Sigma_Exact::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  const Index rpasize = _rpa_solution.omega.size();  // TODO: Rename
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
#pragma omp parallel for
  for (Index m = 0; m < _qptotal; m++) {
    double sigmc = 0.0;
    const Eigen::MatrixXd& res_m = _residues[m];
    for (Index s = 0; s < rpasize; s++) {
      const Eigen::VectorXd res_mm = res_m.col(s).cwiseAbs2();
      double eigenvalue = _rpa_solution.omega(s);
      sigmc += CalcSigmaC(res_mm, eigenvalue, frequencies(m));
    }
    // Multiply with factor 2.0 to sum over both (identical) spin states
    result(m) = 2.0 * sigmc;
  }
  return result;
}

Eigen::MatrixXd Sigma_Exact::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  const Index rpasize = _rpa_solution.omega.size();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
#pragma omp parallel for
  for (Index m = 0; m < _qptotal; m++) {
    const Eigen::MatrixXd& res_m = _residues[m];
    for (Index n = m + 1; n < _qptotal; n++) {
      double sigmc = 0.0;
      const Eigen::MatrixXd& res_n = _residues[n];
      for (Index s = 0; s < rpasize; s++) {
        Eigen::VectorXd res_mn = res_m.col(s).cwiseProduct(res_n.col(s));
        double eigenvalue = _rpa_solution.omega(s);
        double sigmc_m = CalcSigmaC(res_mn, eigenvalue, frequencies(m));
        double sigmc_n = CalcSigmaC(res_mn, eigenvalue, frequencies(n));
        sigmc += sigmc_m + sigmc_n;
      }
      // Multiply with factor 2.0 to sum over both (identical) spin states
      result(m, n) = 0.5 * 2.0 * sigmc;
      result(n, m) = 0.5 * 2.0 * sigmc;
    }
  }
  return result;
}

std::vector<Eigen::MatrixXd> Sigma_Exact::CalcResidues() const {
  const Index lumo = _opt.homo + 1;
  const Index n_occ = lumo - _opt.rpamin;
  const Index n_unocc = _opt.rpamax - _opt.homo;
  const Index rpasize = n_occ * n_unocc;
  const Index qpoffset = _opt.qpmin - _opt.rpamin;
  const Index auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);
  std::vector<Eigen::MatrixXd> residues(_qptotal);
#pragma omp parallel for
  for (Index m = 0; m < _qptotal; m++) {
    const Eigen::MatrixXd Mmn_mT = _Mmn[m + qpoffset].transpose();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_rpatotal, rpasize);
    for (Index v = 0; v < n_occ; v++) {  // Sum over v
      const Eigen::MatrixXd fc =
          _Mmn[v].block(n_occ, 0, n_unocc, auxsize) * Mmn_mT;  // Sum over chi
      res += fc.transpose() * _rpa_solution.XpY.block(vc.I(v, 0), 0, n_unocc,
                                                      rpasize);  // Sum over c
    }
    residues[m] = res;
  }
  return residues;
}

double Sigma_Exact::CalcSigmaC(const Eigen::VectorXd& res_mn, double eigenvalue,
                               double freq) const {
  const double eta = _opt.eta;
  const Index lumo = _opt.homo + 1;
  const Index n_occ = lumo - _opt.rpamin;
  const Index n_unocc = _opt.rpamax - _opt.homo;
  Eigen::ArrayXd temp = -_rpa.getRPAInputEnergies().array() + freq;
  temp.segment(0, n_occ) += eigenvalue;
  temp.segment(n_occ, n_unocc) -= eigenvalue;
  const Eigen::ArrayXd numer = res_mn.array() * temp;
  const Eigen::ArrayXd denom = temp.abs2() + eta * eta;
  return (numer / denom).sum();
}

}  // namespace xtp
}  // namespace votca