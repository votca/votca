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

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <sys/stat.h>
#include <votca/ctp/apolarsite.h>
#include <votca/tools/constants.h>
#include <votca/xtp/qmiter.h>

using boost::format;

namespace votca {
namespace xtp {

void QMMIter::UpdateMPSFromGDMA(std::vector<std::vector<double> > &multipoles,
                                std::vector<ctp::PolarSeg *> &psegs) {

  for (unsigned int i = 0, qac = 0; i < psegs.size(); ++i) {
    ctp::PolarSeg *pseg = psegs[i];
    for (unsigned int j = 0; j < pseg->size(); ++j, ++qac) {
      // Retrieve multipole info of this atom
      std::vector<double> update = multipoles[qac];
      while (update.size() < 9) update.push_back(0.0);
      // Convert e*(a_0)^k to e*(nm)^k where k = rank
      for (int m = 1; m < 4; m++) {
        update[m] *= pow(tools::conv::bohr2nm, 1);
      }
      for (int m = 4; m < 9; m++) {
        update[m] *= pow(tools::conv::bohr2nm, 2);
      }
      ctp::APolarSite *aps = (*pseg)[j];
      aps->setQs(update, 0);
    }
  }
  return;
}

void QMMIter::UpdatePosChrgFromQMAtoms(std::vector<QMAtom *> &qmatoms,
                                       std::vector<ctp::PolarSeg *> &psegs) {

  double dR_RMS = 0.0;
  double dQ_RMS = 0.0;
  double dQ_SUM = 0.0;

  for (unsigned int i = 0, qac = 0; i < psegs.size(); ++i) {
    ctp::PolarSeg *pseg = psegs[i];
    for (unsigned int j = 0; j < pseg->size(); ++j, ++qac) {
      // Retrieve info from QMAtom
      QMAtom *qmatm = qmatoms[qac];
      tools::vec upd_r = qmatm->getPos() * tools::conv::bohr2nm;
      double upd_Q00 = qmatm->getPartialcharge();
      // Compare to previous r, Q00
      ctp::APolarSite *aps = (*pseg)[j];
      tools::vec old_r = aps->getPos();
      double old_Q00 = aps->getQ00();
      double dR = tools::abs(upd_r - old_r);
      double dQ00 = upd_Q00 - old_Q00;

      dR_RMS += dR * dR;
      dQ_RMS += dQ00 * dQ00;
      dQ_SUM += dQ00;

      // Forward updated r, Q00 to APS
      aps->setPos(upd_r);
      aps->setQ00(upd_Q00, 0);
    }
  }

  dR_RMS /= qmatoms.size();
  dQ_RMS /= qmatoms.size();
  dR_RMS = sqrt(dR_RMS);
  dQ_RMS = sqrt(dQ_RMS);

  this->setdRdQ(dR_RMS, dQ_RMS, dQ_SUM);
  return;
}

void QMMIter::setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM) {
  _hasdRdQ = true;
  _dR_RMS = dR_RMS;
  _dQ_RMS = dQ_RMS;
  _dQ_SUM = dQ_SUM;
  return;
}

void QMMIter::setQMSF(double energy_QM, double energy_SF, double energy_GWBSE) {
  _hasQM = true;
  _e_QM = energy_QM;
  _e_SF = energy_SF;
  _hasGWBSE = true;
  _e_GWBSE = energy_GWBSE;

  return;
}

void QMMIter::setE_FM(double ef00, double ef01, double ef02, double ef11,
                      double ef12, double em0, double em1, double em2,
                      double efm) {

  _hasMM = true;
  _ef_00 = ef00;
  _ef_01 = ef01;
  _ef_02 = ef02;
  _ef_11 = ef11;
  _ef_12 = ef12;
  _em_0_ = em0;
  _em_1_ = em1;
  _em_2_ = em2;
  _e_fm_ = efm;
  return;
}

double QMMIter::getMMEnergy() { return _ef_11 + _ef_12 + _em_1_ + _em_2_; }

double QMMIter::getQMMMEnergy() {
  return _e_QM + +_e_GWBSE + _ef_11 + _ef_12 + _em_1_ + _em_2_;
}

}  // namespace xtp
}  // namespace votca
