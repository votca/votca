/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>
#include <votca/xtp/qmiter.h>
#include <sys/stat.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <votca/tools/constants.h>
#include <votca/ctp/apolarsite.h>

using boost::format;

namespace votca {
    namespace xtp {

     
void QMMIter::UpdateMPSFromGDMA(std::vector<std::vector<double> > &multipoles, std::vector< ctp::PolarSeg* > &psegs) {


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

        void QMMIter::UpdatePosChrgFromQMAtoms(std::vector< ctp::QMAtom* > &qmatoms,
                std::vector< ctp::PolarSeg* > &psegs) {

            double AA_to_NM = tools::conv::ang2nm;

            double dR_RMS = 0.0;
            double dQ_RMS = 0.0;
            double dQ_SUM = 0.0;

            for (unsigned int i = 0, qac = 0; i < psegs.size(); ++i) {
                ctp::PolarSeg *pseg = psegs[i];
                for (unsigned int j = 0; j < pseg->size(); ++j, ++qac) {

                    // Retrieve info from ctp::QMAtom
                    ctp::QMAtom *qmatm = qmatoms[qac];
                    vec upd_r = vec(qmatm->x, qmatm->y, qmatm->z);
                    upd_r *= AA_to_NM;
                    double upd_Q00 = qmatm->charge;
                    //cout << "updating charge to " << qmatm->charge << endl;

                    // Compare to previous r, Q00
                    ctp::APolarSite *aps = (*pseg)[j];
                    vec old_r = aps->getPos();
                    double old_Q00 = aps->getQ00();
                    double dR = abs(upd_r - old_r);
                    double dQ00 = upd_Q00 - old_Q00;

                    dR_RMS += dR*dR;
                    dQ_RMS += dQ00*dQ00;
                    dQ_SUM += dQ00;

                    // Forward updated r, Q00 to APS
                    aps->setPos(upd_r);
                    aps->setQ00(upd_Q00, 0);
                }
            }


            // cout << " dR_RMS " << dR_RMS << "dQ RMS " << dQ_RMS << endl;
            dR_RMS /= qmatoms.size();
            dQ_RMS /= qmatoms.size();
            dR_RMS = sqrt(dR_RMS);
            dQ_RMS = sqrt(dQ_RMS);

            // cout << " dR_RMS " << dR_RMS << "dQ RMS " << dQ_RMS << endl;


            this->setdRdQ(dR_RMS, dQ_RMS, dQ_SUM);
            return;
        }

        void QMMIter::GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop, Orbitals &orb,
                bool split_dpl, double dpl_spacing) {

            double AA_to_NM = tools::conv::ang2nm;

            // INNER SHELL QM0
            for (unsigned int i = 0; i < ptop->QM0().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->QM0()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {

                    ctp::APolarSite *aps = (*pseg)[j];
                    vec pos = aps->getPos() / AA_to_NM;
                    double Q = aps->getQ00();
                    string type = "qm";

                    orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, false);

                }
            }

            // MIDDLE SHELL MM1
            for (unsigned int i = 0; i < ptop->MM1().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->MM1()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {

                    ctp::APolarSite *aps = (*pseg)[j];
                    if(aps->getRank()>1){
                        cerr<<"WARNING: Segments contain quadrupoles, votca cannot yet split them into charges. Your result will most likely be wrong."<<endl;
                    }
                    vec pos = aps->getPos() / AA_to_NM;
                    double Q = aps->getQ00();
                    string type = "mm";

                    orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, true);

                    if (split_dpl) {
                        vec tot_dpl = aps->getU1();
                        if (aps->getRank() > 0) {
                            tot_dpl += aps->getQ1();
                        }
                        // Calculate virtual charge positions
                        double a = dpl_spacing; // this is in nm
                        double mag_d = abs(tot_dpl); // this is in e * nm
                        vec dir_d_0 = tot_dpl.normalize();
                        vec dir_d = dir_d_0.normalize();
                        vec A = pos + 0.5 * a * dir_d / AA_to_NM; // converted to AA
                        vec B = pos - 0.5 * a * dir_d / AA_to_NM;
                        double qA = mag_d / a;
                        double qB = -qA;
                        // Zero out if magnitude small [e*nm]
                        if (aps->getIsoP() < 1e-9 || mag_d < 1e-9) {
                            A = aps->getPos() + 0.1 * a * vec(1, 0, 0); // != pos since self-energy may diverge
                            B = aps->getPos() - 0.1 * a * vec(1, 0, 0);
                            qA = 0;
                            qB = 0;
                        }
                        orb.AddAtom("A", A.x(), A.y(), A.z(), qA, true);
                        orb.AddAtom("B", B.x(), B.y(), B.z(), qB, true);
                    }
                }
            }

            // OUTER SHELL MM2
            for (unsigned int i = 0; i < ptop->MM2().size(); ++i) {
                ctp::PolarSeg *pseg = ptop->MM2()[i];
                for (unsigned int j = 0; j < pseg->size(); ++j) {

                    ctp::APolarSite *aps = (*pseg)[j];
                    vec pos = aps->getPos() / AA_to_NM;
                    double Q = aps->getQ00();
                    string type = "mm";

                    orb.AddAtom(aps->getName(), pos.x(), pos.y(), pos.z(), Q, true);
                }
            }
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

        void QMMIter::setE_FM(double ef00, double ef01, double ef02,
                double ef11, double ef12, double em0, double em1, double em2, double efm) {

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

        double QMMIter::getMMEnergy() {

            assert(_hasMM);
            return _ef_11 + _ef_12 + _em_1_ + _em_2_;
        }

        double QMMIter::getQMMMEnergy() {

            assert(_hasQM && _hasMM && _hasGWBSE);
            return _e_QM + + _e_GWBSE + _ef_11 + _ef_12 + _em_1_ + _em_2_;
        }




    }
}
