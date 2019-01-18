/*
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <votca/xtp/bse.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/aomatrix.h>
#include "votca/xtp/qmstate.h"
#include "votca/xtp/vc2index.h"
#include "votca/xtp/populationanalysis.h"
using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void BSE::SetupDirectInteractionOperator() {
    RPA rpa = RPA(_Mmn);
    rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
    rpa.UpdateRPAInputEnergies(_orbitals.MOEnergies(), _Hqp.diagonal(), _opt.qpmin);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(rpa.calculate_epsilon_r(0));
    _Mmn.MultiplyRightWithAuxMatrix(es.eigenvectors());
    _epsilon_0_inv = VectorXfd::Zero(es.eigenvalues().size());
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        if (es.eigenvalues()(i) > 1e-8) {
            _epsilon_0_inv(i) = 1 / es.eigenvalues()(i);
        }
    }
}

void BSE::Solve_triplets() {
    MatrixXfd H = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hd<real_gwbse>(H);
    Add_Hqp<real_gwbse>(H);
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Setup TDA triplet hamiltonian " << flush;
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Solving for first " << _opt.nmax << " eigenvectors" << flush;
    tools::linalg_eigenvalues(H, _bse_triplet_energies, _bse_triplet_coefficients, _opt.nmax);
    return;
}

void BSE::Solve_singlets() {
    if (_opt.useTDA) {
        Solve_singlets_TDA();
    } else {
        Solve_singlets_BTDA();
    }
}

void BSE::Solve_singlets_TDA() {
    MatrixXfd H = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hd<real_gwbse>(H);
    Add_Hqp<real_gwbse>(H);
    Add_Hx<real_gwbse, 2>(H);
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Setup TDA singlet hamiltonian " << flush;
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Solving for first " << _opt.nmax << " eigenvectors" << flush;
    tools::linalg_eigenvalues(H, _bse_singlet_energies, _bse_singlet_coefficients, _opt.nmax);
    return;
}

void BSE::SetupHs() {
    _eh_s = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hd<real_gwbse>(_eh_s);
    Add_Hqp<real_gwbse>(_eh_s);
    Add_Hx<real_gwbse, 2>(_eh_s);
}

void BSE::SetupHt() {
    _eh_t = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hd<real_gwbse>(_eh_t);
    Add_Hqp<real_gwbse>(_eh_t);
}

void BSE::Solve_singlets_BTDA() {

    // For details of the method, see EPL,78(2007)12001,
    // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
    // setup resonant (A) and RARC blocks (B)
    //corresponds to
    // _ApB = (_eh_d +_eh_qp + _eh_d2 + 4.0 * _eh_x);
    // _AmB = (_eh_d +_eh_qp - _eh_d2);
    Eigen::MatrixXd ApB = Eigen::MatrixXd::Zero(_bse_size, _bse_size);
    Add_Hd<double>(ApB);
    Add_Hqp<double>(ApB);

    Eigen::MatrixXd AmB = ApB;
    Add_Hd2<double, -1 > (AmB);

    Add_Hx<double, 4>(ApB);
    Add_Hd2<double, 1>(ApB);
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Setup singlet hamiltonian " << flush;

    // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not positive definite
    //(A-B) is not needed any longer and can be overwritten
    XTP_LOG(logDEBUG, _log) << TimeStamp() << " Trying Cholesky decomposition of KAA-KAB" << flush;
    Eigen::LLT< Eigen::Ref<Eigen::MatrixXd> > L(AmB);

    for (int i = 0; i < AmB.rows(); ++i) {
        for (int j = i + 1; j < AmB.cols(); ++j) {
            AmB(i, j) = 0;
        }
    }

    if (L.info() != 0) {
        XTP_LOG(logDEBUG, _log) << TimeStamp() << " Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller basisset. This can indicate a triplet instability." << flush;
        throw std::runtime_error("Cholesky decompostion failed");
    } else {
        XTP_LOG(logDEBUG, _log) << TimeStamp() << " Cholesky decomposition of KAA-KAB was successful" << flush;
    }
    Eigen::MatrixXd temp = ApB*AmB;
    ApB.noalias() = AmB.transpose() * temp;
    temp.resize(0, 0);
    XTP_LOG(logDEBUG, _log) << TimeStamp() << " Calculated H = L^T(A+B)L " << flush;
    Eigen::VectorXd eigenvalues;
    Eigen::MatrixXd eigenvectors;
    XTP_LOG(logDEBUG, _log)
            << TimeStamp() << " Solving for first " << _opt.nmax << " eigenvectors" << flush;
    bool success_diag = tools::linalg_eigenvalues(ApB, eigenvalues, eigenvectors, _opt.nmax);
    if (!success_diag) {
        XTP_LOG(logDEBUG, _log) << TimeStamp() << " Could not solve problem" << flush;
    } else {
        XTP_LOG(logDEBUG, _log) << TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
    }
    ApB.resize(0, 0);
    eigenvalues = eigenvalues.cwiseSqrt();

#if (GWBSE_DOUBLE)
    _bse_singlet_energies = eigenvalues;
#else
    _bse_singlet_energies = eigenvalues.cast<float>();
#endif
    // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 + 1/sqrt(eps_l)L ] R_l
    //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 - 1/sqrt(eps_l)L ] R_l
    // determine inverse of L^T
    Eigen::MatrixXd LmT = AmB.inverse().transpose();
    int dim = LmT.rows();
    _bse_singlet_energies.resize(_opt.nmax);
    _bse_singlet_coefficients.resize(dim, _opt.nmax); // resonant part (_X_evec)
    _bse_singlet_coefficients_AR.resize(dim, _opt.nmax); // anti-resonant part (_Y_evec)
    for (int level = 0; level < _opt.nmax; level++) {
        double sqrt_eval = std::sqrt(_bse_singlet_energies(level));
        // get l-th reduced EV
#if (GWBSE_DOUBLE)
        _bse_singlet_coefficients.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT + AmB) * eigenvectors.col(level));
        _bse_singlet_coefficients_AR.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT - AmB) * eigenvectors.col(level));
#else
        _bse_singlet_coefficients.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT + AmB) * eigenvectors.col(level)).cast<float>();
        _bse_singlet_coefficients_AR.col(level) = (0.5 / sqrt_eval * (_bse_singlet_energies(level) * LmT - AmB) * eigenvectors.col(level)).cast<float>();
#endif

    }

    return;
}

template <typename T>
void BSE::Add_Hqp(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {

    int offset = _opt.vmin - _opt.qpmin;
    vc2index vc = vc2index(0, 0, _bse_ctotal);
#pragma omp parallel for
    for (int v1 = 0; v1 < _bse_vtotal; v1++) {
        for (int c1 = 0; c1 < _bse_ctotal; c1++) {
            int index_vc = vc.I(v1, c1);
            // v->c
            for (int c2 = 0; c2 < _bse_ctotal; c2++) {
                int index_vc2 = vc.I(v1, c2);
                H(index_vc2, index_vc) += _Hqp(c2 + _bse_vtotal + offset, c1 + _bse_vtotal + offset);
            }
            // c-> v
            for (int v2 = 0; v2 < _bse_vtotal; v2++) {
                int index_vc2 = vc.I(v2, c1);
                H(index_vc2, index_vc) -= _Hqp(v2 + offset, v1 + offset);
            }
        }
    }
    return;
}

template <typename T>
void BSE::Add_Hd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
    int auxsize = _Mmn.auxsize();
    vc2index vc = vc2index(0, 0, _bse_ctotal);
#pragma omp parallel for
    for (int v1 = 0; v1 < _bse_vtotal; v1++) {
        const MatrixXfd Mmn1T = (_Mmn[v1 + _opt.vmin ].block(_opt.vmin, 0, _bse_vtotal, auxsize) * _epsilon_0_inv.asDiagonal()).transpose();
        for (int c1 = 0; c1 < _bse_ctotal; c1++) {
            const MatrixXfd& Mmn2 = _Mmn[c1 + _bse_cmin];
            const MatrixXfd Mmn2xMmn1T = Mmn2.block(_bse_cmin, 0, _bse_ctotal, auxsize) * Mmn1T;
            int i1 = vc.I(v1, c1);
            for (int v2 = 0; v2 < _bse_vtotal; v2++) {
                for (int c2 = 0; c2 < _bse_ctotal; c2++) {
                    int i2 = vc.I(v2, c2);
                    H(i2, i1) -= Mmn2xMmn1T(c2, v2);
                }
            }
        }
    }

    return;
}

template <typename T, int factor>
void BSE::Add_Hd2(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
    int auxsize = _Mmn.auxsize();
    vc2index vc = vc2index(0, 0, _bse_ctotal);
#pragma omp parallel for       
    for (int c1 = 0; c1 < _bse_ctotal; c1++) {
        const MatrixXfd Mmn2T = factor * (_Mmn[c1 + _bse_cmin ].block(_opt.vmin, 0, _bse_vtotal, auxsize) * _epsilon_0_inv.asDiagonal()).transpose();
        for (int v1 = 0; v1 < _bse_vtotal; v1++) {
            const MatrixXfd& Mmn1 = _Mmn[v1 + _opt.vmin];
            MatrixXfd Mmn1xMmn2T = Mmn1.block(_bse_cmin, 0, _bse_ctotal, auxsize) * Mmn2T;
            int i1 = vc.I(v1, c1);
            for (int v2 = 0; v2 < _bse_vtotal; v2++) {
                for (int c2 = 0; c2 < _bse_ctotal; c2++) {
                    int i2 = vc.I(v2, c2);
                    H(i2, i1) -= Mmn1xMmn2T(c2, v2);
                }
            }
        }
    }
    return;
}

template <typename T, int factor>
void BSE::Add_Hx(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
    int auxsize = _Mmn.auxsize();
    vc2index vc = vc2index(0, 0, _bse_ctotal);
#pragma omp parallel for
    for (int v1 = 0; v1 < _bse_vtotal; v1++) {
        const MatrixXfd Mmn1 = factor * (_Mmn[v1 + _opt.vmin].block(_bse_cmin, 0, _bse_ctotal, auxsize)).transpose();
        for (int c1 = 0; c1 < _bse_ctotal; c1++) {
            int i1 = vc.I(v1, c1);
            for (int v2 = 0; v2 < _bse_vtotal; v2++) {
                const MatrixXfd& Mmn2 = _Mmn[v2 + _opt.vmin];
                const VectorXfd Mmnx2 = Mmn2.block(_bse_cmin, 0, _bse_ctotal, auxsize) * Mmn1.col(c1);
                for (int c2 = 0; c2 < _bse_ctotal; c2++) {
                    int i2 = vc.I(v2, c2);
                    H(i2, i1) += Mmnx2(c2);
                }
            }
        }
    }
    return;
}

void BSE::printFragInfo(const std::vector<QMFragment<BSE_Population> > & frags, int state) const {
    for (const QMFragment<BSE_Population>& frag : frags) {
        double dq = frag.value().H[state] - frag.value().E[state];
        double qeff = dq + frag.value().Gs;
        XTP_LOG(logINFO, _log) << format("           Fragment %1$8s%% -- hole: %2$5.1f%%  electron: %3$5.1f%%  dQ: %4$+5.2f  Qeff: %5$+5.2f")
                % frag.name() % (100.0 * frag.value().H[state]) % (100.0 * frag.value().E[state]) % dq % qeff << flush;
    }
    return;
}

void BSE::printWeights(int i_bse, double weight)const{

    vc2index vc = vc2index(_opt.vmin, _bse_cmin, _bse_ctotal);
    if (weight > _opt.min_print_weight) {
        XTP_LOG(logINFO, _log) << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                % (_opt.homo - vc.v(i_bse)) % (vc.c(i_bse) - _opt.homo - 1) % (100.0 * weight) << flush;
    }
    return;
}

void BSE::Analyze_singlets(std::vector<QMFragment<BSE_Population> >& singlets) {

    Interaction act;
    QMStateType singlet = QMStateType(QMStateType::Singlet);
    _orbitals.CalcCoupledTransition_Dipoles();
    std::vector<double> oscs = _orbitals.Oscillatorstrengths();
    if (tools::globals::verbose) {
        act = Analyze_eh_interaction(singlet);
    }
    if (singlets.size() > 0) {
        Lowdin low;
        low.CalcChargeperFragment(singlets, _orbitals, singlet);
    }

    double hrt2ev = tools::conv::hrt2ev;
    XTP_LOG(logINFO, _log) << "  ====== singlet energies (eV) ====== " << flush;
    for (int i = 0; i < _opt.nmax; ++i) {
        double osc = oscs[i];
        if (tools::globals::verbose) {
            XTP_LOG(logINFO, _log) << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f")
                    % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i)))
                    % (hrt2ev * act.qp_contrib(i)) % (hrt2ev * act.exchange_contrib(i)) % (hrt2ev * act.direct_contrib(i)) << flush;
        } else {
            XTP_LOG(logINFO, _log) << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                    % (i + 1) % (hrt2ev * _bse_singlet_energies(i)) % (1240.0 / (hrt2ev * _bse_singlet_energies(i))) << flush;
        }
        const Eigen::Vector3d& trdip = _orbitals.TransitionDipoles()[i];
        XTP_LOG(logINFO, _log) << format("           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f")
                % trdip[0] % trdip[1] % trdip[2] % (trdip.squaredNorm()) % osc << flush;
        for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
            // if contribution is larger than 0.2, print
            double weight = std::pow(_bse_singlet_coefficients(i_bse, i), 2);
            if (!_opt.useTDA) {
                weight -= std::pow(_bse_singlet_coefficients_AR(i_bse, i), 2);
            }
            printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (singlets.size() > 0) {
            printFragInfo(singlets, i);
        }

        XTP_LOG(logINFO, _log) << flush;
    }
    return;
}

void BSE::Analyze_triplets(std::vector<QMFragment<BSE_Population> >& triplets) {

    Interaction act;
    QMStateType triplet = QMStateType(QMStateType::Triplet);
    if (tools::globals::verbose) {
        act = Analyze_eh_interaction(triplet);
    }
    if (triplets.size() > 0) {
        Lowdin low;
        low.CalcChargeperFragment(triplets, _orbitals, triplet);
    }

    XTP_LOG(logINFO, _log) << "  ====== triplet energies (eV) ====== " << flush;
    for (int i = 0; i < _opt.nmax; ++i) {
        if (tools::globals::verbose) {
            XTP_LOG(logINFO, _log) << format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> = %4$+1.4f <K_d> = %5$+1.4f")
                    % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i)))
                    % (tools::conv::hrt2ev * act.qp_contrib(i)) % (tools::conv::hrt2ev * act.direct_contrib(i)) << flush;
        } else {
            XTP_LOG(logINFO, _log) << format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm")
                    % (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) % (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i))) << flush;
        }
        for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
            // if contribution is larger than 0.2, print
            double weight = pow(_bse_triplet_coefficients(i_bse, i), 2);
            printWeights(i_bse, weight);
        }
        // results of fragment population analysis 
        if (triplets.size() > 0) {
            printFragInfo(triplets, i);
        }
        XTP_LOG(logINFO, _log) << format("   ") << flush;
    }
    // storage to orbitals object

    return;
}

Eigen::VectorXd BSE::Analyze_IndividualContribution(const QMStateType& type, const MatrixXfd& H) {
    Eigen::VectorXd contrib = Eigen::VectorXd::Zero(_opt.nmax);
    if (type == QMStateType::Singlet) {
        for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
            MatrixXfd slice_R = _bse_singlet_coefficients.block(0, i_exc, _bse_size, 1);
            contrib(i_exc) = (slice_R.transpose() * H * slice_R).value();
            if (_bse_singlet_coefficients_AR.cols() > 0) {
                MatrixXfd slice_AR = _bse_singlet_coefficients_AR.block(0, i_exc, _bse_size, 1);
                // get anti-resonant contribution from direct Keh
                contrib(i_exc) -= (slice_AR.transpose() * H * slice_AR).value();
            }
        }
    } else if (type == QMStateType::Triplet) {
        for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
            MatrixXfd _slice_R = _bse_triplet_coefficients.block(0, i_exc, _bse_size, 1);
            contrib(i_exc) = (_slice_R.transpose() * H * _slice_R).value();
        }
    } else {
        throw std::runtime_error("BSE::Analyze_eh_interaction:Spin not known!");
    }
    return contrib;
}

BSE::Interaction BSE::Analyze_eh_interaction(const QMStateType& type) {

    Interaction analysis;
    MatrixXfd H = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hqp<real_gwbse>(H);
    analysis.qp_contrib = Analyze_IndividualContribution(type, H);

    H = MatrixXfd::Zero(_bse_size, _bse_size);
    Add_Hd<real_gwbse>(H);
    analysis.direct_contrib = Analyze_IndividualContribution(type, H);
    if (type == QMStateType::Singlet) {
        H = MatrixXfd::Zero(_bse_size, _bse_size);
        Add_Hx<real_gwbse, 2>(H);
        analysis.exchange_contrib = Analyze_IndividualContribution(type, H);
    } else {
        analysis.exchange_contrib = Eigen::VectorXd::Zero(0);
    }

    return analysis;
}

}
};
