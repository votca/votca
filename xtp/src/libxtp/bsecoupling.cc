/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/bse.h"
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/bsecoupling.h"

namespace votca {
namespace xtp {
using namespace std;
using boost::format;
using namespace tools;

void BSECoupling::Initialize(Property& options) {
  string spintype = options.get("spin").as<std::string>();
  if (spintype == "all") {
    doSinglets_ = true;
    doTriplets_ = true;
  } else if (spintype == "triplet") {
    doTriplets_ = true;
    doSinglets_ = false;
  } else if (spintype == "singlet") {
    doSinglets_ = true;
    doTriplets_ = false;
  } else {
    throw std::runtime_error(
        (boost::format(
             "Choice % for type not known. Available singlet,triplet,all") %
         spintype)
            .str());
  }
  output_perturbation_ = options.get("use_perturbation").as<bool>();

  levA_ = options.get("moleculeA.states").as<Index>();
  levB_ = options.get("moleculeB.states").as<Index>();
  occA_ = options.get("moleculeA.occLevels").as<Index>();
  occB_ = options.get("moleculeB.occLevels").as<Index>();
  unoccA_ = options.get("moleculeA.unoccLevels").as<Index>();
  unoccB_ = options.get("moleculeB.unoccLevels").as<Index>();
}

// =============================================================================
// WriteMatrixToProperty
// Write an Eigen matrix as space-separated rows into an XML property node.
// Each row is stored as attribute "row_N". Unit conversion is applied if given.
// =============================================================================
void BSECoupling::WriteMatrixToProperty(tools::Property& prop,
                                        const std::string& name,
                                        const Eigen::MatrixXd& mat,
                                        double conversion) {
  tools::Property& mat_prop = prop.add(name, "");
  mat_prop.setAttribute("rows", mat.rows());
  mat_prop.setAttribute("cols", mat.cols());
  for (Index i = 0; i < mat.rows(); i++) {
    std::string row_str;
    for (Index j = 0; j < mat.cols(); j++) {
      if (j > 0) row_str += " ";
      row_str += (format("%1$1.6e") % (mat(i, j) * conversion)).str();
    }
    mat_prop.setAttribute("row_" + std::to_string(i), row_str);
  }
}

// =============================================================================
// WriteToProperty
// Write the pairwise effective couplings for one (stateA, stateB) pair.
// Backward compatible: attributes unchanged.
// =============================================================================
void BSECoupling::WriteToProperty(Property& summary, const QMState& stateA,
                                  const QMState& stateB) const {
  Property& coupling_summary = summary.add("coupling", "");
  double JAB_pert = 0;
  double JAB_diag = 0;
  if (stateA.Type() == QMStateType::Singlet) {
    JAB_pert =
        getSingletCouplingElement(stateA.StateIdx(), stateB.StateIdx(), 0);
    JAB_diag =
        getSingletCouplingElement(stateA.StateIdx(), stateB.StateIdx(), 1);
  } else if (stateA.Type() == QMStateType::Triplet) {
    JAB_pert =
        getTripletCouplingElement(stateA.StateIdx(), stateB.StateIdx(), 0);
    JAB_diag =
        getTripletCouplingElement(stateA.StateIdx(), stateB.StateIdx(), 1);
  }
  coupling_summary.setAttribute("stateA", stateA.ToString());
  coupling_summary.setAttribute("stateB", stateB.ToString());
  coupling_summary.setAttribute("j_pert",
                                (format("%1$1.6e") % JAB_pert).str());
  coupling_summary.setAttribute("j_diag",
                                (format("%1$1.6e") % JAB_diag).str());
}

// =============================================================================
// Addoutput
// Write all results to the XML property tree.
// Extends the existing pairwise coupling output with:
//   - diagnostics node (xi, PT/RM discrepancy, downfolding_safe flag)
//   - tb_matrices node (raw H and S blocks in FE+CT basis, in eV)
// The existing <coupling> elements are unchanged for backward compatibility.
// =============================================================================
void BSECoupling::Addoutput(Property& type_summary, const Orbitals&,
                            const Orbitals&) const {
  tools::Property& bsecoupling = type_summary.add(Identify(), "");
  string algorithm = "j_diag";
  if (output_perturbation_) {
    algorithm = "j_pert";
  }

  // Lambda to write one spin channel — avoids duplicating singlet/triplet logic
  auto WriteSpinChannel = [&](const std::string& spin_name,
                               const std::array<Eigen::MatrixXd, 2>& JAB,
                               const Eigen::MatrixXd& J_dimer,
                               const Eigen::MatrixXd& S_dimer,
                               const Diagnostics& diag,
                               const Eigen::VectorXd& monA_energies,
                               const Eigen::VectorXd& monB_energies,
                               QMStateType spin_type) {
    Property& spin_summary = bsecoupling.add(spin_name, "");
    spin_summary.setAttribute("algorithm", algorithm);

    // --- existing pairwise effective couplings (backward compatible) ---
    for (Index stateA = 0; stateA < levA_; ++stateA) {
      QMState qmstateA = QMState(spin_type, stateA, false);
      for (Index stateB = 0; stateB < levB_; ++stateB) {
        QMState qmstateB = QMState(spin_type, stateB, false);
        WriteToProperty(spin_summary, qmstateA, qmstateB);
      }
    }

    // --- monomer site energies ---
    // Isolated monomer excitation energies, pairwise-consistent for TB.
    // Environmental corrections should be added on top of these externally.
    {
      Property& mono_prop = spin_summary.add("monomer_energies", "");
      Property& propA = mono_prop.add("fragmentA", "");
      for (Index i = 0; i < monA_energies.size(); i++) {
        Property& e = propA.add("energy", "");
        QMState st(spin_type, i, false);
        e.setAttribute("state", st.ToString());
        e.setAttribute("eV", (format("%1$1.6e") % monA_energies(i)).str());
      }
      Property& propB = mono_prop.add("fragmentB", "");
      for (Index i = 0; i < monB_energies.size(); i++) {
        Property& e = propB.add("energy", "");
        QMState st(spin_type, i, false);
        e.setAttribute("state", st.ToString());
        e.setAttribute("eV", (format("%1$1.6e") % monB_energies(i)).str());
      }
    }

    // --- diagnostics ---
    Property& diag_prop = spin_summary.add("diagnostics", "");
    diag_prop.setAttribute("xi", (format("%1$1.4f") % diag.xi).str());
    diag_prop.setAttribute(
        "pt_rm_discrepancy_eV",
        (format("%1$1.6e") %
         (diag.pt_rm_discrepancy * votca::tools::conv::hrt2ev))
            .str());
    diag_prop.setAttribute("downfolding_safe",
                           diag.downfolding_safe ? "true" : "false");

    // --- raw TB matrices ---
    // Written in eV for H, dimensionless for S.
    // Block layout mirrors the comment in CalcJ_dimer:
    //   [ H_FE_FE  H_FE_CT ]   [ S_FE_FE  S_FE_CT ]
    //   [ H_CT_FE  H_CT_CT ]   [ S_CT_FE  S_CT_CT ]
    Index bse_exc = levA_ + levB_;
    Index ct = J_dimer.rows() - bse_exc;

    Property& tb_prop = spin_summary.add("tb_matrices", "");
    tb_prop.setAttribute("n_FE", bse_exc);
    tb_prop.setAttribute("n_CT", ct);
    tb_prop.setAttribute("n_occA", occA_);
    tb_prop.setAttribute("n_unoccA", unoccA_);
    tb_prop.setAttribute("n_occB", occB_);
    tb_prop.setAttribute("n_unoccB", unoccB_);

    // FE-FE block always present
    WriteMatrixToProperty(tb_prop, "H_FE_FE",
                          J_dimer.topLeftCorner(bse_exc, bse_exc),
                          votca::tools::conv::hrt2ev);
    WriteMatrixToProperty(tb_prop, "S_FE_FE",
                          S_dimer.topLeftCorner(bse_exc, bse_exc));

    if (ct > 0) {
      // FE-CT off-diagonal blocks
      WriteMatrixToProperty(tb_prop, "H_FE_CT",
                            J_dimer.topRightCorner(bse_exc, ct),
                            votca::tools::conv::hrt2ev);
      WriteMatrixToProperty(tb_prop, "S_FE_CT",
                            S_dimer.topRightCorner(bse_exc, ct));

      // CT-CT block
      WriteMatrixToProperty(tb_prop, "H_CT_CT",
                            J_dimer.bottomRightCorner(ct, ct),
                            votca::tools::conv::hrt2ev);
      WriteMatrixToProperty(tb_prop, "S_CT_CT",
                            S_dimer.bottomRightCorner(ct, ct));
    }
  };

  if (doSinglets_) {
    WriteSpinChannel(QMStateType(QMStateType::Singlet).ToLongString(),
                     JAB_singlet, J_dimer_singlet_, S_dimer_singlet_,
                     diag_singlet_,
                     monomerA_energies_singlet_, monomerB_energies_singlet_,
                     QMStateType(QMStateType::Singlet));
  }
  if (doTriplets_) {
    WriteSpinChannel(QMStateType(QMStateType::Triplet).ToLongString(),
                     JAB_triplet, J_dimer_triplet_, S_dimer_triplet_,
                     diag_triplet_,
                     monomerA_energies_triplet_, monomerB_energies_triplet_,
                     QMStateType(QMStateType::Triplet));
  }
}

double BSECoupling::getSingletCouplingElement(Index levelA, Index levelB,
                                              Index methodindex) const {
  return JAB_singlet[methodindex](levelA, levelB + levA_) *
         votca::tools::conv::hrt2ev;
}

double BSECoupling::getTripletCouplingElement(Index levelA, Index levelB,
                                              Index methodindex) const {
  return JAB_triplet[methodindex](levelA, levelB + levA_) *
         votca::tools::conv::hrt2ev;
}

Eigen::MatrixXd BSECoupling::SetupCTStates(Index bseA_vtotal, Index bseB_vtotal,
                                           Index bseAB_vtotal,
                                           Index bseAB_ctotal,
                                           const Eigen::MatrixXd& A_AB,
                                           const Eigen::MatrixXd& B_AB) const {
  Index noAB = occA_ * unoccB_;
  Index noBA = unoccA_ * occB_;
  Index bseAB_size = bseAB_vtotal * bseAB_ctotal;
  Eigen::MatrixXd CTstates = Eigen::MatrixXd::Zero(bseAB_size, noAB + noBA);

  auto A_occ = A_AB.middleCols(bseA_vtotal - occA_, occA_);
  auto A_unocc = A_AB.middleCols(bseA_vtotal, unoccA_);
  auto B_occ = B_AB.middleCols(bseB_vtotal - occB_, occB_);
  auto B_unocc = B_AB.middleCols(bseB_vtotal, unoccB_);

  const Eigen::MatrixXd A_occ_occ = A_occ.topRows(bseAB_vtotal);
  const Eigen::MatrixXd B_unocc_unocc = B_unocc.bottomRows(bseAB_ctotal);

  // notation AB is CT states with A+B-, BA is the counterpart
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   Setting up CT-states" << flush;

#pragma omp parallel for
  for (Index a_occ = 0; a_occ < occA_; a_occ++) {
    for (Index b_unocc = 0; b_unocc < unoccB_; b_unocc++) {
      Index index = a_occ * unoccB_ + b_unocc;
      Eigen::MatrixXd Coeff =
          B_unocc_unocc.col(b_unocc) * A_occ_occ.col(a_occ).transpose();
      CTstates.col(index) =
          Eigen::Map<Eigen::VectorXd>(Coeff.data(), bseAB_size);
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "  " << noBA << " CT states A+B- created" << flush;

  const Eigen::MatrixXd A_unocc_unocc = A_unocc.bottomRows(bseAB_ctotal);
  const Eigen::MatrixXd B_occ_occ = B_occ.topRows(bseAB_vtotal);

#pragma omp parallel for
  for (Index b_occ = 0; b_occ < occB_; b_occ++) {
    for (Index a_unocc = 0; a_unocc < unoccA_; a_unocc++) {
      Index index = b_occ * unoccA_ + a_unocc + noAB;
      Eigen::MatrixXd Coeff =
          A_unocc_unocc.col(a_unocc) * B_occ_occ.col(b_occ).transpose();
      CTstates.col(index) =
          Eigen::Map<Eigen::VectorXd>(Coeff.data(), bseAB_size);
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "  " << noBA << " CT states A-B+ created" << flush;
  return CTstates;
}

Eigen::MatrixXd BSECoupling::ProjectFrenkelExcitons(
    const Eigen::MatrixXd& BSE_Coeffs, const Eigen::MatrixXd& X_AB,
    Index bseX_vtotal, Index bseX_ctotal, Index bseAB_vtotal,
    Index bseAB_ctotal) const {
  Index bseAB_size = bseAB_vtotal * bseAB_ctotal;
  auto X_occ = X_AB.leftCols(bseX_vtotal);
  auto X_unocc = X_AB.rightCols(bseX_ctotal);
  const Eigen::MatrixXd X_occ_occ = X_occ.topRows(bseAB_vtotal);
  const Eigen::MatrixXd X_unocc_unocc = X_unocc.bottomRows(bseAB_ctotal);
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(bseAB_size, BSE_Coeffs.cols());
  // no pragma here because often we will only have one Coeff
  for (Index i = 0; i < BSE_Coeffs.cols(); i++) {
    Eigen::VectorXd coeff = BSE_Coeffs.col(i);
    Eigen::Map<Eigen::MatrixXd> coeffmatrix =
        Eigen::Map<Eigen::MatrixXd>(coeff.data(), bseX_ctotal, bseX_vtotal);
    Eigen::MatrixXd proj =
        X_unocc_unocc * coeffmatrix * X_occ_occ.transpose();
    result.col(i) = Eigen::Map<Eigen::VectorXd>(proj.data(), proj.size());
  }
  return result;
}

int GetSign(double value) {
  if (value < 0) {
    return -1;
  } else if (value > 0) {
    return 1;
  }
  return 0;
}

void BSECoupling::CalculateCouplings(const Orbitals& orbitalsA,
                                     const Orbitals& orbitalsB,
                                     const Orbitals& orbitalsAB) {
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "  Calculating exciton couplings" << flush;
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Using "
                              << OPENMP::getMaxThreads() << " threads" << flush;

  CheckAtomCoordinates(orbitalsA, orbitalsB, orbitalsAB);

  Index basisB = orbitalsB.getBasisSetSize();
  Index basisA = orbitalsA.getBasisSetSize();
  if ((basisA == 0) || (basisB == 0)) {
    throw std::runtime_error("Basis set size is not stored in monomers");
  }

  // get exciton information of molecule A
  Index bseA_cmax = orbitalsA.getBSEcmax();
  Index bseA_cmin = orbitalsA.getBSEcmin();
  Index bseA_vmax = orbitalsA.getBSEvmax();
  Index bseA_vmin = orbitalsA.getBSEvmin();
  Index bseA_vtotal = bseA_vmax - bseA_vmin + 1;
  Index bseA_ctotal = bseA_cmax - bseA_cmin + 1;
  Index bseA_total = bseA_vtotal + bseA_ctotal;
  Index bseA_size = bseA_vtotal * bseA_ctotal;
  Index bseA_singlet_exc = orbitalsA.BSESinglets().eigenvectors().cols();
  Index bseA_triplet_exc = orbitalsA.BSETriplets().eigenvectors().cols();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   molecule A has " << bseA_singlet_exc
      << " singlet excitons with dimension " << bseA_size << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   molecule A has " << bseA_triplet_exc
      << " triplet excitons with dimension " << bseA_size << flush;

  // get exciton information of molecule B
  Index bseB_cmax = orbitalsB.getBSEcmax();
  Index bseB_cmin = orbitalsB.getBSEcmin();
  Index bseB_vmax = orbitalsB.getBSEvmax();
  Index bseB_vmin = orbitalsB.getBSEvmin();
  Index bseB_vtotal = bseB_vmax - bseB_vmin + 1;
  Index bseB_ctotal = bseB_cmax - bseB_cmin + 1;
  Index bseB_total = bseB_vtotal + bseB_ctotal;
  Index bseB_size = bseB_vtotal * bseB_ctotal;
  Index bseB_singlet_exc = orbitalsB.BSESinglets().eigenvectors().cols();
  Index bseB_triplet_exc = orbitalsB.BSETriplets().eigenvectors().cols();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   molecule B has " << bseB_singlet_exc
      << " singlet excitons with dimension " << bseB_size << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   molecule B has " << bseB_triplet_exc
      << " triplet excitons with dimension " << bseB_size << flush;

  if (doSinglets_) {
    if (levA_ > bseA_singlet_exc) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << "  Number of excitons you want is greater than stored for molecule "
             "A. Setting to max number available"
          << flush;
      levA_ = bseA_singlet_exc;
    }
    if (levB_ > bseB_singlet_exc) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << "  Number of excitons you want is greater than stored for molecule "
             "B. Setting to max number available"
          << flush;
      levB_ = bseB_singlet_exc;
    }
  }
  if (doTriplets_) {
    if (levA_ > bseA_triplet_exc) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << "  Number of Frenkel states you want is greater than stored for "
             "molecule A. Setting to max number available"
          << flush;
      levA_ = bseA_triplet_exc;
    }
    if (levB_ > bseB_triplet_exc) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << "  Number of Frenkel states you want is greater than stored for "
             "molecule B. Setting to max number available"
          << flush;
      levB_ = bseB_triplet_exc;
    }
  }
  if (unoccA_ > bseA_ctotal || unoccA_ < 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << "  Number of unoccupied orbitals in molecule A for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    unoccA_ = bseA_ctotal;
  }
  if (unoccB_ > bseB_ctotal || unoccB_ < 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << "  Number of unoccupied orbitals in molecule B for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    unoccB_ = bseB_ctotal;
  }
  if (occA_ > bseA_vtotal || occA_ < 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << "  Number of occupied orbitals in molecule A for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    occA_ = bseA_vtotal;
  }
  if (occB_ > bseB_vtotal || occB_ < 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << "  Number of occupied orbitals in molecule B for CT creation "
           "exceeds number of KS-orbitals in BSE"
        << flush;
    occB_ = bseB_vtotal;
  }

  // get exciton information of pair AB
  Index bseAB_cmax = orbitalsAB.getBSEcmax();
  Index bseAB_cmin = orbitalsAB.getBSEcmin();
  Index bseAB_vmax = orbitalsAB.getBSEvmax();
  Index bseAB_vmin = orbitalsAB.getBSEvmin();
  Index bseAB_vtotal = bseAB_vmax - bseAB_vmin + 1;
  Index bseAB_ctotal = bseAB_cmax - bseAB_cmin + 1;
  Index bseAB_total = bseAB_vtotal + bseAB_ctotal;
  Index bseAB_size = bseAB_vtotal * bseAB_ctotal;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   levels used for BSE of molA: " << bseA_vmin
      << " to " << bseA_cmax << " total: " << bseA_total << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   levels used for BSE of molB: " << bseB_vmin
      << " to " << bseB_cmax << " total: " << bseB_total << flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "   levels used for BSE of dimer AB: " << bseAB_vmin
      << " to " << bseAB_cmax << " total: " << bseAB_total << flush;

  Eigen::MatrixXd MOsA =
      orbitalsA.MOs().eigenvectors().middleCols(bseA_vmin, bseA_total);
  Eigen::MatrixXd MOsB =
      orbitalsB.MOs().eigenvectors().middleCols(bseB_vmin, bseB_total);
  Eigen::MatrixXd MOsAB =
      orbitalsAB.MOs().eigenvectors().middleCols(bseAB_vmin, bseAB_total);

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Calculating overlap matrix for basisset: "
      << orbitalsAB.getDFTbasisName() << flush;

  Eigen::MatrixXd overlap = CalculateOverlapMatrix(orbitalsAB) * MOsAB;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Projecting monomers onto dimer orbitals" << flush;

  Eigen::MatrixXd A_AB = overlap.topRows(basisA).transpose() * MOsA;
  Eigen::MatrixXd B_AB = overlap.bottomRows(basisB).transpose() * MOsB;
  Eigen::VectorXd mag_A = A_AB.colwise().squaredNorm();
  if (mag_A.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer A on dimer is insufficient,mag="
        << mag_A.minCoeff() << flush;
  }
  Eigen::VectorXd mag_B = B_AB.colwise().squaredNorm();
  if (mag_B.any() < 0.95) {
    XTP_LOG(Log::error, *pLog_)
        << "\nWarning: "
        << "Projection of orbitals of monomer B on dimer is insufficient,mag="
        << mag_B.minCoeff() << flush;
  }

  AOBasis dftbasis = orbitalsAB.getDftBasis();
  AOBasis auxbasis = orbitalsAB.getAuxBasis();
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(auxbasis.AOBasisSize(), orbitalsAB.getRPAmin(),
                 orbitalsAB.getGWAmax(), orbitalsAB.getRPAmin(),
                 orbitalsAB.getRPAmax());
  Mmn.Fill(auxbasis, dftbasis, orbitalsAB.MOs().eigenvectors());

  const Eigen::MatrixXd& qpcoeff = orbitalsAB.QPdiag().eigenvectors();
  Eigen::MatrixXd Hqp = qpcoeff *
                        orbitalsAB.QPdiag().eigenvalues().asDiagonal() *
                        qpcoeff.transpose();

  BSE::options opt;
  opt.cmax = orbitalsAB.getBSEcmax();
  opt.homo = orbitalsAB.getHomo();
  opt.qpmin = orbitalsAB.getGWAmin();
  opt.qpmax = orbitalsAB.getGWAmax();
  opt.rpamax = orbitalsAB.getRPAmax();
  opt.rpamin = orbitalsAB.getRPAmin();
  opt.useTDA = true;
  opt.vmin = orbitalsAB.getBSEvmin();
  opt.use_Hqp_offdiag = orbitalsAB.GetFlagUseHqpOffdiag();
  BSE bse(*pLog_, Mmn);
  bse.configure(opt, orbitalsAB.RPAInputEnergies(), Hqp);
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Setup BSE operator" << flush;

  if (doSinglets_) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   Evaluating singlets" << flush;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   Setup Hamiltonian" << flush;

    // Store isolated monomer singlet excitation energies [eV].
    // These are pairwise-consistent site energies for TB assembly —
    // independent of which dimer partner is used.
    monomerA_energies_singlet_ =
        orbitalsA.BSESinglets().eigenvalues().head(levA_) *
        votca::tools::conv::hrt2ev;
    monomerB_energies_singlet_ =
        orbitalsB.BSESinglets().eigenvalues().head(levB_) *
        votca::tools::conv::hrt2ev;

    Eigen::MatrixXd FE_AB = Eigen::MatrixXd::Zero(bseAB_size, levA_ + levB_);
    const Eigen::MatrixXd bseA =
        orbitalsA.BSESinglets().eigenvectors().leftCols(levA_);
    FE_AB.leftCols(levA_) = ProjectFrenkelExcitons(
        bseA, A_AB, bseA_vtotal, bseA_ctotal, bseAB_vtotal, bseAB_ctotal);
    const Eigen::MatrixXd bseB =
        orbitalsB.BSESinglets().eigenvectors().leftCols(levB_);
    FE_AB.rightCols(levB_) = ProjectFrenkelExcitons(
        bseB, B_AB, bseB_vtotal, bseB_ctotal, bseAB_vtotal, bseAB_ctotal);
    Eigen::MatrixXd CTStates = SetupCTStates(
        bseA_vtotal, bseB_vtotal, bseAB_vtotal, bseAB_ctotal, A_AB, B_AB);
    JAB_singlet =
        ProjectExcitons(FE_AB, CTStates, bse.getSingletOperator_TDA(),
                        J_dimer_singlet_, S_dimer_singlet_, diag_singlet_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   calculated singlet couplings " << flush;
  }
  if (doTriplets_) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   Evaluating triplets" << flush;

    // Store isolated monomer triplet excitation energies [eV].
    monomerA_energies_triplet_ =
        orbitalsA.BSETriplets().eigenvalues().head(levA_) *
        votca::tools::conv::hrt2ev;
    monomerB_energies_triplet_ =
        orbitalsB.BSETriplets().eigenvalues().head(levB_) *
        votca::tools::conv::hrt2ev;

    Eigen::MatrixXd FE_AB = Eigen::MatrixXd::Zero(bseAB_size, levA_ + levB_);
    const Eigen::MatrixXd bseA =
        orbitalsA.BSETriplets().eigenvectors().leftCols(levA_);
    FE_AB.leftCols(levA_) = ProjectFrenkelExcitons(
        bseA, A_AB, bseA_vtotal, bseA_ctotal, bseAB_vtotal, bseAB_ctotal);
    const Eigen::MatrixXd bseB =
        orbitalsB.BSETriplets().eigenvectors().leftCols(levB_);
    FE_AB.rightCols(levB_) = ProjectFrenkelExcitons(
        bseB, B_AB, bseB_vtotal, bseB_ctotal, bseAB_vtotal, bseAB_ctotal);
    Eigen::MatrixXd CTStates = SetupCTStates(
        bseA_vtotal, bseB_vtotal, bseAB_vtotal, bseAB_ctotal, A_AB, B_AB);
    JAB_triplet =
        ProjectExcitons(FE_AB, CTStates, bse.getTripletOperator_TDA(),
                        J_dimer_triplet_, S_dimer_triplet_, diag_triplet_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << "   calculated triplet couplings " << flush;
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << "  Done with exciton couplings" << flush;
}

// =============================================================================
// OrthogonalizeCTs
// Merge FE and CT projection vectors into a single matrix without any
// pre-orthogonalization. The non-orthogonality between FE and CT states,
// and within each block, is handled correctly downstream:
//   - For the reduction method: by the joint Lowdin in CalcJ_dimer.
//   - For TB use of J_dimer/S_dimer: by the generalized eigenvalue problem
//     or SRG.
//
// Pre-orthogonalizing CT states against FE states (as done previously) used
// the incorrect projector P = F*F^T, which is only exact when FE_AB columns
// are orthonormal. Since projected monomer FE states are generally not
// orthonormal, this introduced a systematic error — particularly for larger
// CT manifolds (many occ/unocc pairs) where higher orbital pairs have
// non-negligible overlap with the FE space. Removing the pre-orthogonalization
// preserves the physical FE and CT state characters and lets the subsequent
// linear algebra handle the non-orthogonality correctly and consistently.
//
// A linear dependence check is performed on the CT states: if any CT state
// has norm below threshold after projection onto the dimer basis, it indicates
// near-linear dependence with other basis states and is flagged.
// =============================================================================
Eigen::MatrixXd BSECoupling::OrthogonalizeCTs(Eigen::MatrixXd& FE_AB,
                                              Eigen::MatrixXd& CTStates) const {
  Index ct = CTStates.cols();
  Index bse_exc = levA_ + levB_;
  Index bseAB_size = FE_AB.rows();

  // Check for near-linear dependence in CT states via their norms.
  // A very small norm indicates a CT state that is nearly identical to
  // another basis state and would cause numerical problems in the Lowdin step.
  if (ct > 0) {
    Eigen::VectorXd norm = CTStates.colwise().norm();
    Index minstateindex = 0;
    double minnorm = norm.minCoeff(&minstateindex);
    if (minnorm < 0.95) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " WARNING: CT-state " << minstateindex
          << " norm is only " << minnorm
          << " -- near-linear dependence in projection basis" << flush;
    }
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Merging " << ct
        << " CT states into projection (no pre-orthogonalization)" << flush;
  }

  Eigen::MatrixXd projection(bseAB_size, bse_exc + ct);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " merging projections into one vector  " << flush;
  projection.leftCols(bse_exc) = FE_AB;
  FE_AB.resize(0, 0);
  if (ct > 0) {
    projection.rightCols(ct) = CTStates;
  }
  CTStates.resize(0, 0);
  return projection;
}

// =============================================================================
// CalcJ_dimer
// Form J_dimer and S_dimer from the projection, expose them via output
// references for TB use, then Lowdin-orthogonalize to produce J_ortho.
// The block structure of J_dimer and S_dimer is:
//   [ H_FE_FE  H_FE_CT ]   [ S_FE_FE  S_FE_CT ]
//   [ H_CT_FE  H_CT_CT ]   [ S_CT_FE  S_CT_CT ]
// =============================================================================
template <class BSE_OPERATOR>
Eigen::MatrixXd BSECoupling::CalcJ_dimer(BSE_OPERATOR& H,
                                         Eigen::MatrixXd& projection,
                                         Eigen::MatrixXd& J_dimer_out,
                                         Eigen::MatrixXd& S_dimer_out) const {
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << "   Setting up coupling matrix size "
      << projection.cols() << flush;

  // this only works for hermitian/symmetric H so only in TDA
  J_dimer_out = projection.transpose() * (H * projection).eval();

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << "   Setting up overlap matrix size "
      << projection.cols() << flush;
  S_dimer_out = projection.transpose() * projection;

  projection.resize(0, 0);

  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;
  XTP_LOG(Log::debug, *pLog_) << " J_dimer_[Hrt]" << flush;
  XTP_LOG(Log::debug, *pLog_) << J_dimer_out << flush;
  XTP_LOG(Log::debug, *pLog_) << " S_dimer_" << flush;
  XTP_LOG(Log::debug, *pLog_) << S_dimer_out << flush;
  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_dimer_out);
  XTP_LOG(Log::debug, *pLog_)
      << TimeStamp() << "   Smallest value of dimer overlapmatrix is "
      << es.eigenvalues()(0) << flush;

  Eigen::MatrixXd Sm1 = es.operatorInverseSqrt();
  Eigen::MatrixXd J_ortho = Sm1 * J_dimer_out * Sm1;

  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;
  XTP_LOG(Log::debug, *pLog_) << " J_ortho_[Hrt]" << flush;
  XTP_LOG(Log::debug, *pLog_) << J_ortho << flush;
  XTP_LOG(Log::debug, *pLog_) << " S_-1/2" << flush;
  XTP_LOG(Log::debug, *pLog_) << Sm1 << flush;
  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;

  return J_ortho;
}

// =============================================================================
// ComputeDiagnostics
// Assess whether CT downfolding is reliable for this pair.
//   xi:                max |H_FE_CT| / |E_FE - E_CT| over all FE-CT pairs.
//                      Small xi => perturbation theory is safe.
//   pt_rm_discrepancy: max |J_pert - J_diag| over all FE pairs [Hrt].
//                      Large discrepancy flags near-resonant CT states.
//   downfolding_safe:  true if xi < 0.3 and discrepancy < 1e-4 Hrt.
// Note: J_dimer here is the raw pre-Lowdin matrix (in Hrt), not J_ortho.
// =============================================================================
BSECoupling::Diagnostics BSECoupling::ComputeDiagnostics(
    const Eigen::MatrixXd& J_dimer, const Eigen::MatrixXd& J_pert,
    const Eigen::MatrixXd& J_diag) const {
  Index bse_exc = levA_ + levB_;
  Index ct = J_dimer.rows() - bse_exc;
  Diagnostics diag;

  if (ct > 0) {
    for (Index i = 0; i < bse_exc; i++) {
      double Ei = J_dimer(i, i);
      for (Index k = bse_exc; k < bse_exc + ct; k++) {
        double Ek = J_dimer(k, k);
        double dE = std::abs(Ei - Ek);
        double coupling = std::abs(J_dimer(i, k));
        if (dE > 1e-10) {
          diag.xi = std::max(diag.xi, coupling / dE);
        } else {
          // exact degeneracy — downfolding is meaningless
          diag.xi = std::numeric_limits<double>::infinity();
        }
      }
    }
  }

  // PT vs RM discrepancy over all FE pairs
  for (Index i = 0; i < levA_; i++) {
    for (Index j = 0; j < levB_; j++) {
      Index jd = j + levA_;
      double diff = std::abs(J_pert(i, jd) - J_diag(i, jd));
      diag.pt_rm_discrepancy = std::max(diag.pt_rm_discrepancy, diff);
    }
  }

  diag.downfolding_safe =
      std::isfinite(diag.xi) && (diag.xi < 0.3) &&
      (diag.pt_rm_discrepancy < 1e-4);

  return diag;
}

// =============================================================================
// ProjectExcitons
// Top-level driver: orthogonalize CTs, form J_dimer/S_dimer, Lowdin-
// orthogonalize, run perturbation and reduction methods, compute diagnostics.
// Raw J_dimer and S_dimer (pre-Lowdin, in Hrt) and diagnostics are returned
// via output references for storage as member variables.
// =============================================================================
template <class BSE_OPERATOR>
std::array<Eigen::MatrixXd, 2> BSECoupling::ProjectExcitons(
    Eigen::MatrixXd& FE_AB, Eigen::MatrixXd& CTStates, BSE_OPERATOR H,
    Eigen::MatrixXd& J_dimer_out, Eigen::MatrixXd& S_dimer_out,
    Diagnostics& diag_out) const {

  Eigen::MatrixXd projection = OrthogonalizeCTs(FE_AB, CTStates);
  Eigen::MatrixXd J_ortho =
      CalcJ_dimer(H, projection, J_dimer_out, S_dimer_out);

  std::array<Eigen::MatrixXd, 2> J;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << "   Running Perturbation algorithm" << flush;
  J[0] = Perturbation(J_ortho);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << "    Running Projection algorithm" << flush;
  J[1] = Fulldiag(J_ortho);

  diag_out = ComputeDiagnostics(J_dimer_out, J[0], J[1]);

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << "   Diagnostics: xi=" << diag_out.xi
      << " PT/RM discrepancy[Hrt]=" << diag_out.pt_rm_discrepancy
      << " downfolding_safe=" << (diag_out.downfolding_safe ? "true" : "false")
      << flush;

  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;
  XTP_LOG(Log::debug, *pLog_) << "Jeff_pert[Hrt]" << flush;
  XTP_LOG(Log::debug, *pLog_) << J[0] << flush;
  XTP_LOG(Log::debug, *pLog_) << "Jeff_diag[Hrt]" << flush;
  XTP_LOG(Log::debug, *pLog_) << J[1] << flush;
  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;

  return J;
}

// =============================================================================
// Perturbation
// First-order perturbation theory correction for CT-mediated coupling.
// The factor of 0.5 in the correction term comes from the symmetric treatment
// of corrections |delta Phi_A> and |delta Phi_B> in eq 28 of Wehner/Baumeier
// JCTC 2017 — each contributes half the total CT-mediated coupling.
// =============================================================================
Eigen::MatrixXd BSECoupling::Perturbation(
    const Eigen::MatrixXd& J_dimer) const {
  Index bse_exc = levA_ + levB_;
  Index ct = J_dimer.rows() - bse_exc;
  Eigen::MatrixXd J_result = J_dimer;
  if (ct > 0) {
    Eigen::MatrixXd transformation =
        Eigen::MatrixXd::Identity(J_dimer.rows(), J_dimer.cols());
    Eigen::MatrixXd Ct = J_dimer.bottomRightCorner(ct, ct);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ct);
    transformation.bottomRightCorner(ct, ct) = es.eigenvectors();
    Ct.resize(0, 0);

    XTP_LOG(Log::debug, *pLog_) << "FE state hamiltonian" << flush;
    XTP_LOG(Log::debug, *pLog_)
        << J_dimer.topLeftCorner(bse_exc, bse_exc) << flush;
    XTP_LOG(Log::debug, *pLog_) << "eigenvalues of CT states" << flush;
    XTP_LOG(Log::debug, *pLog_) << es.eigenvalues().transpose() << flush;

    J_result = transformation.transpose() * J_dimer * transformation;
    XTP_LOG(Log::debug, *pLog_)
        << "---------------------------------------" << flush;
    XTP_LOG(Log::debug, *pLog_) << " J_ortho_[Hrt] CT-state diag" << flush;
    XTP_LOG(Log::debug, *pLog_) << J_result << flush;
    XTP_LOG(Log::debug, *pLog_)
        << "---------------------------------------" << flush;
  }

  Eigen::MatrixXd Jmatrix = Eigen::MatrixXd::Zero(bse_exc, bse_exc);
  for (Index stateA = 0; stateA < levA_; stateA++) {
    double Ea = J_result(stateA, stateA);
    for (Index stateB = 0; stateB < levB_; stateB++) {
      Index stateBd = stateB + levA_;
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << "   Calculating coupling between exciton A"
          << stateA + 1 << " and exciton B" << stateB + 1 << flush;
      double J = J_result(stateA, stateBd);
      double Eb = J_result(stateBd, stateBd);
      for (Index k = bse_exc; k < (bse_exc + ct); k++) {
        double Eab = J_result(k, k);
        if (std::abs(Eab - Ea) < 0.001) {
          XTP_LOG(Log::error, *pLog_)
              << TimeStamp() << "Energydifference between state A "
              << stateA + 1 << "and CT state " << k + 1 << " is " << Eab - Ea
              << "[Hrt]" << flush;
        }
        if (std::abs(Eab - Eb) < 0.001) {
          XTP_LOG(Log::error, *pLog_)
              << TimeStamp() << "Energydifference between state B "
              << stateB + 1 << "and CT state " << k + 1 << " is " << Eab - Eb
              << "[Hrt]" << flush;
        }
        // Factor 0.5 from symmetric first-order correction: see eq 28,
        // Wehner/Baumeier JCTC 2017. Each of |delta Phi_A> and |delta Phi_B>
        // contributes equally, giving the 1/2 * (1/Delta_A + 1/Delta_B) form.
        J += 0.5 * J_result(k, stateA) * J_result(k, stateBd) *
             (1.0 / (Ea - Eab) + 1.0 / (Eb - Eab));
      }
      Jmatrix(stateA, stateBd) = J;
      Jmatrix(stateBd, stateA) = J;
    }
  }
  return Jmatrix;
}

// =============================================================================
// Fulldiag (reduction method)
// Diagonalizes J_ortho (full FE+CT system), then for each (stateA, stateB)
// pair independently extracts the two eigenvectors most resembling those
// states, projects to the 2x2 FE subspace, Lowdin-orthogonalizes, and
// rotates back to recover the effective coupling. This per-pair 2x2 approach
// is consistent with the reduction method of Wehner/Baumeier JCTC 2017 and
// is robust at CT/FE near-resonance where perturbation theory diverges.
// =============================================================================
Eigen::MatrixXd BSECoupling::Fulldiag(const Eigen::MatrixXd& J_dimer) const {
  Index bse_exc = levA_ + levB_;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J_dimer);
  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;
  XTP_LOG(Log::debug, *pLog_) << "Eigenvectors of J" << flush;
  XTP_LOG(Log::debug, *pLog_) << es.eigenvectors() << flush;
  XTP_LOG(Log::debug, *pLog_) << "J_eigenvalues[Hrt]" << flush;
  XTP_LOG(Log::debug, *pLog_) << es.eigenvalues() << flush;
  XTP_LOG(Log::debug, *pLog_)
      << "---------------------------------------" << flush;

  Eigen::MatrixXd Jmat = Eigen::MatrixXd::Zero(bse_exc, bse_exc);
  for (Index stateA = 0; stateA < levA_; stateA++) {
    for (Index stateB = 0; stateB < levB_; stateB++) {
      Index stateBd = stateB + levA_;
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << "   Calculating coupling between exciton A"
          << stateA + 1 << " and exciton B" << stateB + 1 << flush;

      std::array<int, 2> indexes;
      std::array<int, 2> signs;

      // Find the eigenstates closest to stateA and stateBd by maximum overlap
      es.eigenvectors().row(stateA).cwiseAbs().maxCoeff(&indexes[0]);
      es.eigenvectors().row(stateBd).cwiseAbs().maxCoeff(&indexes[1]);
      if (indexes[0] == indexes[1]) {
        // same eigenvector selected for both: pick next best for stateBd
        Eigen::RowVectorXd stateamplitudes =
            es.eigenvectors().row(stateBd).cwiseAbs();
        stateamplitudes[indexes[1]] = 0.0;
        stateamplitudes.maxCoeff(&indexes[1]);
      }

      signs[0] = GetSign(es.eigenvectors()(stateA, indexes[0]));
      signs[1] = GetSign(es.eigenvectors()(stateBd, indexes[1]));

      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << "   Order is: [Initial state n->nth eigenvalue]"
          << flush;
      XTP_LOG(Log::info, *pLog_) << "    A" << stateA + 1 << ":" << stateA + 1
                                 << "->" << indexes[0] + 1 << " ";
      XTP_LOG(Log::info, *pLog_) << "    B" << stateB + 1 << ":" << stateBd + 1
                                 << "->" << indexes[1] + 1 << " " << flush;

      // Build 2x2 transformation matrix from the two selected eigenvectors
      // (rows stateA and stateBd only) and diagonal eigenvalue matrix
      Eigen::Matrix2d Emat = Eigen::Matrix2d::Zero();
      Eigen::Matrix2d Tmat = Eigen::Matrix2d::Zero();
      for (Index i = 0; i < 2; i++) {
        Index k = indexes[i];
        double sign = signs[i];
        Tmat(0, i) = sign * es.eigenvectors()(stateA, k);
        Tmat(1, i) = sign * es.eigenvectors()(stateBd, k);
        Emat(i, i) = es.eigenvalues()(k);
      }
      Tmat.colwise().normalize();

      if (Tmat.determinant() < 0) {
        XTP_LOG(Log::info, *pLog_)
            << " Reduced state matrix is not in a right handed basis, "
               "multiplying second eigenvector by -1 "
            << flush;
        Tmat.col(1) *= -1;
      }

      XTP_LOG(Log::debug, *pLog_)
          << "---------------------------------------" << flush;
      XTP_LOG(Log::debug, *pLog_) << " T_" << flush;
      XTP_LOG(Log::debug, *pLog_) << Tmat << flush;

      // Lowdin orthogonalize the 2x2 projected subspace
      Eigen::Matrix2d S_small = Tmat * Tmat.transpose();
      XTP_LOG(Log::debug, *pLog_) << "S_small" << flush;
      XTP_LOG(Log::debug, *pLog_) << S_small << flush;

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> ss(S_small);
      Eigen::Matrix2d sm1 = ss.operatorInverseSqrt();
      Emat = sm1 * Emat * sm1;

      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << "   Smallest value of dimer overlapmatrix is "
          << ss.eigenvalues()(0) << flush;
      XTP_LOG(Log::debug, *pLog_) << "S-1/2" << flush;
      XTP_LOG(Log::debug, *pLog_) << sm1 << flush;
      XTP_LOG(Log::debug, *pLog_) << "E_ortho" << flush;
      XTP_LOG(Log::debug, *pLog_) << Emat << flush;

      Tmat = Tmat * sm1;
      XTP_LOG(Log::debug, *pLog_) << "T_ortho" << flush;
      XTP_LOG(Log::debug, *pLog_) << Tmat << flush;
      XTP_LOG(Log::debug, *pLog_)
          << "---------------------------------------" << flush;

      // H_eff = T * E * T^T; off-diagonal is the effective coupling
      Eigen::Matrix2d J_small = Tmat * Emat * Tmat.transpose();
      XTP_LOG(Log::debug, *pLog_) << "T_ortho*E_ortho*T_ortho^T" << flush;
      XTP_LOG(Log::debug, *pLog_) << J_small << flush;

      Jmat(stateA, stateBd) = J_small(0, 1);
      Jmat(stateBd, stateA) = J_small(1, 0);
    }
  }
  return Jmat;
}

}  // namespace xtp
}  // namespace votca
