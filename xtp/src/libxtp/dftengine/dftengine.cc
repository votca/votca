/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/IncrementalFockBuilder.h"
#include "votca/xtp/IndexParser.h"
#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/mmregion.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/pmdecomposition.h"

namespace votca {
namespace xtp {

void DFTEngine::Initialize(tools::Property& options) {

  const std::string key_xtpdft = "xtpdft";
  dftbasis_name_ = options.get(".basisset").as<std::string>();

  if (options.exists(".auxbasisset")) {
    auxbasis_name_ = options.get(".auxbasisset").as<std::string>();
  }

  if (!auxbasis_name_.empty()) {
    screening_eps_ = options.get(key_xtpdft + ".screening_eps").as<double>();
    fock_matrix_reset_ =
        options.get(key_xtpdft + ".fock_matrix_reset").as<Index>();
  }
  if (options.exists(".ecp")) {
    ecp_name_ = options.get(".ecp").as<std::string>();
  }

  initial_guess_ = options.get(".initial_guess").as<std::string>();

  grid_name_ = options.get(key_xtpdft + ".integration_grid").as<std::string>();
  xc_functional_name_ = options.get(".functional").as<std::string>();

  if (options.exists(key_xtpdft + ".externaldensity")) {
    integrate_ext_density_ = true;
    orbfilename_ =
        options.get(key_xtpdft + ".externaldensity.orbfile").as<std::string>();
    gridquality_ = options.get(key_xtpdft + ".externaldensity.gridquality")
                       .as<std::string>();
    state_ =
        options.get(key_xtpdft + ".externaldensity.state").as<std::string>();
  }

  if (options.exists(".externalfield")) {
    integrate_ext_field_ = true;
    extfield_ = options.get(".externalfield").as<Eigen::Vector3d>();
  }

  conv_opt_.Econverged =
      options.get(key_xtpdft + ".convergence.energy").as<double>();
  conv_opt_.error_converged =
      options.get(key_xtpdft + ".convergence.error").as<double>();
  max_iter_ =
      options.get(key_xtpdft + ".convergence.max_iterations").as<Index>();

  std::string method =
      options.get(key_xtpdft + ".convergence.method").as<std::string>();
  if (method == "DIIS") {
    conv_opt_.usediis = true;
  } else if (method == "mixing") {
    conv_opt_.usediis = false;
  }
  if (!conv_opt_.usediis) {
    conv_opt_.histlength = 1;
    conv_opt_.maxout = false;
  }
  conv_opt_.mixingparameter =
      options.get(key_xtpdft + ".convergence.mixing").as<double>();
  conv_opt_.levelshift =
      options.get(key_xtpdft + ".convergence.levelshift").as<double>();
  conv_opt_.levelshiftend =
      options.get(key_xtpdft + ".convergence.levelshift_end").as<double>();
  conv_opt_.maxout =
      options.get(key_xtpdft + ".convergence.DIIS_maxout").as<bool>();
  conv_opt_.histlength =
      options.get(key_xtpdft + ".convergence.DIIS_length").as<Index>();
  conv_opt_.diis_start =
      options.get(key_xtpdft + ".convergence.DIIS_start").as<double>();
  conv_opt_.adiis_start =
      options.get(key_xtpdft + ".convergence.ADIIS_start").as<double>();

  if (options.exists(key_xtpdft + ".dft_in_dft.activeatoms")) {
     active_atoms_as_string = options.get(key_xtpdft + ".dft_in_dft.activeatoms").as<std::string>();
    std::cout << active_atoms_as_string << std::endl;
  }
}

void DFTEngine::PrintMOs(const Eigen::VectorXd& MOEnergies, Log::Level level) {
  XTP_LOG(level, *pLog_) << "  Orbital energies: " << std::flush;
  XTP_LOG(level, *pLog_) << "  index occupation energy(Hartree) " << std::flush;
  for (Index i = 0; i < MOEnergies.size(); i++) {
    Index occupancy = 0;
    if (i < numofelectrons_ / 2) {
      occupancy = 2;
    }
    XTP_LOG(level, *pLog_) << (boost::format(" %1$5d      %2$1d   %3$+1.10f") %
                               i % occupancy % MOEnergies(i))
                                  .str()
                           << std::flush;
  }
  return;
}

void DFTEngine::CalcElDipole(const Orbitals& orb) const {
  QMState state = QMState("n");
  Eigen::Vector3d result = orb.CalcElDipole(state);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Electric Dipole is[e*bohr]:\n\t\t dx=" << result[0]
      << "\n\t\t dy=" << result[1] << "\n\t\t dz=" << result[2] << std::flush;
  return;
}

std::array<Eigen::MatrixXd, 2> DFTEngine::CalcERIs_EXX(
    const Eigen::MatrixXd& MOCoeff, const Eigen::MatrixXd& Dmat,
    double error) const {
  if (!auxbasis_name_.empty()) {
    if (conv_accelerator_.getUseMixing() || MOCoeff.rows() == 0) {
      return ERIs_.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), Dmat);
    } else {
      Eigen::MatrixXd occblock = MOCoeff.leftCols(numofelectrons_ / 2);
      return ERIs_.CalculateERIs_EXX_3c(occblock, Dmat);
    }
  } else {
    return ERIs_.CalculateERIs_EXX_4c(Dmat, error);
  }
}

Eigen::MatrixXd DFTEngine::CalcERIs(const Eigen::MatrixXd& Dmat,
                                    double error) const {
  if (!auxbasis_name_.empty()) {
    return ERIs_.CalculateERIs_3c(Dmat);
  } else {
    return ERIs_.CalculateERIs_4c(Dmat, error);
  }
}

tools::EigenSystem DFTEngine::IndependentElectronGuess(
    const Mat_p_Energy& H0) const {
  return conv_accelerator_.SolveFockmatrix(H0.matrix());
}

tools::EigenSystem DFTEngine::ModelPotentialGuess(
    const Mat_p_Energy& H0, const QMMolecule& mol,
    const Vxc_Potential<Vxc_Grid>& vxcpotential) const {
  Eigen::MatrixXd Dmat = AtomicGuess(mol);
  Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Vxc matrix " << std::flush;

  Eigen::MatrixXd H = H0.matrix() + e_vxc.matrix();

  if (ScaHFX_ > 0) {
    std::array<Eigen::MatrixXd, 2> both =
        CalcERIs_EXX(Eigen::MatrixXd::Zero(0, 0), Dmat, 1e-12);
    H += both[0];
    H += ScaHFX_ * both[1];
  } else {
    H += CalcERIs(Dmat, 1e-12);
  }
  return conv_accelerator_.SolveFockmatrix(H);
}

bool DFTEngine::Evaluate(Orbitals& orb) {
  Prepare(orb);
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());
  tools::EigenSystem MOs;
  MOs.eigenvalues() = Eigen::VectorXd::Zero(H0.cols());
  MOs.eigenvectors() = Eigen::MatrixXd::Zero(H0.rows(), H0.cols());
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  ConfigOrbfile(orb);

  if (initial_guess_ == "orbfile") {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Reading guess from orbitals object/file"
        << std::flush;
    MOs = orb.MOs();
    MOs.eigenvectors() = OrthogonalizeGuess(MOs.eigenvectors());
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Setup Initial Guess using: " << initial_guess_
        << std::flush;
    if (initial_guess_ == "independent") {
      MOs = IndependentElectronGuess(H0);
    } else if (initial_guess_ == "atom") {
      MOs = ModelPotentialGuess(H0, orb.QMAtoms(), vxcpotential);
    } else {
      throw std::runtime_error("Initial guess method not known/implemented");
    }
  }

  Eigen::MatrixXd Dmat = conv_accelerator_.DensityMatrix(MOs);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Guess Matrix gives N=" << std::setprecision(9)
      << Dmat.cwiseProduct(dftAOoverlap_.Matrix()).sum() << " electrons."
      << std::flush;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " STARTING SCF cycle" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << " ----------------------------------------------"
         "----------------------------"
      << std::flush;

  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(Dmat.rows(), Dmat.cols());
  Eigen::MatrixXd K;
  if (ScaHFX_ > 0) {
    K = Eigen::MatrixXd::Zero(Dmat.rows(), Dmat.cols());
  }

  double start_incremental_F_threshold = 1e-4;
  if (!auxbasis_name_.empty()) {
    start_incremental_F_threshold = 0.0;  // Disable if RI is used
  }
  IncrementalFockBuilder incremental_fock(*pLog_, start_incremental_F_threshold,
                                          fock_matrix_reset_);
  incremental_fock.Configure(Dmat);

  for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << max_iter_ << std::flush;
    Mat_p_Energy e_vxc = vxcpotential.IntegrateVXC(Dmat);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled DFT Vxc matrix " << std::flush;

    Eigen::MatrixXd H = H0.matrix() + e_vxc.matrix();
    double Eone = Dmat.cwiseProduct(H0.matrix()).sum();
    double Etwo = e_vxc.energy();
    double exx = 0.0;

    incremental_fock.Start(this_iter, conv_accelerator_.getDIIsError());
    incremental_fock.resetMatrices(J, K, Dmat);
    incremental_fock.UpdateCriteria(conv_accelerator_.getDIIsError(),
                                    this_iter);

    double integral_error =
        std::min(conv_accelerator_.getDIIsError() * 1e-5, 1e-5);
    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> both = CalcERIs_EXX(
          MOs.eigenvectors(), incremental_fock.getDmat_diff(), integral_error);
      J += both[0];
      H += J;
      Etwo += 0.5 * Dmat.cwiseProduct(J).sum();
      K += both[1];
      H += 0.5 * ScaHFX_ * K;
      exx = 0.25 * ScaHFX_ * Dmat.cwiseProduct(K).sum();
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " Filled F+K matrix " << std::flush;
    } else {
      J += CalcERIs(incremental_fock.getDmat_diff(), integral_error);
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " Filled F matrix " << std::flush;
      H += J;
      Etwo += 0.5 * Dmat.cwiseProduct(J).sum();
    }

    Etwo += exx;
    double totenergy = Eone + H0.energy() + Etwo;
    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Single particle energy "
                               << std::setprecision(12) << Eone << std::flush;
    XTP_LOG(Log::info, *pLog_) << TimeStamp() << " Two particle energy "
                               << std::setprecision(12) << Etwo << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << std::setprecision(12) << " Local Exc contribution "
        << e_vxc.energy() << std::flush;
    if (ScaHFX_ > 0) {
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << std::setprecision(12)
          << " Non local Ex contribution " << exx << std::flush;
    }
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Total Energy " << std::setprecision(12) << totenergy
        << std::flush;

    Dmat = conv_accelerator_.Iterate(Dmat, H, MOs, totenergy);
    incremental_fock.UpdateDmats(Dmat, conv_accelerator_.getDIIsError(),
                                 this_iter);

    PrintMOs(MOs.eigenvalues(), Log::info);

    XTP_LOG(Log::info, *pLog_) << "\t\tGAP "
                               << MOs.eigenvalues()(numofelectrons_ / 2) -
                                      MOs.eigenvalues()(numofelectrons_ / 2 - 1)
                               << std::flush;

    if (conv_accelerator_.isConverged()) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Total Energy has converged to "
          << std::setprecision(9) << conv_accelerator_.getDeltaE()
          << "[Ha] after " << this_iter + 1
          << " iterations. DIIS error is converged up to "
          << conv_accelerator_.getDIIsError() << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Final Single Point Energy "
          << std::setprecision(12) << totenergy << " Ha" << std::flush;
      XTP_LOG(Log::error, *pLog_) << TimeStamp() << std::setprecision(12)
                                  << " Final Local Exc contribution "
                                  << e_vxc.energy() << " Ha" << std::flush;
      if (ScaHFX_ > 0) {
        XTP_LOG(Log::error, *pLog_) << TimeStamp() << std::setprecision(12)
                                    << " Final Non Local Ex contribution "
                                    << exx << " Ha" << std::flush;
      }

      PrintMOs(MOs.eigenvalues(), Log::error);
      orb.setQMEnergy(totenergy);
      orb.MOs() = MOs;
      CalcElDipole(orb);
      break;
    } else if (this_iter == max_iter_ - 1) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " DFT calculation has not converged after "
          << max_iter_
          << " iterations. Use more iterations or another convergence "
             "acceleration scheme."
          << std::flush;
      return false;
    }
  }
  return true;
}

bool DFTEngine::EvaluateActiveRegion(Orbitals& orb) {

  // reading in the orbitals of the full DFT calculation
  tools::EigenSystem embeddingMOs;
  embeddingMOs.eigenvalues() = orb.MOs().eigenvalues();
  embeddingMOs.eigenvectors() = orb.MOs().eigenvectors();

  // constructing the full electron density matrix
  const Eigen::MatrixXd FullDensityMatrix = orb.DensityMatrixGroundState();



  // reading the localized orbitals and update the occupied orbitals in MOs
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Passing localized orbitals as the initial guess"
      << std::endl;
  Eigen::MatrixXd LMOs = orb.getPMLocalizedOrbital();

  embeddingMOs.eigenvectors().leftCols(orb.getNumberOfAlphaElectrons()) = LMOs;

  // determine the active and inactive electron densities
  std::vector<Index> activeatoms =
      IndexParser().CreateIndexVector(active_atoms_as_string);

  ActiveDensityMatrix DMAT_A(orb, activeatoms);
  const Eigen::MatrixXd InitialActiveDensityMatrix = DMAT_A.compute_Dmat_A();

        //    std::cout << InitialActiveDensityMatrix  << std::endl;


  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Active density formation done" << std::endl;

  Eigen::MatrixXd InactiveDensityMatrix =
      FullDensityMatrix - InitialActiveDensityMatrix;

  // AOoverlap to calculate number of active/inactive electrons
  AOBasis aobasis = orb.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);

   // std::cout << overlap.Matrix() << std::endl;


  Index all_electrons = static_cast<Index>(std::round(
      FullDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));
  Index active_electrons = static_cast<Index>(std::round(
      InitialActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));
  Index inactive_electrons = static_cast<Index>(std::round(
      InactiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));

  XTP_LOG(Log::error, *pLog_) << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Total electrons: " << all_electrons << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Active electrons: " << active_electrons << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Inactive electrons: " << inactive_electrons
      << std::flush;

  // check for consistency
  if ((active_electrons + inactive_electrons) != all_electrons) {
    XTP_LOG(Log::error, *pLog_) << TimeStamp()
                                << " Sum of active and inactive electrons does "
                                   "not match full number of electrons!"
                                << std::flush;
    return false;
  }

  // setup the DFT engine with the active electrons only
  Prepare(orb, active_electrons);

  // setup one-electron part of the Hamiltonian
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());

  // energy of the one-lectron part of the Hamiltonian
  const double E0_full = FullDensityMatrix.cwiseProduct(H0.matrix()).sum();
  const double E0_initial_active =
      InitialActiveDensityMatrix.cwiseProduct(H0.matrix()).sum();
  const double E_nuc = H0.energy();

  // setup the Vxc integrator
  Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(orb.QMAtoms());
  // get the constant XC contributions for the embedding
  const Mat_p_Energy xc_full = vxcpotential.IntegrateVXC(FullDensityMatrix);
  const Mat_p_Energy xc_initial_active =
      vxcpotential.IntegrateVXC(InitialActiveDensityMatrix);

  // get the constant Coulomb matrices for the embedding
  Eigen::MatrixXd J_full =
      Eigen::MatrixXd::Zero(FullDensityMatrix.rows(), FullDensityMatrix.cols());
  Eigen::MatrixXd J_initial_active =
      Eigen::MatrixXd::Zero(FullDensityMatrix.rows(), FullDensityMatrix.cols());
  Eigen::MatrixXd K_full =
      Eigen::MatrixXd::Zero(FullDensityMatrix.rows(), FullDensityMatrix.cols());
  Eigen::MatrixXd K_initial_active =
      Eigen::MatrixXd::Zero(FullDensityMatrix.rows(), FullDensityMatrix.cols());
  if (ScaHFX_ > 0) {
    std::array<Eigen::MatrixXd, 2> JandK_full =
        CalcERIs_EXX(embeddingMOs.eigenvectors(), FullDensityMatrix, 1e-12);
    std::array<Eigen::MatrixXd, 2> JandK_initial_active =
        CalcERIs_EXX(embeddingMOs.eigenvectors(), InitialActiveDensityMatrix, 1e-12);
    J_full = JandK_full[0];
    K_full = 0.5 * ScaHFX_ * JandK_full[1];
    J_initial_active = JandK_initial_active[0];
    K_initial_active = 0.5 * ScaHFX_ * JandK_initial_active[1];
  } else {
    J_full = CalcERIs(FullDensityMatrix, 1e-12);
    J_initial_active = CalcERIs(InitialActiveDensityMatrix, 1e-12);
  }
  // get the Hartree energies
  const double E_Hartree_full =
      0.5 * FullDensityMatrix.cwiseProduct(J_full + K_full).sum();
  const double E_Hartree_initial_active =
      0.5 * InitialActiveDensityMatrix
                .cwiseProduct(J_initial_active + K_initial_active)
                .sum();

  // Energy of the full reference system
  const double Total_E_full =
      E0_full + E_Hartree_full + xc_full.energy() + E_nuc;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " Reference total energy full electron system: " << Total_E_full
      << " Ha" << std::flush;

  // projection parameter, to be made an option
  double mu = 10000;
  const Eigen::MatrixXd Level_Shift_Operator =
      mu * overlap.Matrix() * InactiveDensityMatrix * overlap.Matrix();

  // XC and Hartree contribution to the external embedding potential/energy
  const Eigen::MatrixXd v_embedding = J_full + K_full + xc_full.matrix() -
                                      J_initial_active - K_initial_active -
                                      xc_initial_active.matrix();

  const double constant_embedding_energy = Total_E_full - E0_initial_active -
                                           E_Hartree_initial_active -
                                           xc_initial_active.energy();

  XTP_LOG(Log::error, *pLog_) << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Constant energy embedding terms: " << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " 1) E_DFT[full]: " << Total_E_full << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " 2) E_0[initial active]: " << E0_initial_active
      << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " 3) E_H[initial active]: " << E_Hartree_initial_active
      << std::flush;
  XTP_LOG(Log::info, *pLog_) << TimeStamp() << " 4) E_XC[initial active]: "
                             << xc_initial_active.energy() << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Total 1 -2 -3 -4: " << constant_embedding_energy
      << std::flush;

  // Self-consistent calculation for the active electrons in the embedding
  // potential
  Eigen::MatrixXd ActiveDensityMatrix = InitialActiveDensityMatrix;
  for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Iteration " << this_iter + 1
                                << " of " << max_iter_ << std::flush;

    // get the XC contribution for the updated active density matrix
    Mat_p_Energy xc_active = vxcpotential.IntegrateVXC(ActiveDensityMatrix);
    double E_xc_active = xc_active.energy();

    // get the Hartree contribution for the updated active density matrix
    Eigen::MatrixXd J_active = Eigen::MatrixXd::Zero(
        ActiveDensityMatrix.rows(), ActiveDensityMatrix.cols());
    Eigen::MatrixXd K_active = Eigen::MatrixXd::Zero(
        ActiveDensityMatrix.rows(), ActiveDensityMatrix.cols());
    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> JandK_active =
          CalcERIs_EXX(embeddingMOs.eigenvectors(), ActiveDensityMatrix, 1e-12);
      J_active = JandK_active[0];
      K_active = 0.5 * ScaHFX_ * JandK_active[1];

    } else {
      J_active = CalcERIs(ActiveDensityMatrix, 1e-12);
    }
    double E_Hartree_active =
        0.5 * ActiveDensityMatrix.cwiseProduct(J_active + K_active).sum();

    // update the active Hamiltonian
    Eigen::MatrixXd H_active = H0.matrix() + v_embedding +
                               Level_Shift_Operator + xc_active.matrix() +
                               J_active + K_active;

    const double E0_active =
        ActiveDensityMatrix.cwiseProduct(H0.matrix()).sum();
    const double E_embedding =
        (ActiveDensityMatrix - InitialActiveDensityMatrix)
            .cwiseProduct(v_embedding)
            .sum();
    const double E_levelshift =
        ActiveDensityMatrix.cwiseProduct(Level_Shift_Operator).sum();

    double TotalEnergy = E0_active + E_Hartree_active + E_xc_active +
                         constant_embedding_energy + E_levelshift + E_embedding;

    XTP_LOG(Log::error, *pLog_) << std::flush;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Total Energy after embedding: " << TotalEnergy
        << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "   I) E_0[active]: " << E0_active << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "  II) E_H[active]: " << E_Hartree_active
        << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " III) E_xc[active]: " << E_xc_active << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "  IV) E_ext_emb :" << constant_embedding_energy
        << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "   V) E_corr_emb :" << E_embedding << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "  VI) E_levelshift :" << E_levelshift << std::flush;
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << "  Total I + II + III + IV + V +VI :" << TotalEnergy
        << std::flush;

    // get a new density matrix in the active region
    std::cout << "\n"
              << "Active electrons in = "
              << ActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
              << std::endl;
    ActiveDensityMatrix = conv_accelerator_.Iterate(ActiveDensityMatrix,
                                                    H_active, embeddingMOs, TotalEnergy);
    std::cout << "\n"
              << "Active electrons out= "
              << ActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
              << std::endl;

    PrintMOs(embeddingMOs.eigenvalues(), Log::info);

    if (conv_accelerator_.isConverged()) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Total embedding energy has converged to "
          << std::setprecision(9) << conv_accelerator_.getDeltaE()
          << "[Ha] after " << this_iter + 1
          << " iterations. DIIS error is converged up to "
          << conv_accelerator_.getDIIsError() << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Final Single Point Energy of embedding "
          << std::setprecision(12) << TotalEnergy << " Ha" << std::flush;
      PrintMOs(embeddingMOs.eigenvalues(), Log::error);
      orb.setEmbeddedMOs(embeddingMOs);
      CalcElDipole(orb);
      break;
    } else if (this_iter == max_iter_ - 1) {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " DFT calculation for active region has not converged after "
          << max_iter_
          << " iterations. Use more iterations or another convergence "
             "acceleration scheme."
          << std::flush;
      return false;
    }
  }

  return true;
}  // namespace xtp

Mat_p_Energy DFTEngine::SetupH0(const QMMolecule& mol) const {

  AOKinetic dftAOkinetic;

  dftAOkinetic.Fill(dftbasis_);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Kinetic energy matrix ." << std::flush;

  AOMultipole dftAOESP;
  dftAOESP.FillPotential(dftbasis_, mol);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT nuclear potential matrix." << std::flush;

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Constructed independent particle hamiltonian "
      << std::flush;
  double E0 = NuclearRepulsion(mol);
  XTP_LOG(Log::error, *pLog_) << TimeStamp() << " Nuclear Repulsion Energy is "
                              << std::setprecision(9) << E0 << std::flush;

  if (!ecp_name_.empty()) {
    AOECP dftAOECP;
    dftAOECP.FillPotential(dftbasis_, ecp_);
    H0 += dftAOECP.Matrix();
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled DFT ECP matrix" << std::flush;
  }

  if (externalsites_ != nullptr) {
    XTP_LOG(Log::error, *pLog_) << TimeStamp() << " " << externalsites_->size()
                                << " External sites" << std::flush;
    bool has_quadrupoles = std::any_of(
        externalsites_->begin(), externalsites_->end(),
        [](const std::unique_ptr<StaticSite>& s) { return s->getRank() == 2; });
    std::string header =
        " Name      Coordinates[a0]     charge[e]         dipole[e*a0]    ";
    if (has_quadrupoles) {
      header += "              quadrupole[e*a0^2]";
    }
    XTP_LOG(Log::error, *pLog_) << header << std::flush;
    Index limit = 50;
    Index counter = 0;
    for (const std::unique_ptr<StaticSite>& site : *externalsites_) {
      if (counter == limit) {
        break;
      }
      std::string output =
          (boost::format("  %1$s"
                         "   %2$+1.4f %3$+1.4f %4$+1.4f"
                         "   %5$+1.4f") %
           site->getElement() % site->getPos()[0] % site->getPos()[1] %
           site->getPos()[2] % site->getCharge())
              .str();
      const Eigen::Vector3d& dipole = site->getDipole();
      output += (boost::format("   %1$+1.4f %2$+1.4f %3$+1.4f") % dipole[0] %
                 dipole[1] % dipole[2])
                    .str();
      if (site->getRank() > 1) {
        Eigen::VectorXd quadrupole = site->Q().tail<5>();
        output +=
            (boost::format("   %1$+1.4f %2$+1.4f %3$+1.4f %4$+1.4f %5$+1.4f") %
             quadrupole[0] % quadrupole[1] % quadrupole[2] % quadrupole[3] %
             quadrupole[4])
                .str();
      }
      XTP_LOG(Log::error, *pLog_) << output << std::flush;
      counter++;
    }
    if (counter == limit) {
      XTP_LOG(Log::error, *pLog_)
          << "              ... (" << externalsites_->size() - limit
          << " sites not displayed)\n"
          << std::flush;
    }

    Mat_p_Energy ext_multipoles =
        IntegrateExternalMultipoles(mol, *externalsites_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Nuclei-external site interaction energy "
        << std::setprecision(9) << ext_multipoles.energy() << std::flush;
    E0 += ext_multipoles.energy();
    H0 += ext_multipoles.matrix();
  }

  if (integrate_ext_density_) {
    Orbitals extdensity;
    extdensity.ReadFromCpt(orbfilename_);
    Mat_p_Energy extdensity_result = IntegrateExternalDensity(mol, extdensity);
    E0 += extdensity_result.energy();
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Nuclei-external density interaction energy "
        << std::setprecision(9) << extdensity_result.energy() << std::flush;
    H0 += extdensity_result.matrix();
  }

  if (integrate_ext_field_) {

    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Integrating external electric field with F[Hrt]="
        << extfield_.transpose() << std::flush;
    H0 += IntegrateExternalField(mol);
  }

  return Mat_p_Energy(E0, H0);
}

void DFTEngine::SetupInvariantMatrices() {

  dftAOoverlap_.Fill(dftbasis_);

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Filled DFT Overlap matrix." << std::flush;

  conv_opt_.numberofelectrons = numofelectrons_;
  conv_accelerator_.Configure(conv_opt_);
  conv_accelerator_.setLogger(pLog_);
  conv_accelerator_.setOverlap(dftAOoverlap_, 1e-8);
  conv_accelerator_.PrintConfigOptions();

  if (!auxbasis_name_.empty()) {
    // prepare invariant part of electron repulsion integrals
    ERIs_.Initialize(dftbasis_, auxbasis_);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Inverted AUX Coulomb matrix, removed "
        << ERIs_.Removedfunctions() << " functions from aux basis"
        << std::flush;
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Setup invariant parts of Electron Repulsion integrals "
        << std::flush;
  } else {
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculating 4c diagonals. " << std::flush;
    ERIs_.Initialize_4c(dftbasis_);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Calculated 4c diagonals. " << std::flush;
  }

  return;
}

Eigen::MatrixXd DFTEngine::RunAtomicDFT_unrestricted(
    const QMAtom& uniqueAtom) const {
  bool with_ecp = !ecp_name_.empty();
  if (uniqueAtom.getElement() == "H" || uniqueAtom.getElement() == "He") {
    with_ecp = false;
  }

  QMMolecule atom = QMMolecule("individual_atom", 0);
  atom.push_back(uniqueAtom);

  BasisSet basisset;
  basisset.Load(dftbasis_name_);
  AOBasis dftbasis;
  dftbasis.Fill(basisset, atom);
  Vxc_Grid grid;
  grid.GridSetup(grid_name_, atom, dftbasis);
  Vxc_Potential<Vxc_Grid> gridIntegration(grid);
  gridIntegration.setXCfunctional(xc_functional_name_);

  ECPAOBasis ecp;
  if (with_ecp) {
    ECPBasisSet ecps;
    ecps.Load(ecp_name_);
    ecp.Fill(ecps, atom);
  }

  Index numofelectrons = uniqueAtom.getNuccharge();
  Index alpha_e = 0;
  Index beta_e = 0;

  if ((numofelectrons % 2) != 0) {
    alpha_e = numofelectrons / 2 + numofelectrons % 2;
    beta_e = numofelectrons / 2;
  } else {
    alpha_e = numofelectrons / 2;
    beta_e = alpha_e;
  }

  AOOverlap dftAOoverlap;
  AOKinetic dftAOkinetic;
  AOMultipole dftAOESP;
  AOECP dftAOECP;
  ERIs ERIs_atom;

  // DFT AOOverlap matrix

  dftAOoverlap.Fill(dftbasis);
  dftAOkinetic.Fill(dftbasis);

  dftAOESP.FillPotential(dftbasis, atom);
  ERIs_atom.Initialize_4c(dftbasis);

  ConvergenceAcc Convergence_alpha;
  ConvergenceAcc Convergence_beta;
  ConvergenceAcc::options opt_alpha = conv_opt_;
  opt_alpha.mode = ConvergenceAcc::KSmode::open;
  opt_alpha.histlength = 20;
  opt_alpha.levelshift = 0.1;
  opt_alpha.levelshiftend = 0.0;
  opt_alpha.usediis = true;
  opt_alpha.adiis_start = 0.0;
  opt_alpha.diis_start = 0.0;
  opt_alpha.numberofelectrons = alpha_e;

  ConvergenceAcc::options opt_beta = opt_alpha;
  opt_beta.numberofelectrons = beta_e;

  Logger log;
  Convergence_alpha.Configure(opt_alpha);
  Convergence_alpha.setLogger(&log);
  Convergence_alpha.setOverlap(dftAOoverlap, 1e-8);
  Convergence_beta.Configure(opt_beta);
  Convergence_beta.setLogger(&log);
  Convergence_beta.setOverlap(dftAOoverlap, 1e-8);
  /**** Construct initial density  ****/

  Eigen::MatrixXd H0 = dftAOkinetic.Matrix() + dftAOESP.Matrix();
  if (with_ecp) {
    dftAOECP.FillPotential(dftbasis, ecp);
    H0 += dftAOECP.Matrix();
  }
  tools::EigenSystem MOs_alpha = Convergence_alpha.SolveFockmatrix(H0);

  Eigen::MatrixXd dftAOdmat_alpha = Convergence_alpha.DensityMatrix(MOs_alpha);
  if (uniqueAtom.getElement() == "H") {
    return dftAOdmat_alpha;
  }
  tools::EigenSystem MOs_beta = Convergence_beta.SolveFockmatrix(H0);
  Eigen::MatrixXd dftAOdmat_beta = Convergence_beta.DensityMatrix(MOs_beta);

  Index maxiter = 80;
  for (Index this_iter = 0; this_iter < maxiter; this_iter++) {

    Eigen::MatrixXd H_alpha = H0;
    Eigen::MatrixXd H_beta = H_alpha;

    double E_two_alpha = 0.0;
    double E_two_beta = 0.0;

    double integral_error = std::min(1e-5 * 0.5 *
                                         (Convergence_alpha.getDIIsError() +
                                          Convergence_beta.getDIIsError()),
                                     1e-5);
    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> both_alpha =
          ERIs_atom.CalculateERIs_EXX_4c(dftAOdmat_alpha, integral_error);

      std::array<Eigen::MatrixXd, 2> both_beta =
          ERIs_atom.CalculateERIs_EXX_4c(dftAOdmat_beta, integral_error);
      Eigen::MatrixXd Hartree = both_alpha[0] + both_beta[0];
      E_two_alpha += Hartree.cwiseProduct(dftAOdmat_alpha).sum();
      E_two_beta += Hartree.cwiseProduct(dftAOdmat_beta).sum();
      E_two_alpha += 0.5 * both_alpha[1].cwiseProduct(dftAOdmat_alpha).sum();
      E_two_beta += 0.5 * both_beta[1].cwiseProduct(dftAOdmat_beta).sum();
      H_alpha += Hartree + ScaHFX_ * both_alpha[1];
      H_beta += Hartree + ScaHFX_ * both_beta[1];

    } else {
      Eigen::MatrixXd Hartree = ERIs_atom.CalculateERIs_4c(
          dftAOdmat_alpha + dftAOdmat_beta, integral_error);
      E_two_alpha += Hartree.cwiseProduct(dftAOdmat_alpha).sum();
      E_two_beta += Hartree.cwiseProduct(dftAOdmat_beta).sum();
      H_alpha += Hartree;
      H_beta += Hartree;
    }

    Mat_p_Energy e_vxc_alpha = gridIntegration.IntegrateVXC(dftAOdmat_alpha);
    H_alpha += e_vxc_alpha.matrix();
    E_two_alpha += e_vxc_alpha.energy();

    Mat_p_Energy e_vxc_beta = gridIntegration.IntegrateVXC(dftAOdmat_beta);
    H_beta += e_vxc_beta.matrix();
    E_two_beta += e_vxc_beta.energy();

    double E_one_alpha = dftAOdmat_alpha.cwiseProduct(H0).sum();
    double E_one_beta = dftAOdmat_beta.cwiseProduct(H0).sum();
    double E_alpha = E_one_alpha + E_two_alpha;
    double E_beta = E_one_beta + E_two_beta;
    double totenergy = E_alpha + E_beta;
    // evolve alpha
    dftAOdmat_alpha =
        Convergence_alpha.Iterate(dftAOdmat_alpha, H_alpha, MOs_alpha, E_alpha);
    // evolve beta
    dftAOdmat_beta =
        Convergence_beta.Iterate(dftAOdmat_beta, H_beta, MOs_beta, E_beta);

    XTP_LOG(Log::debug, *pLog_)
        << TimeStamp() << " Iter " << this_iter << " of " << maxiter << " Etot "
        << totenergy << " diise_a " << Convergence_alpha.getDIIsError()
        << " diise_b " << Convergence_beta.getDIIsError() << "\n\t\t a_gap "
        << MOs_alpha.eigenvalues()(alpha_e) -
               MOs_alpha.eigenvalues()(alpha_e - 1)
        << " b_gap "
        << MOs_beta.eigenvalues()(beta_e) - MOs_beta.eigenvalues()(beta_e - 1)
        << " Nalpha="
        << dftAOoverlap.Matrix().cwiseProduct(dftAOdmat_alpha).sum()
        << " Nbeta=" << dftAOoverlap.Matrix().cwiseProduct(dftAOdmat_beta).sum()
        << std::flush;

    bool converged =
        Convergence_alpha.isConverged() && Convergence_beta.isConverged();
    if (converged || this_iter == maxiter - 1) {

      if (converged) {
        XTP_LOG(Log::info, *pLog_)
            << TimeStamp() << " Converged after " << this_iter + 1
            << " iterations" << std::flush;
      } else {
        XTP_LOG(Log::info, *pLog_)
            << TimeStamp() << " Not converged after " << this_iter + 1
            << " iterations. Unconverged density.\n\t\t\t"
            << " DIIsError_alpha=" << Convergence_alpha.getDIIsError()
            << " DIIsError_beta=" << Convergence_beta.getDIIsError()
            << std::flush;
      }
      break;
    }
  }
  Eigen::MatrixXd avgmatrix =
      SphericalAverageShells(dftAOdmat_alpha + dftAOdmat_beta, dftbasis);
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Atomic density Matrix for " << uniqueAtom.getElement()
      << " gives N=" << std::setprecision(9)
      << avgmatrix.cwiseProduct(dftAOoverlap.Matrix()).sum() << " electrons."
      << std::flush;
  return avgmatrix;
}

Eigen::MatrixXd DFTEngine::AtomicGuess(const QMMolecule& mol) const {

  std::vector<std::string> elements = mol.FindUniqueElements();
  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Scanning molecule of size " << mol.size()
      << " for unique elements" << std::flush;
  QMMolecule uniqueelements = QMMolecule("uniqueelements", 0);
  for (auto element : elements) {
    uniqueelements.push_back(QMAtom(0, element, Eigen::Vector3d::Zero()));
  }

  XTP_LOG(Log::info, *pLog_) << TimeStamp() << " " << uniqueelements.size()
                             << " unique elements found" << std::flush;
  std::vector<Eigen::MatrixXd> uniqueatom_guesses;
  for (QMAtom& unique_atom : uniqueelements) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Calculating atom density for "
        << unique_atom.getElement() << std::flush;
    Eigen::MatrixXd dmat_unrestricted = RunAtomicDFT_unrestricted(unique_atom);
    uniqueatom_guesses.push_back(dmat_unrestricted);
  }

  Eigen::MatrixXd guess =
      Eigen::MatrixXd::Zero(dftbasis_.AOBasisSize(), dftbasis_.AOBasisSize());
  Index start = 0;
  for (const QMAtom& atom : mol) {
    Index index = 0;
    for (; index < uniqueelements.size(); index++) {
      if (atom.getElement() == uniqueelements[index].getElement()) {
        break;
      }
    }
    Eigen::MatrixXd& dmat_unrestricted = uniqueatom_guesses[index];
    guess.block(start, start, dmat_unrestricted.rows(),
                dmat_unrestricted.cols()) = dmat_unrestricted;
    start += dmat_unrestricted.rows();
  }

  return guess;
}

void DFTEngine::ConfigOrbfile(Orbitals& orb) {
  if (initial_guess_ == "orbfile") {

    if (orb.hasDFTbasisName()) {
      if (orb.getDFTbasisName() != dftbasis_name_) {
        throw std::runtime_error(
            (boost::format("Basisset Name in guess orb file "
                           "and in dftengine option file differ %1% vs %2%") %
             orb.getDFTbasisName() % dftbasis_name_)
                .str());
      }
    } else {
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " WARNING: "
             "Orbital file has no basisset information,"
             "using it as a guess might work or not for calculation with "
          << dftbasis_name_ << std::flush;
    }
  }
  orb.setXCFunctionalName(xc_functional_name_);
  orb.setXCGrid(grid_name_);
  orb.setScaHFX(ScaHFX_);
  if (!ecp_name_.empty()) {
    orb.setECPName(ecp_name_);
  }
  if (!auxbasis_name_.empty()) {
    orb.SetupAuxBasis(auxbasis_name_);
  }

  if (initial_guess_ == "orbfile") {
    if (orb.hasECPName() || !ecp_name_.empty()) {
      if (orb.getECPName() != ecp_name_) {
        throw std::runtime_error(
            (boost::format("ECPs in orb file: %1% and options %2% differ") %
             orb.getECPName() % ecp_name_)
                .str());
      }
    }
    if (orb.getNumberOfAlphaElectrons() != numofelectrons_ / 2) {
      throw std::runtime_error(
          (boost::format("Number of electron in guess orb file: %1% and in "
                         "dftengine: %2% differ.") %
           orb.getNumberOfAlphaElectrons() % (numofelectrons_ / 2))
              .str());
    }
    if (orb.getBasisSetSize() != dftbasis_.AOBasisSize()) {
      throw std::runtime_error(
          (boost::format("Number of levels in guess orb file: "
                         "%1% and in dftengine: %2% differ.") %
           orb.getBasisSetSize() % dftbasis_.AOBasisSize())
              .str());
    }
  } else {
    orb.setNumberOfAlphaElectrons(numofelectrons_ / 2);
    orb.setNumberOfOccupiedLevels(numofelectrons_ / 2);
  }
  return;
}

void DFTEngine::Prepare(Orbitals& orb, Index numofelectrons) {
  QMMolecule& mol = orb.QMAtoms();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Using " << OPENMP::getMaxThreads() << " threads"
      << std::flush;

  if (XTP_HAS_MKL_OVERLOAD()) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using MKL overload for Eigen " << std::flush;
  } else {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp()
        << " Using native Eigen implementation, no BLAS overload "
        << std::flush;
  }

  XTP_LOG(Log::error, *pLog_) << " Molecule Coordinates [A] " << std::flush;
  for (const QMAtom& atom : mol) {
    const Eigen::Vector3d pos = atom.getPos() * tools::conv::bohr2ang;
    std::string output = (boost::format("  %1$s"
                                        "   %2$+1.4f %3$+1.4f %4$+1.4f") %
                          atom.getElement() % pos[0] % pos[1] % pos[2])
                             .str();

    XTP_LOG(Log::error, *pLog_) << output << std::flush;
  }

  orb.SetupDftBasis(dftbasis_name_);
  dftbasis_ = orb.getDftBasis();

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Loaded DFT Basis Set " << dftbasis_name_ << " with "
      << dftbasis_.AOBasisSize() << " functions" << std::flush;

  if (!auxbasis_name_.empty()) {
    BasisSet auxbasisset;
    auxbasisset.Load(auxbasis_name_);
    auxbasis_.Fill(auxbasisset, mol);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Loaded AUX Basis Set " << auxbasis_name_ << " with "
        << auxbasis_.AOBasisSize() << " functions" << std::flush;
  }
  if (!ecp_name_.empty()) {
    ECPBasisSet ecpbasisset;
    ecpbasisset.Load(ecp_name_);
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Loaded ECP library " << ecp_name_ << std::flush;

    std::vector<std::string> results = ecp_.Fill(ecpbasisset, mol);
    XTP_LOG(Log::info, *pLog_)
        << TimeStamp() << " Filled ECP Basis" << std::flush;
    if (results.size() > 0) {
      std::string message = "";
      for (const std::string& element : results) {
        message += " " + element;
      }
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Found no ECPs for elements" << message
          << std::flush;
    }
  }

  if (numofelectrons < 0) {
    for (const QMAtom& atom : mol) {
      numofelectrons_ += atom.getNuccharge();
    }
  } else {
    numofelectrons_ = numofelectrons;
  }

  // here number of electrons is actually the total number, everywhere else in
  // votca it is just alpha_electrons
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Total number of electrons: " << numofelectrons_
      << std::flush;

  SetupInvariantMatrices();
  return;
}

Vxc_Potential<Vxc_Grid> DFTEngine::SetupVxc(const QMMolecule& mol) {
  ScaHFX_ = Vxc_Potential<Vxc_Grid>::getExactExchange(xc_functional_name_);
  if (ScaHFX_ > 0) {
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " Using hybrid functional with alpha=" << ScaHFX_
        << std::flush;
  }
  Vxc_Grid grid;
  grid.GridSetup(grid_name_, mol, dftbasis_);
  Vxc_Potential<Vxc_Grid> vxc(grid);
  vxc.setXCfunctional(xc_functional_name_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Setup numerical integration grid " << grid_name_
      << " for vxc functional " << xc_functional_name_ << std::flush;
  XTP_LOG(Log::info, *pLog_)
      << "\t\t "
      << " with " << grid.getGridSize() << " points"
      << " divided into " << grid.getBoxesSize() << " boxes" << std::flush;
  return vxc;
}

double DFTEngine::NuclearRepulsion(const QMMolecule& mol) const {
  double E_nucnuc = 0.0;

  for (Index i = 0; i < mol.size(); i++) {
    const Eigen::Vector3d& r1 = mol[i].getPos();
    double charge1 = double(mol[i].getNuccharge());
    for (Index j = 0; j < i; j++) {
      const Eigen::Vector3d& r2 = mol[j].getPos();
      double charge2 = double(mol[j].getNuccharge());
      E_nucnuc += charge1 * charge2 / (r1 - r2).norm();
    }
  }
  return E_nucnuc;
}

// spherically average the density matrix belonging to two shells
Eigen::MatrixXd DFTEngine::SphericalAverageShells(
    const Eigen::MatrixXd& dmat, const AOBasis& dftbasis) const {
  Eigen::MatrixXd avdmat = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  for (const AOShell& shellrow : dftbasis) {
    Index size_row = shellrow.getNumFunc();
    Index start_row = shellrow.getStartIndex();
    for (const AOShell& shellcol : dftbasis) {
      Index size_col = shellcol.getNumFunc();
      Index start_col = shellcol.getStartIndex();
      Eigen::MatrixXd shelldmat =
          dmat.block(start_row, start_col, size_row, size_col);
      if (shellrow.getL() == shellcol.getL()) {
        double diagavg = shelldmat.diagonal().sum() / double(shelldmat.rows());
        Index offdiagelements =
            shelldmat.rows() * shelldmat.cols() - shelldmat.cols();
        double offdiagavg = (shelldmat.sum() - shelldmat.diagonal().sum()) /
                            double(offdiagelements);
        avdmat.block(start_row, start_col, size_row, size_col).array() =
            offdiagavg;
        avdmat.block(start_row, start_col, size_row, size_col)
            .diagonal()
            .array() = diagavg;
      } else {
        double avg = shelldmat.sum() / double(shelldmat.size());
        avdmat.block(start_row, start_col, size_row, size_col).array() = avg;
      }
    }
  }
  return avdmat;
}

double DFTEngine::ExternalRepulsion(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  if (multipoles.size() == 0) {
    return 0;
  }

  double E_ext = 0;
  eeInteractor interactor;
  for (const QMAtom& atom : mol) {
    StaticSite nucleus = StaticSite(atom, double(atom.getNuccharge()));
    for (const std::unique_ptr<StaticSite>& site : *externalsites_) {
      if ((site->getPos() - nucleus.getPos()).norm() < 1e-7) {
        XTP_LOG(Log::error, *pLog_) << TimeStamp()
                                    << " External site sits on nucleus, "
                                       "interaction between them is ignored."
                                    << std::flush;
        continue;
      }
      E_ext += interactor.CalcStaticEnergy_site(*site, nucleus);
    }
  }
  return E_ext;
}

Eigen::MatrixXd DFTEngine::IntegrateExternalField(const QMMolecule& mol) const {

  AODipole dipole;
  dipole.setCenter(mol.getPos());
  dipole.Fill(dftbasis_);
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(dipole.Dimension(), dipole.Dimension());
  for (Index i = 0; i < 3; i++) {
    result -= dipole.Matrix()[i] * extfield_[i];
  }
  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalMultipoles(
    const QMMolecule& mol,
    const std::vector<std::unique_ptr<StaticSite> >& multipoles) const {

  Mat_p_Energy result(dftbasis_.AOBasisSize(), dftbasis_.AOBasisSize());
  AOMultipole dftAOESP;

  dftAOESP.FillPotential(dftbasis_, multipoles);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Filled DFT external multipole potential matrix"
      << std::flush;
  result.matrix() = dftAOESP.Matrix();
  result.energy() = ExternalRepulsion(mol, multipoles);

  return result;
}

Mat_p_Energy DFTEngine::IntegrateExternalDensity(
    const QMMolecule& mol, const Orbitals& extdensity) const {
  BasisSet basis;
  basis.Load(extdensity.getDFTbasisName());
  AOBasis aobasis;
  aobasis.Fill(basis, extdensity.QMAtoms());
  Vxc_Grid grid;
  grid.GridSetup(gridquality_, extdensity.QMAtoms(), aobasis);
  DensityIntegration<Vxc_Grid> numint(grid);
  Eigen::MatrixXd dmat = extdensity.DensityMatrixFull(state_);

  numint.IntegrateDensity(dmat);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated external density" << std::flush;
  Eigen::MatrixXd e_contrib = numint.IntegratePotential(dftbasis_);
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated potential from electron density"
      << std::flush;
  AOMultipole esp;
  esp.FillPotential(dftbasis_, extdensity.QMAtoms());

  double nuc_energy = 0.0;
  for (const QMAtom& atom : mol) {
    nuc_energy +=
        numint.IntegratePotential(atom.getPos()) * double(atom.getNuccharge());
    for (const QMAtom& extatom : extdensity.QMAtoms()) {
      const double dist = (atom.getPos() - extatom.getPos()).norm();
      nuc_energy +=
          double(atom.getNuccharge()) * double(extatom.getNuccharge()) / dist;
    }
  }
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Calculated potential from nuclei" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Electrostatic: " << nuc_energy << std::flush;
  return Mat_p_Energy(nuc_energy, e_contrib + esp.Matrix());
}

Eigen::MatrixXd DFTEngine::OrthogonalizeGuess(
    const Eigen::MatrixXd& GuessMOs) const {
  Eigen::MatrixXd nonortho =
      GuessMOs.transpose() * dftAOoverlap_.Matrix() * GuessMOs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(nonortho);
  Eigen::MatrixXd result = GuessMOs * es.operatorInverseSqrt();
  return result;
}

}  // namespace xtp
}  // namespace votca
