/*
 *            Copyright 2009-2023 The VOTCA Development Team
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
/*
 *References- (1) A Simple, Exact Density-Functional-Theory
 *Embedding Scheme. Frederick R. Manby, Martina Stella, Jason D. Goodpaster, and
 *Thomas F. Miller Journal of Chemical Theory and Computation 2012 8 (8),
 *2564-2568 DOI: 10.1021/ct300544e
 *(2) Taylor A. Barnes, Jason D. Goodpaster, Frederick R. Manby, and Thomas F.
 *Miller III , "Accurate basis set truncation for wavefunction embedding", J.
 *Chem. Phys. 139, 024103 (2013) https://doi.org/10.1063/1.4811112
 *(3) Simon J. Bennie, Martina Stella, Thomas F. Miller III, and Frederick R.
 *Manby , "Accelerating wavefunction in density-functional-theory embedding by
 *truncating the active basis set", J. Chem. Phys. 143, 024105 (2015)
 *https://doi.org/10.1063/1.4923367
 *(3) Projection-Based Wavefunction-in-DFT Embedding
 *Sebastian J. R. Lee, Matthew Welborn, Frederick R. Manby, and Thomas F. Miller
 *III Accounts of Chemical Research 2019 52 (5), 1359-1368
 *DOI: 10.1021/acs.accounts.8b00672
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
#include "votca/xtp/pmlocalization.h"

namespace votca {
namespace xtp {

bool DFTEngine::EvaluateActiveRegion(Orbitals& orb) {

  // reading in the orbitals of the full DFT calculation
  tools::EigenSystem embeddingMOs = orb.MOs();

  // constructing the full electron density matrix
  const Eigen::MatrixXd FullDensityMatrix = orb.DensityMatrixGroundState();
  XTP_LOG(Log::error, *pLog_) << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Starting Projection based dft-in-dft embedding"
      << std::flush;

  // reading the localized orbitals and update the occupied orbitals in MOs
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Passing localized orbitals as the initial guess"
      << std::flush;
  Eigen::MatrixXd LMOs = orb.getLMOs();

  embeddingMOs.eigenvectors().leftCols(orb.getNumberOfAlphaElectrons()) = LMOs;

  // determine the active and inactive electron densities
  std::vector<Index> activeatoms =
      IndexParser().CreateIndexVector(active_atoms_as_string_);

  XTP_LOG(Log::error, *pLog_)
      << "Indices of active atoms selected are: " << active_atoms_as_string_
      << std::flush;
  ActiveDensityMatrix DMAT_A(orb, activeatoms, active_threshold_);
  const Eigen::MatrixXd InitialActiveDensityMatrix = DMAT_A.compute_Dmat_A()[0];
  Eigen::MatrixXd InitialInactiveMOs = DMAT_A.compute_Dmat_A()[2];

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Active density formation done" << std::flush;

  Eigen::MatrixXd InactiveDensityMatrix =
      FullDensityMatrix - InitialActiveDensityMatrix;

  orb.setInactiveDensity(InactiveDensityMatrix);

  // AOoverlap to calculate number of active/inactive electrons
  AOBasis aobasis = orb.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);

  Index all_electrons = static_cast<Index>(
      std::round(FullDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));
  active_electrons_ = static_cast<Index>(std::round(
      InitialActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));
  Index inactive_electrons = static_cast<Index>(
      std::round(InactiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()));

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Total electrons: " << all_electrons << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Active electrons: " << active_electrons_
      << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " Inactive electrons: " << inactive_electrons
      << std::flush;

  // check for consistency
  if ((active_electrons_ + inactive_electrons) != all_electrons) {
    throw std::runtime_error(
        " Sum of active and inactive electrons does "
        "not match full number of electrons!");
    return false;
  }

  // setup the DFT engine with the active electrons only
  Prepare(orb, active_electrons_);

  // setup one-electron part of the Hamiltonian
  Mat_p_Energy H0 = SetupH0(orb.QMAtoms());

  // energy of the one-electron part of the Hamiltonian
  const double E0_full = FullDensityMatrix.cwiseProduct(H0.matrix()).sum();
  const double E0_initial_active =
      InitialActiveDensityMatrix.cwiseProduct(H0.matrix()).sum();
  E_nuc_ = H0.energy();

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
    std::array<Eigen::MatrixXd, 2> JandK_initial_active = CalcERIs_EXX(
        embeddingMOs.eigenvectors(), InitialActiveDensityMatrix, 1e-12);
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
  Total_E_full_ = E0_full + E_Hartree_full + xc_full.energy() + E_nuc_;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " Reference total energy full electron system: " << Total_E_full_
      << " Ha" << std::flush;

  // projection parameter, to be made an option
  Eigen::MatrixXd ProjectionOperator =
      overlap.Matrix() * InactiveDensityMatrix * overlap.Matrix();

  // XC and Hartree contribution to the external embedding potential/energy
  const Eigen::MatrixXd v_embedding = J_full + K_full + xc_full.matrix() -
                                      J_initial_active - K_initial_active -
                                      xc_initial_active.matrix();

  const double constant_embedding_energy = Total_E_full_ - E0_initial_active -
                                           E_Hartree_initial_active -
                                           xc_initial_active.energy();

  XTP_LOG(Log::info, *pLog_)
      << TimeStamp() << " Constant energy embedding terms: " << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp() << " 1) E_DFT[full]: " << Total_E_full_ << std::flush;
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

  if (truncate_) {  // Truncation starts here
    orb.setCalculationType("Truncated");
    TruncateBasis(orb, activeatoms, H0, InitialActiveDensityMatrix, v_embedding,
                  InitialInactiveMOs);
    active_and_border_atoms_ = activeatoms;
  }
  // SCF loop if you don't truncate active region
  else {
    orb.setCalculationType("Embedded_noTrunc");
    Eigen::MatrixXd ActiveDensityMatrix = InitialActiveDensityMatrix;
    for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
      XTP_LOG(Log::error, *pLog_) << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Iteration " << this_iter + 1 << " of "
          << max_iter_ << std::flush;

      // get the XC contribution for the updated active density matrix
      Mat_p_Energy xc_active = vxcpotential.IntegrateVXC(ActiveDensityMatrix);
      double E_xc_active = xc_active.energy();

      // get the Hartree contribution for the updated active density matrix
      Eigen::MatrixXd J_active = Eigen::MatrixXd::Zero(
          ActiveDensityMatrix.rows(), ActiveDensityMatrix.cols());
      Eigen::MatrixXd K_active = Eigen::MatrixXd::Zero(
          ActiveDensityMatrix.rows(), ActiveDensityMatrix.cols());
      if (ScaHFX_ > 0) {
        std::array<Eigen::MatrixXd, 2> JandK_active = CalcERIs_EXX(
            embeddingMOs.eigenvectors(), ActiveDensityMatrix, 1e-12);
        J_active = JandK_active[0];
        K_active = 0.5 * ScaHFX_ * JandK_active[1];

      } else {
        J_active = CalcERIs(ActiveDensityMatrix, 1e-12);
      }
      double E_Hartree_active =
          0.5 * ActiveDensityMatrix.cwiseProduct(J_active + K_active).sum();

      // update the active Hamiltonian
      const Eigen::MatrixXd Level_Shift_Operator =
          levelshift_ * ProjectionOperator;
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
                           constant_embedding_energy + E_levelshift +
                           E_embedding;

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
          << TimeStamp() << "  VI) E_levelshift :" << E_levelshift
          << std::flush;
      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << "  Total I + II + III + IV + V +VI :" << TotalEnergy
          << std::flush;

      // get a new density matrix in the active region
      XTP_LOG(Log::info, *pLog_)
          << std::endl
          << "Active electrons in = "
          << ActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
          << std::flush;
      ActiveDensityMatrix = conv_accelerator_.Iterate(
          ActiveDensityMatrix, H_active, embeddingMOs, TotalEnergy);
      XTP_LOG(Log::info, *pLog_)
          << std::endl
          << "Active electrons out= "
          << ActiveDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
          << std::flush;

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
        XTP_LOG(Log::error, *pLog_)
            << "The difference between the reference energy and the embedding "
               "energy is "
            << Total_E_full_ - TotalEnergy << " Ha" << std::flush;

        if (abs(Total_E_full_ - TotalEnergy) > 1e-3) {
          XTP_LOG(Log::error, *pLog_)
              << "Warning!! The difference is greater than 1e-03 Ha"
              << std::flush;
        }

        PrintMOs(embeddingMOs.eigenvalues(), Log::error);
        orb.setEmbeddedMOs(embeddingMOs);
        orb.setNumofActiveElectrons(active_electrons_);
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
  }
  return true;
}

bool DFTEngine::EvaluateTruncatedActiveRegion(Orbitals& trunc_orb) {
  if (truncate_) {
    activemol_.ReorderAtomIDs();
    trunc_orb.QMAtoms() = activemol_;
    Prepare(trunc_orb, active_electrons_);
    Vxc_Potential<Vxc_Grid> vxcpotential = SetupVxc(trunc_orb.QMAtoms());
    ConfigOrbfile(trunc_orb);
    AOBasis aobasis = trunc_orb.getDftBasis();
    AOOverlap overlap;
    overlap.Fill(aobasis);

    const double E0_initial_truncated =
        InitialActiveDmat_trunc_.cwiseProduct(H0_trunc_).sum();
    const Mat_p_Energy xc_initial_truncated =
        vxcpotential.IntegrateVXC(InitialActiveDmat_trunc_);
    Eigen::MatrixXd J_initial_truncated = Eigen::MatrixXd::Zero(
        InitialActiveDmat_trunc_.rows(), InitialActiveDmat_trunc_.cols());
    Eigen::MatrixXd K_initial_truncated = Eigen::MatrixXd::Zero(
        InitialActiveDmat_trunc_.rows(), InitialActiveDmat_trunc_.cols());
    tools::EigenSystem MOs_trunc;
    MOs_trunc.eigenvalues() = Eigen::VectorXd::Zero(H0_trunc_.cols());
    MOs_trunc.eigenvectors() =
        Eigen::MatrixXd::Zero(H0_trunc_.rows(), H0_trunc_.cols());
    if (ScaHFX_ > 0) {
      std::array<Eigen::MatrixXd, 2> JandK_initial_truncated = CalcERIs_EXX(
          MOs_trunc.eigenvectors(), InitialActiveDmat_trunc_, 1e-12);
      J_initial_truncated = JandK_initial_truncated[0];
      K_initial_truncated = 0.5 * ScaHFX_ * JandK_initial_truncated[1];
    } else {
      J_initial_truncated = CalcERIs(InitialActiveDmat_trunc_, 1e-12);
    }
    // get the Hartree energies
    const double E_Hartree_initial_truncated =
        0.5 * InitialActiveDmat_trunc_
                  .cwiseProduct(J_initial_truncated + K_initial_truncated)
                  .sum();
    const double Initial_truncated_energy = E0_initial_truncated +
                                            E_Hartree_initial_truncated +
                                            xc_initial_truncated.energy();
    XTP_LOG(Log::error, *pLog_)
        << "Initial truncated energy = " << Initial_truncated_energy
        << std::flush;

    Eigen::MatrixXd PurifiedActiveDmat_trunc =
        McWeenyPurification(InitialActiveDmat_trunc_, overlap);
    Eigen::MatrixXd TruncatedDensityMatrix = PurifiedActiveDmat_trunc;

    for (Index this_iter = 0; this_iter < max_iter_; this_iter++) {
      XTP_LOG(Log::error, *pLog_) << std::flush;
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp() << " Iteration " << this_iter + 1 << " of "
          << max_iter_ << std::flush;

      // get the XC contribution for the updated truncated density matrix
      Mat_p_Energy xc_truncated =
          vxcpotential.IntegrateVXC(TruncatedDensityMatrix);
      double E_xc_truncated = xc_truncated.energy();

      // get the Hartree contribution for the updated truncated density matrix
      Eigen::MatrixXd J_truncated = Eigen::MatrixXd::Zero(
          TruncatedDensityMatrix.rows(), TruncatedDensityMatrix.cols());
      Eigen::MatrixXd K_truncated = Eigen::MatrixXd::Zero(
          TruncatedDensityMatrix.rows(), TruncatedDensityMatrix.cols());
      if (ScaHFX_ > 0) {
        std::array<Eigen::MatrixXd, 2> JandK_truncated = CalcERIs_EXX(
            MOs_trunc.eigenvectors(), TruncatedDensityMatrix, 1e-12);
        J_truncated = JandK_truncated[0];
        K_truncated = 0.5 * ScaHFX_ * JandK_truncated[1];

      } else {
        J_truncated = CalcERIs(TruncatedDensityMatrix, 1e-12);
      }
      double E_Hartree_truncated =
          0.5 *
          TruncatedDensityMatrix.cwiseProduct(J_truncated + K_truncated).sum();
      // update the truncated Hamiltonian
      Eigen::MatrixXd H_truncated =
          H0_trunc_ + xc_truncated.matrix() + J_truncated + K_truncated;

      double TruncatedTotalEnergy =
          TruncatedDensityMatrix.cwiseProduct(H0_trunc_).sum() +
          E_Hartree_truncated + E_xc_truncated + Total_E_full_ -
          Initial_truncated_energy +
          ((TruncatedDensityMatrix - InitialActiveDmat_trunc_)
               .cwiseProduct(v_embedding_trunc_)
               .sum());
      // get the new truncated density matrix
      XTP_LOG(Log::info, *pLog_)
          << "Electrons in = "
          << TruncatedDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
          << std::flush;
      TruncatedDensityMatrix = conv_accelerator_.Iterate(
          TruncatedDensityMatrix, H_truncated, MOs_trunc, TruncatedTotalEnergy);
      XTP_LOG(Log::info, *pLog_)
          << "Electrons out = "
          << TruncatedDensityMatrix.cwiseProduct(overlap.Matrix()).sum()
          << std::flush;

      XTP_LOG(Log::info, *pLog_)
          << TimeStamp() << " E_trunc (Ha) = "
          << TruncatedDensityMatrix.cwiseProduct(H0_trunc_).sum() +
                 E_Hartree_truncated + E_xc_truncated
          << "\n \t \t \t"
          << " E_fullDFT (Ha) = " << Total_E_full_ << "\n \t \t \t"
          << " E_initial_trunc(Ha) = " << Initial_truncated_energy
          << "\n \t \t \t"
          << " E_embedding_correction(Ha) = "
          << ((TruncatedDensityMatrix - InitialActiveDmat_trunc_)
                  .cwiseProduct(v_embedding_trunc_)
                  .sum())
          << std::flush;

      XTP_LOG(Log::error, *pLog_)
          << " Truncated Energy of the system: E_trunc + E_fullDFT - "
          << "\n \t \t"
          << " E_initial_trunc + E_embedding_correction = "
          << TruncatedTotalEnergy << std::flush;

      PrintMOs(MOs_trunc.eigenvalues(), Log::info);

      if (conv_accelerator_.isConverged()) {
        XTP_LOG(Log::error, *pLog_)
            << TimeStamp() << " Total embedding energy has converged to "
            << std::setprecision(9) << conv_accelerator_.getDeltaE()
            << "[Ha] after " << this_iter + 1
            << " iterations. DIIS error is converged up to "
            << conv_accelerator_.getDIIsError() << std::flush;
        XTP_LOG(Log::error, *pLog_)
            << TimeStamp() << " Final Single Point Energy of embedding "
            << std::setprecision(12) << TruncatedTotalEnergy << " Ha"
            << std::flush;

        PrintMOs(MOs_trunc.eigenvalues(), Log::error);
        trunc_orb.setEmbeddedMOs(MOs_trunc);
        trunc_orb.setNumofActiveElectrons(active_electrons_);
        break;
      }
    }
    TruncMOsFullBasis(trunc_orb, active_and_border_atoms_, numfuncpatom_);
  }
  return true;
}

void DFTEngine::TruncateBasis(Orbitals& orb, std::vector<Index>& activeatoms,
                              Mat_p_Energy& H0,
                              Eigen::MatrixXd InitialActiveDensityMatrix,
                              Eigen::MatrixXd v_embedding,
                              Eigen::MatrixXd InitialInactiveMOs) {
  AOBasis aobasis;
  aobasis = orb.getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);
  // Mulliken population per basis function on every active atom
  Eigen::VectorXd DiagonalofDmatA = InitialActiveDensityMatrix.diagonal();
  Eigen::VectorXd DiagonalofOverlap = overlap.Matrix().diagonal();
  Eigen::VectorXd MnP = DiagonalofDmatA.cwiseProduct(DiagonalofOverlap);

  // Get a vector containing the number of basis functions per atom
  const std::vector<Index>& numfuncpatom = aobasis.getFuncPerAtom();
  numfuncpatom_ = numfuncpatom;
  Index numofactivebasisfunction = 0;
  // Store start indices. Will be used later
  std::vector<Index> start_indices, start_indices_activemolecule;
  Index start_idx = 0, start_idx_activemolecule = 0;
  std::vector<Index> borderatoms;
  for (Index atom_num = 0; atom_num < orb.QMAtoms().size(); atom_num++) {
    start_indices.push_back(start_idx);
    /* Condition for basis fn on an atom to be counted: either in active
     region or MnP of any function > threshold (border atoms) */
    bool partOfActive =
        (std::find(activeatoms.begin(), activeatoms.end(),
                   orb.QMAtoms()[atom_num].getId()) != activeatoms.end());
    if (partOfActive == true) {
      activemol_.push_back(orb.QMAtoms()[atom_num]);
      /* if active append the atom to the new molecule. Also increment active
       basis function size and start indices */
      start_indices_activemolecule.push_back(start_idx_activemolecule);
      start_idx_activemolecule += numfuncpatom[atom_num];
      numofactivebasisfunction += numfuncpatom[atom_num];
    } else {
      std::vector<const AOShell*> inactiveshell =
          aobasis.getShellsofAtom(orb.QMAtoms()[atom_num].getId());
      bool loop_break = false;
      for (const AOShell* shell : inactiveshell) {
        for (Index shell_fn_no = shell->getStartIndex();
             shell_fn_no < shell->getStartIndex() + shell->getNumFunc();
             shell_fn_no++) {
          if (MnP[shell_fn_no] > truncation_threshold_) {
            activeatoms.push_back(atom_num);
            borderatoms.push_back(atom_num);  // push this index to border
                                              // atoms
            activemol_.push_back(orb.QMAtoms()[atom_num]);
            loop_break = true;  // if any function in the whole atom satisfies
                                // we break loop to check next atom
            start_indices_activemolecule.push_back(start_idx_activemolecule);
            start_idx_activemolecule += numfuncpatom[atom_num];
            numofactivebasisfunction += numfuncpatom[atom_num];
            break;
          }
        }
        if (loop_break) break;
      }
    }
    start_idx += numfuncpatom[atom_num];
  }

  // Sort atoms before you cut the Hamiltonian
  sort(activeatoms.begin(), activeatoms.end());
  sort(borderatoms.begin(), borderatoms.end());
  Eigen::MatrixXd InactiveDensityMatrix = orb.getInactiveDensity();
  XTP_LOG(Log::error, *pLog_) << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << "Active + Border Molecule Size = " << activeatoms.size() << "\n \t \t "
      << "Border Molecule Size = " << borderatoms.size() << std::flush;
  Eigen::MatrixXd ProjectionOperator =
      overlap.Matrix() * InactiveDensityMatrix * overlap.Matrix();

  Eigen::MatrixXd H_embedding =
      H0.matrix() + v_embedding + levelshift_ * ProjectionOperator;

  if (borderatoms.size() != 0) {
    Eigen::MatrixXd BorderMOs;
    for (Index lmo_index = 0; lmo_index < InitialInactiveMOs.cols();
         lmo_index++) {
      double mullikenpop_lmo_borderatoms = 0;
      for (Index borderatom : borderatoms) {
        Index start = start_indices[borderatom];
        Index size = numfuncpatom[borderatom];
        mullikenpop_lmo_borderatoms +=
            (InitialInactiveMOs.col(lmo_index) *
             InitialInactiveMOs.col(lmo_index).transpose() * overlap.Matrix())
                .diagonal()
                .segment(start, size)
                .sum();
      }
      /*If more than half of a MO contributes on a border atom include that in
       * the Border MOs list. Not made an user option. If you want to fix
       * further for density leak you can reduce this number*/
      if (mullikenpop_lmo_borderatoms > 0.25) {
        BorderMOs.conservativeResize(InitialInactiveMOs.rows(),
                                     BorderMOs.cols() + 1);
        BorderMOs.col(BorderMOs.cols() - 1) = InitialInactiveMOs.col(lmo_index);
      }
    }

    Eigen::MatrixXd BorderDmat = 2 * BorderMOs * BorderMOs.transpose();
    Eigen::MatrixXd BorderProjectionOperator =
        overlap.Matrix() * BorderDmat * overlap.Matrix();
    Eigen::MatrixXd DistantProjectionOperator =
        ProjectionOperator - BorderProjectionOperator;
    // 1e+2 is the suitable projection value for border MOs: see citations
    H_embedding +=
        1e+2 * BorderProjectionOperator +
        levelshift_ * (DistantProjectionOperator - ProjectionOperator);
  }

  // from here it is time to truncate Hamiltonian H0
  H0_trunc_ =
      Eigen::MatrixXd::Zero(numofactivebasisfunction, numofactivebasisfunction);
  InitialActiveDmat_trunc_ =
      Eigen::MatrixXd::Zero(numofactivebasisfunction, numofactivebasisfunction);
  v_embedding_trunc_ =
      Eigen::MatrixXd::Zero(numofactivebasisfunction, numofactivebasisfunction);
  // Defined needed matrices with basisset_size * basisset_size

  for (Index activeatom1_idx = 0; activeatom1_idx < Index(activeatoms.size());
       activeatom1_idx++) {
    Index activeatom1 = activeatoms[activeatom1_idx];
    for (Index activeatom2_idx = 0; activeatom2_idx < Index(activeatoms.size());
         activeatom2_idx++) {
      Index activeatom2 = activeatoms[activeatom2_idx];
      // Fill the H0, Dmat and v_embedding at the right indices
      H0_trunc_.block(start_indices_activemolecule[activeatom1_idx],
                      start_indices_activemolecule[activeatom2_idx],
                      numfuncpatom[activeatom1], numfuncpatom[activeatom2]) =
          H_embedding.block(
              start_indices[activeatom1], start_indices[activeatom2],
              numfuncpatom[activeatom1], numfuncpatom[activeatom2]);
      InitialActiveDmat_trunc_.block(
          start_indices_activemolecule[activeatom1_idx],
          start_indices_activemolecule[activeatom2_idx],
          numfuncpatom[activeatom1], numfuncpatom[activeatom2]) =
          InitialActiveDensityMatrix.block(
              start_indices[activeatom1], start_indices[activeatom2],
              numfuncpatom[activeatom1], numfuncpatom[activeatom2]);
      v_embedding_trunc_.block(start_indices_activemolecule[activeatom1_idx],
                               start_indices_activemolecule[activeatom2_idx],
                               numfuncpatom[activeatom1],
                               numfuncpatom[activeatom2]) =
          v_embedding.block(
              start_indices[activeatom1], start_indices[activeatom2],
              numfuncpatom[activeatom1], numfuncpatom[activeatom2]);
    }
  }
}

Eigen::MatrixXd DFTEngine::McWeenyPurification(Eigen::MatrixXd& Dmat_in,
                                               AOOverlap& overlap) {
  Eigen::MatrixXd Ssqrt = overlap.Sqrt();
  Eigen::MatrixXd InvSsqrt = overlap.Pseudo_InvSqrt(1e-8);
  Eigen::MatrixXd ModifiedDmat = 0.5 * Ssqrt * Dmat_in * Ssqrt;
  for (Index iter = 0; iter < 100; iter++) {
    Eigen::MatrixXd Dmat_new = (3 * ModifiedDmat * ModifiedDmat) -
                               (2 * ModifiedDmat * ModifiedDmat * ModifiedDmat);
    double IdempotencyError =
        ((Dmat_new * Dmat_new - Dmat_new) * (Dmat_new * Dmat_new - Dmat_new))
            .trace();
    XTP_LOG(Log::info, *pLog_)
        << "Idempotency Error: " << IdempotencyError << std::flush;

    ModifiedDmat = Dmat_new;
    if (IdempotencyError < 1e-20) break;
  }
  Eigen::MatrixXd Dmat_out = 2 * InvSsqrt * ModifiedDmat * InvSsqrt;
  return 0.5 * (Dmat_out + Dmat_out.transpose());
}

void DFTEngine::TruncMOsFullBasis(Orbitals& orb, std::vector<Index> activeatoms,
                                  std::vector<Index> numfuncpatom) {
  Eigen::MatrixXd expandtruncorb = orb.getEmbeddedMOs().eigenvectors();
  Index start_index = 0;
  Index numofactualactoms = numfuncpatom.size();
  for (Index atomindex = 0; atomindex < numofactualactoms; atomindex++) {
    bool partOfActive = (std::find(activeatoms.begin(), activeatoms.end(),
                                   atomindex) != activeatoms.end());
    if (partOfActive == false) {
      expandtruncorb =
          InsertZeroRows(expandtruncorb, start_index, numfuncpatom[atomindex]);
    }
    start_index += numfuncpatom[atomindex];
  }
  orb.setTruncMOsFullBasis(expandtruncorb);
}

Eigen::MatrixXd DFTEngine::InsertZeroRows(Eigen::MatrixXd MOsMatrix,
                                          Index startidx, Index numofzerorows) {
  Eigen::MatrixXd FinalMatrix =
      Eigen::MatrixXd::Zero(MOsMatrix.rows() + numofzerorows, MOsMatrix.cols());
  FinalMatrix.topRows(startidx) = MOsMatrix.topRows(startidx);
  FinalMatrix.bottomRows(MOsMatrix.rows() - startidx) =
      MOsMatrix.bottomRows(MOsMatrix.rows() - startidx);
  return FinalMatrix;
}
}  // namespace xtp
}  // namespace votca
