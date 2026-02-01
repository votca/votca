
/*
 * Copyright 2009-2023 The VOTCA Development Team
 * (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_ORBITALS_H
#define VOTCA_XTP_ORBITALS_H

// VOTCA includes
#include <votca/tools/globals.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "aobasis.h"
#include "checkpoint.h"
#include "classicalsegment.h"
#include "eigen.h"
#include "qmmolecule.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

/**
 * \brief container for molecular orbitals
 *
 * The Orbitals class stores orbital id, energy, MO coefficients, basis set
 *
 */
class Orbitals {
 public:
  Orbitals();

  bool hasBasisSetSize() const {
    return (dftbasis_.AOBasisSize() > 0) ? true : false;
  }

  void setEmbeddedMOs(tools::EigenSystem &system) { mos_embedding_ = system; }

  const tools::EigenSystem &getEmbeddedMOs() const { return mos_embedding_; }

  void setTruncMOsFullBasis(const Eigen::MatrixXd &expandedMOs) {
    expandedMOs_ = expandedMOs;
  }

  const Eigen::MatrixXd getTruncMOsFullBasis() const { return expandedMOs_; }

  Index getBasisSetSize() const { return dftbasis_.AOBasisSize(); }

  Index getLumo() const { return occupied_levels_; }

  Index getHomo() const { return occupied_levels_ - 1; }

  // access to DFT number of levels, new, tested
  bool hasNumberOfLevels() const {
    return ((occupied_levels_ > 0) ? true : false);
  }
  bool hasNumberOfLevelsBeta() const {
    return ((occupied_levels_beta_ > 0) ? true : false);
  }

  void setNumberOfOccupiedLevels(Index occupied_levels) {
    occupied_levels_ = occupied_levels;
  }
  void setNumberOfOccupiedLevelsBeta(Index occupied_levels_beta) {
    occupied_levels_beta_ = occupied_levels_beta;
  }

  // access to DFT number of electrons, new, tested
  bool hasNumberOfAlphaElectrons() const {
    return (number_alpha_electrons_ > 0) ? true : false;
  }
  bool hasNumberOfBetaElectrons() const {
    return (number_beta_electrons_ > 0) ? true : false;
  }

  Index getNumberOfAlphaElectrons() const { return number_alpha_electrons_; };
  Index getNumberOfBetaElectrons() const { return number_beta_electrons_; };

  void setNumberOfAlphaElectrons(Index electrons) {
    number_alpha_electrons_ = electrons;
  }
  void setNumberOfBetaElectrons(Index electrons) {
    number_beta_electrons_ = electrons;
  }
  bool hasECPName() const { return (ECP_ != "") ? true : false; }

  const std::string &getECPName() const { return ECP_; };

  void setECPName(const std::string &ECP) { ECP_ = ECP; };

  // access to QM package name, new, tested

  bool hasQMpackage() const { return (!qm_package_.empty()); }

  const std::string &getQMpackage() const { return qm_package_; }

  void setQMpackage(const std::string &qmpackage) { qm_package_ = qmpackage; }

  // access to DFT molecular orbital energies, new, tested
  bool hasMOs() const { return (mos_.eigenvalues().size() > 0) ? true : false; }
  bool hasBetaMOs() const {
    return (mos_beta_.eigenvalues().size() > 0) ? true : false;
  }

  const tools::EigenSystem &MOs() const { return mos_; }
  tools::EigenSystem &MOs() { return mos_; }

  const Eigen::MatrixXd &Occupations() const { return occupations_; }
  Eigen::MatrixXd &Occupations() { return occupations_; }

  const tools::EigenSystem &MOs_beta() const { return mos_beta_; }
  tools::EigenSystem &MOs_beta() { return mos_beta_; }

  // determine (pseudo-)degeneracy of a DFT molecular orbital
  std::vector<Index> CheckDegeneracy(Index level,
                                     double energy_difference) const;

  Index NumberofStates(QMStateType type) const {
    switch (type.Type()) {
      case QMStateType::Singlet:
        return Index(BSE_singlet_.eigenvalues().size());
        break;
      case QMStateType::Triplet:
        return Index(BSE_triplet_.eigenvalues().size());
        break;
      case QMStateType::KSstate:
        return Index(mos_.eigenvalues().size());
        break;
      case QMStateType::PQPstate:
        return Index(QPpert_energies_.size());
        break;
      case QMStateType::DQPstate:
        return Index(QPdiag_.eigenvalues().size());
        break;
      default:
        return 1;
        break;
    }
  }

  void setCalculationType(std::string CalcType) { CalcType_ = CalcType; }
  std::string getCalculationType() const { return CalcType_; }

  void setChargeAndSpin(Index charge, Index spin) {
    total_charge_ = charge;
    total_spin_ = spin;
  }

  Index getSpin() const { return total_spin_; }
  Index getCharge() const { return total_charge_; }
  bool isOpenShell() const { return (total_spin_ > 1) ? true : false; }

  bool hasQMAtoms() const { return (atoms_.size() > 0) ? true : false; }

  const QMMolecule &QMAtoms() const { return atoms_; }

  QMMolecule &QMAtoms() { return atoms_; }

  void updateAtomPostion(Index atom_index, Eigen::Vector3d new_position) {
    this->QMAtoms()[atom_index].setPos(new_position);
    dftbasis_.UpdateShellPositions(this->QMAtoms());
    auxbasis_.UpdateShellPositions(this->QMAtoms());
  }

  void setXCFunctionalName(std::string functionalname) {
    functionalname_ = functionalname;
  }
  const std::string &getXCFunctionalName() const { return functionalname_; }

  void setXCGrid(std::string grid) { grid_quality_ = grid; }
  const std::string &getXCGrid() const { return grid_quality_; }

  // access to QM total energy, new, tested
  bool hasQMEnergy() const { return (qm_energy_ != 0.0) ? true : false; }

  double getDFTTotalEnergy() const { return qm_energy_; }

  void setQMEnergy(double qmenergy) { qm_energy_ = qmenergy; }

  // access to DFT basis set name
  bool hasDFTbasisName() const {
    return (!dftbasis_.Name().empty()) ? true : false;
  }

  const std::string &getDFTbasisName() const { return dftbasis_.Name(); }

  void SetupDftBasis(std::string basis_name);

  void SetupAuxBasis(std::string aux_basis_name);

  const AOBasis &getDftBasis() const {
    if (dftbasis_.AOBasisSize() == 0) {
      throw std::runtime_error(
          "Requested the DFT basis, but no basis is present. Make sure "
          "SetupDftBasis is called.");
    } else {
      return dftbasis_;
    }
  }

  const AOBasis &getAuxBasis() const {
    if (auxbasis_.AOBasisSize() == 0) {
      throw std::runtime_error(
          "Requested the Aux basis, but no basis is present. Make sure "
          "SetupAuxBasis is called.");
    } else {
      return auxbasis_;
    }
  }

  /*
   * ======= GW-BSE related functions =======
   */

  // access to auxiliary basis set name

  bool hasAuxbasisName() const {
    return (!auxbasis_.Name().empty()) ? true : false;
  }

  const std::string getAuxbasisName() const { return auxbasis_.Name(); }

  // access to list of indices used in GWA

  bool hasGWAindices() const { return (qpmax_ > 0) ? true : false; }

  void setGWindices(Index qpmin, Index qpmax) {
    qpmin_ = qpmin;
    qpmax_ = qpmax;
  }

  Index getGWAmin() const { return qpmin_; }

  Index getGWAmax() const { return qpmax_; }

  // access to list of indices used in RPA

  bool hasRPAindices() const { return (rpamax_ > 0) ? true : false; }

  void setRPAindices(Index rpamin, Index rpamax) {
    rpamin_ = rpamin;
    rpamax_ = rpamax;
  }

  Index getRPAmin() const { return rpamin_; }

  Index getRPAmax() const { return rpamax_; }

  // access to list of indices used in BSE

  void setTDAApprox(bool usedTDA) { useTDA_ = usedTDA; }
  bool getTDAApprox() const { return useTDA_; }

  bool hasBSEindices() const { return (bse_cmax_ > 0) ? true : false; }

  void setBSEindices(Index vmin, Index cmax) {
    bse_vmin_ = vmin;
    bse_vmax_ = getHomo();
    bse_cmin_ = getLumo();
    bse_cmax_ = cmax;
    bse_vtotal_ = bse_vmax_ - bse_vmin_ + 1;
    bse_ctotal_ = bse_cmax_ - bse_cmin_ + 1;
    bse_size_ = bse_vtotal_ * bse_ctotal_;
    return;
  }

  Index getBSEvmin() const { return bse_vmin_; }

  Index getBSEvmax() const { return bse_vmax_; }

  Index getBSEcmin() const { return bse_cmin_; }

  Index getBSEcmax() const { return bse_cmax_; }

  double getScaHFX() const { return ScaHFX_; }

  void setScaHFX(double ScaHFX) { ScaHFX_ = ScaHFX; }

  // access to perturbative QP energies
  bool hasRPAInputEnergies() const { return (rpa_inputenergies_.size() > 0); }

  const Eigen::VectorXd &RPAInputEnergies() const { return rpa_inputenergies_; }

  Eigen::VectorXd &RPAInputEnergies() { return rpa_inputenergies_; }

  // access to RPA input energies energies
  bool hasQPpert() const {
    return (QPpert_energies_.size() > 0) ? true : false;
  }

  const Eigen::VectorXd &QPpertEnergies() const { return QPpert_energies_; }

  Eigen::VectorXd &QPpertEnergies() { return QPpert_energies_; }

  // access to diagonalized QP energies and wavefunctions

  bool hasQPdiag() const {
    return (QPdiag_.eigenvalues().size() > 0) ? true : false;
  }
  const tools::EigenSystem &QPdiag() const { return QPdiag_; }
  tools::EigenSystem &QPdiag() { return QPdiag_; }

  bool hasBSETriplets() const {
    return (BSE_triplet_.eigenvectors().cols() > 0) ? true : false;
  }

  const tools::EigenSystem &BSETriplets() const { return BSE_triplet_; }

  tools::EigenSystem &BSETriplets() { return BSE_triplet_; }

  // access to singlet energies and wave function coefficients

  bool hasBSESinglets() const {
    return (BSE_singlet_.eigenvectors().cols() > 0) ? true : false;
  }

  const tools::EigenSystem &BSESinglets() const { return BSE_singlet_; }

  tools::EigenSystem &BSESinglets() { return BSE_singlet_; }

  // access to BSE energies with dynamical screening
  bool hasBSESinglets_dynamic() const {
    return (BSE_singlet_energies_dynamic_.size() > 0) ? true : false;
  }

  const Eigen::VectorXd &BSESinglets_dynamic() const {
    return BSE_singlet_energies_dynamic_;
  }

  Eigen::VectorXd &BSESinglets_dynamic() {
    return BSE_singlet_energies_dynamic_;
  }

  bool hasBSETriplets_dynamic() const {
    return (BSE_triplet_energies_dynamic_.size() > 0) ? true : false;
  }

  const Eigen::VectorXd &BSETriplets_dynamic() const {
    return BSE_triplet_energies_dynamic_;
  }

  Eigen::VectorXd &BSETriplets_dynamic() {
    return BSE_triplet_energies_dynamic_;
  }

  // access to transition dipole moments

  bool hasTransitionDipoles() const {
    return (transition_dipoles_.size() > 0) ? true : false;
  }

  const std::vector<Eigen::Vector3d> &TransitionDipoles() const {
    return transition_dipoles_;
  }

  Eigen::VectorXd Oscillatorstrengths() const;

  Eigen::Vector3d CalcElDipole(const QMState &state) const;

  // Calculates full electron density for state or transition density, if you
  // want to calculate only the density contribution of hole or electron use
  // DensityMatrixExcitedState
  Eigen::MatrixXd DensityMatrixFull(const QMState &state) const;
  Eigen::MatrixXd DensityMatrixWithoutGS(const QMState &state) const;

  // functions for calculating density matrices
  Eigen::MatrixXd DensityMatrixGroundState() const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState(
      const QMState &state) const;
  Eigen::MatrixXd DensityMatrixQuasiParticle(const QMState &state) const;
  Eigen::MatrixXd DensityMatrixKSstate(const QMState &state) const;
  Eigen::MatrixXd CalculateQParticleAORepresentation() const;
  double getTotalStateEnergy(const QMState &state) const;    // Hartree
  double getExcitedStateEnergy(const QMState &state) const;  // Hartree

  void OrderMOsbyEnergy();

  void PrepareDimerGuess(const Orbitals &orbitalsA, const Orbitals &orbitalsB);

  void CalcCoupledTransition_Dipoles();

  void WriteToCpt(const std::string &filename) const;

  void ReadFromCpt(const std::string &filename);

  void WriteToCpt(CheckpointWriter &w) const;
  void WriteBasisSetsToCpt(CheckpointWriter &w) const;
  void ReadFromCpt(CheckpointReader &r);
  void ReadBasisSetsFromCpt(CheckpointReader &r);

  bool GetFlagUseHqpOffdiag() const { return use_Hqp_offdiag_; };
  void SetFlagUseHqpOffdiag(bool flag) { use_Hqp_offdiag_ = flag; };

  const Eigen::MatrixXd &getLMOs() const { return lmos_; };
  void setLMOs(const Eigen::MatrixXd &matrix) { lmos_ = matrix; }

  const Eigen::VectorXd &getLMOs_energies() const { return lmos_energies_; };
  void setLMOs_energies(const Eigen::VectorXd &energies) {
    lmos_energies_ = energies;
  }

  Index getNumOfActiveElectrons() { return active_electrons_; }
  void setNumofActiveElectrons(const Index active_electrons) {
    active_electrons_ = active_electrons;
  }

  const Eigen::MatrixXd &getInactiveDensity() const { return inactivedensity_; }
  void setInactiveDensity(Eigen::MatrixXd inactivedensity) {
    inactivedensity_ = inactivedensity;
  }

 private:
  std::array<Eigen::MatrixXd, 3> CalcFreeTransition_Dipoles() const;

  // returns indeces of a re-sorted vector of energies from lowest to highest
  std::vector<Index> SortEnergies();

  void WriteToCpt(CheckpointFile &f) const;

  void ReadFromCpt(CheckpointFile &f);
  Eigen::MatrixXd TransitionDensityMatrix(const QMState &state) const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState_R(
      const QMState &state) const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState_AR(
      const QMState &state) const;
  Eigen::MatrixXd CalcAuxMat_cc(const Eigen::VectorXd &coeffs) const;
  Eigen::MatrixXd CalcAuxMat_vv(const Eigen::VectorXd &coeffs) const;

  Index occupied_levels_;
  Index occupied_levels_beta_ = 0;
  Index number_alpha_electrons_;
  Index number_beta_electrons_;
  std::string ECP_ = "";
  bool useTDA_;

  std::string CalcType_ = "NoEmbedding";

  tools::EigenSystem mos_;
  tools::EigenSystem mos_beta_;
  Eigen::MatrixXd occupations_;

  tools::EigenSystem mos_embedding_;

  Eigen::MatrixXd lmos_;
  Eigen::VectorXd lmos_energies_;
  Index active_electrons_;
  Eigen::MatrixXd inactivedensity_;
  Eigen::MatrixXd expandedMOs_;

  QMMolecule atoms_;

  AOBasis dftbasis_;
  AOBasis auxbasis_;

  double qm_energy_ = 0;

  Index total_charge_;
  Index total_spin_;

  // new variables for GW-BSE storage
  Index rpamin_ = 0;
  Index rpamax_ = 0;

  Index qpmin_ = 0;
  Index qpmax_ = 0;
  Index bse_vmin_ = 0;
  Index bse_vmax_ = 0;
  Index bse_cmin_ = 0;
  Index bse_cmax_ = 0;
  Index bse_size_ = 0;
  Index bse_vtotal_ = 0;
  Index bse_ctotal_ = 0;

  double ScaHFX_ = 0;

  std::string functionalname_ = "";
  std::string grid_quality_ = "";
  std::string qm_package_ = "";

  Eigen::VectorXd rpa_inputenergies_;
  // perturbative quasiparticle energies
  Eigen::VectorXd QPpert_energies_;

  // quasiparticle energies and coefficients after diagonalization
  tools::EigenSystem QPdiag_;

  tools::EigenSystem BSE_singlet_;
  std::vector<Eigen::Vector3d> transition_dipoles_;
  tools::EigenSystem BSE_triplet_;

  // singlet and triplet energies after perturbative dynamical screening
  Eigen::VectorXd BSE_singlet_energies_dynamic_;
  Eigen::VectorXd BSE_triplet_energies_dynamic_;

  bool use_Hqp_offdiag_ = true;

  // Version 2: adds BSE energies after perturbative dynamical screening
  // Version 3: changed shell ordering
  // Version 4: added vxc grid quality
  // Version 5: added the dft and aux basisset
  // Version 6: added spin in dft
  static constexpr int orbitals_version() { return 6; }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORBITALS_H
