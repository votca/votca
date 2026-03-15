
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
  /// Construct an empty orbital container with default metadata.
  Orbitals();

  /// Report whether a DFT AO basis has already been attached.
  bool hasBasisSetSize() const {
    return (dftbasis_.AOBasisSize() > 0) ? true : false;
  }

  /// Store molecular orbitals obtained from an embedding calculation.
  void setEmbeddedMOs(tools::EigenSystem &system) { mos_embedding_ = system; }

  /// Return molecular orbitals obtained from an embedding calculation.
  const tools::EigenSystem &getEmbeddedMOs() const { return mos_embedding_; }

  /// Store truncated active-region orbitals expanded back to the full AO basis.
  void setTruncMOsFullBasis(const Eigen::MatrixXd &expandedMOs) {
    expandedMOs_ = expandedMOs;
  }

  /// Return truncated active-region orbitals represented in the full AO basis.
  const Eigen::MatrixXd getTruncMOsFullBasis() const { return expandedMOs_; }

  /// Return the number of AO basis functions in the DFT basis.
  Index getBasisSetSize() const { return dftbasis_.AOBasisSize(); }

  /// Return the alpha-spin LUMO index inferred from the stored electron counts.
  Index getLumoAlpha() const {
    if (number_alpha_electrons_ > 0) {
      return number_alpha_electrons_;
    }
    return occupied_levels_;
  }

  /// Return the alpha-spin HOMO index.
  Index getHomoAlpha() const { return getLumoAlpha() - 1; }

  /// Return the beta-spin LUMO index, including legacy restricted closed-shell fallback logic.
  Index getLumoBeta() const {
    if (number_beta_electrons_ > 0) {
      return number_beta_electrons_;
    }

    // Legacy restricted closed-shell fallback
    if (!hasUnrestrictedOrbitals() && total_spin_ == 1) {
      return occupied_levels_;
    }

    return occupied_levels_beta_;
  }

  /// Return the beta-spin HOMO index.
  Index getHomoBeta() const { return getLumoBeta() - 1; }

  // Current generic convention: use alpha frontier orbitals
  /// Return the default LUMO index used by spin-agnostic callers.
  Index getLumo() const { return getLumoAlpha(); }
  /// Return the default HOMO index used by spin-agnostic callers.
  Index getHomo() const { return getHomoAlpha(); }

  // access to DFT number of levels, new, tested
  /// Report whether the number of occupied spatial orbitals has been set.
  bool hasNumberOfLevels() const {
    return ((occupied_levels_ > 0) ? true : false);
  }
  /// Report whether the beta-spin occupied-orbital count has been set explicitly.
  bool hasNumberOfLevelsBeta() const {
    return ((occupied_levels_beta_ > 0) ? true : false);
  }

  /// Store the number of occupied spatial orbitals and update closed-shell electron counts when applicable.
  void setNumberOfOccupiedLevels(Index occupied_levels) {
    occupied_levels_ = occupied_levels;

    // Backward compatibility for legacy restricted/singlet workflows:
    // many callers only set occupied_levels_ and expect a closed-shell density.
    if (!hasUnrestrictedOrbitals() && total_spin_ == 1) {
      number_alpha_electrons_ = occupied_levels;
      number_beta_electrons_ = occupied_levels;
    }
  }

  /// Store the number of occupied beta-spin orbitals.
  void setNumberOfOccupiedLevelsBeta(Index occupied_levels_beta) {
    occupied_levels_beta_ = occupied_levels_beta;
  }

  /// Store the total number of alpha electrons.
  void setNumberOfAlphaElectrons(Index electrons) {
    number_alpha_electrons_ = electrons;
  }

  /// Store the total number of beta electrons.
  void setNumberOfBetaElectrons(Index electrons) {
    number_beta_electrons_ = electrons;
  }

  // access to DFT number of electrons, new, tested
  /// Report whether the alpha-electron count has been set explicitly.
  bool hasNumberOfAlphaElectrons() const {
    return (number_alpha_electrons_ > 0) ? true : false;
  }
  /// Report whether the beta-electron count has been set explicitly.
  bool hasNumberOfBetaElectrons() const {
    return (number_beta_electrons_ > 0) ? true : false;
  }

  /// Return the stored number of alpha electrons.
  Index getNumberOfAlphaElectrons() const { return number_alpha_electrons_; };
  /// Return the stored number of beta electrons.
  Index getNumberOfBetaElectrons() const { return number_beta_electrons_; };


  /// Report whether an effective core potential label has been stored.
  bool hasECPName() const { return (ECP_ != "") ? true : false; }

  /// Return the effective core potential label.
  const std::string &getECPName() const { return ECP_; };

  /// Store the effective core potential label.
  void setECPName(const std::string &ECP) { ECP_ = ECP; };

  // access to QM package name, new, tested

  /// Report whether the originating QM package name has been stored.
  bool hasQMpackage() const { return (!qm_package_.empty()); }

  /// Return the stored QM package name.
  const std::string &getQMpackage() const { return qm_package_; }

  /// Store the name of the QM package that produced these orbitals.
  void setQMpackage(const std::string &qmpackage) { qm_package_ = qmpackage; }

  // access to DFT molecular orbital energies, new, tested
  /// Report whether alpha/restricted molecular orbitals are available.
  bool hasMOs() const { return (mos_.eigenvalues().size() > 0) ? true : false; }
  /// Report whether beta-spin molecular orbitals are available.
  bool hasBetaMOs() const {
    return (mos_beta_.eigenvalues().size() > 0) ? true : false;
  }

  /// Return read-only access to alpha/restricted molecular orbitals.
  const tools::EigenSystem &MOs() const { return mos_; }
  /// Return writable access to alpha/restricted molecular orbitals.
  tools::EigenSystem &MOs() { return mos_; }

  /// Return the stored fractional occupation matrix.
  const Eigen::MatrixXd &Occupations() const { return occupations_; }
  /// Return writable access to the fractional occupation matrix.
  Eigen::MatrixXd &Occupations() { return occupations_; }

  /// Return read-only access to beta-spin molecular orbitals.
  const tools::EigenSystem &MOs_beta() const { return mos_beta_; }
  /// Return writable access to beta-spin molecular orbitals.
  tools::EigenSystem &MOs_beta() { return mos_beta_; }

  // determine (pseudo-)degeneracy of a DFT molecular orbital
  /// Return all orbitals whose energies are degenerate with the requested level within a tolerance.
  std::vector<Index> CheckDegeneracy(Index level,
                                     double energy_difference) const;

  /// Return the number of states available for the requested state family.
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

  /// Store the calculation-type tag used by downstream workflows.
  void setCalculationType(std::string CalcType) { CalcType_ = CalcType; }
  /// Return the stored calculation-type tag.
  std::string getCalculationType() const { return CalcType_; }

  /// Store the total charge and spin multiplicity associated with the orbital set.
  void setChargeAndSpin(Index charge, Index spin) {
    total_charge_ = charge;
    total_spin_ = spin;
  }

  /// Return the stored spin multiplicity.
  Index getSpin() const { return total_spin_; }
  /// Return the stored total charge.
  Index getCharge() const { return total_charge_; }
  /// Report whether the stored state corresponds to an open-shell system.
  bool isOpenShell() const { return (total_spin_ > 1) ? true : false; }

  /// Report whether a molecular geometry has been stored.
  bool hasQMAtoms() const { return (atoms_.size() > 0) ? true : false; }

  /// Return read-only access to the molecular geometry.
  const QMMolecule &QMAtoms() const { return atoms_; }

  /// Return writable access to the molecular geometry.
  QMMolecule &QMAtoms() { return atoms_; }

  /// Update one atomic position and keep attached AO basis shells synchronized.
  void updateAtomPostion(Index atom_index, Eigen::Vector3d new_position) {
    this->QMAtoms()[atom_index].setPos(new_position);
    dftbasis_.UpdateShellPositions(this->QMAtoms());
    auxbasis_.UpdateShellPositions(this->QMAtoms());
  }

  /// Store the exchange-correlation functional label associated with the orbitals.
  void setXCFunctionalName(std::string functionalname) {
    functionalname_ = functionalname;
  }
  /// Return the exchange-correlation functional label.
  const std::string &getXCFunctionalName() const { return functionalname_; }

  /// Store the numerical XC grid quality label.
  void setXCGrid(std::string grid) { grid_quality_ = grid; }
  /// Return the numerical XC grid quality label.
  const std::string &getXCGrid() const { return grid_quality_; }

  // access to QM total energy, new, tested
  /// Report whether a total electronic energy has been stored.
  bool hasQMEnergy() const { return (qm_energy_ != 0.0) ? true : false; }

  /// Return the stored total DFT energy.
  double getDFTTotalEnergy() const { return qm_energy_; }

  /// Store the total DFT energy.
  void setQMEnergy(double qmenergy) { qm_energy_ = qmenergy; }

  // access to DFT basis set name
  /// Report whether a DFT basis-set name has been stored.
  bool hasDFTbasisName() const {
    return (!dftbasis_.Name().empty()) ? true : false;
  }

  /// Return the DFT basis-set name.
  const std::string &getDFTbasisName() const { return dftbasis_.Name(); }

  /// Build and attach the DFT AO basis from the stored molecular geometry.
  void SetupDftBasis(std::string basis_name);

  /// Build and attach the auxiliary AO basis from the stored molecular geometry.
  void SetupAuxBasis(std::string aux_basis_name);

  /// Return the DFT AO basis, throwing if it has not been initialized.
  const AOBasis &getDftBasis() const {
    if (dftbasis_.AOBasisSize() == 0) {
      throw std::runtime_error(
          "Requested the DFT basis, but no basis is present. Make sure "
          "SetupDftBasis is called.");
    } else {
      return dftbasis_;
    }
  }

  /// Return the auxiliary AO basis, throwing if it has not been initialized.
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

  /// Report whether an auxiliary basis-set name has been stored.
  bool hasAuxbasisName() const {
    return (!auxbasis_.Name().empty()) ? true : false;
  }

  /// Return the auxiliary basis-set name.
  const std::string getAuxbasisName() const { return auxbasis_.Name(); }

  // access to list of indices used in GWA

  /// Report whether the GW quasiparticle window has been defined.
  bool hasGWAindices() const { return (qpmax_ > 0) ? true : false; }

  /// Store the orbital window used in GW calculations.
  void setGWindices(Index qpmin, Index qpmax) {
    qpmin_ = qpmin;
    qpmax_ = qpmax;
  }

  /// Return the lower GW orbital index.
  Index getGWAmin() const { return qpmin_; }

  /// Return the upper GW orbital index.
  Index getGWAmax() const { return qpmax_; }

  // access to list of indices used in RPA

  /// Report whether the RPA window has been defined.
  bool hasRPAindices() const { return (rpamax_ > 0) ? true : false; }

  /// Store the orbital window used in RPA calculations.
  void setRPAindices(Index rpamin, Index rpamax) {
    rpamin_ = rpamin;
    rpamax_ = rpamax;
  }

  /// Return the lower RPA orbital index.
  Index getRPAmin() const { return rpamin_; }

  /// Return the upper RPA orbital index.
  Index getRPAmax() const { return rpamax_; }

  // access to list of indices used in BSE

  /// Enable or disable the Tamm-Dancoff approximation flag.
  void setTDAApprox(bool usedTDA) { useTDA_ = usedTDA; }
  /// Return whether the Tamm-Dancoff approximation is enabled.
  bool getTDAApprox() const { return useTDA_; }

  /// Report whether the BSE valence/conduction window has been defined.
  bool hasBSEindices() const { return (bse_cmax_ > 0) ? true : false; }

  /// Define the BSE excitation window from a valence lower bound and conduction upper bound.
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

  /// Return the lowest valence orbital included in BSE.
  Index getBSEvmin() const { return bse_vmin_; }

  /// Return the highest valence orbital included in BSE.
  Index getBSEvmax() const { return bse_vmax_; }

  /// Return the lowest conduction orbital included in BSE.
  Index getBSEcmin() const { return bse_cmin_; }

  /// Return the highest conduction orbital included in BSE.
  Index getBSEcmax() const { return bse_cmax_; }

  /// Return the fraction of exact exchange associated with the functional.
  double getScaHFX() const { return ScaHFX_; }

  /// Store the fraction of exact exchange associated with the functional.
  void setScaHFX(double ScaHFX) { ScaHFX_ = ScaHFX; }

  // access to perturbative QP energies
  /// Report whether RPA input energies are available.
  bool hasRPAInputEnergies() const { return (rpa_inputenergies_.size() > 0); }

  /// Return read-only access to the RPA input energies.
  const Eigen::VectorXd &RPAInputEnergies() const { return rpa_inputenergies_; }

  /// Return writable access to the RPA input energies.
  Eigen::VectorXd &RPAInputEnergies() { return rpa_inputenergies_; }

  // access to RPA input energies energies
  /// Report whether perturbative quasiparticle energies are available.
  bool hasQPpert() const {
    return (QPpert_energies_.size() > 0) ? true : false;
  }

  /// Return read-only access to perturbative quasiparticle energies.
  const Eigen::VectorXd &QPpertEnergies() const { return QPpert_energies_; }

  /// Return writable access to perturbative quasiparticle energies.
  Eigen::VectorXd &QPpertEnergies() { return QPpert_energies_; }

  // access to diagonalized QP energies and wavefunctions

  /// Report whether diagonalized quasiparticle energies and orbitals are available.
  bool hasQPdiag() const {
    return (QPdiag_.eigenvalues().size() > 0) ? true : false;
  }
  /// Return read-only access to the diagonalized quasiparticle representation.
  const tools::EigenSystem &QPdiag() const { return QPdiag_; }
  /// Return writable access to the diagonalized quasiparticle representation.
  tools::EigenSystem &QPdiag() { return QPdiag_; }

  /// Report whether triplet BSE eigenpairs are available.
  bool hasBSETriplets() const {
    return (BSE_triplet_.eigenvectors().cols() > 0) ? true : false;
  }

  /// Return read-only access to triplet BSE eigenpairs.
  const tools::EigenSystem &BSETriplets() const { return BSE_triplet_; }

  /// Return writable access to triplet BSE eigenpairs.
  tools::EigenSystem &BSETriplets() { return BSE_triplet_; }

  // access to singlet energies and wave function coefficients

  /// Report whether singlet BSE eigenpairs are available.
  bool hasBSESinglets() const {
    return (BSE_singlet_.eigenvectors().cols() > 0) ? true : false;
  }

  /// Return read-only access to singlet BSE eigenpairs.
  const tools::EigenSystem &BSESinglets() const { return BSE_singlet_; }

  /// Return writable access to singlet BSE eigenpairs.
  tools::EigenSystem &BSESinglets() { return BSE_singlet_; }

  // access to BSE energies with dynamical screening
  /// Report whether dynamically screened singlet BSE energies are available.
  bool hasBSESinglets_dynamic() const {
    return (BSE_singlet_energies_dynamic_.size() > 0) ? true : false;
  }

  /// Return dynamically screened singlet BSE energies.
  const Eigen::VectorXd &BSESinglets_dynamic() const {
    return BSE_singlet_energies_dynamic_;
  }

  /// Return writable access to dynamically screened singlet BSE energies.
  Eigen::VectorXd &BSESinglets_dynamic() {
    return BSE_singlet_energies_dynamic_;
  }

  /// Report whether dynamically screened triplet BSE energies are available.
  bool hasBSETriplets_dynamic() const {
    return (BSE_triplet_energies_dynamic_.size() > 0) ? true : false;
  }

  /// Return dynamically screened triplet BSE energies.
  const Eigen::VectorXd &BSETriplets_dynamic() const {
    return BSE_triplet_energies_dynamic_;
  }

  /// Return writable access to dynamically screened triplet BSE energies.
  Eigen::VectorXd &BSETriplets_dynamic() {
    return BSE_triplet_energies_dynamic_;
  }

  // access to transition dipole moments

  /// Report whether transition dipole moments have been computed.
  bool hasTransitionDipoles() const {
    return (transition_dipoles_.size() > 0) ? true : false;
  }

  /// Return the stored transition dipole moments.
  const std::vector<Eigen::Vector3d> &TransitionDipoles() const {
    return transition_dipoles_;
  }

  /// Compute oscillator strengths from the stored excitation energies and transition dipoles.
  Eigen::VectorXd Oscillatorstrengths() const;

  /// Compute the electronic dipole moment associated with a state density.
  Eigen::Vector3d CalcElDipole(const QMState &state) const;

  // Calculates full electron density for state or transition density, if you
  // want to calculate only the density contribution of hole or electron use
  // DensityMatrixExcitedState
  /// Build the full AO density matrix for an excited or transition state including the ground-state part.
  Eigen::MatrixXd DensityMatrixFull(const QMState &state) const;
  /// Build the excited- or transition-state density contribution without the ground-state density.
  Eigen::MatrixXd DensityMatrixWithoutGS(const QMState &state) const;

  // functions for calculating density matrices
  /// Build the ground-state AO density matrix from the stored orbitals.
  Eigen::MatrixXd DensityMatrixGroundState() const;
  /// Build separate hole and electron AO densities for an excited state.
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState(
      const QMState &state) const;
  /// Build the AO density matrix for a quasiparticle state.
  Eigen::MatrixXd DensityMatrixQuasiParticle(const QMState &state) const;
  /// Build the AO density matrix for a single KS orbital state.
  Eigen::MatrixXd DensityMatrixKSstate(const QMState &state) const;
  /// Transform quasiparticle amplitudes into the AO representation.
  Eigen::MatrixXd CalculateQParticleAORepresentation() const;
  /// Return the absolute total energy of the requested state in Hartree.
  double getTotalStateEnergy(const QMState &state) const;    // Hartree
  /// Return the excitation energy of the requested state in Hartree.
  double getExcitedStateEnergy(const QMState &state) const;  // Hartree

  /// Resort molecular orbitals and associated quantities by increasing orbital energy.
  void OrderMOsbyEnergy();

  /// Build a dimer starting guess by combining orbital information from two fragments.
  void PrepareDimerGuess(const Orbitals &orbitalsA, const Orbitals &orbitalsB);

  /// Compute transition dipoles for coupled excited states from the stored BSE data.
  void CalcCoupledTransition_Dipoles();

  /// Write the orbital container to a checkpoint file on disk.
  void WriteToCpt(const std::string &filename) const;

  /// Read the orbital container from a checkpoint file on disk.
  void ReadFromCpt(const std::string &filename);

  /// Serialize the orbital container into an already-open checkpoint writer.
  void WriteToCpt(CheckpointWriter &w) const;
  /// Serialize the attached AO basis sets into an already-open checkpoint writer.
  void WriteBasisSetsToCpt(CheckpointWriter &w) const;
  /// Deserialize the orbital container from an already-open checkpoint reader.
  void ReadFromCpt(CheckpointReader &r);
  /// Deserialize attached AO basis sets from an already-open checkpoint reader.
  void ReadBasisSetsFromCpt(CheckpointReader &r);

  /// Return whether off-diagonal quasiparticle Hamiltonian elements should be used.
  bool GetFlagUseHqpOffdiag() const { return use_Hqp_offdiag_; };
  /// Set whether off-diagonal quasiparticle Hamiltonian elements should be used.
  void SetFlagUseHqpOffdiag(bool flag) { use_Hqp_offdiag_ = flag; };

  /// Return localized molecular orbitals, if available.
  const Eigen::MatrixXd &getLMOs() const { return lmos_; };
  /// Store localized molecular orbitals.
  void setLMOs(const Eigen::MatrixXd &matrix) { lmos_ = matrix; }

  /// Return the energies associated with localized molecular orbitals.
  const Eigen::VectorXd &getLMOs_energies() const { return lmos_energies_; };
  /// Store the energies associated with localized molecular orbitals.
  void setLMOs_energies(const Eigen::VectorXd &energies) {
    lmos_energies_ = energies;
  }

  /// Return the number of electrons assigned to the active region.
  Index getNumOfActiveElectrons() { return active_electrons_; }
  /// Store the number of electrons assigned to the active region.
  void setNumofActiveElectrons(const Index active_electrons) {
    active_electrons_ = active_electrons;
  }

  /// Return the inactive-region density matrix used in embedding workflows.
  const Eigen::MatrixXd &getInactiveDensity() const { return inactivedensity_; }
  /// Store the inactive-region density matrix used in embedding workflows.
  void setInactiveDensity(Eigen::MatrixXd inactivedensity) {
    inactivedensity_ = inactivedensity;
  }

  /************************************
  * Extensions spin DFT
  *************************************/
/// Report whether separate beta-spin orbitals are present.
bool hasUnrestrictedOrbitals() const { return hasBetaMOs(); }

/// Report whether the orbitals represent a restricted open-shell reference.
bool isRestrictedOpenShell() const {
  return total_spin_ > 1 && !hasUnrestrictedOrbitals();
}

/// Report whether the orbitals represent an unrestricted open-shell reference.
bool isUnrestrictedOpenShell() const {
  return total_spin_ > 1 && hasUnrestrictedOrbitals();
}
 /// Build separate alpha and beta ground-state density matrices from the stored orbitals.
 std::array<Eigen::MatrixXd, 2> DensityMatrixGroundStateSpinResolved() const;

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

  Index occupied_levels_ =0 ;
  Index occupied_levels_beta_ = 0;
  Index number_alpha_electrons_ = 0;
  Index number_beta_electrons_ = 0;
  std::string ECP_ = "";
  bool useTDA_;

  std::string CalcType_ = "NoEmbedding";

  tools::EigenSystem mos_;
  tools::EigenSystem mos_beta_;
  Eigen::MatrixXd occupations_;

  tools::EigenSystem mos_embedding_;

  Eigen::MatrixXd lmos_;
  Eigen::VectorXd lmos_energies_;
  Index active_electrons_ = 0;
  Eigen::MatrixXd inactivedensity_;
  Eigen::MatrixXd expandedMOs_;

  QMMolecule atoms_;

  AOBasis dftbasis_;
  AOBasis auxbasis_;

  double qm_energy_ = 0;

  Index total_charge_ =0;
  Index total_spin_ =1;

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

  bool use_Hqp_offdiag_ = false;

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
