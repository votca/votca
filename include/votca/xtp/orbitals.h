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

#pragma once
#ifndef VOTCA_XTP_ORBITALS_H
#define VOTCA_XTP_ORBITALS_H

#include "aobasis.h"
#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/qmstate.h>

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

  bool hasBasisSetSize() const { return (_basis_set_size > 0) ? true : false; }

  int getBasisSetSize() const { return _basis_set_size; }

  void setBasisSetSize(int basis_set_size) { _basis_set_size = basis_set_size; }

  int getLumo() const { return _occupied_levels; }

  int getHomo() const { return _occupied_levels - 1; }
  // access to DFT number of levels, new, tested

  bool hasNumberOfLevels() const {
    return ((_occupied_levels > 0) ? true : false);
  }

  void setNumberOfOccupiedLevels(int occupied_levels) {
    _occupied_levels = occupied_levels;
  }

  // access to DFT number of electrons, new, tested

  bool hasNumberOfAlphaElectrons() const {
    return (_number_alpha_electrons > 0) ? true : false;
  }

  int getNumberOfAlphaElectrons() const { return _number_alpha_electrons; };

  void setNumberOfAlphaElectrons(int electrons) {
    _number_alpha_electrons = electrons;
  }

  bool hasECPName() const { return (_ECP != "") ? true : false; }

  const std::string &getECPName() const { return _ECP; };

  void setECPName(const std::string &ECP) { _ECP = ECP; };

  // access to QM package name, new, tested

  bool hasQMpackage() const { return (!_qm_package.empty()); }

  const std::string &getQMpackage() const { return _qm_package; }

  void setQMpackage(const std::string &qmpackage) { _qm_package = qmpackage; }

  // access to DFT molecular orbital energies, new, tested
  bool hasMOs() const { return (_mos.eigenvalues().size() > 0) ? true : false; }

  const tools::EigenSystem &MOs() const { return _mos; }
  tools::EigenSystem &MOs() { return _mos; }

  // access to DFT molecular orbital energy of a specific level
  double getMOEnergy(int level) const {
    if (level < _mos.eigenvalues().size()) {
      return _mos.eigenvalues()[level];
    } else {
      throw std::runtime_error("Level index is outside array range");
    }
    return 0;
  }

  // determine (pseudo-)degeneracy of a DFT molecular orbital
  std::vector<int> CheckDegeneracy(int level, double energy_difference) const;

  int NumberofStates(QMStateType type) const {
    switch (type.Type()) {
      case QMStateType::Singlet:
        return _BSE_singlet.eigenvalues().size();
        break;
      case QMStateType::Triplet:
        return _BSE_triplet.eigenvalues().size();
        break;
      case QMStateType::KSstate:
        return _mos.eigenvalues().size();
        break;
      case QMStateType::PQPstate:
        return _QPpert_energies.rows();
        break;
      case QMStateType::DQPstate:
        return _QPdiag.eigenvalues().size();
        break;
      default:
        return 1;
        break;
    }
  }

  bool hasQMAtoms() const { return (_atoms.size() > 0) ? true : false; }

  const QMMolecule &QMAtoms() const { return _atoms; }

  QMMolecule &QMAtoms() { return _atoms; }

  // access to classical self-energy in MM environment, new, tested
  bool hasSelfEnergy() const { return (_self_energy != 0.0) ? true : false; }

  double getSelfEnergy() const { return _self_energy; }

  void setSelfEnergy(double selfenergy) { _self_energy = selfenergy; }

  void setXCFunctionalName(std::string functionalname) {
    _functionalname = functionalname;
  }
  const std::string &getXCFunctionalName() const { return _functionalname; }

  // access to QM total energy, new, tested
  bool hasQMEnergy() const { return (_qm_energy != 0.0) ? true : false; }

  double getDFTTotalEnergy() const { return _qm_energy; }

  void setQMEnergy(double qmenergy) { _qm_energy = qmenergy; }

  // access to DFT basis set name

  bool hasDFTbasisName() const { return (!_dftbasis.empty()) ? true : false; }

  void setDFTbasisName(const std::string basis) { _dftbasis = basis; }

  const std::string &getDFTbasisName() const { return _dftbasis; }

  AOBasis SetupDftBasis() const { return SetupBasis<true>(); }
  AOBasis SetupAuxBasis() const { return SetupBasis<false>(); }

  /*
   *  ======= GW-BSE related functions =======
   */

  // access to auxiliary basis set name

  bool hasAuxbasisName() const { return (!_auxbasis.empty()) ? true : false; }

  void setAuxbasisName(std::string basis) { _auxbasis = basis; }

  const std::string &getAuxbasisName() const { return _auxbasis; }

  // access to list of indices used in GWA

  bool hasGWAindices() const { return (_qpmax > 0) ? true : false; }

  void setGWindices(int qpmin, int qpmax) {
    _qpmin = qpmin;
    _qpmax = qpmax;
  }

  int getGWAmin() const { return _qpmin; }

  int getGWAmax() const { return _qpmax; }

  // access to list of indices used in RPA

  bool hasRPAindices() const { return (_rpamax > 0) ? true : false; }

  void setRPAindices(int rpamin, int rpamax) {
    _rpamin = rpamin;
    _rpamax = rpamax;
  }

  int getRPAmin() const { return _rpamin; }

  int getRPAmax() const { return _rpamax; }

  // access to list of indices used in BSE

  void setTDAApprox(bool usedTDA) { _useTDA = usedTDA; }
  bool getTDAApprox() const { return _useTDA; }

  bool hasBSEindices() const { return (_bse_cmax > 0) ? true : false; }

  void setBSEindices(int vmin, int cmax) {
    _bse_vmin = vmin;
    _bse_vmax = getHomo();
    _bse_cmin = getLumo();
    _bse_cmax = cmax;
    _bse_vtotal = _bse_vmax - _bse_vmin + 1;
    _bse_ctotal = _bse_cmax - _bse_cmin + 1;
    _bse_size = _bse_vtotal * _bse_ctotal;
    return;
  }

  int getBSEvmin() const { return _bse_vmin; }

  int getBSEvmax() const { return _bse_vmax; }

  int getBSEcmin() const { return _bse_cmin; }

  int getBSEcmax() const { return _bse_cmax; }

  double getScaHFX() const { return _ScaHFX; }

  void setScaHFX(double ScaHFX) { _ScaHFX = ScaHFX; }

  // access to perturbative QP energies
  bool hasQPpert() const {
    return (_QPpert_energies.size() > 0) ? true : false;
  }

  const Eigen::MatrixXd &QPpertEnergies() const { return _QPpert_energies; }

  Eigen::MatrixXd &QPpertEnergies() { return _QPpert_energies; }

  // access to diagonalized QP energies and wavefunctions

  bool hasQPdiag() const {
    return (_QPdiag.eigenvalues().size() > 0) ? true : false;
  }
  const tools::EigenSystem &QPdiag() const { return _QPdiag; }
  tools::EigenSystem &QPdiag() { return _QPdiag; }

  bool hasBSETriplets() const {
    return (_BSE_triplet.eigenvectors().cols() > 0) ? true : false;
  }

  const tools::EigenSystem &BSETriplets() const { return _BSE_triplet; }

  tools::EigenSystem &BSETriplets() { return _BSE_triplet; }

  // access to singlet energies and wave function coefficients

  bool hasBSESinglets() const {
    return (_BSE_singlet.eigenvectors().cols() > 0) ? true : false;
  }

  const tools::EigenSystem &BSESinglets() const { return _BSE_singlet; }

  tools::EigenSystem &BSESinglets() { return _BSE_singlet; }

  // access to transition dipole moments

  bool hasTransitionDipoles() const {
    return (_transition_dipoles.size() > 0) ? true : false;
  }

  const std::vector<Eigen::Vector3d> &TransitionDipoles() const {
    return _transition_dipoles;
  }

  Eigen::VectorXd Oscillatorstrengths() const;

  Eigen::Vector3d CalcElDipole(const QMState &state) const;

  // Calculates full electron density for state or transition density, if you
  // want to calculate only the density contribution of hole or electron use
  // DensityMatrixExcitedState
  Eigen::MatrixXd DensityMatrixFull(const QMState &state) const;

  // functions for calculating density matrices
  Eigen::MatrixXd DensityMatrixGroundState() const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState(
      const QMState &state) const;
  Eigen::MatrixXd DensityMatrixQuasiParticle(const QMState &state) const;
  Eigen::MatrixXd CalculateQParticleAORepresentation() const;
  double getTotalStateEnergy(const QMState &state) const;    // Hartree
  double getExcitedStateEnergy(const QMState &state) const;  // Hartree

  void OrderMOsbyEnergy();

  void PrepareDimerGuess(const Orbitals &orbitalsA, const Orbitals &orbitalsB);

  void CalcCoupledTransition_Dipoles();

  void WriteToCpt(const std::string &filename) const;

  void ReadFromCpt(const std::string &filename);

  void WriteToCpt(CheckpointWriter w) const;
  void ReadFromCpt(CheckpointReader parent);

 private:
  std::array<Eigen::MatrixXd, 3> CalcFreeTransition_Dipoles() const;

  // returns indeces of a re-sorted vector of energies from lowest to highest
  std::vector<int> SortEnergies();

  template <bool dftbasis>
  AOBasis SetupBasis() const {
    BasisSet bs;
    if (dftbasis) {
      bs.Load(this->getDFTbasisName());
    } else {
      bs.Load(this->getAuxbasisName());
    }
    AOBasis basis;
    basis.Fill(bs, this->QMAtoms());
    return basis;
  }

  void WriteToCpt(CheckpointFile f) const;

  void ReadFromCpt(CheckpointFile f);
  Eigen::MatrixXd TransitionDensityMatrix(const QMState &state) const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState_R(
      const QMState &state) const;
  std::array<Eigen::MatrixXd, 2> DensityMatrixExcitedState_AR(
      const QMState &state) const;
  Eigen::MatrixXd CalcAuxMat_cc(const Eigen::VectorXd &coeffs) const;
  Eigen::MatrixXd CalcAuxMat_vv(const Eigen::VectorXd &coeffs) const;

  int _basis_set_size;
  int _occupied_levels;
  int _number_alpha_electrons;
  std::string _ECP = "";
  bool _useTDA;

  tools::EigenSystem _mos;

  QMMolecule _atoms;

  double _qm_energy = 0;
  double _self_energy = 0;

  // new variables for GW-BSE storage
  int _rpamin = 0;
  int _rpamax = 0;

  int _qpmin = 0;
  int _qpmax = 0;

  int _bse_vmin = 0;
  int _bse_vmax = 0;
  int _bse_cmin = 0;
  int _bse_cmax = 0;
  int _bse_size = 0;
  int _bse_vtotal = 0;
  int _bse_ctotal = 0;

  double _ScaHFX = 0;

  std::string _dftbasis = "";
  std::string _auxbasis = "";

  std::string _functionalname = "";
  std::string _qm_package = "";

  // perturbative quasiparticle energies
  Eigen::MatrixXd _QPpert_energies;

  // quasiparticle energies and coefficients after diagonalization

  tools::EigenSystem _QPdiag;

  tools::EigenSystem _BSE_singlet;
  std::vector<Eigen::Vector3d> _transition_dipoles;
  tools::EigenSystem _BSE_triplet;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORBITALS_H
