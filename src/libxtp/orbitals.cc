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

#include "votca/xtp/orbitals.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aomatrix3d.h"
#include "votca/xtp/qmstate.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdio.h>
#include <votca/xtp/vc2index.h>
#include <votca/xtp/version.h>

namespace votca {
namespace xtp {

Orbitals::Orbitals() : _atoms("", 0) { ; }

/**
 *
 * @param _energy_difference [ev] Two levels are degenerate if their energy is
 * smaller than this value
 * @return vector with indices off all orbitals degenerate to this including
 * itself
 */
std::vector<int> Orbitals::CheckDegeneracy(int level,
                                           double energy_difference) const {

  std::vector<int> result;
  if (level > _mos.eigenvalues().size()) {
    throw std::runtime_error(
        "Level for degeneracy is higher than maximum level");
  }
  double MOEnergyLevel = _mos.eigenvalues()(level);

  for (int i = 0; i < _mos.eigenvalues().size(); ++i) {
    if (std::abs(_mos.eigenvalues()(i) - MOEnergyLevel) < energy_difference) {
      result.push_back(i);
    }
  }

  if (result.empty()) {
    result.push_back(level);
  }
  return result;
}

std::vector<int> Orbitals::SortEnergies() {
  std::vector<int> index = std::vector<int>(_mos.eigenvalues().size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [this](int i1, int i2) {
    return this->MOs().eigenvalues()[i1] < this->MOs().eigenvalues()[i2];
  });
  return index;
}

Eigen::MatrixXd Orbitals::DensityMatrixFull(const QMState& state) const {
  if (state.isTransition()) {
    return this->TransitionDensityMatrix(state);
  }
  Eigen::MatrixXd result = this->DensityMatrixGroundState();
  if (state.Type().isExciton()) {
    std::array<Eigen::MatrixXd, 2> DMAT = DensityMatrixExcitedState(state);
    result = result - DMAT[0] + DMAT[1];  // Ground state + hole_contribution +
                                          // electron contribution
  } else if (state.Type() == QMStateType::DQPstate) {
    Eigen::MatrixXd DMATQP = DensityMatrixQuasiParticle(state);
    if (state.Index() > getHomo()) {
      result += DMATQP;
    } else {
      result -= DMATQP;
    }
  } else if (state.Type() != QMStateType::Gstate) {
    throw std::runtime_error(
        "DensityMatrixFull does not yet implement QMStateType:" +
        state.Type().ToLongString());
  }
  return result;
}

// Determine ground state density matrix

Eigen::MatrixXd Orbitals::DensityMatrixGroundState() const {
  if (!hasMOs()) {
    throw std::runtime_error("Orbitals file does not contain MO coefficients");
  }
  Eigen::MatrixXd occstates = _mos.eigenvectors().leftCols(_occupied_levels);
  Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
  return dmatGS;
}

Eigen::MatrixXd Orbitals::CalculateQParticleAORepresentation() const {
  if (!hasQPdiag()) {
    throw std::runtime_error("Orbitals file does not contain QP coefficients");
  }
  return _mos.eigenvectors().block(0, _qpmin, _mos.eigenvectors().rows(),
                                   _qpmax - _qpmin + 1) *
         _QPdiag.eigenvectors();
}

// Determine QuasiParticle Density Matrix

Eigen::MatrixXd Orbitals::DensityMatrixQuasiParticle(
    const QMState& state) const {
  if (state.Type() != QMStateType::DQPstate) {
    throw std::runtime_error("State:" + state.ToString() +
                             " is not a quasiparticle state");
  }
  Eigen::MatrixXd lambda = CalculateQParticleAORepresentation();
  Eigen::MatrixXd dmatQP = lambda.col(state.Index() - _qpmin) *
                           lambda.col(state.Index() - _qpmin).transpose();
  return dmatQP;
}

Eigen::Vector3d Orbitals::CalcElDipole(const QMState& state) const {
  Eigen::Vector3d nuclei_dip = Eigen::Vector3d::Zero();
  if (!state.isTransition()) {
    for (const QMAtom& atom : _atoms) {
      nuclei_dip += (atom.getPos() - _atoms.getPos()) * atom.getNuccharge();
    }
  }
  AOBasis basis = SetupDftBasis();
  AODipole dipole;
  dipole.setCenter(_atoms.getPos());
  dipole.Fill(basis);

  Eigen::MatrixXd dmat = this->DensityMatrixFull(state);
  Eigen::Vector3d electronic_dip;
  for (int i = 0; i < 3; ++i) {
    electronic_dip(i) = dmat.cwiseProduct(dipole.Matrix()[i]).sum();
  }
  return nuclei_dip - electronic_dip;
}

Eigen::MatrixXd Orbitals::TransitionDensityMatrix(const QMState& state) const {
  if (state.Type() != QMStateType::Singlet) {
    throw std::runtime_error(
        "Spin type not known for transition density matrix. Available only for "
        "singlet");
  }
  const Eigen::MatrixXd& BSECoefs = _BSE_singlet.eigenvectors();
  if (BSECoefs.cols() < state.Index() + 1 || BSECoefs.rows() < 2) {
    throw std::runtime_error("Orbitals object has no information about state:" +
                             state.ToString());
  }

  // The Transition dipole is sqrt2 bigger because of the spin, the excited
  // state is a linear combination of 2 slater determinants, where either alpha
  // or beta spin electron is excited

  /*Trying to implement D_{alpha,beta}=
   * sqrt2*sum_{i}^{occ}sum_{j}^{virt}{BSEcoef(i,j)*MOcoef(alpha,i)*MOcoef(beta,j)}
   */
  // c stands for conduction band and thus virtual orbitals
  // v stand for valence band and thus occupied orbitals

  Eigen::VectorXd coeffs = BSECoefs.col(state.Index());

  if (!_useTDA) {
    coeffs += _BSE_singlet.eigenvectors2().col(state.Index());
  }
  coeffs *= std::sqrt(2.0);
  auto occlevels = _mos.eigenvectors().block(
      0, _bse_vmin, _mos.eigenvectors().rows(), _bse_vtotal);
  auto virtlevels = _mos.eigenvectors().block(
      0, _bse_cmin, _mos.eigenvectors().rows(), _bse_ctotal);
  Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), _bse_ctotal,
                                        _bse_vtotal);

  return occlevels * mat.transpose() * virtlevels.transpose();
}

std::array<Eigen::MatrixXd, 2> Orbitals::DensityMatrixExcitedState(
    const QMState& state) const {
  std::array<Eigen::MatrixXd, 2> dmat = DensityMatrixExcitedState_R(state);
  if (!_useTDA) {
    std::array<Eigen::MatrixXd, 2> dmat_AR =
        DensityMatrixExcitedState_AR(state);
    dmat[0] -= dmat_AR[0];
    dmat[1] -= dmat_AR[1];
  }
  return dmat;
}

// Excited state density matrix

std::array<Eigen::MatrixXd, 2> Orbitals::DensityMatrixExcitedState_R(
    const QMState& state) const {
  if (!state.Type().isExciton()) {
    throw std::runtime_error(
        "Spin type not known for density matrix. Available are singlet and "
        "triplet");
  }

  const Eigen::MatrixXd& BSECoefs = (state.Type() == QMStateType::Singlet)
                                        ? _BSE_singlet.eigenvectors()
                                        : _BSE_triplet.eigenvectors();
  if (BSECoefs.cols() < state.Index() + 1 || BSECoefs.rows() < 2) {
    throw std::runtime_error("Orbitals object has no information about state:" +
                             state.ToString());
  }
  /******
   *
   *    Density matrix for GW-BSE based excitations
   *
   *    - electron contribution
   *      D_ab = \sum{vc} \sum{c'} A_{vc}A_{vc'} mo_a(c)mo_b(c')
   *
   *    - hole contribution
   *      D_ab = \sum{vc} \sum{v'} A_{vc}A_{v'c} mo_a(v)mo_b(v')
   *
   */

  Eigen::VectorXd coeffs = BSECoefs.col(state.Index());

  std::array<Eigen::MatrixXd, 2> dmatEX;
  // hole part as matrix products
  Eigen::MatrixXd occlevels = _mos.eigenvectors().block(
      0, _bse_vmin, _mos.eigenvectors().rows(), _bse_vtotal);
  dmatEX[0] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

  // electron part as matrix products
  Eigen::MatrixXd virtlevels = _mos.eigenvectors().block(
      0, _bse_cmin, _mos.eigenvectors().rows(), _bse_ctotal);
  dmatEX[1] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();

  return dmatEX;
}

Eigen::MatrixXd Orbitals::CalcAuxMat_vv(const Eigen::VectorXd& coeffs) const {
  const Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), _bse_ctotal,
                                              _bse_vtotal);
  return mat.transpose() * mat;
}

Eigen::MatrixXd Orbitals::CalcAuxMat_cc(const Eigen::VectorXd& coeffs) const {
  const Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), _bse_ctotal,
                                              _bse_vtotal);
  return mat * mat.transpose();
}

std::array<Eigen::MatrixXd, 2> Orbitals::DensityMatrixExcitedState_AR(
    const QMState& state) const {

  if (!state.Type().isExciton()) {
    throw std::runtime_error(
        "Spin type not known for density matrix. Available are singlet and "
        "triplet");
  }

  const Eigen::MatrixXd& BSECoefs_AR = (state.Type() == QMStateType::Singlet)
                                           ? _BSE_singlet.eigenvectors2()
                                           : _BSE_triplet.eigenvectors2();
  if (BSECoefs_AR.cols() < state.Index() + 1 || BSECoefs_AR.rows() < 2) {
    throw std::runtime_error("Orbitals object has no information about state:" +
                             state.ToString());
  }
  /******
   *
   *    Density matrix for GW-BSE based excitations
   *
   *    - electron contribution
   *      D_ab = \sum{vc} \sum{v'} B_{vc}B_{v'c} mo_a(v)mo_b(v')
   *
   *    - hole contribution
   *      D_ab = \sum{vc} \sum{c'} B_{vc}B_{vc'} mo_a(c)mo_b(c')
   *
   *
   *   more efficient:
   *
   *   - electron contribution
   *      D_ab = \sum{v} \sum{v'} mo_a(v)mo_b(v') [ \sum{c} B_{vc}B_{v'c} ]
   *           = \sum{v} \sum{v'} mo_a(v)mo_b(v') B_{vv'}
   *
   *   - hole contribution
   *      D_ab = \sum{c} \sum{c'} mo_a(c)mo_b(c') [ \sum{v} B_{vc}B_{vc'} ]
   *           = \sum{c} \sum{c'} mo_a(c)mo_b(c') B_{cc'}
   *
   */

  Eigen::VectorXd coeffs = BSECoefs_AR.col(state.Index());

  std::array<Eigen::MatrixXd, 2> dmatAR;
  Eigen::MatrixXd virtlevels = _mos.eigenvectors().block(
      0, _bse_cmin, _mos.eigenvectors().rows(), _bse_ctotal);
  dmatAR[0] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();
  // electron part as matrix products
  Eigen::MatrixXd occlevels = _mos.eigenvectors().block(
      0, _bse_vmin, _mos.eigenvectors().rows(), _bse_vtotal);
  dmatAR[1] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

  return dmatAR;
}

Eigen::VectorXd Orbitals::Oscillatorstrengths() const {

  long size = long(_transition_dipoles.size());
  if (size > _BSE_singlet.eigenvalues().size()) {
    size = _BSE_singlet.eigenvalues().size();
  }
  Eigen::VectorXd oscs = Eigen::VectorXd::Zero(size);
  for (long i = 0; i < size; ++i) {
    oscs(i) = _transition_dipoles[i].squaredNorm() * 2.0 / 3.0 *
              (_BSE_singlet.eigenvalues()(i));
  }
  return oscs;
}

double Orbitals::getTotalStateEnergy(const QMState& state) const {
  double total_energy = getDFTTotalEnergy();
  if (state.Type() == QMStateType::Gstate) {
    return total_energy;
  }
  total_energy += getExcitedStateEnergy(state);
  return total_energy;
}

double Orbitals::getExcitedStateEnergy(const QMState& state) const {

  double omega = 0.0;
  if (state.isTransition()) {
    throw std::runtime_error(
        "Total Energy does not exist for transition state");
  }

  if (state.Type() == QMStateType::Singlet) {
    if (_BSE_singlet.eigenvalues().size() < state.Index() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    omega = _BSE_singlet.eigenvalues()[state.Index()];
  } else if (state.Type() == QMStateType::Triplet) {
    if (_BSE_triplet.eigenvalues().size() < state.Index() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    omega = _BSE_triplet.eigenvalues()[state.Index()];
  } else if (state.Type() == QMStateType::DQPstate) {
    if (_QPdiag.eigenvalues().size() < state.Index() + 1 - getGWAmin()) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return _QPdiag.eigenvalues()[state.Index() - getGWAmin()];
  } else if (state.Type() == QMStateType::KSstate) {
    if (_mos.eigenvalues().size() < state.Index() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return _mos.eigenvalues()[state.Index()];
  } else if (state.Type() == QMStateType::PQPstate) {
    if (this->_QPpert_energies.rows() < state.Index() + 1 - getGWAmin()) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return _QPpert_energies(state.Index() - getGWAmin(), 3);
  } else {
    throw std::runtime_error(
        "GetTotalEnergy only knows states:singlet,triplet,KS,DQP,PQP");
  }
  return omega;  //  e.g. hartree
}

std::array<Eigen::MatrixXd, 3> Orbitals::CalcFreeTransition_Dipoles() const {
  const Eigen::MatrixXd& dft_orbitals = _mos.eigenvectors();
  AOBasis basis = SetupDftBasis();
  // Testing electric dipole AOMatrix
  AODipole dft_dipole;
  dft_dipole.Fill(basis);

  // now transition dipole elements for free interlevel transitions
  std::array<Eigen::MatrixXd, 3> interlevel_dipoles;

  Eigen::MatrixXd empty =
      dft_orbitals.block(0, _bse_cmin, basis.AOBasisSize(), _bse_ctotal);
  Eigen::MatrixXd occ =
      dft_orbitals.block(0, _bse_vmin, basis.AOBasisSize(), _bse_vtotal);
  for (int i = 0; i < 3; i++) {
    interlevel_dipoles[i] = empty.transpose() * dft_dipole.Matrix()[i] * occ;
  }
  return interlevel_dipoles;
}

void Orbitals::CalcCoupledTransition_Dipoles() {
  std::array<Eigen::MatrixXd, 3> interlevel_dipoles =
      CalcFreeTransition_Dipoles();
  long numofstates = _BSE_singlet.eigenvalues().size();
  _transition_dipoles.resize(0);
  _transition_dipoles.reserve(numofstates);
  const double sqrt2 = sqrt(2.0);
  for (long i_exc = 0; i_exc < numofstates; i_exc++) {

    Eigen::VectorXd coeffs = _BSE_singlet.eigenvectors().col(i_exc);
    if (!_useTDA) {
      coeffs += _BSE_singlet.eigenvectors2().col(i_exc);
    }
    Eigen::Map<Eigen::MatrixXd> mat(coeffs.data(), _bse_ctotal, _bse_vtotal);
    Eigen::Vector3d tdipole = Eigen::Vector3d::Zero();
    for (int i = 0; i < 3; i++) {
      tdipole[i] = mat.cwiseProduct(interlevel_dipoles[i]).sum();
    }
    // The Transition dipole is sqrt2 bigger because of the spin, the
    // excited state is a linear combination of 2 slater determinants,
    // where either alpha or beta spin electron is excited
    _transition_dipoles.push_back(-sqrt2 * tdipole);  //- because electrons are
                                                      // negative
  }
}

void Orbitals::OrderMOsbyEnergy() {
  std::vector<int> sort_index = SortEnergies();
  tools::EigenSystem MO_copy = _mos;
  long size = _mos.eigenvalues().size();
  for (long i = 0; i < size; ++i) {
    _mos.eigenvalues()(i) = MO_copy.eigenvalues()(sort_index[i]);
  }
  for (long i = 0; i < size; ++i) {
    _mos.eigenvectors().col(i) = MO_copy.eigenvectors().col(sort_index[i]);
  }
}

/**
 * \brief Guess for a dimer based on monomer orbitals
 *
 * Given two monomer orbitals (A and B) constructs a guess for dimer
 * orbitals: | A 0 | and energies: [EA, EB]
 *           | 0 B |
 */
void Orbitals::PrepareDimerGuess(const Orbitals& orbitalsA,
                                 const Orbitals& orbitalsB) {

  // constructing the direct product orbA x orbB
  int basisA = orbitalsA.getBasisSetSize();
  int basisB = orbitalsB.getBasisSetSize();

  int electronsA = orbitalsA.getNumberOfAlphaElectrons();
  int electronsB = orbitalsB.getNumberOfAlphaElectrons();

  _mos.eigenvectors() = Eigen::MatrixXd::Zero(basisA + basisB, basisA + basisB);

  // AxB = | A 0 |  //   A = [EA, EB]  //
  //       | 0 B |  //                 //
  if (orbitalsA.getDFTbasisName() != orbitalsB.getDFTbasisName()) {
    throw std::runtime_error("Basissets of Orbitals A and B differ " +
                             orbitalsA.getDFTbasisName() + ":" +
                             orbitalsB.getDFTbasisName());
  }
  this->setDFTbasisName(orbitalsA.getDFTbasisName());
  if (orbitalsA.getECPName() != orbitalsB.getECPName()) {
    throw std::runtime_error("ECPs of Orbitals A and B differ " +
                             orbitalsA.getECPName() + ":" +
                             orbitalsB.getECPName());
  }
  this->setECPName(orbitalsA.getECPName());
  this->setBasisSetSize(basisA + basisB);
  this->setNumberOfOccupiedLevels(electronsA + electronsB);
  this->setNumberOfAlphaElectrons(electronsA + electronsB);

  _mos.eigenvectors().topLeftCorner(basisA, basisA) =
      orbitalsA.MOs().eigenvectors();
  _mos.eigenvectors().bottomRightCorner(basisB, basisB) =
      orbitalsB.MOs().eigenvectors();

  _mos.eigenvalues().resize(basisA + basisB);

  _mos.eigenvalues().head(basisA) = orbitalsA.MOs().eigenvalues();
  _mos.eigenvalues().tail(basisB) = orbitalsB.MOs().eigenvalues();

  OrderMOsbyEnergy();

  return;
}

void Orbitals::WriteToCpt(const std::string& filename) const {
  CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
  WriteToCpt(cpf);
}

void Orbitals::WriteToCpt(CheckpointFile f) const {
  WriteToCpt(f.getWriter("/QMdata"));
}

void Orbitals::WriteToCpt(CheckpointWriter w) const {
  w(XtpVersionStr(), "Version");
  w(_basis_set_size, "basis_set_size");
  w(_occupied_levels, "occupied_levels");
  w(_number_alpha_electrons, "number_alpha_electrons");

  w(_mos, "mos");

  CheckpointWriter molgroup = w.openChild("qmmolecule");
  _atoms.WriteToCpt(molgroup);

  w(_qm_energy, "qm_energy");
  w(_qm_package, "qm_package");

  w(_dftbasis, "dftbasis");
  w(_auxbasis, "auxbasis");

  w(_rpamin, "rpamin");
  w(_rpamax, "rpamax");
  w(_qpmin, "qpmin");
  w(_qpmax, "qpmax");
  w(_bse_vmin, "bse_vmin");
  w(_bse_cmax, "bse_cmax");
  w(_functionalname, "XCFunctional");
  w(_ScaHFX, "ScaHFX");

  w(_useTDA, "useTDA");
  w(_ECP, "ECP");

  w(_QPpert_energies, "QPpert_energies");

  w(_QPdiag, "QPdiag");

  w(_BSE_singlet, "BSE_singlet");

  w(_transition_dipoles, "transition_dipoles");

  w(_BSE_triplet, "BSE_triplet");
}

void Orbitals::ReadFromCpt(const std::string& filename) {
  CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
  ReadFromCpt(cpf);
}

void Orbitals::ReadFromCpt(CheckpointFile f) {
  ReadFromCpt(f.getReader("/QMdata"));
}

void Orbitals::ReadFromCpt(CheckpointReader r) {
  r(_basis_set_size, "basis_set_size");
  r(_occupied_levels, "occupied_levels");
  r(_number_alpha_electrons, "number_alpha_electrons");

  r(_mos, "mos");

  // Read qmatoms
  CheckpointReader molgroup = r.openChild("qmmolecule");
  _atoms.ReadFromCpt(molgroup);

  r(_qm_energy, "qm_energy");
  r(_qm_package, "qm_package");

  r(_dftbasis, "dftbasis");
  r(_auxbasis, "auxbasis");

  r(_rpamin, "rpamin");
  r(_rpamax, "rpamax");
  r(_qpmin, "qpmin");
  r(_qpmax, "qpmax");
  r(_bse_vmin, "bse_vmin");
  r(_bse_cmax, "bse_cmax");
  setBSEindices(_bse_vmin, _bse_cmax);
  try {
    r(_functionalname, "XCFunctional");
  } catch (std::runtime_error& e) {
    ;
  }
  r(_ScaHFX, "ScaHFX");
  r(_useTDA, "useTDA");
  r(_ECP, "ECP");

  r(_QPpert_energies, "QPpert_energies");
  r(_QPdiag, "QPdiag");

  r(_BSE_singlet, "BSE_singlet");

  r(_transition_dipoles, "transition_dipoles");

  r(_BSE_triplet, "BSE_triplet");
}
}  // namespace xtp
}  // namespace votca
