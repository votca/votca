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

// Standard includes
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <numeric>

// Local VOTCA includes
#include "votca/tools/version.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

Orbitals::Orbitals() : atoms_("", 0) { ; }

/**
 *
 * @param  level Index of the level that is to be checked for degeneracy
 * @param  energy_difference [ev] Two levels are degenerate if their energy is
 * smaller than this value
 * @return vector with indices off all orbitals degenerate to this including
 * itself
 */
std::vector<Index> Orbitals::CheckDegeneracy(Index level,
                                             double energy_difference) const {

  std::vector<Index> result;
  if (level > mos_.eigenvalues().size()) {
    throw std::runtime_error(
        "Level for degeneracy is higher than maximum level");
  }
  double MOEnergyLevel = mos_.eigenvalues()(level);

  for (Index i = 0; i < mos_.eigenvalues().size(); ++i) {
    if (std::abs(mos_.eigenvalues()(i) - MOEnergyLevel) < energy_difference) {
      result.push_back(i);
    }
  }

  if (result.empty()) {
    result.push_back(level);
  }
  return result;
}

std::vector<Index> Orbitals::SortEnergies() {
  std::vector<Index> index = std::vector<Index>(mos_.eigenvalues().size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [this](Index i1, Index i2) {
    return this->MOs().eigenvalues()[i1] < this->MOs().eigenvalues()[i2];
  });
  return index;
}

/**
 * SetupDftBasis constructs the dft basis, to do this the overlap integral needs
 * to be evaluated with libint. Hence libint should be initialized for it to
 * work.
 */
void Orbitals::SetupDftBasis(std::string basis_name) {
  if (this->QMAtoms().size() == 0) {
    throw std::runtime_error("Can't setup AOBasis without atoms");
  }
  BasisSet bs;
  bs.Load(basis_name);
  dftbasis_.Fill(bs, this->QMAtoms());
}

void Orbitals::SetupAuxBasis(std::string aux_basis_name) {
  if (this->QMAtoms().size() == 0) {
    throw std::runtime_error("Can't setup Aux AOBasis without atoms");
  }
  BasisSet bs;
  bs.Load(aux_basis_name);
  auxbasis_.Fill(bs, this->QMAtoms());
}

/*
 * Returns the density matrix relative to the ground state, for the full density
 * use DensityMatrixFull
 */
Eigen::MatrixXd Orbitals::DensityMatrixWithoutGS(const QMState& state) const {
  if (state.Type().isExciton()) {
    std::array<Eigen::MatrixXd, 2> DMAT = DensityMatrixExcitedState(state);
    return DMAT[1] - DMAT[0];
  } else if (state.Type().isKSState() || state.Type().isPQPState()) {
    return DensityMatrixKSstate(state);
  } else if (state.Type() == QMStateType::DQPstate) {
    Eigen::MatrixXd DMATQP = DensityMatrixQuasiParticle(state);
    if (state.StateIdx() > getHomo()) {
      return DMATQP;
    } else {
      return -DMATQP;
    }
  } else {
    throw std::runtime_error(
        "DensityMatrixWithoutGS does not yet implement QMStateType:" +
        state.Type().ToLongString());
  }
}

/*
 * Returns the density matrix with the ground state density, for the partial
 * density relative to the ground state use DensityMatrixWithoutGS
 */
Eigen::MatrixXd Orbitals::DensityMatrixFull(const QMState& state) const {
  if (state.isTransition()) {
    return this->TransitionDensityMatrix(state);
  }
  Eigen::MatrixXd result = this->DensityMatrixGroundState();;
  if (getCalculationType() != "") {
    result += getInactiveDensity();
  }
  std::cout << "DID YOU GET HERE? " << std::endl;
  for (Index atom = 0; atom < QMAtoms().size(); atom++) {
    std::cout << atom << ": " << QMAtoms()[atom] << std::endl;
  }

  AOBasis aobasis = getDftBasis();
  AOOverlap overlap;
  overlap.Fill(aobasis);
  double elec = result.cwiseProduct(overlap.Matrix()).sum();
  if (getCalculationType() == "Truncated") {
    double inact = getInactiveDensity().cwiseProduct(overlap.Matrix()).sum();
    std::cout << "Electrons after this point are: " << elec << " & " << inact;
  }
  else{
    std::cout << "Electrons after this point are: " << elec;
  }
  

  if (state.Type().isExciton()) {
    std::array<Eigen::MatrixXd, 2> DMAT = DensityMatrixExcitedState(state);
    result = result - DMAT[0] + DMAT[1];  // Ground state + hole_contribution +
                                          // electron contribution
  } else if (state.Type() == QMStateType::DQPstate) {
    Eigen::MatrixXd DMATQP = DensityMatrixQuasiParticle(state);
    if (state.StateIdx() > getHomo()) {
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
  Eigen::MatrixXd occstates = mos_.eigenvectors().leftCols(occupied_levels_);
  Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
  return dmatGS;
}

// Eigen::MatrixXd Orbitals::EmbDensityMatrixGroundState() const {
//   Eigen::MatrixXd occstates =
//       mos_embedding_.eigenvectors().leftCols(active_electrons_ / 2);
//   Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
//   return dmatGS;
// }

// Eigen::MatrixXd Orbitals::TruncDensityMatrixGroundState() const {
//   Eigen::MatrixXd occstates = expandedMOs_.leftCols(active_electrons_ / 2);
//   Eigen::MatrixXd dmatGS = 2.0 * occstates * occstates.transpose();
//   return dmatGS;
// }

// Density matrix for a single KS orbital
Eigen::MatrixXd Orbitals::DensityMatrixKSstate(const QMState& state) const {
  if (!hasMOs()) {
    throw std::runtime_error("Orbitals file does not contain MO coefficients");
  }
  if (state.Type() != QMStateType::KSstate &&
      state.Type() != QMStateType::PQPstate) {
    throw std::runtime_error("State:" + state.ToString() +
                             " is not a Kohn Sham state");
  }
  Eigen::VectorXd KSstate = mos_.eigenvectors().col(state.StateIdx());
  Eigen::MatrixXd dmatKS = KSstate * KSstate.transpose();
  return dmatKS;
}

Eigen::MatrixXd Orbitals::CalculateQParticleAORepresentation() const {
  if (!hasQPdiag()) {
    throw std::runtime_error("Orbitals file does not contain QP coefficients");
  }
  return mos_.eigenvectors().middleCols(qpmin_, qpmax_ - qpmin_ + 1) *
         QPdiag_.eigenvectors();
}

// Determine QuasiParticle Density Matrix
Eigen::MatrixXd Orbitals::DensityMatrixQuasiParticle(
    const QMState& state) const {
  if (state.Type() != QMStateType::DQPstate) {
    throw std::runtime_error("State:" + state.ToString() +
                             " is not a quasiparticle state");
  }
  Eigen::MatrixXd lambda = CalculateQParticleAORepresentation();
  Eigen::MatrixXd dmatQP = lambda.col(state.StateIdx() - qpmin_) *
                           lambda.col(state.StateIdx() - qpmin_).transpose();
  return dmatQP;
}

Eigen::Vector3d Orbitals::CalcElDipole(const QMState& state) const {
  Eigen::Vector3d nuclei_dip = Eigen::Vector3d::Zero();
  if (!state.isTransition()) {
    for (const QMAtom& atom : atoms_) {
      nuclei_dip += (atom.getPos() - atoms_.getPos()) * atom.getNuccharge();
    }
  }
  AOBasis basis = getDftBasis();
  AODipole dipole;
  dipole.setCenter(atoms_.getPos());
  dipole.Fill(basis);

  Eigen::MatrixXd dmat = this->DensityMatrixFull(state);
  Eigen::Vector3d electronic_dip;
  for (Index i = 0; i < 3; ++i) {
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
  const Eigen::MatrixXd& BSECoefs = BSE_singlet_.eigenvectors();
  if (BSECoefs.cols() < state.StateIdx() + 1 || BSECoefs.rows() < 2) {
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

  Eigen::VectorXd coeffs = BSECoefs.col(state.StateIdx());

  if (!useTDA_) {
    coeffs += BSE_singlet_.eigenvectors2().col(state.StateIdx());
  }
  coeffs *= std::sqrt(2.0);
  auto occlevels = mos_.eigenvectors().middleCols(bse_vmin_, bse_vtotal_);
  auto virtlevels = mos_.eigenvectors().middleCols(bse_cmin_, bse_ctotal_);
  Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), bse_ctotal_,
                                        bse_vtotal_);

  return occlevels * mat.transpose() * virtlevels.transpose();
}

std::array<Eigen::MatrixXd, 2> Orbitals::DensityMatrixExcitedState(
    const QMState& state) const {
  std::array<Eigen::MatrixXd, 2> dmat = DensityMatrixExcitedState_R(state);
  if (!useTDA_) {
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
                                        ? BSE_singlet_.eigenvectors()
                                        : BSE_triplet_.eigenvectors();
  if (BSECoefs.cols() < state.StateIdx() + 1 || BSECoefs.rows() < 2) {
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

  Eigen::VectorXd coeffs = BSECoefs.col(state.StateIdx());

  std::array<Eigen::MatrixXd, 2> dmatEX;
  // hole part as matrix products
  Eigen::MatrixXd occlevels =
      mos_.eigenvectors().middleCols(bse_vmin_, bse_vtotal_);
  dmatEX[0] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

  // electron part as matrix products
  Eigen::MatrixXd virtlevels =
      mos_.eigenvectors().middleCols(bse_cmin_, bse_ctotal_);
  dmatEX[1] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();

  return dmatEX;
}

Eigen::MatrixXd Orbitals::CalcAuxMat_vv(const Eigen::VectorXd& coeffs) const {
  const Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), bse_ctotal_,
                                              bse_vtotal_);
  return mat.transpose() * mat;
}

Eigen::MatrixXd Orbitals::CalcAuxMat_cc(const Eigen::VectorXd& coeffs) const {
  const Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), bse_ctotal_,
                                              bse_vtotal_);
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
                                           ? BSE_singlet_.eigenvectors2()
                                           : BSE_triplet_.eigenvectors2();
  if (BSECoefs_AR.cols() < state.StateIdx() + 1 || BSECoefs_AR.rows() < 2) {
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

  Eigen::VectorXd coeffs = BSECoefs_AR.col(state.StateIdx());

  std::array<Eigen::MatrixXd, 2> dmatAR;
  Eigen::MatrixXd virtlevels =
      mos_.eigenvectors().middleCols(bse_cmin_, bse_ctotal_);
  dmatAR[0] = virtlevels * CalcAuxMat_cc(coeffs) * virtlevels.transpose();
  // electron part as matrix products
  Eigen::MatrixXd occlevels =
      mos_.eigenvectors().middleCols(bse_vmin_, bse_vtotal_);
  dmatAR[1] = occlevels * CalcAuxMat_vv(coeffs) * occlevels.transpose();

  return dmatAR;
}

Eigen::VectorXd Orbitals::Oscillatorstrengths() const {

  Index size = Index(transition_dipoles_.size());
  if (size > BSE_singlet_.eigenvalues().size()) {
    size = BSE_singlet_.eigenvalues().size();
  }
  Eigen::VectorXd oscs = Eigen::VectorXd::Zero(size);
  for (Index i = 0; i < size; ++i) {
    oscs(i) = transition_dipoles_[i].squaredNorm() * 2.0 / 3.0 *
              (BSE_singlet_.eigenvalues()(i));
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
    if (BSE_singlet_.eigenvalues().size() < state.StateIdx() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    omega = BSE_singlet_.eigenvalues()[state.StateIdx()];
  } else if (state.Type() == QMStateType::Triplet) {
    if (BSE_triplet_.eigenvalues().size() < state.StateIdx() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    omega = BSE_triplet_.eigenvalues()[state.StateIdx()];
  } else if (state.Type() == QMStateType::DQPstate) {
    if (QPdiag_.eigenvalues().size() < state.StateIdx() + 1 - getGWAmin()) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return QPdiag_.eigenvalues()[state.StateIdx() - getGWAmin()];
  } else if (state.Type() == QMStateType::KSstate) {
    if (mos_.eigenvalues().size() < state.StateIdx() + 1) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return mos_.eigenvalues()[state.StateIdx()];
  } else if (state.Type() == QMStateType::PQPstate) {
    if (this->QPpert_energies_.rows() < state.StateIdx() + 1 - getGWAmin()) {
      throw std::runtime_error("Orbitals::getTotalEnergy You want " +
                               state.ToString() +
                               " which has not been calculated");
    }
    return QPpert_energies_(state.StateIdx() - getGWAmin(), 3);
  } else if (state.Type() == QMStateType::LMOstate) {
    if (lmos_energies_.size() < state.StateIdx() + 1) {
      throw std::runtime_error(
          "Orbitals::getTotalEnergy You want " + state.ToString() +
          " which is a LMO for virtual orbitals. Not implemented.");
    }
    return lmos_energies_(state.StateIdx());
  } else {
    throw std::runtime_error(
        "GetTotalEnergy only knows states:singlet,triplet,KS,DQP,PQP,LMOs");
  }
  return omega;  //  e.g. hartree
}

std::array<Eigen::MatrixXd, 3> Orbitals::CalcFreeTransition_Dipoles() const {
  const Eigen::MatrixXd& dft_orbitals = mos_.eigenvectors();
  AOBasis basis = getDftBasis();
  // Testing electric dipole AOMatrix
  AODipole dft_dipole;
  dft_dipole.Fill(basis);

  // now transition dipole elements for free interlevel transitions
  std::array<Eigen::MatrixXd, 3> interlevel_dipoles;

  Eigen::MatrixXd empty = dft_orbitals.middleCols(bse_cmin_, bse_ctotal_);
  Eigen::MatrixXd occ = dft_orbitals.middleCols(bse_vmin_, bse_vtotal_);
  for (Index i = 0; i < 3; i++) {
    interlevel_dipoles[i] = empty.transpose() * dft_dipole.Matrix()[i] * occ;
  }
  return interlevel_dipoles;
}

void Orbitals::CalcCoupledTransition_Dipoles() {
  std::array<Eigen::MatrixXd, 3> interlevel_dipoles =
      CalcFreeTransition_Dipoles();
  Index numofstates = BSE_singlet_.eigenvalues().size();
  transition_dipoles_.resize(0);
  transition_dipoles_.reserve(numofstates);
  const double sqrt2 = std::sqrt(2.0);
  for (Index i_exc = 0; i_exc < numofstates; i_exc++) {

    Eigen::VectorXd coeffs = BSE_singlet_.eigenvectors().col(i_exc);
    if (!useTDA_) {
      coeffs += BSE_singlet_.eigenvectors2().col(i_exc);
    }
    Eigen::Map<Eigen::MatrixXd> mat(coeffs.data(), bse_ctotal_, bse_vtotal_);
    Eigen::Vector3d tdipole = Eigen::Vector3d::Zero();
    for (Index i = 0; i < 3; i++) {
      tdipole[i] = mat.cwiseProduct(interlevel_dipoles[i]).sum();
    }
    // The Transition dipole is sqrt2 bigger because of the spin, the
    // excited state is a linear combination of 2 slater determinants,
    // where either alpha or beta spin electron is excited
    transition_dipoles_.push_back(-sqrt2 * tdipole);  //- because electrons are
                                                      // negatively charged
  }
}

void Orbitals::OrderMOsbyEnergy() {
  std::vector<Index> sort_index = SortEnergies();
  tools::EigenSystem MO_copy = mos_;
  Index size = mos_.eigenvalues().size();
  for (Index i = 0; i < size; ++i) {
    mos_.eigenvalues()(i) = MO_copy.eigenvalues()(sort_index[i]);
  }
  for (Index i = 0; i < size; ++i) {
    mos_.eigenvectors().col(i) = MO_copy.eigenvectors().col(sort_index[i]);
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
  Index basisA = orbitalsA.getBasisSetSize();
  Index basisB = orbitalsB.getBasisSetSize();

  Index electronsA = orbitalsA.getNumberOfAlphaElectrons();
  Index electronsB = orbitalsB.getNumberOfAlphaElectrons();

  mos_.eigenvectors() = Eigen::MatrixXd::Zero(basisA + basisB, basisA + basisB);

  // AxB = | A 0 |  //   A = [EA, EB]  //
  //       | 0 B |  //                 //
  if (orbitalsA.getDFTbasisName() != orbitalsB.getDFTbasisName()) {
    throw std::runtime_error("Basissets of Orbitals A and B differ " +
                             orbitalsA.getDFTbasisName() + ":" +
                             orbitalsB.getDFTbasisName());
  }
  this->SetupDftBasis(orbitalsA.getDFTbasisName());
  if (orbitalsA.getECPName() != orbitalsB.getECPName()) {
    throw std::runtime_error("ECPs of Orbitals A and B differ " +
                             orbitalsA.getECPName() + ":" +
                             orbitalsB.getECPName());
  }
  this->setECPName(orbitalsA.getECPName());
  this->setNumberOfOccupiedLevels(electronsA + electronsB);
  this->setNumberOfAlphaElectrons(electronsA + electronsB);

  mos_.eigenvectors().topLeftCorner(basisA, basisA) =
      orbitalsA.MOs().eigenvectors();
  mos_.eigenvectors().bottomRightCorner(basisB, basisB) =
      orbitalsB.MOs().eigenvectors();

  mos_.eigenvalues().resize(basisA + basisB);

  mos_.eigenvalues().head(basisA) = orbitalsA.MOs().eigenvalues();
  mos_.eigenvalues().tail(basisB) = orbitalsB.MOs().eigenvalues();

  OrderMOsbyEnergy();

  return;
}

void Orbitals::WriteToCpt(const std::string& filename) const {
  CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
  WriteToCpt(cpf);
}

void Orbitals::WriteToCpt(CheckpointFile f) const {
  CheckpointWriter writer = f.getWriter("/QMdata");
  WriteToCpt(writer);
  WriteBasisSetsToCpt(writer);
}

void Orbitals::WriteBasisSetsToCpt(CheckpointWriter w) const {
  CheckpointWriter dftWriter = w.openChild("dft");
  dftbasis_.WriteToCpt(dftWriter);
  CheckpointWriter auxWriter = w.openChild("aux");
  auxbasis_.WriteToCpt(auxWriter);
}

void Orbitals::WriteToCpt(CheckpointWriter w) const {
  w(votca::tools::ToolsVersionStr(), "XTPVersion");
  w(orbitals_version(), "version");
  w(occupied_levels_, "occupied_levels");
  w(number_alpha_electrons_, "number_alpha_electrons");

  w(mos_, "mos");
  w(active_electrons_, "active_electrons");
  w(mos_embedding_, "mos_embedding");
  w(lmos_, "LMOs");
  w(lmos_energies_, "LMOs_energies");
  w(inactivedensity_, "inactivedensity");
  w(expandedMOs_, "TruncMOsFullBasis");

  CheckpointWriter molgroup = w.openChild("qmmolecule");
  atoms_.WriteToCpt(molgroup);

  w(qm_energy_, "qm_energy");
  w(qm_package_, "qm_package");

  w(rpamin_, "rpamin");
  w(rpamax_, "rpamax");
  w(qpmin_, "qpmin");
  w(qpmax_, "qpmax");
  w(bse_vmin_, "bse_vmin");
  w(bse_cmax_, "bse_cmax");
  w(functionalname_, "XCFunctional");
  w(grid_quality_, "XC_grid_quality");
  w(ScaHFX_, "ScaHFX");

  w(useTDA_, "useTDA");
  w(ECP_, "ECP");

  w(rpa_inputenergies_, "RPA_inputenergies");
  w(QPpert_energies_, "QPpert_energies");

  w(QPdiag_, "QPdiag");

  w(BSE_singlet_, "BSE_singlet");

  w(transition_dipoles_, "transition_dipoles");

  w(BSE_triplet_, "BSE_triplet");

  w(use_Hqp_offdiag_, "use_Hqp_offdiag");

  w(BSE_singlet_energies_dynamic_, "BSE_singlet_dynamic");

  w(BSE_triplet_energies_dynamic_, "BSE_triplet_dynamic");

  w(CalcType_, "CalcType");
}

void Orbitals::ReadFromCpt(const std::string& filename) {
  CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
  ReadFromCpt(cpf);
}

void Orbitals::ReadFromCpt(CheckpointFile f) {
  CheckpointReader reader = f.getReader("/QMdata");
  ReadFromCpt(reader);
  ReadBasisSetsFromCpt(reader);
}

void Orbitals::ReadBasisSetsFromCpt(CheckpointReader r) {
  CheckpointReader dftReader = r.openChild("dft");
  dftbasis_.ReadFromCpt(dftReader);
  CheckpointReader auxReader = r.openChild("aux");
  auxbasis_.ReadFromCpt(auxReader);
}

void Orbitals::ReadFromCpt(CheckpointReader r) {
  r(occupied_levels_, "occupied_levels");
  r(number_alpha_electrons_, "number_alpha_electrons");
  int version;
  r(version, "version");
  // Read qmatoms
  CheckpointReader molgroup = r.openChild("qmmolecule");
  atoms_.ReadFromCpt(molgroup);

  r(qm_energy_, "qm_energy");
  r(qm_package_, "qm_package");
  try {
    r(lmos_, "LMOs");
    r(lmos_energies_, "LMOs_energies");
  } catch (std::runtime_error& e) {
    ;
  }

  r(version, "version");
  r(mos_, "mos");
  r(mos_embedding_, "mos_embedding");
  r(active_electrons_, "active_electrons");
  r(inactivedensity_, "inactivedensity");
  r(CalcType_, "CalcType");
  r(expandedMOs_, "TruncMOsFullBasis");

  if (version < 3) {
    // clang-format off
    std::array<Index, 49> votcaOrder_old = {
        0,                             // s
        0, -1, 1,                      // p
        0, -1, 1, -2, 2,               // d
        0, -1, 1, -2, 2, -3, 3,        // f
        0, -1, 1, -2, 2, -3, 3, -4, 4,  // g
        0, -1, 1, -2, 2, -3, 3, -4, 4,-5,5,  // h
        0, -1, 1, -2, 2, -3, 3, -4, 4,-5,5,-6,6  // i
    };
    // clang-format on

    std::array<Index, 49> multiplier;
    multiplier.fill(1);
    OrbReorder ord(votcaOrder_old, multiplier);
    ord.reorderOrbitals(mos_.eigenvectors(), this->getDftBasis());
  }

  if (version < 5) {  // we need to construct the basissets, NB. can only be
                      // done after reading the atoms.
    std::string dft_basis_name;
    std::string aux_basis_name;
    r(dft_basis_name, "dftbasis");
    r(aux_basis_name, "auxbasis");
    this->SetupDftBasis(dft_basis_name);
    this->SetupAuxBasis(aux_basis_name);
  }

  r(rpamin_, "rpamin");
  r(rpamax_, "rpamax");
  r(qpmin_, "qpmin");
  r(qpmax_, "qpmax");
  r(bse_vmin_, "bse_vmin");
  r(bse_cmax_, "bse_cmax");
  setBSEindices(bse_vmin_, bse_cmax_);
  try {
    r(functionalname_, "XCFunctional");
    r(grid_quality_, "XC_grid_quality");
  } catch (std::runtime_error& e) {
    grid_quality_ = "medium";
  }
  r(ScaHFX_, "ScaHFX");
  r(useTDA_, "useTDA");
  r(ECP_, "ECP");

  r(rpa_inputenergies_, "RPA_inputenergies");
  r(QPpert_energies_, "QPpert_energies");
  r(QPdiag_, "QPdiag");

  r(BSE_singlet_, "BSE_singlet");

  r(transition_dipoles_, "transition_dipoles");

  r(BSE_triplet_, "BSE_triplet");

  r(use_Hqp_offdiag_, "use_Hqp_offdiag");

  if (version > 1) {
    r(BSE_singlet_energies_dynamic_, "BSE_singlet_dynamic");

    r(BSE_triplet_energies_dynamic_, "BSE_triplet_dynamic");
  }
}
}  // namespace xtp
}  // namespace votca
