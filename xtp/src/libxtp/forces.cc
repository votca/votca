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
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/forces.h"
#include "votca/xtp/statetracker.h"

namespace votca {
namespace xtp {

using std::flush;
void Forces::Initialize(tools::Property& options) {
  force_method_ = options.get(".method").as<std::string>();

  displacement_ = options.get(".displacement").as<double>();  // Angstrom
  displacement_ *= tools::conv::ang2bohr;

  remove_total_force_ = options.get(".CoMforce_removal").as<bool>();
}

void Forces::Calculate(const Orbitals& orbitals) {

  Index natoms = orbitals.QMAtoms().size();
  forces_ = Eigen::MatrixX3d::Zero(natoms, 3);

  for (Index atom_index = 0; atom_index < natoms; atom_index++) {

    XTP_LOG(Log::debug, *pLog_)
        << "FORCES--DEBUG working on atom " << atom_index << flush;
    Eigen::Vector3d atom_force = Eigen::Vector3d::Zero();
    // Calculate Force on this atom
    if (force_method_ == "forward") {
      atom_force = NumForceForward(orbitals, atom_index);
    }
    if (force_method_ == "central") {
      atom_force = NumForceCentral(orbitals, atom_index);
    }
    forces_.row(atom_index) = atom_force.transpose();
  }
  if (remove_total_force_) {
    RemoveTotalForce();
  }
  return;
}

void Forces::Report() const {

  XTP_LOG(Log::error, *pLog_)
      << (boost::format(" ---- FORCES (Hartree/Bohr)   ")).str() << flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("      %1$s differences   ") % force_method_).str()
      << flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format("      displacement %1$1.4f Angstrom   ") %
          (displacement_ * tools::conv::bohr2ang))
             .str()
      << flush;
  XTP_LOG(Log::error, *pLog_)
      << (boost::format(" Atom\t x\t  y\t  z ")).str() << flush;

  for (Index i = 0; i < forces_.rows(); i++) {
    XTP_LOG(Log::error, *pLog_)
        << (boost::format("%1$4d    %2$+1.4f  %3$+1.4f  %4$+1.4f") % i %
            forces_(i, 0) % forces_(i, 1) % forces_(i, 2))
               .str()
        << flush;
  }
  return;
}

Eigen::Vector3d Forces::NumForceForward(Orbitals orbitals, Index atom_index) {
  Eigen::Vector3d force = Eigen::Vector3d::Zero();
  // get this atoms's current coordinates
  double energy_center =
      orbitals.getTotalStateEnergy(tracker_.CalcState(orbitals));
  const Eigen::Vector3d current_pos = orbitals.QMAtoms()[atom_index].getPos();
  for (Index i_cart = 0; i_cart < 3; i_cart++) {
    Eigen::Vector3d displacement_vec = Eigen::Vector3d::Zero();
    displacement_vec[i_cart] = displacement_;
    // update the coordinate
    Eigen::Vector3d pos_displaced = current_pos + displacement_vec;
    orbitals.updateAtomPostion(atom_index, pos_displaced);
    gwbse_engine_.ExcitationEnergies(orbitals);
    double energy_displaced =
        orbitals.getTotalStateEnergy(tracker_.CalcState(orbitals));
    force(i_cart) = (energy_center - energy_displaced) / displacement_;
    orbitals.updateAtomPostion(
        atom_index, current_pos);  // restore original coordinate into segment
  }                                // Cartesian directions
  return force;
}

Eigen::Vector3d Forces::NumForceCentral(Orbitals orbitals, Index atom_index) {
  Eigen::Vector3d force = Eigen::Vector3d::Zero();
  const Eigen::Vector3d current_pos = orbitals.QMAtoms()[atom_index].getPos();
  for (Index i_cart = 0; i_cart < 3; i_cart++) {
    XTP_LOG(Log::debug, *pLog_)
        << "FORCES--DEBUG           Cartesian component " << i_cart << flush;
    Eigen::Vector3d displacement_vec = Eigen::Vector3d::Zero();
    displacement_vec[i_cart] = displacement_;
    // update the coordinate
    Eigen::Vector3d pos_displaced = current_pos + displacement_vec;
    orbitals.updateAtomPostion(atom_index, pos_displaced);
    gwbse_engine_.ExcitationEnergies(orbitals);
    double energy_displaced_plus =
        orbitals.getTotalStateEnergy(tracker_.CalcState(orbitals));
    // update the coordinate
    pos_displaced = current_pos - displacement_vec;
    orbitals.updateAtomPostion(atom_index, pos_displaced);
    gwbse_engine_.ExcitationEnergies(orbitals);
    double energy_displaced_minus =
        orbitals.getTotalStateEnergy(tracker_.CalcState(orbitals));
    force(i_cart) =
        0.5 * (energy_displaced_minus - energy_displaced_plus) / displacement_;
    orbitals.updateAtomPostion(
        atom_index, current_pos);  // restore original coordinate into orbital
  }
  return force;
}

void Forces::RemoveTotalForce() {
  Eigen::Vector3d avgtotal_force =
      forces_.colwise().sum() / double(forces_.rows());
  for (Index i_atom = 0; i_atom < forces_.rows(); i_atom++) {
    forces_.row(i_atom) -= avgtotal_force;
  }
  return;
}

}  // namespace xtp
}  // namespace votca
