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

// Standard includes
#include <cmath>
#include <stdexcept>
#include <string>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/ewald_potential.h"

namespace votca {
namespace xtp {
template <class Grid>
Ewald_Potential<Grid>::~Ewald_Potential() {}



template <class Grid>
Mat_p_Energy Ewald_Potential<Grid>::IntegrateEwald(Index basissize) const {

  Mat_p_Energy vewald = Mat_p_Energy(basissize,basissize);

#pragma omp parallel for schedule(guided) reduction(+ : vewald)
  for (Index i = 0; i < grid_.getBoxesSize(); ++i) {
    const GridBox& box = grid_[i];
    if (!box.Matrixsize()) {
      continue;
    }
    double Eewald_box = 0.0;

    Eigen::MatrixXd Vewald_here =
        Eigen::MatrixXd::Zero(box.Matrixsize(),box.Matrixsize());
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    const std::vector<double>& pot = box.getPotentialValues();

    // The potential values are copied in from another grid, so a
    // mismatch in box structure would otherwise read past the end of a
    // vector rather than fail.
    if (Index(pot.size()) != box.size()) {
      throw std::runtime_error(
          "Ewald_Potential::IntegrateEwald: box " + std::to_string(i) +
          " has " + std::to_string(box.size()) + " grid points but " +
          std::to_string(pot.size()) +
          " potential values. The grid the potential was evaluated on is "
          "not the grid being integrated over.");
    }

    // iterate over gridpoints
    for (Index p = 0; p < box.size(); p++) {
      AOShell::AOValues ao = box.CalcAOValues(points[p]);
      const double weight = weights[p] * pot[p];

      // ABSOLUTE VALUE, and not for tidiness. The equivalent guard in
      // Vxc_Potential reads rho*weight < 1e-20, where rho is a DENSITY
      // and the product cannot be negative, so it means "negligibly
      // small". Here the integrand carries a POTENTIAL, which is
      // negative over roughly half of space around a neutral periodic
      // background -- the unsigned test discarded every one of those
      // points, keeping only the half that does not cancel.
      if (std::abs(weight) < 1.e-20) {
        continue;  // skip the rest, if integrand is very small
      }

      Eewald_box += weight;
      Vewald_here.noalias() += weight * ao.values * ao.values.transpose();
    }
    box.AddtoBigMatrix(vewald.matrix(), Vewald_here);
    vewald.energy() += Eewald_box;
  }

  // The accumulated Eewald_box is sum_p w_p phi_p, i.e. the integral of
  // the potential over the grid -- not an energy of anything. The energy
  // of the electron density in this potential is the contraction of the
  // matrix below with the density matrix, which is the caller's to form.
  // Returned as zero rather than as a plausible-looking number nobody
  // should use; the nuclear share has no route through here at all and
  // is supplied separately.
  //
  // ELECTRON CHARGE -- the minus sign. What the loop accumulated is the
  // plain potential integral, +<chi_mu|phi|chi_nu>. Electrons carry
  // charge -1, so their energy in this potential is -\int rho phi, and
  // the caller adds this matrix straight into H0, where Tr(D H0) must
  // already BE that energy. AOMultipole::FillPotential does exactly this
  // for the nuclei and for the external multipoles -- `aopotential_ -=
  // Fill(aobasis)` in both overloads -- and this is the same conversion
  // for the same reason.
  //
  // Nearly invisible if wrong. For a NEUTRAL QM region the nuclear term
  // +sum_A Z_A phi(R_A) and the electronic term -\int rho phi differ only
  // through the shape of the density, not its total charge, so they very
  // nearly cancel. Dropping the sign turns that cancellation into a
  // doubling. Measured on a neutral methane in a rank-0 methane
  // background, against an otherwise identical job with no ewaldregion:
  // nuclear -6.689 meV, electronic -6.817 meV, total -13.506 meV, where
  // the correct total is +0.128 meV -- the right order for a
  // near-spherical neutral molecule in a slowly varying periodic
  // potential, and consistent with the classical channel's own
  // -0.46 meV per segment.
  //
  // No existing test could catch it. A zeroed background has phi = 0, so
  // both halves vanish whatever the sign. EwaldRegion::PotentialAt is
  // pinned against lattice sums and unit probes, none of which build an
  // AO matrix. And ApplyFieldTo, which carries the entire validation
  // against legacy, never reaches this file.
  return Mat_p_Energy(0.0, -vewald.matrix());
}

template class Ewald_Potential<Vxc_Grid>;

}  // namespace xtp
}  // namespace votca
