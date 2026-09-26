/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ewald_potential_test

// Standard includes
#include <cmath>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/ewald_potential.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "xtp_libint2.h"

using namespace votca;
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(ewald_potential_test)

namespace {

AOBasis CreateBasis(const QMMolecule& mol) {
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/vxc_grid/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

void FillPotential(Vxc_Grid& grid, double value) {
  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    GridBox& box = grid[i];
    box.getPotentialValues().assign(box.getGridPoints().size(), value);
  }
}

}  // namespace

// A CONSTANT potential is the one case where the answer is known in
// closed form without doing the integral: with phi = c everywhere,
//
//   V_uv = -integral chi_u(r) chi_v(r) c dr = -c * S_uv
//
// so the integrator must reproduce MINUS c times the AO overlap matrix.
//
// The minus is the electron charge, and it is the whole point of this
// case. What comes back goes straight into H0, where Tr(D H0) must
// already be an ENERGY, and the electrons that D counts carry charge -1.
// AOMultipole::FillPotential applies the same conversion for the nuclei
// and the external multipoles (`aopotential_ -= Fill(aobasis)`).
//
// This case previously asserted +c * S, i.e. it encoded the missing
// sign. In a job that error is nearly invisible: for a neutral QM
// region the nuclear term +sum_A Z_A phi(R_A) and the electronic term
// -\int rho phi differ only through the shape of the density, not its
// total charge, so they very nearly cancel -- and flipping one turns
// the cancellation into a doubling. Measured on a neutral methane in a
// rank-0 methane background against an identical job with no
// ewaldregion: -13.506 meV where the answer is +0.128 meV, the former
// being almost exactly twice the -6.689 meV nuclear half.
//
// The decisive assertion is the pair, not either half. phi = -1 must
// give exactly minus what phi = +1 gives -- same arithmetic, one sign --
// which needs no guess about how well a "medium" grid integrates an
// overlap. It is also what pins down the bug this case was written for:
// the small-integrand guard was copied from Vxc_Potential, where it
// reads rho*weight < 1e-20 and rho is a DENSITY, hence non-negative.
// Applied to weight*phi it discards every point where the potential is
// negative -- roughly half of space around a neutral background -- so
// the phi = -1 matrix came back all zeros while phi = +1 looked fine.
BOOST_AUTO_TEST_CASE(constant_potential_gives_minus_the_overlap_matrix) {
  libint2::initialize();
  QMMolecule mol("none", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_grid/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  AOOverlap overlap;
  overlap.Fill(aobasis);
  const Eigen::MatrixXd& S = overlap.Matrix();

  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);

  FillPotential(grid, 1.0);
  const Eigen::MatrixXd plus = Ewald_Potential<Vxc_Grid>(grid)
                                   .IntegrateEwald(aobasis.AOBasisSize())
                                   .matrix();

  FillPotential(grid, -1.0);
  const Eigen::MatrixXd minus = Ewald_Potential<Vxc_Grid>(grid)
                                    .IntegrateEwald(aobasis.AOBasisSize())
                                    .matrix();

  BOOST_REQUIRE_EQUAL(plus.rows(), S.rows());
  BOOST_REQUIRE_GT(plus.cwiseAbs().maxCoeff(), 1e-6);

  // Exact: nothing about the grid enters a sign flip.
  BOOST_CHECK_SMALL((minus + plus).cwiseAbs().maxCoeff(),
                    1e-12 * plus.cwiseAbs().maxCoeff());

  // The physical direction, stated bluntly so that no future rewrite of
  // the comparison below can quietly restore the wrong convention: S is
  // positive definite, so a POSITIVE potential everywhere must give a
  // matrix whose contraction with any density is NEGATIVE. Electrons are
  // attracted into a region of positive potential.
  BOOST_CHECK_LT(plus.trace(), 0.0);

  // Loose, and only an anchor that the right object is being integrated:
  // a "medium" grid does not reproduce an overlap matrix to machine
  // precision, and this case is not a statement about grid quality.
  BOOST_CHECK_SMALL((plus + S).cwiseAbs().maxCoeff(),
                    1e-4 * S.cwiseAbs().maxCoeff());

  libint2::finalize();
}

// A potential that changes sign across the molecule, which is what the
// periodic background actually produces. Integrating phi and -phi must
// again differ only in sign; under the discarded-negative-points bug the
// two are unrelated, because each run drops a different half of the grid.
BOOST_AUTO_TEST_CASE(sign_changing_potential_is_antisymmetric) {
  libint2::initialize();
  QMMolecule mol("none", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_grid/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);

  auto fill_with_z = [&](double scale) {
    for (Index i = 0; i < grid.getBoxesSize(); ++i) {
      GridBox& box = grid[i];
      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      std::vector<double>& values = box.getPotentialValues();
      values.resize(points.size());
      for (std::size_t p = 0; p < points.size(); ++p) {
        values[p] = scale * points[p].z();
      }
    }
  };

  fill_with_z(1.0);
  const Eigen::MatrixXd plus = Ewald_Potential<Vxc_Grid>(grid)
                                   .IntegrateEwald(aobasis.AOBasisSize())
                                   .matrix();

  fill_with_z(-1.0);
  const Eigen::MatrixXd minus = Ewald_Potential<Vxc_Grid>(grid)
                                    .IntegrateEwald(aobasis.AOBasisSize())
                                    .matrix();

  BOOST_REQUIRE_GT(plus.cwiseAbs().maxCoeff(), 1e-6);
  BOOST_CHECK_SMALL((minus + plus).cwiseAbs().maxCoeff(),
                    1e-12 * plus.cwiseAbs().maxCoeff());

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
