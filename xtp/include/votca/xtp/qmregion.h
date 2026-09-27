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

#pragma once
#ifndef VOTCA_XTP_QMREGION_H
#define VOTCA_XTP_QMREGION_H

// Local VOTCA includes
#include "hist.h"
#include "orbitals.h"
#include "qmpackagefactory.h"
#include "region.h"
#include "statetracker.h"
#include "vxc_grid.h"

/**
 * \brief defines a qm region and runs dft and gwbse calculations
 *
 *
 *
 */

namespace votca {
namespace xtp {

class PolarRegion;
class StaticRegion;
class EwaldRegion;
class QMRegion : public Region {

 public:
  QMRegion(Index id, Logger& log, std::string workdir)
      : Region(id, log), workdir_(workdir) {};
  ~QMRegion() override = default;

  void Initialize(const tools::Property& prop) override;

  bool Converged() const override;

  void Evaluate(std::vector<std::unique_ptr<Region> >& regions) override;

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

  void ApplyQMFieldToPolarSegments(std::vector<PolarSegment>& segments) const;

  // Builds the integration grid the background's potential is sampled on.
  // Called lazily by InteractwithEwaldRegion rather than from Initialize,
  // for two reasons. A job with no EwaldRegion must not pay for a grid it
  // never uses -- and, less obviously, must not have ewald_grid_ready_ set,
  // because Evaluate takes that flag as the signal to reject any qmpackage
  // other than xtp. Building it eagerly would break every existing
  // qmmm job that runs orca.
  //
  // Takes no options: the grid name and basis come from dftoptions_, which
  // Initialize has already stored. They MUST be the ones DFTEngine uses --
  // see the comment in the definition.
  void PrepareEwaldPotentialGrid();

  Index size() const override { return size_; }

  void WritePDB(csg::PDBWriter& writer) const override;

  std::string identify() const override { return "qmregion"; }

  void push_back(const QMMolecule& mol);

  void Reset() override;

  double charge() const override;
  double Etotal() const override { return E_hist_.back(); }

  // Coordinates only, by value. That is not a convenience -- see
  // ewaldgrid_ below for why the grid object itself must not leave this
  // class.
  std::vector<Eigen::Vector3d> copyEwaldGrid();

  /**
   * \brief Embedded GW-BSE, screened by the polar regions (see
   * EnvironmentScreening, GWBSE::setScreeningEnvironment).
   *
   * With environment_screening on, the inter-region loop converges the
   * GROUND state: the QM density polarizes the polar regions, and their
   * induced dipoles act back on it -- the embedding H0. The excitation's
   * own polarization of the environment is not iterated; it enters the
   * screened interaction of one GW-BSE run instead, which this does, once
   * the loop has converged. The polar regions respond through their Thole
   * operator; those listed as shell regions through alpha/epsilon without
   * coupling. Updates the region energy with the state energy found.
   */
  bool EnvironmentScreeningEnabled() const { return screening_; }
  void EvaluateScreenedGWBSE(std::vector<std::unique_ptr<Region> >& regions);

 protected:
  void AppendResult(tools::Property& prop) const override;
  double InteractwithQMRegion(const QMRegion& region) override;
  double InteractwithPolarRegion(const PolarRegion& region) override;
  double InteractwithStaticRegion(const StaticRegion& region) override;
  double InteractwithEwaldRegion(const EwaldRegion& region) override;

 private:
  void AddNucleiFields(std::vector<PolarSegment>& segments,
                       const StaticSegment& seg) const;

  Index size_ = 0;
  Orbitals orb_;

  QMState initstate_;
  std::string workdir_ = "";
  std::unique_ptr<QMPackage> qmpackage_ = nullptr;

  std::string grid_accuracy_for_ext_interaction_ = "medium";

  hist<double> E_hist_;
  hist<Eigen::MatrixXd> Dmat_hist_;

  // convergence options
  double DeltaD_ = 5e-5;
  double DeltaE_ = 5e-5;
  double DeltaDmax_ = 5e-5;

  bool do_gwbse_ = false;

  // environment screening (embedded GW-BSE)
  bool screening_ = false;
  bool screening_include_kreac_ = true;
  std::vector<Index> screening_shell_regions_;
  double screening_shell_dielectric_ = 4.0;
  // DFT-in-DFT rewrite of orb_ and closed-shell check, before any GW-BSE
  void PrepareOrbitalsForGWBSE();
  // state energy of the tracked state, added to the DFT total energy
  double StateEnergy(const QMState& state) const;
  bool do_localize_ = false;
  bool do_dft_in_dft_ = false;

  tools::Property dftoptions_;
  tools::Property gwbseoptions_;
  tools::Property localize_options_;

  StateTracker statetracker_;

  // for QMEwald. The periodic background reaches the Hamiltonian as a
  // potential sampled on this grid, which the DFT engine integrates
  // against the density.
  //
  // There used to be a second route, handing the engine multipole
  // moments and k-vectors to build its own AO matrices from, gated by a
  // second flag. Nothing ever set it up -- the legacy jobcalculator that
  // once did was removed -- so it has been deleted along with the
  // machinery behind it (DFTEngine's IntegrateEwald* helpers, the
  // AOEwald* matrices, the ewaldcontainer types). Two things are worth
  // recording about that, since the decision was to remove code that had
  // been deliberately kept:
  //
  //  - It was blocked anyway. A rank-1 (induced dipole) source needs
  //    operator-centre derivatives from libint2, which is why
  //    AOEwaldRealSpaceDipoles was never instantiated even from the dead
  //    path -- the real-space route split every dipole into a pair of
  //    point charges instead.
  //  - It would not have been the faster route in any case. It replaces
  //    a loop over grid points (linear in QM size) with one over shell
  //    pairs (quadratic), so it wins only for small QM regions -- the
  //    opposite of what it was being kept for.
  //
  // Recoverable from the git history if either of those ever changes.
  //
  // ONLY ITS POINTS AND VALUES ARE VALID. PrepareEwaldPotentialGrid
  // builds this grid against a BasisSet, an AOBasis and a QMMolecule
  // that are all locals of that function, and GridBox::
  // FindSignificantShells stores raw `const AOShell*` into the basis it
  // is handed (gridbox.cc, addShell(&store)). Those shells die with the
  // function, so from the moment PrepareEwaldPotentialGrid returns this
  // object holds dangling pointers.
  //
  // What remains safe is everything that does not follow them:
  // getGridpoints, getPotentialValues, getBoxesSize, GridBox::size.
  // CalcAOValues, Matrixsize, AddtoBigMatrix and anything else touching
  // significant_shells is undefined behaviour. That is why the grid is
  // never integrated on here or in QMPackage, and why DFTEngine rebuilds
  // its own from grid_name_ and copies only the values across -- see
  // dftengine.cc.
  //
  // Nothing dereferences them today. The public Vxc_Grid& accessor that
  // used to sit next to copyEwaldGrid() was removed because it handed
  // this object out with no way to know that, and had no callers.
  // Giving the basis a longer life would remove the hazard, but it buys
  // nothing on its own: the transport is by value either way, and
  // DFTEngine has to integrate against the AO ordering of its own
  // dftbasis_ regardless.
  Vxc_Grid ewaldgrid_;
  bool ewald_grid_ready_ = false;
  // Whether the background's potential has already been laid down on that
  // grid. Separate from ewald_grid_ready_ because the grid is built once
  // and the potential is evaluated once, but for different reasons: the
  // grid because geometry and basis are fixed, the potential because the
  // background is frozen. See InteractwithEwaldRegion.
  bool ewald_potential_evaluated_ = false;
  // sum_A Z_A phi(R_A). The grid carries the potential the ELECTRONS
  // feel; the nuclei sit in the same potential and have no other way in.
  double ewald_nuclear_energy_ = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMREGION_H
