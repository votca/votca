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
#include <votca/xtp/ewaldcontainer.h>

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

  std::vector<Eigen::Vector3d> copyEwaldGrid();

  Vxc_Grid& getEwaldGrid() { return ewaldgrid_; };

  // =========== EWALD MOMENTS SETTER AND ACCESS ==========
  // +++++++++++ BACKGROUND +++++++++++++++++++++++++++++++
  void setEwaldBackground(ewaldcontainer::PotentialData* bg) {
    ewald_background_ = bg;
    ewald_moments_ready_ = true;
  }

  ewaldcontainer::PotentialData& ewaldBackground() {
    assert(ewald_background_ != nullptr);
    return *ewald_background_;
  }

  const ewaldcontainer::PotentialData& ewaldBackground() const {
    assert(ewald_background_ != nullptr);
    return *ewald_background_;
  }

  // +++++++++++ FOREGROUND CORRECTION ++++++++++++++++++++++
  void setEwaldForegroundCorrection(ewaldcontainer::PotentialData* fg_corr) {
    ewald_foreground_correction_ = fg_corr;
  }

  ewaldcontainer::PotentialData& ewaldForegroundCorrection() {
    assert(ewald_foreground_correction_ != nullptr);
    return *ewald_foreground_correction_;
  }

  const ewaldcontainer::PotentialData& ewaldForegroundCorrection() const {
    assert(ewald_foreground_correction_ != nullptr);
    return *ewald_foreground_correction_;
  }

  // +++++++++++ SHAPE CORRECTION +++++++++++++++++++++++++++
  void setEwaldShapeCorrection(ewaldcontainer::PotentialData* shape_corr) {
    ewald_shape_correction_ = shape_corr;
  }

  ewaldcontainer::PotentialData& ewaldShapeCorrection() {
    assert(ewald_shape_correction_ != nullptr);
    return *ewald_shape_correction_;
  }

  const ewaldcontainer::PotentialData& ewaldShapeCorrection() const {
    assert(ewald_shape_correction_ != nullptr);
    return *ewald_shape_correction_;
  }

  // +++++++++++ MM1 REGION +++++++++++++++++++++++++++
  void setEwaldMM1(ewaldcontainer::PotentialData* mm1) { ewald_mm1_ = mm1; }

  ewaldcontainer::PotentialData& ewaldMM1() {
    assert(ewald_mm1_ != nullptr);
    return *ewald_mm1_;
  }

  const ewaldcontainer::PotentialData& ewaldMM1() const {
    assert(ewald_mm1_ != nullptr);
    return *ewald_mm1_;
  }

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
  bool do_localize_ = false;
  bool do_dft_in_dft_ = false;

  tools::Property dftoptions_;
  tools::Property gwbseoptions_;
  tools::Property localize_options_;

  StateTracker statetracker_;

  // for QMEwald
  //
  // TWO INDEPENDENT ROUTES, TWO FLAGS. The periodic background can reach
  // the Hamiltonian either as a potential sampled on this grid, which the
  // DFT engine integrates against the density, or as multipoles and
  // k-vectors the engine builds its own AO matrices from. They were
  // sharing one is_qmewald_, so preparing the grid alone also sent
  // Evaluate down the moments path and into ewaldBackground()'s
  // assert(ewald_background_ != nullptr).
  Vxc_Grid ewaldgrid_;
  bool ewald_grid_ready_ = false;
  // Whether the background's potential has already been laid down on that
  // grid. Separate from ewald_grid_ready_ because the grid is built once
  // and the potential is evaluated once, but for different reasons: the
  // grid because geometry and basis are fixed, the potential because the
  // background is frozen. See InteractwithEwaldRegion.
  bool ewald_potential_evaluated_ = false;
  // ORPHANED, deliberately kept. Set only by setEwaldBackground, whose
  // only remaining caller is XTPDFT's own pass-through -- the code that
  // originally drove it was the legacy xtp/ewald/ jobcalculator, which
  // has been removed. So nothing in a current job can make this true,
  // and the analytic QM-coupling machinery it gates (DFTEngine's
  // IntegrateEwaldRealSpaceMultipoles / IntegrateEwaldReciprocalSpace,
  // the foreground and shape corrections, the AOEwald* matrices and
  // aoplanewave) compiles but is unreachable.
  //
  // Kept because that machinery is most of an analytic alternative to
  // the grid route -- see EwaldRegion::PotentialAt -- which is worth
  // having once libint2 can do the operator-centre derivatives a rank-1
  // source needs. Reviving it means feeding these moments from
  // EwaldRegion rather than from the calculator that used to.
  //
  // Said here explicitly so nobody has to work out from scratch why a
  // whole code path never fires.
  bool ewald_moments_ready_ = false;
  // sum_A Z_A phi(R_A). The grid carries the potential the ELECTRONS
  // feel; the nuclei sit in the same potential and have no other way in.
  double ewald_nuclear_energy_ = 0.0;
  ewaldcontainer::PotentialData* ewald_background_ = nullptr;
  ewaldcontainer::PotentialData* ewald_foreground_correction_ = nullptr;
  ewaldcontainer::PotentialData* ewald_shape_correction_ = nullptr;
  ewaldcontainer::PotentialData* ewald_mm1_ = nullptr;
  // ewaldcontainer::PotentialData* ewald_qm0_ = nullptr;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMREGION_H
