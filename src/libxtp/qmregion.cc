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

// Local VOTCA includes
#include "votca/xtp/qmregion.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/gwbse.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/staticregion.h"
#include "votca/xtp/vxc_grid.h"

namespace votca {
namespace xtp {

void QMRegion::Initialize(const tools::Property& prop) {
  if (this->id_ != 0) {
    throw std::runtime_error(
        this->identify() +
        " must always be region 0. Currently only one qm region is possible.");
  }

  initstate_ = prop.get("state").as<QMState>();
  if (initstate_.Type() == QMStateType::Hole ||
      initstate_.Type() == QMStateType::Electron) {
    throw std::runtime_error(
        "Charged QM Regions are not implemented currently");
  }
  if (initstate_.Type().isExciton() || initstate_.Type().isGWState()) {

    do_gwbse_ = true;
    gwbseoptions_ = prop.get("gwbse");
    if (prop.exists("statetracker")) {
      tools::Property filter = prop.get("statetracker");
      statetracker_.setLogger(&log_);
      statetracker_.Initialize(filter);
      statetracker_.setInitialState(initstate_);
      statetracker_.PrintInfo();
    } else {
      throw std::runtime_error(
          "No statetracker options for excited states found");
    }
  }

  grid_accuracy_for_ext_interaction_ =
      prop.get("grid_for_potential").as<std::string>();
  DeltaE_ = prop.get("tolerance_energy").as<double>();
  DeltaD_ = prop.get("tolerance_density").as<double>();

  dftoptions_ = prop.get("dftpackage");
}

bool QMRegion::Converged() const {
  if (!E_hist_.filled()) {
    return false;
  }

  double Echange = E_hist_.getDiff();
  double Dchange =
      Dmat_hist_.getDiff().norm() / double(Dmat_hist_.back().cols());
  double Dmax = Dmat_hist_.getDiff().cwiseAbs().maxCoeff();
  std::string info = "not converged";
  bool converged = false;
  if (Dchange < DeltaD_ && Dmax < DeltaD_ && std::abs(Echange) < DeltaE_) {
    info = "converged";
    converged = true;
  }
  XTP_LOG(Log::error, log_)
      << " Region:" << this->identify() << " " << this->getId() << " is "
      << info << " deltaE=" << Echange << " RMS Dmat=" << Dchange
      << " MaxDmat=" << Dmax << std::flush;
  return converged;
}

void QMRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {

  std::vector<double> interact_energies = ApplyInfluenceOfOtherRegions(regions);
  double e_ext =
      std::accumulate(interact_energies.begin(), interact_energies.end(), 0.0);
  XTP_LOG(Log::info, log_)
      << TimeStamp()
      << " Calculated interaction potentials with other regions. E[hrt]= "
      << e_ext << std::flush;
  XTP_LOG(Log::info, log_) << "Writing inputs" << std::flush;
  qmpackage_->setRunDir(workdir_);
  qmpackage_->WriteInputFile(orb_);
  XTP_LOG(Log::error, log_) << "Running DFT calculation" << std::flush;
  bool run_success = qmpackage_->Run();
  if (!run_success) {
    info_ = false;
    errormsg_ = "DFT-run failed";
    return;
  }

  bool Logfile_parse = qmpackage_->ParseLogFile(orb_);
  if (!Logfile_parse) {
    info_ = false;
    errormsg_ = "Parsing DFT logfile failed.";
    return;
  }
  bool Orbfile_parse = qmpackage_->ParseMOsFile(orb_);
  if (!Orbfile_parse) {
    info_ = false;
    errormsg_ = "Parsing DFT orbfile failed.";
    return;
  }
  QMState state = QMState("groundstate");
  double energy = orb_.getDFTTotalEnergy();
  if (do_gwbse_) {
    GWBSE gwbse(orb_);
    gwbse.setLogger(&log_);
    gwbse.Initialize(gwbseoptions_);
    gwbse.Evaluate();
    state = statetracker_.CalcStateAndUpdate(orb_);
    if (state.Type().isExciton()) {
      energy += orb_.getExcitedStateEnergy(state);
    } else {
      // if unoccupied, add QP level energy
      if (state.StateIdx() > orb_.getHomo()) {
        energy += orb_.getExcitedStateEnergy(state);
      } else {
        // if unoccupied, subtract QP level energy
        energy -= orb_.getExcitedStateEnergy(state);
      }
    }
  }
  E_hist_.push_back(energy);
  Dmat_hist_.push_back(orb_.DensityMatrixFull(state));
  return;
}

void QMRegion::push_back(const QMMolecule& mol) {
  if (orb_.QMAtoms().size() == 0) {
    orb_.QMAtoms() = mol;
  } else {
    orb_.QMAtoms().AddContainer(mol);
  }
  size_++;
}

double QMRegion::charge() const {
  double charge = 0.0;
  if (!do_gwbse_) {
    Index nuccharge = 0;
    for (const QMAtom& a : orb_.QMAtoms()) {
      nuccharge += a.getNuccharge();
    }

    Index electrons = orb_.getNumberOfAlphaElectrons() * 2;
    charge = double(nuccharge - electrons);
  } else {
    QMState state = statetracker_.InitialState();
    if (state.Type().isExciton()) {
      charge = 0.0;
    } else if (state.Type().isSingleParticleState()) {
      if (state.StateIdx() <= orb_.getHomo()) {
        charge = +1.0;
      } else {
        charge = -1.0;
      }
    }
  }
  return charge;
}

void QMRegion::AppendResult(tools::Property& prop) const {
  prop.add("E_total", std::to_string(E_hist_.back() * tools::conv::hrt2ev));
  if (do_gwbse_) {
    prop.add("Initial_State", statetracker_.InitialState().ToString());
    prop.add("Final_State", statetracker_.CalcState(orb_).ToString());
  }
}

void QMRegion::Reset() {

  std::string dft_package_name = dftoptions_.get("name").as<std::string>();
  qmpackage_ = std::unique_ptr<QMPackage>(
      QMPackageFactory::QMPackages().Create(dft_package_name));
  qmpackage_->setLog(&log_);
  qmpackage_->Initialize(dftoptions_);
  Index charge = 0;
  if (initstate_.Type() == QMStateType::Electron) {
    charge = -1;
  } else if (initstate_.Type() == QMStateType::Hole) {
    charge = +1;
  }
  qmpackage_->setCharge(charge);
  return;
}
double QMRegion::InteractwithQMRegion(const QMRegion&) {
  throw std::runtime_error(
      "QMRegion-QMRegion interaction is not implemented yet.");
  return 0.0;
}
double QMRegion::InteractwithPolarRegion(const PolarRegion& region) {
  qmpackage_->AddRegion(region);
  return 0.0;
}
double QMRegion::InteractwithStaticRegion(const StaticRegion& region) {
  qmpackage_->AddRegion(region);
  return 0.0;
}

void QMRegion::WritePDB(csg::PDBWriter& writer) const {
  writer.WriteContainer(orb_.QMAtoms());
}

void QMRegion::AddNucleiFields(std::vector<PolarSegment>& segments,
                               const StaticSegment& seg) const {
  eeInteractor e;
#pragma omp parallel for
  for (Index i = 0; i < Index(segments.size()); ++i) {
    e.ApplyStaticField<StaticSegment, Estatic::noE_V>(seg, segments[i]);
  }
}

void QMRegion::ApplyQMFieldToPolarSegments(
    std::vector<PolarSegment>& segments) const {

  Vxc_Grid grid;
  AOBasis basis =
      orb_.SetupDftBasis();  // grid needs a basis in scope all the time
  grid.GridSetup(grid_accuracy_for_ext_interaction_, orb_.QMAtoms(), basis);
  DensityIntegration<Vxc_Grid> numint(grid);

  QMState state = QMState("groundstate");
  if (do_gwbse_) {
    state = statetracker_.CalcState(orb_);
  }
  Eigen::MatrixXd dmat = orb_.DensityMatrixFull(state);
  double Ngrid = numint.IntegrateDensity(dmat);
  AOOverlap overlap;
  overlap.Fill(basis);
  double N_comp = dmat.cwiseProduct(overlap.Matrix()).sum();
  if (std::abs(Ngrid - N_comp) > 0.001) {
    XTP_LOG(Log::error, log_) << "=======================" << std::flush;
    XTP_LOG(Log::error, log_)
        << "WARNING: Calculated Densities at Numerical Grid, Number of "
           "electrons "
        << Ngrid << " is far away from the the real value " << N_comp
        << ", you should increase the accuracy of the integration grid."
        << std::flush;
  }
#pragma omp parallel for
  for (Index i = 0; i < Index(segments.size()); ++i) {
    PolarSegment& seg = segments[i];
    for (PolarSite& site : seg) {
      site.V_noE() += numint.IntegrateField(site.getPos());
    }
  }

  StaticSegment seg(orb_.QMAtoms().getType(), orb_.QMAtoms().getId());
  for (const QMAtom& atom : orb_.QMAtoms()) {
    seg.push_back(StaticSite(atom, double(atom.getNuccharge())));
  }
  AddNucleiFields(segments, seg);
}

void QMRegion::WriteToCpt(CheckpointWriter& w) const {
  w(id_, "id");
  w(identify(), "type");
  w(size_, "size");
  w(do_gwbse_, "GWBSE");
  w(initstate_.ToString(), "initial_state");
  w(grid_accuracy_for_ext_interaction_, "ext_grid");
  CheckpointWriter v = w.openChild("orbitals");
  orb_.WriteToCpt(v);

  CheckpointWriter v2 = w.openChild("E-hist");
  E_hist_.WriteToCpt(v2);

  CheckpointWriter v3 = w.openChild("D-hist");
  Dmat_hist_.WriteToCpt(v3);
  if (do_gwbse_) {
    CheckpointWriter v4 = w.openChild("statefilter");
    statetracker_.WriteToCpt(v4);
  }
}

void QMRegion::ReadFromCpt(CheckpointReader& r) {
  r(id_, "id");
  r(size_, "size");
  r(do_gwbse_, "GWBSE");
  std::string state;
  r(state, "initial_state");
  initstate_.FromString(state);
  r(grid_accuracy_for_ext_interaction_, "ext_grid");
  CheckpointReader rr = r.openChild("orbitals");
  orb_.ReadFromCpt(rr);

  CheckpointReader rr2 = r.openChild("E-hist");
  E_hist_.ReadFromCpt(rr2);

  CheckpointReader rr3 = r.openChild("D-hist");
  Dmat_hist_.ReadFromCpt(rr3);
  if (do_gwbse_) {
    CheckpointReader rr4 = r.openChild("statefilter");
    statetracker_.ReadFromCpt(rr4);
  }
}

}  // namespace xtp
}  // namespace votca
