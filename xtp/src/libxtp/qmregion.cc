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

// Local VOTCA includes
#include "votca/xtp/qmregion.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/density_integration.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/ewaldregion.h"
#include "votca/xtp/gwbse.h"
#include "votca/xtp/pmlocalization.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/staticregion.h"
#include "votca/xtp/vxc_grid.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace votca {
namespace xtp {

void QMRegion::Initialize(const tools::Property& prop) {
  if (this->id_ != 0) {
    throw std::runtime_error(
        this->identify() +
        " must always be region 0. Currently only one qm region is possible.");
  }

  initstate_ = prop.get("state").as<QMState>();
  /*if (initstate_.Type() == QMStateType::Hole ||
      initstate_.Type() == QMStateType::Electron) {
    throw std::runtime_error(
        "Charged QM Regions are not implemented currently");
  }*/
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

  if (initstate_.Type().isKSState()) {

    if (prop.exists("statetracker")) {
      tools::Property filter = prop.get("statetracker");
      statetracker_.setLogger(&log_);
      statetracker_.Initialize(filter);
      statetracker_.setInitialState(initstate_);
      statetracker_.PrintInfo();
    } else {
      throw std::runtime_error("No statetracker options for KS states found");
    }
  }

  grid_accuracy_for_ext_interaction_ =
      prop.get("grid_for_potential").as<std::string>();
  DeltaE_ = prop.get("tolerance_energy").as<double>();
  DeltaD_ = prop.get("tolerance_density").as<double>();
  DeltaDmax_ = prop.get("tolerance_density_max").as<double>();

  dftoptions_ = prop.get("dftpackage");
  localize_options_ = prop.get("localize");

  if (prop.get("dftpackage.xtpdft.dft_in_dft.activeatoms")
          .as<std::string>()
          .size() != 0) {
    do_localize_ = true;
    do_dft_in_dft_ = true;
  }

  if (prop.exists("environment_screening")) {
    const tools::Property env = prop.get("environment_screening");
    screening_ = env.ifExistsReturnElseReturnDefault<bool>("enabled", true);
  }
  if (screening_) {
    const tools::Property env = prop.get("environment_screening");
    // Excitations are what the screening is for. A charged DFT run would
    // put the extra charge's polarization into the ground-state loop and
    // again into W; a plain ground-state or KS-level run has no GW-BSE to
    // screen.
    if (initstate_.Type() == QMStateType::Electron ||
        initstate_.Type() == QMStateType::Hole) {
      throw std::runtime_error(
          "environment_screening with state " + initstate_.ToString() +
          ": the charged DFT state would be polarized in the QM/MM loop and "
          "screened again in GW. Run the neutral molecule and ask for the "
          "quasiparticle level instead: state pqpN (or dqpN) with N the "
          "HOMO or LUMO index.");
    }
    if (initstate_.Type() == QMStateType::ExcitonUKS) {
      throw std::runtime_error(
          "environment_screening is implemented for closed-shell GW-BSE "
          "only; the unrestricted exciton state " +
          initstate_.ToString() + " is not supported.");
    }
    if (!do_gwbse_) {
      throw std::runtime_error(
          "environment_screening needs a GW-BSE state (an exciton or a "
          "quasiparticle level); " +
          initstate_.ToString() + " has no GW-BSE calculation to screen.");
    }
    screening_include_kreac_ =
        env.ifExistsReturnElseReturnDefault<bool>("include_kreac", true);
    screening_shell_dielectric_ =
        env.ifExistsReturnElseReturnDefault<double>("shell_dielectric", 4.0);
    screening_site_width_ =
        env.ifExistsReturnElseReturnDefault<double>("site_width", 0.5);
    if (!(screening_site_width_ >= 0.0)) {
      throw std::runtime_error(
          "environment_screening: site_width must be >= 0 (0 for point "
          "sites), got " +
          std::to_string(screening_site_width_) + ".");
    }
    std::string shells =
        env.ifExistsReturnElseReturnDefault<std::string>("shell_regions", "");
    std::replace(shells.begin(), shells.end(), ',', ' ');
    std::istringstream in(shells);
    Index id;
    while (in >> id) {
      screening_shell_regions_.push_back(id);
    }
    XTP_LOG(Log::error, log_)
        << " Environment screening on: ground-state QM/MM loop, then one "
           "screened GW-BSE run. K_reac "
        << (screening_include_kreac_ ? "included" : "excluded")
        << ", shell dielectric " << screening_shell_dielectric_
        << ", site width " << screening_site_width_ << " alpha^(1/3)"
        << std::flush;
  }
}

void QMRegion::PrepareOrbitalsForGWBSE() {
  if (orb_.isOpenShell()) {
    throw std::runtime_error("GWBSE not implemented for open-shell systems");
  }
  if (do_dft_in_dft_) {
    Index active_electrons = orb_.getNumOfActiveElectrons();

    if (active_electrons % 2 != 0) {
      throw std::runtime_error(
          "DFT-in-DFT embedded GWBSE currently supports only closed-shell "
          "active spaces with an even number of electrons");
    }

    orb_.MOs() = orb_.getEmbeddedMOs();

    const Index active_occ = active_electrons / 2;

    // Restricted closed-shell rewrite of the embedded orbital object
    orb_.setNumberOfOccupiedLevels(active_occ);
    orb_.setNumberOfOccupiedLevelsBeta(active_occ);
    orb_.setNumberOfAlphaElectrons(active_occ);
    orb_.setNumberOfBetaElectrons(active_occ);
    orb_.setChargeAndSpin(orb_.getCharge(), 1);
  }
}

double QMRegion::StateEnergy(const QMState& state) const {
  if (state.Type().isExciton()) {
    return orb_.getExcitedStateEnergy(state);
  }
  // quasiparticle: adding an electron to an unoccupied level, removing one
  // from an occupied level
  if (state.StateIdx() > orb_.getHomo()) {
    return orb_.getExcitedStateEnergy(state);
  }
  return -orb_.getExcitedStateEnergy(state);
}

ScreeningEnvironment QMRegion::BuildScreeningEnvironment(
    const std::vector<std::unique_ptr<Region> >& regions) const {
  ScreeningEnvironment env;
  env.shell_dielectric = screening_shell_dielectric_;
  env.include_kreac = screening_include_kreac_;
  env.site_width = screening_site_width_;
  bool have_damp = false;
  PolarRegion Polardummy(0, log_);
  for (const std::unique_ptr<Region>& reg : regions) {
    if (reg->identify() != Polardummy.identify()) {
      continue;
    }
    const PolarRegion& polar = dynamic_cast<const PolarRegion&>(*reg);
    const bool shell =
        std::find(screening_shell_regions_.begin(),
                  screening_shell_regions_.end(),
                  polar.getId()) != screening_shell_regions_.end();
    if (shell) {
      env.shell_segments.insert(env.shell_segments.end(), polar.begin(),
                                polar.end());
      continue;
    }
    if (have_damp && polar.ExpDamp() != env.exp_damp) {
      throw std::runtime_error(
          "environment_screening: the explicit polar regions respond as one "
          "coupled Thole system and must share exp_damp; found " +
          std::to_string(env.exp_damp) + " and " +
          std::to_string(polar.ExpDamp()) + ".");
    }
    env.exp_damp = polar.ExpDamp();
    have_damp = true;
    env.explicit_segments.insert(env.explicit_segments.end(), polar.begin(),
                                 polar.end());
  }
  for (Index id : screening_shell_regions_) {
    const bool is_polar = std::any_of(
        regions.begin(), regions.end(),
        [&](const std::unique_ptr<Region>& reg) {
          return reg->getId() == id && reg->identify() == Polardummy.identify();
        });
    if (!is_polar) {
      throw std::runtime_error("environment_screening: shell region " +
                               std::to_string(id) +
                               " is not a polar region of this job.");
    }
  }
  if (env.empty()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " WARNING: environment_screening is on but there is no polar "
           "region; GW-BSE runs unscreened."
        << std::flush;
  }
  return env;
}

std::string QMRegion::ScreeningAuxBasisName() const {
  // The auxiliary basis GWBSE::Initialize will pick: the DFT run's own
  // (dftpackage.auxbasisset), else gwbse.auxbasisset, else aux-<basis>.
  if (dftoptions_.exists("auxbasisset") &&
      !dftoptions_.get("auxbasisset").as<std::string>().empty()) {
    return dftoptions_.get("auxbasisset").as<std::string>();
  }
  if (gwbseoptions_.exists("auxbasisset") &&
      !gwbseoptions_.get("auxbasisset").as<std::string>().empty()) {
    return gwbseoptions_.get("auxbasisset").as<std::string>();
  }
  return "aux-" + dftoptions_.get("basisset").as<std::string>();
}

void QMRegion::CheckScreeningEnvironment(
    const std::vector<std::unique_ptr<Region> >& regions) const {
  if (!screening_) {
    return;
  }
  // R depends on the QM atoms, the auxiliary basis on them and the polar
  // sites -- none of which the QM/MM loop changes -- so it is checked
  // here, before the loop, rather than after the ground state has
  // converged. Failing now costs seconds; failing in EvaluateScreenedGWBSE
  // costs the whole loop.
  const ScreeningEnvironment env = BuildScreeningEnvironment(regions);
  if (env.empty()) {
    return;
  }
  const std::string auxname = ScreeningAuxBasisName();
  BasisSet bs;
  bs.Load(auxname);
  AOBasis auxbasis;
  auxbasis.Fill(bs, orb_.QMAtoms());
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Checking the screening environment (" << auxname
      << ", " << auxbasis.AOBasisSize() << " functions)" << std::flush;
  const ScreeningCheck check = EnvironmentScreening::Check(
      auxbasis, orb_.QMAtoms(), env, EnvironmentScreening::Metric(auxbasis));
  if (!check.ok()) {
    throw std::runtime_error(
        "environment_screening: 1 + R is not positive definite (lowest "
        "eigenvalue of R " +
        std::to_string(check.lowest) +
        "): the environment would screen some charge fluctuation of the QM "
        "region by more than the fluctuation itself. Diagnosis:\n" +
        check.report);
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Screening environment is stable: lowest eigenvalue "
      << "of R " << check.lowest << " (must be > -1)" << std::flush;
  XTP_LOG(Log::info, log_) << check.report << std::flush;
}

void QMRegion::EvaluateScreenedGWBSE(
    std::vector<std::unique_ptr<Region> >& regions) {
  if (!screening_) {
    return;
  }
  if (orb_.getCalculationType() == "Truncated") {
    throw std::runtime_error(
        "environment_screening is not implemented for truncated active "
        "regions.");
  }

  const ScreeningEnvironment env = BuildScreeningEnvironment(regions);

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Running environment-screened GW-BSE" << std::flush;
  const double e_dft = orb_.getDFTTotalEnergy();
  PrepareOrbitalsForGWBSE();
  GWBSE gwbse(orb_);
  gwbse.setLogger(&log_);
  gwbse.Initialize(gwbseoptions_);
  gwbse.setScreeningEnvironment(env);
  gwbse.Evaluate();
  const QMState state = statetracker_.CalcStateAndUpdate(orb_);
  const double energy = e_dft + StateEnergy(state);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Screened state " << state.ToString() << ": "
      << std::setprecision(10) << StateEnergy(state) * tools::conv::hrt2ev
      << " eV, region energy " << energy << " Hartree" << std::flush;
  E_hist_.push_back(energy);
}

// helper function to hand the grid changable over to Ewald
std::vector<Eigen::Vector3d> QMRegion::copyEwaldGrid() {

  std::vector<const Eigen::Vector3d*> src = ewaldgrid_.getGridpoints();
  std::vector<Eigen::Vector3d> dst(src.size());
  std::transform(src.begin(), src.end(), dst.begin(),
                 [](const Eigen::Vector3d* v) { return *v; });
  return dst;
}

void QMRegion::PrepareEwaldPotentialGrid() {

  // purpose: prepare a grid integration object that can
  // - hand over the grid to Ewald
  // - after Ewald calculates the potential at the grid points, be put into
  // xtpdft
  //
  // The grid name is xtpdft's own integration_grid, NOT grid_for_potential.
  // That is a constraint, not a preference: what comes back from Ewald is a
  // potential VALUE at each quadrature point, and DFTEngine forms <rho|phi>
  // on the grid it rebuilds from grid_name_, which is this same key. A
  // different quality here gives a different number of boxes and a
  // different number of points per box, and the size checks on the
  // DFTEngine side would throw. grid_for_potential governs the opposite
  // direction -- the QM density's influence outward on classical regions,
  // in ApplyQMFieldToPolarSegments -- and has no business here.
  if (orb_.QMAtoms().size() == 0) {
    throw std::runtime_error(
        "QMRegion::PrepareEwaldPotentialGrid: the region holds no atoms. The "
        "grid is built over the QM molecule, so this must run after the "
        "region has been filled.");
  }
  std::string dftbasis_name = dftoptions_.get("basisset").as<std::string>();
  std::string grid_name =
      dftoptions_.get("xtpdft.integration_grid").as<std::string>();
  QMMolecule mol = orb_.QMAtoms();
  BasisSet bs;
  bs.Load(dftbasis_name);
  AOBasis dftbasis;
  dftbasis.Fill(bs, mol);
  ewaldgrid_.GridSetup(grid_name, mol, dftbasis);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Constructed Grid for Integration of Ewald Potential"
      << std::flush;
  // Vxc_Potential<Vxc_Grid> vxc(grid);
  ewald_grid_ready_ = true;
  return;
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
  if (Dchange < DeltaD_ && Dmax < DeltaDmax_ && std::abs(Echange) < DeltaE_) {
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

  Index crg = 0;
  if (initstate_.Type() == QMStateType::Electron) {
    crg = -1;
  } else if (initstate_.Type() == QMStateType::Hole) {
    crg = +1;
  }
  qmpackage_->setCharge(crg);
  qmpackage_->setRunDir(workdir_);
  qmpackage_->WriteInputFile(orb_);

  // Attach the periodic background to xtpdft as a potential sampled on
  // the integration grid. This used to be one of two routes -- the other
  // handed over multipole moments and k-vectors for the engine to build
  // its own AO matrices from -- but nothing ever set that one up, so it
  // has been removed. See the git history of this file if it is ever
  // wanted back.
  if (ewald_grid_ready_) {
    if (qmpackage_->getPackageName() != "xtp") {
      throw std::runtime_error("QMEwald can only run with XTP as qmpackage.");
    }
    qmpackage_->setEwaldgrid(ewaldgrid_);
    qmpackage_->setEwaldNuclearEnergy(ewald_nuclear_energy_);
  }

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
  if (do_localize_) {
    if (orb_.isOpenShell()) {
      throw std::runtime_error(
          "MO localization not implemented for open-shell systems");
    }
    PMLocalization pml(log_, localize_options_);
    pml.computePML(orb_);
  }
  QMMolecule originalmol = orb_.QMAtoms();
  if (do_dft_in_dft_) {
    if (orb_.isOpenShell()) {
      throw std::runtime_error(
          "DFT-in-DFT not implemented for open-shell systems");
    }
    // this only works with XTPDFT, so locally override global qmpackage_
    std::unique_ptr<QMPackage> xtpdft =
        std::unique_ptr<QMPackage>(QMPackageFactory().Create("xtp"));
    xtpdft->setLog(&log_);
    xtpdft->Initialize(dftoptions_);
    xtpdft->setRunDir(workdir_);
    Index charge = 0;
    if (initstate_.Type() == QMStateType::Electron) {
      charge = -1;
    } else if (initstate_.Type() == QMStateType::Hole) {
      charge = +1;
    }
    xtpdft->setCharge(charge);

    std::vector<double> energies = std::vector<double>(regions.size(), 0.0);
    for (std::unique_ptr<Region>& reg : regions) {
      Index id = reg->getId();
      if (id == this->getId()) {
        continue;
      }

      StaticRegion Staticdummy(0, log_);
      PolarRegion Polardummy(0, log_);
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " Evaluating interaction between "
          << this->identify() << " " << this->getId() << " and "
          << reg->identify() << " " << reg->getId() << std::flush;
      if (reg->identify() == Staticdummy.identify()) {
        StaticRegion* staticregion = dynamic_cast<StaticRegion*>(reg.get());
        xtpdft->AddRegion(*staticregion);
      } else if (reg->identify() == Polardummy.identify()) {
        PolarRegion* polarregion = dynamic_cast<PolarRegion*>(reg.get());
        xtpdft->AddRegion(*polarregion);
      } else {
        throw std::runtime_error(
            "Interaction of regions with types:" + this->identify() + " and " +
            reg->identify() + " not implemented");
      }
    }

    e_ext = std::accumulate(interact_energies.begin(), interact_energies.end(),
                            0.0);
    XTP_LOG(Log::info, log_)
        << TimeStamp()
        << " Calculated interaction potentials with other regions. E[hrt]= "
        << e_ext << std::flush;
    XTP_LOG(Log::info, log_) << "Writing inputs" << std::flush;
    xtpdft->WriteInputFile(orb_);
    bool active_run = xtpdft->RunActiveRegion();
    if (!active_run) {
      throw std::runtime_error("\n DFT in DFT embedding failed. Stopping!");
    }
    Logfile_parse = xtpdft->ParseLogFile(orb_);
    if (!Logfile_parse) {
      throw std::runtime_error("\n Parsing DFT logfile failed. Stopping!");
    }
  }

  QMState state = QMState("groundstate");
  double energy = orb_.getDFTTotalEnergy();

  if (initstate_.Type().isKSState()) {
    state = statetracker_.CalcStateAndUpdate(orb_);
    // if unoccupied, add QP level energy
    if (state.StateIdx() > orb_.getHomo()) {
      energy += orb_.getExcitedStateEnergy(state);
    } else {
      // if unoccupied, subtract QP level energy
      energy -= orb_.getExcitedStateEnergy(state);
    }
  }

  // With environment screening the loop converges the ground state only;
  // GW-BSE runs once afterwards (EvaluateScreenedGWBSE).
  if (do_gwbse_ && !screening_) {
    PrepareOrbitalsForGWBSE();
    GWBSE gwbse(orb_);
    gwbse.setLogger(&log_);
    gwbse.Initialize(gwbseoptions_);
    gwbse.Evaluate();
    state = statetracker_.CalcStateAndUpdate(orb_);
    energy += StateEnergy(state);
  }
  // If truncation was enabled then rewrite full basis/aux-basis, MOs in full
  // basis and full QMAtoms
  if (orb_.getCalculationType() == "Truncated") {
    orb_.MOs().eigenvectors() = orb_.getTruncMOsFullBasis();
    orb_.QMAtoms().clearAtoms();
    orb_.QMAtoms() = originalmol;
    orb_.SetupDftBasis(orb_.getDftBasis().Name());
    orb_.SetupAuxBasis(orb_.getAuxBasis().Name());
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

    Index electrons = orb_.getNumberOfAlphaElectrons();
    if (orb_.isOpenShell()) {
      electrons += orb_.getNumberOfBetaElectrons();
    } else {
      electrons *= 2;
    }

    charge = double(nuccharge - electrons);
  } else {
    QMState state = statetracker_.InitialState();
    if (state.Type().isExciton()) {
      charge = 0.0;
    } else if (state.Type().isSingleParticleState() ||
               state.Type().isKSState()) {
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
  qmpackage_ =
      std::unique_ptr<QMPackage>(QMPackageFactory().Create(dft_package_name));
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
      orb_.getDftBasis();  // grid needs a basis in scope all the time
  grid.GridSetup(grid_accuracy_for_ext_interaction_, orb_.QMAtoms(), basis);
  DensityIntegration<Vxc_Grid> numint(grid);

  // With environment screening the environment is polarized by the
  // ground state only: the excitation's own polarization is in W.
  QMState state = QMState("groundstate");
  if ((do_gwbse_ && !screening_) || initstate_.Type().isKSState()) {
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
  CheckpointWriter v = w.openChild("QMdata");
  orb_.WriteToCpt(v);
  orb_.WriteBasisSetsToCpt(v);

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
  CheckpointReader rr = r.openChild("QMdata");
  orb_.ReadFromCpt(rr);
  orb_.ReadBasisSetsFromCpt(rr);

  CheckpointReader rr2 = r.openChild("E-hist");
  E_hist_.ReadFromCpt(rr2);

  CheckpointReader rr3 = r.openChild("D-hist");
  Dmat_hist_.ReadFromCpt(rr3);
  if (do_gwbse_) {
    CheckpointReader rr4 = r.openChild("statefilter");
    statetracker_.ReadFromCpt(rr4);
  }
}

double QMRegion::InteractwithEwaldRegion(const EwaldRegion& region) {
  // First contact with an EwaldRegion is what decides this region takes its
  // environment as a potential on a grid, so the grid is built here rather
  // than in Initialize. See the declaration for why eagerly would be wrong.
  if (!ewald_grid_ready_) {
    PrepareEwaldPotentialGrid();
  }

  // ONCE PER JOB, not once per iteration. phi is the potential of a FROZEN
  // background evaluated at FIXED points: the registry is converged before
  // the job starts and nothing re-polarizes it, the foreground registration
  // is made once by JobTopology, and the grid follows the QM geometry,
  // which does not move inside the inter-region SCF loop. The QM density
  // enters nowhere in it. So re-evaluating would give bit-identical numbers
  // at the cost of a full lattice sum over every quadrature point, every
  // iteration. The values stay where they were written -- in ewaldgrid_'s
  // own boxes and in ewald_nuclear_energy_, both members of this region
  // that Reset() does not touch.
  //
  // This stops being true the moment the background is allowed to respond
  // to the foreground. If that ever changes, this is the guard to remove.
  if (ewald_potential_evaluated_) {
    return 0.0;
  }

  // ELECTRONS. The potential on the integration grid, which the DFT
  // engine integrates against the density into H0.
  const std::vector<Eigen::Vector3d> points = copyEwaldGrid();
  const Eigen::VectorXd phi = region.PotentialAt(points);

  // Vxc_Grid::getGridpoints() walks grid_boxes_ in order and each box's
  // own points in order, so the flat vector maps back by giving each box
  // its share in turn. Checked rather than assumed: a mispairing here is
  // a wrong Hamiltonian with nothing to notice it by.
  Index offset = 0;
  for (Index i = 0; i < ewaldgrid_.getBoxesSize(); ++i) {
    GridBox& box = ewaldgrid_[i];
    const Index n_box = Index(box.getGridPoints().size());
    std::vector<double>& values = box.getPotentialValues();
    values.resize(std::size_t(n_box));
    for (Index p = 0; p < n_box; ++p) {
      values[std::size_t(p)] = phi[offset + p];
    }
    offset += n_box;
  }
  if (offset != phi.size()) {
    throw std::runtime_error(
        "QMRegion::InteractwithEwaldRegion: the grid holds " +
        std::to_string(offset) + " points across its boxes but " +
        std::to_string(phi.size()) +
        " potential values were evaluated. getGridpoints() and the box "
        "walk above disagree about the point ordering.");
  }

  // NUCLEI. The grid potential reaches the electron density only; the
  // nuclei sit in the same potential and have no other route in. Handed
  // to the engine as a scalar to add to E0, NOT returned from here --
  // QMRegion::Evaluate accumulates the Interactwith* return values into
  // e_ext, logs it, and drops it on the floor.
  std::vector<Eigen::Vector3d> nuclei;
  nuclei.reserve(std::size_t(orb_.QMAtoms().size()));
  for (const QMAtom& atom : orb_.QMAtoms()) {
    nuclei.push_back(atom.getPos());
  }
  const Eigen::VectorXd phi_nuclei = region.PotentialAt(nuclei);
  ewald_nuclear_energy_ = 0.0;
  Index a = 0;
  for (const QMAtom& atom : orb_.QMAtoms()) {
    ewald_nuclear_energy_ += double(atom.getNuccharge()) * phi_nuclei[a];
    ++a;
  }

  ewald_potential_evaluated_ = true;

  // 0.0, matching this class's own InteractwithPolarRegion and
  // InteractwithStaticRegion: the environment enters the Hamiltonian, and
  // its energy comes out in the DFT total rather than here.
  return 0.0;
}

}  // namespace xtp
}  // namespace votca
