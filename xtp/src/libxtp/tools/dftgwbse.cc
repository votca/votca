/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// VOTCA includes
#include <algorithm>
#include <sstream>
#include <votca/tools/constants.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/dftengine.h"
#include "votca/xtp/gwbseengine.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/staticregion.h"

// Local private VOTCA includes
#include "dftgwbse.h"
#include "votca/xtp/excitonfollow.h"

namespace votca {
namespace xtp {

void DftGwBse::ParseOptions(const tools::Property& options) {

  // molecule coordinates
  xyzfile_ = job_name_ + ".xyz";

  // options for dft package
  package_options_ = options.get(".dftpackage");

  // GWBSEENGINE options
  gwbseengine_options_ = options;

  // lets get the archive file name from the xyz file name
  archive_file_ = job_name_ + ".orb";

  // XML OUTPUT
  xml_output_ = job_name_ + "_summary.xml";

  if (options.exists(".excited_state_forces")) {
    excited_state_forces_ =
        options.get(".excited_state_forces").as<std::string>();
  }
  {
    const std::string tasks = options.get(".tasks").as<std::string>();
    tools::Tokenizer tok(tasks, " ,");
    const std::vector<std::string> list = tok.ToVector();
    const bool has =
        std::find(list.begin(), list.end(), "gradients") != list.end();
    if (has && list.size() > 1) {
      throw std::runtime_error(
          "tasks: gradients cannot be combined with other tasks (with dft "
          "and gwbse the gradients are computed in the same run anyway)");
    }
    gradients_only_ = has;
  }
  if (options.exists(".qp_gradients")) {
    qp_gradients_ = options.get(".qp_gradients").as<std::string>();
  }
  if (options.exists(".excited_state_forces_follow")) {
    excited_state_forces_follow_ =
        options.get(".excited_state_forces_follow").as<std::string>();
  }
  if (options.exists(".excited_state_forces_qp_response")) {
    excited_state_forces_qp_response_ =
        options.get(".excited_state_forces_qp_response").as<bool>();
  }

  if (options.exists(".mpsfile")) {
    mpsfile_ = options.get(".mpsfile").as<std::string>();
  }

  // check if guess is requested
  if (options.exists(".guess")) {
    guess_file_ = options.get(".guess").as<std::string>();
  }

  // check if a separate MO-only guess is requested (see this class's
  // own header comment on moguess_file_ for why this is a distinct
  // option from .guess, not a variant of it)
  if (options.exists(".moguess")) {
    moguess_file_ = options.get(".moguess").as<std::string>();
  }

  // register all QM packages
  QMPackageFactory{};
}

bool DftGwBse::Run() {

  log_.setReportLevel(Log::current_level);

  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  // tasks = gradients: no DFT or GW-BSE, the requested gradients
  // (qp_gradients, excited_state_forces) are computed from the data of an
  // earlier run in the orb file and added to it
  if (gradients_only_) {
    if (excited_state_forces_.empty() && qp_gradients_.empty()) {
      throw std::runtime_error(
          "tasks gradients: set qp_gradients and/or excited_state_forces");
    }
    XTP_LOG(Log::error, log_)
        << "Gradients from the GW-BSE data in " << archive_file_ << std::flush;
    orbitals.ReadFromCpt(archive_file_);
    CheckStoredRun(orbitals);
    if (!excited_state_forces_.empty()) {
      ComputeExcitedStateForces(orbitals);
    }
    if (!qp_gradients_.empty()) {
      ComputeQPGradients(orbitals);
    }
    XTP_LOG(Log::error, log_)
        << "Saving data to " << archive_file_ << std::flush;
    orbitals.WriteToCpt(archive_file_);
    return true;
  }

  if (!guess_file_.empty()) {
    XTP_LOG(Log::error, log_)
        << "Reading guess from " << guess_file_ << std::flush;
    orbitals.ReadFromCpt(guess_file_);
  } else {
    XTP_LOG(Log::error, log_)
        << "Reading structure from " << xyzfile_ << std::flush;
    orbitals.QMAtoms().LoadFromFile(xyzfile_);
  }

  // Warm-start the MOs from a PREVIOUS geometry step's own, separate
  // .orb file, while keeping the geometry just loaded above (the
  // CURRENT step's own xyzfile_, or guess_file_'s geometry if that
  // path was used instead) -- deliberately independent of guess_file_
  // itself, which ties geometry and MOs together from the same file
  // and therefore cannot represent this "new geometry, warm-started
  // MOs" case a geometry optimization needs between steps. Loaded
  // into a separate, temporary Orbitals object so only its own
  // MOs()/MOs_beta() get copied in -- moguess_orbitals's own geometry
  // (the PREVIOUS step's) is deliberately discarded, never touching
  // orbitals.QMAtoms() itself.
  if (!moguess_file_.empty()) {
    XTP_LOG(Log::error, log_) << "Reading MO guess (geometry unchanged) from "
                              << moguess_file_ << std::flush;
    Orbitals moguess_orbitals;
    moguess_orbitals.ReadFromCpt(moguess_file_);
    orbitals.MOs() = moguess_orbitals.MOs();
    if (moguess_orbitals.hasBetaMOs()) {
      orbitals.MOs_beta() = moguess_orbitals.MOs_beta();
    }
    // Defensive, likely-redundant given charge/spin (and therefore the
    // electron counts) should not change between geometry-optimization
    // steps -- DFTEngine::Evaluate() itself may set these independently
    // from package_options_ regardless of what orbitals already holds
    // here. Kept as a safety net rather than removed outright, since
    // this has not been confirmed either way by tracing through
    // DFTEngine's own guess-path logic.
    orbitals.setNumberOfAlphaElectrons(
        moguess_orbitals.getNumberOfAlphaElectrons());
    orbitals.setNumberOfBetaElectrons(
        moguess_orbitals.getNumberOfBetaElectrons());
  }

  std::unique_ptr<QMPackage> qmpackage =
      std::unique_ptr<QMPackage>(QMPackageFactory().Create(
          package_options_.get("name").as<std::string>()));
  qmpackage->setLog(&log_);
  qmpackage->Initialize(package_options_);
  qmpackage->setRunDir(".");

  if (!mpsfile_.empty()) {
    StaticRegion region(0, log_);
    StaticSegment seg = StaticSegment("", 0);
    seg.LoadFromFile(mpsfile_);
    region.push_back(seg);
    qmpackage->AddRegion(region);
  }

  GWBSEEngine gwbse_engine;
  gwbse_engine.setLog(&log_);
  gwbse_engine.setQMPackage(qmpackage.get());
  gwbse_engine.Initialize(gwbseengine_options_, archive_file_);

  QMMolecule fullMol = orbitals.QMAtoms();
  gwbse_engine.ExcitationEnergies(orbitals);
  // If truncation was enabled then rewrite full basis/aux-basis, MOs in full
  // basis and full QMAtoms
  if (orbitals.getCalculationType() == "Truncated") {
    orbitals.QMAtoms().clearAtoms();
    orbitals.QMAtoms() = fullMol;
    orbitals.MOs().eigenvectors() = orbitals.getTruncMOsFullBasis();
    orbitals.SetupDftBasis(orbitals.getDftBasis().Name());
    if (orbitals.hasAuxbasisName()) {
      orbitals.SetupAuxBasis(orbitals.getAuxBasis().Name());
    }
  }

  if (!excited_state_forces_.empty()) {
    ComputeExcitedStateForces(orbitals);
  }
  if (!qp_gradients_.empty()) {
    ComputeQPGradients(orbitals);
  }

  XTP_LOG(Log::error, log_) << "Saving data to " << archive_file_ << std::flush;
  orbitals.WriteToCpt(archive_file_);

  tools::Property summary = gwbse_engine.ReportSummary();
  if (summary.exists("output")) {  // only do gwbse summary output if we
                                   // actually did gwbse
    XTP_LOG(Log::error, log_)
        << "Writing output to " << xml_output_ << std::flush;
    std::ofstream ofout(xml_output_, std::ofstream::out);
    ofout << (summary.get("output"));
    ofout.close();
  }
  return true;
}

// Root following: the state of the requested type with the largest exciton
// overlap (ExcitonOverlaps) with the state whose forces are stored in the
// orb file of an earlier geometry (excited_state_forces_follow).
QMState DftGwBse::FollowState(const Orbitals& orbitals) {
  Orbitals ref;
  ref.ReadFromCpt(excited_state_forces_follow_);
  if (ref.getExcitedStateForcesState().empty()) {
    throw std::runtime_error(
        "excited_state_forces_follow: " + excited_state_forces_follow_ +
        " holds no excited-state forces");
  }
  const QMState ref_state(ref.getExcitedStateForcesState());
  const Eigen::VectorXd ov =
      ExcitonOverlaps(ref, ref_state, orbitals, ref_state.Type());
  Index best = 0;
  ov.maxCoeff(&best);
  double second = 0.0;
  for (Index i = 0; i < ov.size(); ++i) {
    if (i != best) {
      second = std::max(second, ov(i));
    }
  }
  const QMState state(ref_state.Type(), best, false);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Root following: " << ref_state.ToString() << " of "
      << excited_state_forces_follow_ << " is " << state.ToString()
      << " here (overlap " << ov(best) << ", next largest " << second << ")"
      << std::flush;
  if (ov(best) < 0.5 || second > 0.8 * ov(best)) {
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " WARNING: root following is ambiguous (states mix strongly or "
           "the geometry step is large)"
        << std::flush;
  }
  return state;
}

// Forces -dE/dR of one BSE state: the ground-state forces of the SCF
// (xtpdft.compute_forces) minus the analytic gradient of the excitation
// energy with fixed QP corrections (DFTEngine::BSEExcitationGradient).
// Response of the QP corrections in the gradients (BSE states and QP
// levels): the options of the GW run that the gradient has to reproduce.
void DftGwBse::SetupQPResponse(DFTEngine::BSEGradientInput& input,
                               const Orbitals& orbitals) {
  // response of the QP corrections: implemented for G0W0 and evGW with
  // the plasmon-pole and the exact self-energy; otherwise (cda) the
  // corrections are kept fixed
  const std::string mode =
      gwbseengine_options_.get("gwbse.gw.mode").as<std::string>();
  const std::string integrator =
      gwbseengine_options_.get("gwbse.gw.sigma_integrator").as<std::string>();
  if (excited_state_forces_qp_response_ && (mode == "G0W0" || mode == "evGW") &&
      (integrator == "ppm" || integrator == "exact")) {
    input.qp_response = true;
    input.sigma_integrator = integrator;
    input.qp_selfconsistent = (mode == "evGW");
    if (input.qp_selfconsistent) {
      input.out_of_window_shift =
          gwbseengine_options_.get("gwbse.gw.out_of_window_shift")
              .as<std::string>();
      if (input.out_of_window_shift == "max" &&
          (orbitals.getGWAmax() < orbitals.getRPAmax() ||
           orbitals.getGWAmin() > orbitals.getRPAmin())) {
        // max of the window corrections: the energy has kinks where the
        // level with the largest correction changes; the forces are those
        // of the branch at the geometry, and optimizations can stall at
        // such a kink (formaldehyde S1, PBE0/def2-TZVPP)
        XTP_LOG(Log::error, log_)
            << TimeStamp()
            << " WARNING: evGW with gw.out_of_window_shift max: the energy "
               "is not differentiable where the level with the largest "
               "correction changes; use nearest for geometry optimizations"
            << std::flush;
      }
      input.out_of_window_fit_fraction =
          gwbseengine_options_.get("gwbse.gw.out_of_window_fit_fraction")
              .as<double>();
    }
    input.qpmin = orbitals.getGWAmin();
    input.qpmax = orbitals.getGWAmax();
    input.qp_energies = orbitals.QPpertEnergies();
    input.sigma_eta = gwbseengine_options_.get("gwbse.gw.eta").as<double>();
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Including the response of the QP corrections"
        << std::flush;
  } else {
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " QP corrections kept fixed (excited_state_forces_qp_response "
           "false, or not G0W0/evGW with sigma_integrator ppm or exact)"
        << std::flush;
  }
}

// "homo-2:lumo+3", "lumo", "12:20" (absolute, 0-based MO indices)
std::vector<Index> DftGwBse::ParseLevels(const std::string& text, Index homo) {
  auto parse_one = [homo](std::string token) {
    token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
                token.end());
    std::transform(token.begin(), token.end(), token.begin(), ::tolower);
    Index base = 0;
    std::string rest = token;
    if (token.rfind("homo", 0) == 0) {
      base = homo;
      rest = token.substr(4);
    } else if (token.rfind("lumo", 0) == 0) {
      base = homo + 1;
      rest = token.substr(4);
    }
    if (rest.empty() && rest.size() != token.size()) {
      return base;
    }
    std::size_t used = 0;
    Index offset = 0;
    try {
      offset = Index(std::stol(rest, &used));
    } catch (std::exception&) {
      used = std::string::npos;
    }
    if (used != rest.size() || rest.empty() ||
        (rest.size() != token.size() && rest[0] != '+' && rest[0] != '-')) {
      throw std::runtime_error("qp_gradients: cannot read level '" + token +
                               "' (e.g. homo-2, lumo+3 or an MO index)");
    }
    return base + offset;
  };
  std::vector<Index> levels;
  std::stringstream list(text);
  std::string item;
  while (std::getline(list, item, ',')) {
    const std::size_t colon = item.find(':');
    const Index first = parse_one(item.substr(0, colon));
    const Index last = (colon == std::string::npos)
                           ? first
                           : parse_one(item.substr(colon + 1));
    if (last < first) {
      throw std::runtime_error("qp_gradients: empty range '" + item + "'");
    }
    for (Index level = first; level <= last; ++level) {
      if (std::find(levels.begin(), levels.end(), level) == levels.end()) {
        levels.push_back(level);
      }
    }
  }
  std::sort(levels.begin(), levels.end());
  return levels;
}

void DftGwBse::ComputeQPGradients(Orbitals& orbitals) {
  if (package_options_.get("name").as<std::string>() != "xtp") {
    throw std::runtime_error("qp_gradients: needs the xtp DFT package");
  }
  if (!mpsfile_.empty() || orbitals.getCalculationType() == "Truncated") {
    throw std::runtime_error(
        "qp_gradients: embedding and truncated bases are not supported");
  }
  const std::vector<Index> requested =
      ParseLevels(qp_gradients_, orbitals.getHomo());
  // levels whose gradients are already in the orb file (tasks gradients on
  // a file of an earlier run) are kept and not recomputed
  const std::vector<Index> stored = orbitals.getQPGradientLevels();
  const Eigen::MatrixXd stored_gradients = orbitals.getQPGradients();
  std::vector<Index> levels;
  for (Index level : requested) {
    if (std::find(stored.begin(), stored.end(), level) == stored.end()) {
      levels.push_back(level);
    }
  }
  if (levels.empty()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " QP gradients of all requested levels are "
        << "already stored" << std::flush;
    return;
  }
  if (levels.size() < requested.size()) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " QP gradients of "
        << requested.size() - levels.size()
        << " requested level(s) are already stored" << std::flush;
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Analytic gradients of the QP energies of "
      << levels.size() << " level(s), " << levels.front() << " to "
      << levels.back() << " (HOMO is " << orbitals.getHomo() << ")"
      << std::flush;
  DFTEngine::BSEGradientInput input =
      DFTEngine::QPGradientInputFromOrbitals(orbitals, levels);
  SetupQPResponse(input, orbitals);
  DFTEngine dft;
  dft.setLogger(&log_);
  dft.Initialize(package_options_);
  const DFTEngine::BSEGradientResult result =
      dft.BSEExcitationGradient(orbitals, input);
  if (input.qp_response) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " QP energies rebuilt for the gradient: max "
        << "deviation " << result.qp_rebuild_error << " Ha" << std::flush;
    if (result.qp_rebuild_error > 1e-4) {
      throw std::runtime_error(
          "qp_gradients: the QP energies rebuilt for the gradient differ "
          "from the GW ones by " +
          std::to_string(result.qp_rebuild_error) + " Ha");
    }
  }
  const Index natoms = orbitals.QMAtoms().size();
  // stored and new levels, sorted
  std::vector<std::pair<Index, Eigen::MatrixXd>> merged;
  for (std::size_t l = 0; l < stored.size(); ++l) {
    merged.emplace_back(stored[l],
                        stored_gradients.middleRows(Index(l) * natoms, natoms));
  }
  for (std::size_t l = 0; l < levels.size(); ++l) {
    merged.emplace_back(levels[l], result.qp_gradients[l]);
  }
  std::sort(merged.begin(), merged.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  Eigen::MatrixXd all(Index(merged.size()) * natoms, 3);
  std::vector<Index> all_levels;
  for (std::size_t l = 0; l < merged.size(); ++l) {
    all.middleRows(Index(l) * natoms, natoms) = merged[l].second;
    all_levels.push_back(merged[l].first);
  }
  orbitals.setQPGradients(all, all_levels);
}

// tasks = gradients: the options must be those of the run that wrote the
// orb file. What is stored there is compared; a different GW setup shows
// up in the rebuilt QP and excitation energies of the gradient.
void DftGwBse::CheckStoredRun(const Orbitals& orbitals) const {
  auto lower = [](std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), ::tolower);
    return text;
  };
  auto check = [&](const std::string& what, const std::string& stored,
                   const std::string& option) {
    if (lower(stored) != lower(option)) {
      throw std::runtime_error("tasks gradients: " + what + " of " +
                               archive_file_ + " is " + stored +
                               ", the options say " + option +
                               "; use the options of the run that wrote "
                               "the file");
    }
  };
  if (!orbitals.hasQMAtoms() || orbitals.MOs().eigenvalues().size() == 0 ||
      orbitals.RPAInputEnergies().size() == 0) {
    throw std::runtime_error("tasks gradients: " + archive_file_ +
                             " has no GW data");
  }
  check("the basis set", orbitals.getDftBasis().Name(),
        package_options_.get("basisset").as<std::string>());
  check("the functional", orbitals.getXCFunctionalName(),
        package_options_.get("functional").as<std::string>());
  if (orbitals.hasAuxbasisName()) {
    check("the auxiliary basis set", orbitals.getAuxBasis().Name(),
          package_options_.get("auxbasisset").as<std::string>());
  }
}

void DftGwBse::ComputeExcitedStateForces(Orbitals& orbitals) {
  QMState state(excited_state_forces_);
  if (!excited_state_forces_follow_.empty()) {
    state = FollowState(orbitals);
  }
  const bool singlet = state.Type() == QMStateType::Singlet;
  if (!singlet && state.Type() != QMStateType::Triplet) {
    throw std::runtime_error(
        "excited_state_forces: only singlet or triplet BSE states (e.g. s1)");
  }
  if (package_options_.get("name").as<std::string>() != "xtp") {
    throw std::runtime_error("excited_state_forces: needs the xtp DFT package");
  }
  if (!mpsfile_.empty() || orbitals.getCalculationType() == "Truncated") {
    throw std::runtime_error(
        "excited_state_forces: embedding and truncated bases are not "
        "supported");
  }
  if (!orbitals.hasForces()) {
    throw std::runtime_error(
        "excited_state_forces: needs the ground-state forces, set "
        "dftpackage.xtpdft.compute_forces to true");
  }
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Analytic forces of BSE state " << state.ToString()
      << (orbitals.getTDAApprox() ? " (TDA)" : " (full BSE)") << std::flush;

  DFTEngine::BSEGradientInput input = DFTEngine::BSEGradientInputFromOrbitals(
      orbitals, singlet, state.StateIdx());
  SetupQPResponse(input, orbitals);
  DFTEngine dft;
  dft.setLogger(&log_);
  dft.Initialize(package_options_);
  const DFTEngine::BSEGradientResult result =
      dft.BSEExcitationGradient(orbitals, input);

  const double omega =
      (singlet ? orbitals.BSESinglets() : orbitals.BSETriplets())
          .eigenvalues()(state.StateIdx());
  // the gradient belongs to the Rayleigh quotient of the stored exciton,
  // which differs from the Davidson eigenvalue at the level of the
  // Davidson residual (bse.davidson.tolerance; for full BSE to first
  // order); a different GW-BSE setup gives deviations of 1e-3 Ha and more
  const double omega_error = std::abs(result.omega - omega);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Excitation energy rebuilt for the gradient: "
      << "deviation from the BSE eigenvalue " << omega_error << " Ha"
      << std::flush;
  if (omega_error > 1e-4) {
    throw std::runtime_error(
        "excited_state_forces: the excitation energy rebuilt for the "
        "gradient (" +
        std::to_string(result.omega) +
        " Ha) differs from the BSE eigenvalue (" + std::to_string(omega) +
        " Ha); the GW-BSE setup (e.g. environment screening, dynamical "
        "screening of the BSE kernel) is not covered by the gradient");
  }
  if (omega_error > 1e-6) {
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " WARNING: the BSE eigenvalue is converged only to about "
        << omega_error
        << " Ha; tighten bse.davidson.tolerance (strict or lapack) for "
           "energies consistent with the forces"
        << std::flush;
  }
  if (input.qp_response) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " QP energies rebuilt for the gradient: max "
        << "deviation " << result.qp_rebuild_error << " Ha" << std::flush;
  }
  if (input.qp_selfconsistent) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " RPA energies rebuilt for the gradient: max "
        << "deviation " << result.rpa_rebuild_error << " Ha" << std::flush;
  }
  // the QP solver and the interpolated pole sums of GW are converged to
  // about 1e-6 Ha; a larger deviation means a different setup
  if (input.qp_response && result.qp_rebuild_error > 1e-4) {
    throw std::runtime_error(
        "excited_state_forces: the QP energies rebuilt for the gradient "
        "differ from the GW ones by " +
        std::to_string(result.qp_rebuild_error) + " Ha");
  }
  const Eigen::MatrixXd forces = orbitals.getForces() - result.gradient;
  orbitals.setExcitedStateForces(forces, state.ToString());
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Excited-state forces (Hartree/Bohr) of "
      << state.ToString() << ", Omega = " << omega << " Ha:\n"
      << forces << std::flush;
}

}  // namespace xtp
}  // namespace votca
