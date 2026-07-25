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
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/gwbseengine.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/staticregion.h"

// Local private VOTCA includes
#include "dftgwbse.h"

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
    XTP_LOG(Log::error, log_)
        << "Reading MO guess (geometry unchanged) from " << moguess_file_
        << std::flush;
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

}  // namespace xtp
}  // namespace votca
