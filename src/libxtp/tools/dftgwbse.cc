/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#include "votca/xtp/geometry_optimization.h"
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
  xyzfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".molecule", job_name_ + ".xyz");

  // job tasks
  do_optimize_ = options.get(".optimize").as<bool>();

  // options for dft package
  package_options_ = options.get(".dftpackage");
  package_options_.add("job_name", job_name_);
  package_ = package_options_.get("package.name").as<std::string>();

  // set the basis sets and functional in DFT package
  package_options_.get("package").add(
      "basisset", options.get("basisset").as<std::string>());
  package_options_.get("package").add(
      "auxbasisset", options.get("auxbasisset").as<std::string>());
  package_options_.get("package").add(
      "functional", options.get("functional").as<std::string>());

  // GWBSEENGINE options
  gwbseengine_options_ = options.get(".gwbse_engine");

  // set the basis sets and functional in GWBSE
  gwbseengine_options_.get("gwbse_options.gwbse")
      .add("basisset", options.get("basisset").as<std::string>());
  gwbseengine_options_.get("gwbse_options.gwbse")
      .add("auxbasisset", options.get("auxbasisset").as<std::string>());
  gwbseengine_options_.get("gwbse_options.gwbse.vxc")
      .add("functional", options.get("functional").as<std::string>());

  // lets get the archive file name from the xyz file name
  archive_file_ = job_name_ + ".orb";

  // XML OUTPUT
  xml_output_ = job_name_ + "_summary.xml";

  // check for MPS file with external multipoles for embedding
  do_external_ = options.get("use_mpsfile").as<bool>();
  if (do_external_) {
    mpsfile_ = options.get(".mpsfile").as<std::string>();
  }

  // check if guess is requested
  do_guess_ = options.get("use_guess").as<bool>();
  if (do_guess_) {
    guess_file_ = options.get(".guess").as<std::string>();
  }

  // if optimization is chosen, get options for geometry_optimizer
  if (do_optimize_) {
    geoopt_options_ = options.get(".geometry_optimization");
  }

  // register all QM packages
  QMPackageFactory::RegisterAll();
}

bool DftGwBse::Run() {

  log_.setReportLevel(Log::current_level);

  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (do_guess_) {
    XTP_LOG(Log::error, log_)
        << "Reading guess from " << guess_file_ << std::flush;
    orbitals.ReadFromCpt(guess_file_);
  } else {
    XTP_LOG(Log::error, log_)
        << "Reading structure from " << xyzfile_ << std::flush;
    orbitals.QMAtoms().LoadFromFile(xyzfile_);
  }

  std::unique_ptr<QMPackage> qmpackage = std::unique_ptr<QMPackage>(
      QMPackageFactory::QMPackages().Create(package_));
  qmpackage->setLog(&log_);
  qmpackage->Initialize(package_options_);
  qmpackage->setRunDir(".");

  if (do_external_) {
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

  if (do_optimize_) {
    GeometryOptimization geoopt(gwbse_engine, orbitals);
    geoopt.setLog(&log_);
    geoopt.Initialize(geoopt_options_);
    geoopt.Evaluate();
  } else {
    gwbse_engine.ExcitationEnergies(orbitals);
  }

  XTP_LOG(Log::error, log_) << "Saving data to " << archive_file_ << std::flush;
  orbitals.WriteToCpt(archive_file_);

  tools::Property summary = gwbse_engine.ReportSummary();
  if (summary.exists("output")) {  // only do gwbse summary output if we
                                   // actually did gwbse
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1,
                                        "");
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
