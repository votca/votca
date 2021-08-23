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

  // if optimization is chosen, get options for geometry_optimizer
  if (options.exists(".geometry_optimization")) {
    do_optimize_ = true;
    geoopt_options_ = options.get(".geometry_optimization");
  }
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

  std::unique_ptr<QMPackage> qmpackage =QMPackageFactory().Create(
          package_options_.get("name").as<std::string>());
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
