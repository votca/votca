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
  _xyzfile = options.ifExistsReturnElseReturnDefault<std::string>(
      ".molecule", _job_name + ".xyz");

  // job tasks
  _do_optimize = options.get(".optimize").as<bool>();

  // options for dft package
  _package_options = options.get(".dftpackage");
  _package_options.add("job_name", _job_name);
  _package = _package_options.get("package.name").as<std::string>();

  // set the basis sets and functional in DFT package
  _package_options.get("package").add(
      "basisset", options.get("basisset").as<std::string>());
  _package_options.get("package").add(
      "auxbasisset", options.get("auxbasisset").as<std::string>());
  _package_options.get("package").add(
      "functional", options.get("functional").as<std::string>());

  // GWBSEENGINE options
  _gwbseengine_options = options.get(".gwbse_engine");

  // set the basis sets and functional in GWBSE
  _gwbseengine_options.get("gwbse_options.gwbse")
      .add("basisset", options.get("basisset").as<std::string>());
  _gwbseengine_options.get("gwbse_options.gwbse")
      .add("auxbasisset", options.get("auxbasisset").as<std::string>());
  _gwbseengine_options.get("gwbse_options.gwbse.vxc")
      .add("functional", options.get("functional").as<std::string>());

  // lets get the archive file name from the xyz file name
  _archive_file = _job_name + ".orb";

  // XML OUTPUT
  _xml_output = _job_name + "_summary.xml";

  // check for MPS file with external multipoles for embedding
  _do_external = options.get("use_mpsfile").as<bool>();
  if (_do_external) {
    _mpsfile = options.get(".mpsfile").as<std::string>();
  }

  // check if guess is requested
  _do_guess = options.get("use_guess").as<bool>();
  if (_do_guess) {
    _guess_file = options.get(".guess").as<std::string>();
  }

  // if optimization is chosen, get options for geometry_optimizer
  if (_do_optimize) {
    _geoopt_options = options.get(".geometry_optimization");
  }

  // register all QM packages
  QMPackageFactory::RegisterAll();
}

bool DftGwBse::Run() {

  _log.setReportLevel(Log::current_level);

  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (_do_guess) {
    XTP_LOG(Log::error, _log)
        << "Reading guess from " << _guess_file << std::flush;
    orbitals.ReadFromCpt(_guess_file);
  } else {
    XTP_LOG(Log::error, _log)
        << "Reading structure from " << _xyzfile << std::flush;
    orbitals.QMAtoms().LoadFromFile(_xyzfile);
  }

  std::unique_ptr<QMPackage> qmpackage =
      std::unique_ptr<QMPackage>(QMPackages().Create(_package));
  qmpackage->setLog(&_log);
  qmpackage->Initialize(_package_options);
  qmpackage->setRunDir(".");

  if (_do_external) {
    StaticRegion region(0, _log);
    StaticSegment seg = StaticSegment("", 0);
    seg.LoadFromFile(_mpsfile);
    region.push_back(seg);
    qmpackage->AddRegion(region);
  }

  GWBSEEngine gwbse_engine;
  gwbse_engine.setLog(&_log);
  gwbse_engine.setQMPackage(qmpackage.get());
  gwbse_engine.Initialize(_gwbseengine_options, _archive_file);

  if (_do_optimize) {
    GeometryOptimization geoopt(gwbse_engine, orbitals);
    geoopt.setLog(&_log);
    geoopt.Initialize(_geoopt_options);
    geoopt.Evaluate();
  } else {
    gwbse_engine.ExcitationEnergies(orbitals);
  }

  XTP_LOG(Log::error, _log) << "Saving data to " << _archive_file << std::flush;
  orbitals.WriteToCpt(_archive_file);

  tools::Property summary = gwbse_engine.ReportSummary();
  if (summary.exists("output")) {  // only do gwbse summary output if we
                                   // actually did gwbse
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1,
                                        "");
    XTP_LOG(Log::error, _log)
        << "Writing output to " << _xml_output << std::flush;
    std::ofstream ofout(_xml_output, std::ofstream::out);
    ofout << (summary.get("output"));
    ofout.close();
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
