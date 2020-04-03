/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include "dftgwbse.h"
#include <votca/tools/constants.h>
#include <votca/tools/filesystem.h>
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/staticregion.h>

using namespace std;

namespace votca {
namespace xtp {

void DftGwBse::Initialize(tools::Property& options) {

  // get pre-defined default options from VOTCASHARE/xtp/xml/dftgwbse.xml
  LoadDefaults("xtp");
  // update options with user specified input
  UpdateWithUserOptions(options);

  // molecule coordinates
  _xyzfile = _options.ifExistsReturnElseThrowRuntimeError<string>(".molecule");

  // job tasks
  std::vector<string> choices = {"optimize", "energy"};
  string mode = _options.ifExistsAndinListReturnElseThrowRuntimeError<string>(
      ".mode", choices);

  // options for dft package
  if (_options.exists("dftpackage")) {
    _package_options = _options.get(".dftpackage");
    _package = _package_options.get("package.name").as<string>();
  } else {
    _package_options.add("package", "");
    _package = "xtp";
  }

  // set the basis sets and functional in DFT package
  _package_options.get("package").add("basisset",
                                      _options.get("basisset").as<string>());
  _package_options.get("package").add("auxbasisset",
                                      _options.get("auxbasisset").as<string>());
  _package_options.get("package").add("functional",
                                      _options.get("functional").as<string>());

  // GWBSEENGINE options
  if (_options.exists("gwbse_engine")) {
    _gwbseengine_options = _options.get(".gwbse_engine");
  } else {
    _package_options.add("gwbse_engine", "");
  }

  // set the basis sets and functional in GWBSE
  _gwbseengine_options.get("gwbse_options.gwbse")
      .add("basisset", _options.get("basisset").as<string>());
  _gwbseengine_options.get("gwbse_options.gwbse")
      .add("auxbasisset", _options.get("auxbasisset").as<string>());
  _gwbseengine_options.get("gwbse_options.gwbse.vxc")
      .add("functional", _options.get("functional").as<string>());

  // lets get the archive file name from the xyz file name
  _archive_file = tools::filesystem::GetFileBase(_xyzfile) + ".orb";

  // XML OUTPUT
  _xml_output = tools::filesystem::GetFileBase(_xyzfile) + "_summary.xml";

  // checking for additional requests
  _do_optimize = false;
  _do_external = false;
  _do_guess = false;

  // check for MPS file with external multipoles for embedding
  if (_options.exists(".mpsfile")) {
    _do_external = true;
    _mpsfile = options.get(".mpsfile").as<string>();
  }

  // check if guess is requested
  if (_options.exists(".guess")) {
    _do_guess = true;
    _guess_file = _options.get(".guess").as<string>();
  }

  // if optimization is chosen, get options for geometry_optimizer
  if (mode == "optimize") {
    _do_optimize = true;
    _geoopt_options = _options.get(".geometry_optimization");
  }

  // register all QM packages
  QMPackageFactory::RegisterAll();
}

bool DftGwBse::Evaluate() {
  OPENMP::setMaxThreads(_nThreads);

  _log.setReportLevel(Log::current_level);

  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (_do_guess) {
    XTP_LOG(Log::error, _log) << "Reading guess from " << _guess_file << flush;
    orbitals.ReadFromCpt(_guess_file);
  } else {
    XTP_LOG(Log::error, _log) << "Reading structure from " << _xyzfile << flush;
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

  XTP_LOG(Log::error, _log) << "Saving data to " << _archive_file << flush;
  orbitals.WriteToCpt(_archive_file);

  tools::Property summary = gwbse_engine.ReportSummary();
  if (summary.exists("output")) {  // only do gwbse summary output if we
                                   // actually did gwbse
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1,
                                        "");
    XTP_LOG(Log::error, _log) << "Writing output to " << _xml_output << flush;
    std::ofstream ofout(_xml_output, std::ofstream::out);
    ofout << (summary.get("output"));
    ofout.close();
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
