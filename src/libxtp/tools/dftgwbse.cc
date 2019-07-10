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
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/staticregion.h>

using namespace std;

namespace votca {
namespace xtp {

void DftGwBse::Initialize(tools::Property& options) {

  _do_optimize = false;
  std::string key = "options." + Identify();

  if (options.exists(key + ".mpsfile")) {
    _do_external = true;
    _mpsfile = options.get(key + ".mpsfile").as<string>();
  } else {
    _do_external = false;
  }

  if (options.exists(key + ".guess")) {
    _do_guess = true;
    _guess_file = options.get(key + ".guess").as<string>();

  } else {
    _do_guess = false;
  }

  _archive_file = options.ifExistsReturnElseReturnDefault<string>(
      key + ".archive", "system.orb");
  _reporting = options.ifExistsReturnElseReturnDefault<string>(
      key + ".reporting", "default");

  // job tasks
  std::vector<string> choices = {"optimize", "energy"};
  string mode = options.ifExistsAndinListReturnElseThrowRuntimeError<string>(
      key + ".mode", choices);
  if (mode == "optimize") _do_optimize = true;

  // GWBSEENGINE options
  _gwbseengine_options = options.get(key + ".gwbse_engine");

  // options for dft package
  string _package_xml = options.get(key + ".dftpackage").as<string>();
  load_property_from_xml(_package_options, _package_xml);
  _package = _package_options.get("package.name").as<string>();

  // MOLECULE properties
  _xyzfile =
      options.ifExistsReturnElseThrowRuntimeError<string>(key + ".molecule");

  // XML OUTPUT
  _xml_output = options.ifExistsReturnElseReturnDefault<string>(
      key + ".output", "dftgwbse.out.xml");

  // if optimization is chosen, get options for geometry_optimizer
  if (_do_optimize)
    _geoopt_options = options.get(key + ".geometry_optimization");

  // register all QM packages (Gaussian, NWCHEM, etc)
  QMPackageFactory::RegisterAll();
}

bool DftGwBse::Evaluate() {
  OPENMP::setMaxThreads(_nThreads);
  if (_reporting == "silent")
    _log.setReportLevel(logERROR);  // only output ERRORS, GEOOPT info, and
                                    // excited state info for trial geometry
  if (_reporting == "noisy")
    _log.setReportLevel(logDEBUG);  // OUTPUT ALL THE THINGS
  if (_reporting == "default") _log.setReportLevel(logINFO);  //

  _log.setMultithreading(true);
  _log.setPreface(logINFO, "\n... ...");
  _log.setPreface(logERROR, "\n... ...");
  _log.setPreface(logWARNING, "\n... ...");
  _log.setPreface(logDEBUG, "\n... ...");

  // Get orbitals object
  Orbitals orbitals;

  if (_do_guess) {
    XTP_LOG(logDEBUG, _log) << "Reading guess from " << _guess_file << flush;
    orbitals.ReadFromCpt(_guess_file);
  } else {
    XTP_LOG(logDEBUG, _log) << "Reading structure from " << _xyzfile << flush;
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

  XTP_LOG(logDEBUG, _log) << "Saving data to " << _archive_file << flush;
  orbitals.WriteToCpt(_archive_file);

  tools::Property summary = gwbse_engine.ReportSummary();
  if (summary.exists("output")) {  // only do gwbse summary output if we
                                   // actually did gwbse
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1,
                                        "");
    XTP_LOG(logDEBUG, _log) << "Writing output to " << _xml_output << flush;
    std::ofstream ofout(_xml_output.c_str(), std::ofstream::out);
    ofout << (summary.get("output"));
    ofout.close();
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
