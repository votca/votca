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

#pragma once
#ifndef VOTCA_XTP_COUPLING_H
#define VOTCA_XTP_COUPLING_H

// Local VOTCA includes
#include "votca/xtp/dftcoupling.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmpackagefactory.h"

namespace votca {
namespace xtp {

class Coupling final : public QMTool {
 public:
  Coupling() = default;
  ~Coupling() = default;

  std::string Identify() { return "coupling"; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Run();

 private:
  std::string _MOsA, _MOsB, _MOsAB;
  std::string _logA, _logB, _logAB;

  std::string _package;
  tools::Property _package_options;
  tools::Property _dftcoupling_options;

  std::string _output_file;

  Logger _log;
};

void Coupling::ParseOptions(const tools::Property &options) {

  _MOsA = options.get(".moleculeA.orbitals").as<std::string>();
  _MOsB = options.get(".moleculeB.orbitals").as<std::string>();
  _MOsAB = options.get(".dimerAB.orbitals").as<std::string>();

  _logA = options.get(".moleculeA.log").as<std::string>();
  _logB = options.get(".moleculeB.log").as<std::string>();
  _logAB = options.get(".dimerAB.log").as<std::string>();

  _output_file = options.ifExistsReturnElseReturnDefault<std::string>(
      "output", _job_name + "_coupling.xml");

  _package_options = options.get(".dftpackage");
  _package = _package_options.get("package.name").as<std::string>();
  _dftcoupling_options = options.get(".dftcoupling_options");
}

bool Coupling::Run() {

  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);

  _log.setCommonPreface("\n... ...");

  // get the corresponding object from the QMPackageFactory
  QMPackageFactory factory;
  std::unique_ptr<QMPackage> qmpackage = factory.Create(_package);
  qmpackage->setLog(&_log);
  qmpackage->Initialize(_package_options);
  qmpackage->setRunDir(".");
  Orbitals orbitalsA, orbitalsB, orbitalsAB;

  qmpackage->setLogFileName(_logA);
  bool parse_logA_status = qmpackage->ParseLogFile(orbitalsA);
  if (!parse_logA_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read log of molecule A" << std::flush;
  }

  qmpackage->setLogFileName(_logB);
  bool parse_logB_status = qmpackage->ParseLogFile(orbitalsB);
  if (!parse_logB_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read log of molecule B" << std::flush;
  }

  qmpackage->setLogFileName(_logAB);
  bool parse_logAB_status = qmpackage->ParseLogFile(orbitalsAB);
  if (!parse_logAB_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read log of molecule AB" << std::flush;
  }

  qmpackage->setMOsFileName(_MOsA);
  bool parse_orbitalsA_status = qmpackage->ParseMOsFile(orbitalsA);
  if (!parse_orbitalsA_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read orbitals of molecule A" << std::flush;
  }

  qmpackage->setMOsFileName(_MOsB);
  bool parse_orbitalsB_status = qmpackage->ParseMOsFile(orbitalsB);
  if (!parse_orbitalsB_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read orbitals of molecule B" << std::flush;
  }

  qmpackage->setMOsFileName(_MOsAB);
  bool parse_orbitalsAB_status = qmpackage->ParseMOsFile(orbitalsAB);
  if (!parse_orbitalsAB_status) {
    XTP_LOG(Log::error, _log)
        << "Failed to read orbitals of dimer AB" << std::flush;
  }

  DFTcoupling dftcoupling;
  dftcoupling.setLogger(&_log);
  dftcoupling.Initialize(_dftcoupling_options);

  dftcoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
  std::cout << _log;

  // output the results
  tools::Property summary;
  tools::Property &job_output = summary.add("output", "");
  tools::Property &pair_summary = job_output.add("pair", "");
  dftcoupling.Addoutput(pair_summary, orbitalsA, orbitalsB);

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");

  std::ofstream ofs(_output_file, std::ofstream::out);
  ofs << job_output;
  ofs.close();

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_COUPLING_H
