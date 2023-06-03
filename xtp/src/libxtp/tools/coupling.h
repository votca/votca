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

  std::string Identify() const { return "coupling"; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Run();

 private:
  std::string MOsA_, MOsB_, MOsAB_;
  std::string logA_, logB_, logAB_;

  tools::Property package_options_;
  tools::Property dftcoupling_options_;

  std::string output_file_;

  Logger log_;
};

void Coupling::ParseOptions(const tools::Property &options) {

  MOsA_ = options.get(".moleculeA.orbitals").as<std::string>();
  MOsB_ = options.get(".moleculeB.orbitals").as<std::string>();
  MOsAB_ = options.get(".dimerAB.orbitals").as<std::string>();

  logA_ = options.get(".moleculeA.log").as<std::string>();
  logB_ = options.get(".moleculeB.log").as<std::string>();
  logAB_ = options.get(".dimerAB.log").as<std::string>();

  output_file_ = options.ifExistsReturnElseReturnDefault<std::string>(
      "output", job_name_ + " coupling_.xml");

  package_options_ = options.get(".dftpackage");
  dftcoupling_options_ = options.get(".dftcoupling_options");

  QMPackageFactory{};
}

bool Coupling::Run() {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  // get the corresponding object from the QMPackageFactory
  std::unique_ptr<QMPackage> qmpackage =
      std::unique_ptr<QMPackage>(QMPackageFactory().Create(
          package_options_.get("name").as<std::string>()));
  qmpackage->setLog(&log_);
  qmpackage->Initialize(package_options_);
  qmpackage->setRunDir(".");
  Orbitals orbitalsA, orbitalsB, orbitalsAB;

  qmpackage->setLogFileName(logA_);
  bool parse_logA_status = qmpackage->ParseLogFile(orbitalsA);
  if (!parse_logA_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read log of molecule A" << std::flush;
  }

  qmpackage->setLogFileName(logB_);
  bool parse_logB_status = qmpackage->ParseLogFile(orbitalsB);
  if (!parse_logB_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read log of molecule B" << std::flush;
  }

  qmpackage->setLogFileName(logAB_);
  bool parse_logAB_status = qmpackage->ParseLogFile(orbitalsAB);
  if (!parse_logAB_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read log of molecule AB" << std::flush;
  }

  qmpackage->setMOsFileName(MOsA_);
  bool parse_orbitalsA_status = qmpackage->ParseMOsFile(orbitalsA);
  if (!parse_orbitalsA_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read orbitals of molecule A" << std::flush;
  }

  qmpackage->setMOsFileName(MOsB_);
  bool parse_orbitalsB_status = qmpackage->ParseMOsFile(orbitalsB);
  if (!parse_orbitalsB_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read orbitals of molecule B" << std::flush;
  }

  qmpackage->setMOsFileName(MOsAB_);
  bool parse_orbitalsAB_status = qmpackage->ParseMOsFile(orbitalsAB);
  if (!parse_orbitalsAB_status) {
    XTP_LOG(Log::error, log_)
        << "Failed to read orbitals of dimer AB" << std::flush;
  }

  DFTcoupling dftcoupling;
  dftcoupling.setLogger(&log_);
  dftcoupling.Initialize(dftcoupling_options_);

  dftcoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
  std::cout << log_;

  // output the results
  tools::Property summary;
  tools::Property &job_output = summary.add("output", "");
  tools::Property &pair_summary = job_output.add("pair", "");
  dftcoupling.Addoutput(pair_summary, orbitalsA, orbitalsB);
  std::ofstream ofs(output_file_, std::ofstream::out);
  ofs << job_output;
  ofs.close();

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_COUPLING_H
