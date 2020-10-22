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
#ifndef VOTCA_XTP_LOG2MPS_H
#define VOTCA_XTP_LOG2MPS_H

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class Log2Mps final : public QMTool {
 public:
  Log2Mps() = default;
  ~Log2Mps() = default;

  std::string Identify() { return "log2mps"; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Run();

 private:
  std::string _package;
  std::string _logfile;
  std::string _mpsfile;
};

void Log2Mps::ParseOptions(const tools::Property &options) {

  QMPackageFactory::RegisterAll();

  _package = options.get(".package").as<std::string>();

  if (_package == "xtp") {
    throw std::runtime_error(
        "XTP has no log file. For xtp package just run the partialcharges tool "
        "on you .orb file");
  }
  _logfile = options.ifExistsReturnElseReturnDefault<std::string>(
      ".logfile", _job_name + ".log");

  _mpsfile = options.ifExistsReturnElseReturnDefault<std::string>(
      ".mpsfile", _job_name + ".mps");

  std::cout << "\n... ... " << _logfile << " => " << _mpsfile << "\n";
}

bool Log2Mps::Run() {

  // Logger (required for QM package, so we can just as well use it)
  Logger log;
  log.setCommonPreface("\n... ...");
  log.setReportLevel(Log::current_level);
  log.setMultithreading(true);

  // Set-up QM package
  XTP_LOG(Log::error, log) << "Using package <" << _package << ">"
                           << std::flush;

  std::unique_ptr<QMPackage> qmpack =
      std::unique_ptr<QMPackage>(QMPackages().Create(_package));
  qmpack->setLog(&log);
  qmpack->setRunDir(".");
  qmpack->setLogFileName(_logfile);

  // Create orbitals, fill with life & extract QM atoms

  StaticSegment atoms = qmpack->GetCharges();

  // Sanity checks, total charge

  if (atoms.size() < 1) {
    throw std::runtime_error("ERROR No charges extracted from " + _logfile);
  }

  double Q = atoms.CalcTotalQ();
  XTP_LOG(Log::error, log) << atoms.size()
                           << " QM atoms, total charge Q = " << Q << std::flush;

  std::string tag =
      "::LOG2MPS " + (boost::format("(log-file='%1$s' : %2$d QM atoms)") %
                      _logfile % atoms.size())
                         .str();
  atoms.WriteMPS(_mpsfile, tag);
  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LOG2MPS_H
