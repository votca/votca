/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/format.hpp>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/qmtool.h>

namespace votca {
namespace xtp {

class Log2Mps : public QMTool {
 public:
  Log2Mps(){};
  ~Log2Mps(){};

  std::string Identify() { return "log2mps"; }

  void Initialize(tools::Property &options);
  bool Evaluate();

 private:
  std::string _package;
  std::string _logfile;
  std::string _mpsfile;
};

void Log2Mps::Initialize(tools::Property &opt) {

  QMPackageFactory::RegisterAll();

  std::string key = "options.log2mps";
  _package = opt.get(key + ".package").as<std::string>();

  if (_package == "xtp") {
    throw std::runtime_error(
        "XTP has no log file. For xtp package just run the partialcharges tool "
        "on you .orb file");
  }
  _logfile = opt.get(key + ".logfile").as<std::string>();

  _mpsfile = (opt.exists(key + ".mpsfile"))
                 ? opt.get(key + ".mpsfile").as<std::string>()
                 : "";
  if (_mpsfile == "")
    _mpsfile = _logfile.substr(0, _logfile.size() - 4) + ".mps";

  std::cout << std::endl
            << "... ... " << _logfile << " => " << _mpsfile << std::flush;
}

bool Log2Mps::Evaluate() {

  // Logger (required for QM package, so we can just as well use it)
  Logger log;
  log.setPreface(logINFO, "\n... ...");
  log.setPreface(logDEBUG, "\n... ...");
  log.setReportLevel(logDEBUG);
  log.setMultithreading(true);

  // Set-up QM package
  XTP_LOG_SAVE(logINFO, log)
      << "Using package <" << _package << ">" << std::flush;

  std::unique_ptr<QMPackage> qmpack =
      std::unique_ptr<QMPackage>(QMPackages().Create(_package));
  qmpack->setLog(&log);
  qmpack->setRunDir(".");
  qmpack->setLogFileName(_logfile);

  // Create orbitals, fill with life & extract QM atoms

  StaticSegment atoms = qmpack->GetCharges();

  // Sanity checks, total charge

  if (atoms.size() < 1) {
    std::cout << "\nERROR No charges extracted from " << _logfile
              << ". Abort.\n"
              << std::flush;
    throw std::runtime_error("(see above, input or parsing error)");
  }

  double Q = atoms.CalcTotalQ();
  XTP_LOG_SAVE(logINFO, log)
      << atoms.size() << " QM atoms, total charge Q = " << Q << std::flush;

  std::string tag =
      "::LOG2MPS " + (boost::format("(log-file='%1$s' : %2$d QM atoms)") %
                      _logfile % atoms.size())
                         .str();
  atoms.WriteMPS(_mpsfile, tag);
  return true;
}

}  // namespace xtp
}  // namespace votca

#endif
