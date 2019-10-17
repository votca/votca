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
#ifndef _VOTCA_XTP_PARTIALCHARGES_H
#define _VOTCA_XTP_PARTIALCHARGES_H

#include <boost/filesystem.hpp>
#include <stdio.h>
#include <votca/xtp/esp2multipole.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {

class Partialcharges : public QMTool {
 public:
  Partialcharges() = default;
  ;
  ~Partialcharges() override = default;
  ;

  std::string Identify() override { return "partialcharges"; }

  void Initialize(tools::Property& options) override;
  bool Evaluate() override;

 private:
  std::string _orbfile;
  std::string _output_file;
  tools::Property _esp_options;

  Logger _log;
};

void Partialcharges::Initialize(tools::Property& options) {

  std::string key = "options." + Identify();
  _orbfile = options.get(key + ".input").as<std::string>();
  _output_file = options.get(key + ".output").as<std::string>();
  std::string _esp2multipole_xml =
      options.get(key + ".esp_options").as<std::string>();
  _esp_options.LoadFromXML(_esp2multipole_xml);
}

bool Partialcharges::Evaluate() {
  OPENMP::setMaxThreads(_nThreads);
  _log.setReportLevel(logDEBUG);
  _log.setMultithreading(true);

  _log.setPreface(logINFO, "\n... ...");
  _log.setPreface(logERROR, "\n... ...");
  _log.setPreface(logWARNING, "\n... ...");
  _log.setPreface(logDEBUG, "\n... ...");
  Orbitals orbitals;
  XTP_LOG_SAVE(logDEBUG, _log)
      << " Loading QM data from " << _orbfile << std::flush;
  orbitals.ReadFromCpt(_orbfile);
  XTP_LOG_SAVE(logDEBUG, _log)
      << "Loaded QM data from " << _orbfile << std::flush;
  Esp2multipole esp2multipole = Esp2multipole(_log);
  esp2multipole.Initialize(_esp_options);
  StaticSegment seg = esp2multipole.Extractingcharges(orbitals);
  seg.WriteMPS(_output_file, esp2multipole.GetStateString());

  XTP_LOG_SAVE(logDEBUG, _log)
      << "Written charges to " << _output_file << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif
