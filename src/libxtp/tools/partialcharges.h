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

#ifndef _VOTCA_XTP_PARTIALCHARGES_H
#define _VOTCA_XTP_PARTIALCHARGES_H

#include <boost/filesystem.hpp>
#include <stdio.h>
#include <votca/ctp/logger.h>
#include <votca/xtp/esp2multipole.h>

namespace votca {
namespace xtp {
using namespace std;

class Partialcharges : public ctp::QMTool {
 public:
  Partialcharges(){};
  ~Partialcharges(){};

  string Identify() { return "partialcharges"; }

  void Initialize(Property *options);
  bool Evaluate();
  // two access functions for egwbse interface

 private:
  string _orbfile;
  string _output_file;
  Property _esp_options;

  ctp::Logger _log;
};

void Partialcharges::Initialize(Property *options) {

  // update options with the VOTCASHARE defaults
  UpdateWithDefaults(options, "xtp");
  string key = "options." + Identify();

  _orbfile = options->get(key + ".input").as<string>();
  _output_file = options->get(key + ".output").as<string>();
  string _esp2multipole_xml = options->get(key + ".esp_options").as<string>();
  load_property_from_xml(_esp_options, _esp2multipole_xml.c_str());
  // get the path to the shared folders with xml files
  char *votca_share = getenv("VOTCASHARE");
  if (votca_share == NULL)
    throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
}

bool Partialcharges::Evaluate() {

  _log.setReportLevel(ctp::logDEBUG);
  _log.setMultithreading(true);

  _log.setPreface(ctp::logINFO, "\n... ...");
  _log.setPreface(ctp::logERROR, "\n... ...");
  _log.setPreface(ctp::logWARNING, "\n... ...");
  _log.setPreface(ctp::logDEBUG, "\n... ...");

  CTP_LOG(ctp::logDEBUG, _log)
      << "Converting serialized QM data in " << _orbfile << flush;

  Orbitals _orbitals;
  // load the QM data from serialized orbitals object

  CTP_LOG(ctp::logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
  _orbitals.ReadFromCpt(_orbfile);

  Esp2multipole esp2multipole = Esp2multipole(&_log);
  esp2multipole.Initialize(_esp_options);
  esp2multipole.Extractingcharges(_orbitals);

  esp2multipole.WritetoFile(_output_file);

  CTP_LOG(ctp::logDEBUG, _log)
      << "Written charges to " << _output_file << flush;

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif