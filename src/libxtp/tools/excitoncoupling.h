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
#ifndef _VOTCA_XTP_EXCITONCOUPLINGH_H
#define _VOTCA_XTP_EXCITONCOUPLINGH_H

#include <votca/xtp/logger.h>
#include <votca/xtp/qmtool.h>

#include <stdio.h>
#include <votca/tools/constants.h>
#include <votca/xtp/bsecoupling.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/eeinteractor.h>

#include <votca/xtp/qmpackagefactory.h>

namespace votca {
namespace xtp {

class ExcitonCoupling : public QMTool {
 public:
  std::string Identify() { return "excitoncoupling"; }

  void Initialize(tools::Property& options);
  bool Evaluate();

 private:
  std::string _orbA, _orbB, _orbAB;

  tools::Property _coupling_options;
  std::string _output_file;
  bool _classical;
  std::string _mpsA;
  std::string _mpsB;
  Logger _log;
};

void ExcitonCoupling::Initialize(tools::Property& options) {

  std::string key = "options." + Identify();
  _classical = false;
  if (options.exists(key + ".classical")) {
    _classical = options.get(key + ".classical").as<bool>();

  } else {
    _classical = false;
  }

  if (!_classical) {

    std::string coupling_xml =
        options.get(key + ".bsecoupling_options").as<std::string>();
    load_property_from_xml(_coupling_options, coupling_xml);

    _orbA = options.get(key + ".orbitalsA").as<std::string>();
    _orbB = options.get(key + ".orbitalsB").as<std::string>();
    _orbAB = options.get(key + ".orbitalsAB").as<std::string>();

  } else {
    _mpsA = options.get(key + ".mpsA").as<std::string>();
    _mpsB = options.get(key + ".mpsB").as<std::string>();
  }
  _output_file = options.get(key + ".output").as<std::string>();

  // get the path to the shared folders with xml files
  char* votca_share = getenv("VOTCASHARE");
  if (votca_share == NULL)
    throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
}

bool ExcitonCoupling::Evaluate() {
  OPENMP::setMaxThreads(_nThreads);
  _log.setReportLevel(logDEBUG);
  _log.setMultithreading(true);

  _log.setPreface(logINFO, "\n... ...");
  _log.setPreface(logERROR, "\n... ...");
  _log.setPreface(logWARNING, "\n... ...");
  _log.setPreface(logDEBUG, "\n... ...");
  tools::Property summary;
  tools::Property& job_output = summary.add("output", "");
  // get the corresponding object from the QMPackageFactory
  if (!_classical) {
    Orbitals orbitalsA, orbitalsB, orbitalsAB;
    // load the QM data from serialized orbitals objects

    XTP_LOG_SAVE(logDEBUG, _log)
        << " Loading QM data for molecule A from " << _orbA << std::flush;
    orbitalsA.ReadFromCpt(_orbA);

    XTP_LOG_SAVE(logDEBUG, _log)
        << " Loading QM data for molecule B from " << _orbB << std::flush;
    orbitalsB.ReadFromCpt(_orbB);

    XTP_LOG_SAVE(logDEBUG, _log)
        << " Loading QM data for dimer AB from " << _orbAB << std::flush;
    orbitalsAB.ReadFromCpt(_orbAB);

    BSECoupling bsecoupling;
    bsecoupling.setLogger(&_log);
    bsecoupling.Initialize(_coupling_options);

    bsecoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
    std::cout << _log;

    tools::Property& pair_summary = job_output.add("pair", "");
    tools::Property& type_summary = pair_summary.add("type", "");
    bsecoupling.Addoutput(type_summary, orbitalsA, orbitalsB);

  }

  else if (_classical) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << "Calculating electronic coupling using classical transition charges."
        << _orbB << std::flush;
    PolarSegment seg1 = PolarSegment("A", 0);
    PolarSegment seg2 = PolarSegment("B", 1);
    seg1.LoadFromFile(_mpsA);
    seg2.LoadFromFile(_mpsB);
    eeInteractor ee;
    double J = ee.CalcStaticEnergy(seg1, seg2);

    tools::Property& pair_summary = job_output.add("pair", "");
    pair_summary.setAttribute("idA", 1);
    pair_summary.setAttribute("idB", 2);
    pair_summary.setAttribute("typeA", _mpsA);
    pair_summary.setAttribute("typeB", _mpsB);
    tools::Property& coupling_summary = pair_summary.add("Coupling", "");
    coupling_summary.setAttribute("jABstatic", J);
  }

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");

  std::ofstream ofs(_output_file, std::ofstream::out);
  ofs << job_output;
  ofs.close();
  return true;
}

}  // namespace xtp
}  // namespace votca

#endif
