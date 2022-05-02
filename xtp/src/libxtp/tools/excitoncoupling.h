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
#ifndef VOTCA_XTP_EXCITONCOUPLING_H
#define VOTCA_XTP_EXCITONCOUPLING_H

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/bsecoupling.h"
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class ExcitonCoupling final : public QMTool {
 public:
  std::string Identify() const { return "excitoncoupling"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string orbA_, orbB_, orbAB_;

  tools::Property coupling_options_;
  std::string output_file_;
  bool classical_;
  std::string mpsA_;
  std::string mpsB_;
  Logger log_;
};

void ExcitonCoupling::ParseOptions(const tools::Property& options) {

  classical_ = options.get(".use_classical").as<bool>();

  if (!classical_) {

    coupling_options_.get(".bsecoupling_options");

    orbA_ = options.get(".orbitalsA").as<std::string>();
    orbB_ = options.get(".orbitalsB").as<std::string>();
    orbAB_ = options.get(".orbitalsAB").as<std::string>();

  } else {
    mpsA_ = options.get(".mpsA").as<std::string>();
    mpsB_ = options.get(".mpsB").as<std::string>();
  }
  output_file_ = job_name_ + "_excitoncoupling.xml";
}

bool ExcitonCoupling::Run() {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");
  tools::Property summary;
  tools::Property& job_output = summary.add("output", "");
  // get the corresponding object from the QMPackageFactory
  if (!classical_) {
    Orbitals orbitalsA, orbitalsB, orbitalsAB;
    // load the QM data from serialized orbitals objects

    XTP_LOG(Log::error, log_)
        << " Loading QM data for molecule A from " << orbA_ << std::flush;
    orbitalsA.ReadFromCpt(orbA_);

    XTP_LOG(Log::error, log_)
        << " Loading QM data for molecule B from " << orbB_ << std::flush;
    orbitalsB.ReadFromCpt(orbB_);

    XTP_LOG(Log::error, log_)
        << " Loading QM data for dimer AB from " << orbAB_ << std::flush;
    orbitalsAB.ReadFromCpt(orbAB_);

    BSECoupling bsecoupling;
    bsecoupling.setLogger(&log_);
    bsecoupling.Initialize(coupling_options_);

    bsecoupling.CalculateCouplings(orbitalsA, orbitalsB, orbitalsAB);
    std::cout << log_;

    tools::Property& pair_summary = job_output.add("pair", "");
    tools::Property& type_summary = pair_summary.add("type", "");
    bsecoupling.Addoutput(type_summary, orbitalsA, orbitalsB);

  }

  else if (classical_) {
    XTP_LOG(Log::error, log_)
        << "Calculating electronic coupling using classical transition charges."
        << orbB_ << std::flush;
    PolarSegment seg1 = PolarSegment("A", 0);
    PolarSegment seg2 = PolarSegment("B", 1);
    seg1.LoadFromFile(mpsA_);
    seg2.LoadFromFile(mpsB_);
    eeInteractor ee;
    double J = ee.CalcStaticEnergy(seg1, seg2);

    tools::Property& pair_summary = job_output.add("pair", "");
    pair_summary.setAttribute("idA", 1);
    pair_summary.setAttribute("idB", 2);
    pair_summary.setAttribute("typeA", mpsA_);
    pair_summary.setAttribute("typeB", mpsB_);
    tools::Property& coupling_summary = pair_summary.add("Coupling", "");
    coupling_summary.setAttribute("jABstatic", J);
  }

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");

  std::ofstream ofs(output_file_, std::ofstream::out);
  ofs << job_output;
  ofs.close();
  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EXCITONCOUPLING_H
