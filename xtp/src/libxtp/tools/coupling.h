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

  std::string Identify() const { return "coupling"; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Run();

 private:
  std::string orbA_, orbB_, orbAB_;

  tools::Property package_options_;
  tools::Property dftcoupling_options_;

  std::string output_file_;

  Logger log_;
};

void Coupling::ParseOptions(const tools::Property &options) {

  orbA_ = options.get(".moleculeA").as<std::string>();
  orbB_ = options.get(".moleculeB").as<std::string>();
  orbAB_ = options.get(".dimerAB").as<std::string>();

  output_file_ = options.ifExistsReturnElseReturnDefault<std::string>(
      "output", job_name_ + " coupling_.xml");

  package_options_ = options.get(".dftpackage");
  dftcoupling_options_ = options.get(".dftcoupling_options");

  QMPackageFactory::RegisterAll();
}

bool Coupling::Run() {

  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  Orbitals orbitalsA, orbitalsB, orbitalsAB;
  orbitalsA.ReadFromCpt(orbA_);
  orbitalsB.ReadFromCpt(orbB_);
  orbitalsAB.ReadFromCpt(orbAB_);  

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
