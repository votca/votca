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

#ifndef VOTCA_XTP_QMPACKAGE_H
#define VOTCA_XTP_QMPACKAGE_H

#include <memory>
// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "aobasis.h"
#include "classicalsegment.h"
#include "logger.h"
#include "staticsite.h"
#include "votca/xtp/orbreorder.h"

namespace votca {
namespace xtp {

class Orbitals;

class QMPackage {
 public:
  virtual ~QMPackage() = default;

  virtual std::string getPackageName() const = 0;

  void Initialize(const tools::Property& options);

  /// writes a coordinate file WITHOUT taking into account PBCs
  virtual bool WriteInputFile(const Orbitals& orbitals) = 0;

  bool Run();
  virtual bool ParseLogFile(Orbitals& orbitals) = 0;

  virtual bool ParseMOsFile(Orbitals& orbitals) = 0;

  virtual void CleanUp() = 0;

  template <class MMRegion>
  void AddRegion(const MMRegion& mmregion) {

    using Segmenttype = typename MMRegion::SegmentType;
    using Sitetype = typename Segmenttype::Atom_Type;
    for (const Segmenttype& segment : mmregion) {
      for (const Sitetype& site : segment) {
        externalsites_.push_back(std::make_unique<Sitetype>(site));
      }
    }
    WriteChargeOption();
  }

  void setRunDir(const std::string& run_dir) { run_dir_ = run_dir; }

  void setInputFileName(const std::string& input_file_name) {
    input_file_name_ = input_file_name;
  }

  void setLogFileName(const std::string& log_file_name) {
    log_file_name_ = log_file_name;
  }

  void setMOsFileName(const std::string& mo_file) { mo_file_name_ = mo_file; }

  void setLog(Logger* pLog) { pLog_ = pLog; }

  void setCharge(Index charge) {
    charge_ = charge;
    spin_ = std::abs(charge) + 1;
  }

  bool GuessRequested() const {
    return options_.get("initial_guess").as<std::string>() == "orbfile";
  }

  virtual StaticSegment GetCharges() const = 0;

  virtual Eigen::Matrix3d GetPolarizability() const = 0;

  std::string getLogFile() const { return log_file_name_; };

  std::string getMOFile() const { return mo_file_name_; };

 protected:
  virtual void ParseSpecificOptions(const tools::Property& options) = 0;
  struct MinimalMMCharge {
    MinimalMMCharge(const Eigen::Vector3d& pos, double q) : pos_(pos), q_(q){};
    Eigen::Vector3d pos_;
    double q_;
  };

  virtual bool RunDFT() = 0;
  virtual void WriteChargeOption() = 0;
  std::vector<MinimalMMCharge> SplitMultipoles(const StaticSite& site) const;
  void ReorderOutput(Orbitals& orbitals) const;
  Eigen::MatrixXd ReorderMOsBack(const Orbitals& orbitals) const;
  bool isLinker(std::string name, std::vector<std::string> linker_names);

  std::vector<std::string> GetLineAndSplit(std::ifstream& input_file,
                                           const std::string separators) const;

  // ShellReorder() and ShellMulitplier() specify the order for each
  // QMPackage. Some codes also use different normalisation conditions which
  // lead to other signs for some of the entries, which can be changed via the
  // multipliers.
  virtual const std::array<Index, 49>& ShellMulitplier() const = 0;
  virtual const std::array<Index, 49>& ShellReorder() const = 0;

  Index charge_;
  Index spin_;  // 2S+1mem
  std::string basisset_name_;
  std::string cleanup_ = "";
  std::string input_file_name_;
  std::string log_file_name_;
  std::string mo_file_name_;
  std::string run_dir_;
  std::string scratch_dir_;
  std::string shell_file_name_;
  tools::Property options_;

  Logger* pLog_;

  std::vector<std::unique_ptr<StaticSite> > externalsites_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMPACKAGE_H
