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

#include "votca/xtp/aobasis.h"
#include <votca/tools/property.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/settings.h>
#include <votca/xtp/staticsite.h>

namespace votca {
namespace xtp {

class Orbitals;

// ========================================================================== //
// QMPackage base class for wrappers of ORCA
// ========================================================================== //

class QMPackage {
 public:
  virtual ~QMPackage() = default;

  virtual std::string getPackageName() const = 0;

  virtual void Initialize(tools::Property& options) = 0;

  /// writes a coordinate file WITHOUT taking into account PBCs
  virtual bool WriteInputFile(const Orbitals& orbitals) = 0;

  virtual bool Run() = 0;

  virtual bool ParseLogFile(Orbitals& orbitals) = 0;

  virtual bool ParseMOsFile(Orbitals& orbitals) = 0;

  virtual void CleanUp() = 0;

  template <class MMRegion>
  void AddRegion(const MMRegion& mmregion) {

    using Segmenttype =
        typename std::iterator_traits<typename MMRegion::iterator>::value_type;
    using Sitetype = typename std::iterator_traits<
        typename Segmenttype::iterator>::value_type;
    for (const Segmenttype& segment : mmregion) {
      for (const Sitetype& site : segment) {
        _externalsites.push_back(
            std::unique_ptr<StaticSite>(new Sitetype(site)));
      }
    }
    if (_settings.get<bool>("write_charges")) {
      WriteChargeOption();
    }
  }

  void setRunDir(const std::string& run_dir) { _run_dir = run_dir; }

  void setInputFileName(const std::string& input_file_name) {
    _input_file_name = input_file_name;
  }

  void setLogFileName(const std::string& log_file_name) {
    _log_file_name = log_file_name;
  }

  void setMOsFileName(const std::string& mo_file) { _mo_file_name = mo_file; }

  void setLog(Logger* pLog) { _pLog = pLog; }

  void setCharge(Index charge) {
    _charge = charge;
    _spin = std::abs(charge) + 1;
  }

  bool GuessRequested() const { return _settings.get<bool>("read_guess"); }

  virtual StaticSegment GetCharges() const = 0;

  virtual Eigen::Matrix3d GetPolarizability() const = 0;

 protected:
  struct MinimalMMCharge {
    MinimalMMCharge(const Eigen::Vector3d& pos, double q) : _pos(pos), _q(q){};
    Eigen::Vector3d _pos;
    double _q;
  };

  void ParseCommonOptions(tools::Property& options);
  std::string FindDefaultsFile() const;

  virtual void WriteChargeOption() = 0;
  std::vector<MinimalMMCharge> SplitMultipoles(const StaticSite& site) const;
  void ReorderOutput(Orbitals& orbitals) const;
  Eigen::MatrixXd ReorderMOsBack(const Orbitals& orbitals) const;
  bool isLinker(std::string name, std::vector<std::string> linker_names);

  std::vector<std::string> GetLineAndSplit(std::ifstream& input_file,
                                           const std::string separators) const;

  void ReorderMOsToXTP(Eigen::MatrixXd& v, const AOBasis& basis) const;
  void ReorderMOsToNative(Eigen::MatrixXd& v, const AOBasis& basis) const;

  void ReorderMOs(Eigen::MatrixXd& v, const std::vector<Index>& order) const;
  void MultiplyMOs(Eigen::MatrixXd& v,
                   const std::vector<Index>& multiplier) const;

  std::vector<Index> getMultiplierVector(const AOBasis& basis) const;
  std::vector<Index> getMultiplierShell(const AOShell& shell) const;

  std::vector<Index> getReorderVector(const AOBasis& basis) const;
  std::vector<Index> getReorderShell(const AOShell& shell) const;
  std::vector<Index> invertOrder(const std::vector<Index>& order) const;

  virtual const std::array<Index, 25>& ShellMulitplier() const = 0;
  virtual const std::array<Index, 25>& ShellReorder() const = 0;

  Settings _settings{"package"};

  Index _charge;
  Index _spin;  // 2S+1mem
  std::string _basisset_name;
  std::string _cleanup = "";
  std::string _input_file_name;
  std::string _log_file_name;
  std::string _mo_file_name;
  std::string _options = "";
  std::string _run_dir;
  std::string _scratch_dir;
  std::string _shell_file_name;

  Logger* _pLog;

  std::vector<std::unique_ptr<StaticSite> > _externalsites;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMPACKAGE_H
