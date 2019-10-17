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
#ifndef VOTCA_XTP_QM_PACKAGE_H
#define VOTCA_XTP_QM_PACKAGE_H

#include "votca/xtp/aobasis.h"
#include <votca/tools/property.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/staticsite.h>
namespace votca {
namespace xtp {

class Orbitals;

// ========================================================================== //
// QMPackage base class for wrappers of ORCA, GAUSSIAN, NWCHEM etc //
// ========================================================================== //

class QMPackage {
 public:
  virtual ~QMPackage() = default;
  ;

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
    if (!_write_charges) {
      _write_charges = true;
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

  void setCharge(double charge) {
    _charge = charge;
    _spin = std::abs(charge) + 1;
  }

  bool GuessRequested() const { return _write_guess; }

  virtual StaticSegment GetCharges() const = 0;

  virtual Eigen::Matrix3d GetPolarizability() const = 0;

 protected:
  struct MinimalMMCharge {
    MinimalMMCharge(const Eigen::Vector3d& pos, double q) : _pos(pos), _q(q){};
    Eigen::Vector3d _pos;
    double _q;
  };

  void ParseCommonOptions(tools::Property& options);

  virtual void WriteChargeOption() = 0;
  std::vector<MinimalMMCharge> SplitMultipoles(const StaticSite& site) const;
  void ReorderOutput(Orbitals& orbitals) const;
  Eigen::MatrixXd ReorderMOsBack(const Orbitals& orbitals) const;
  bool isLinker(std::string name, std::vector<std::string> linker_names);

  std::vector<std::string> GetLineAndSplit(std::ifstream& input_file,
                                           const std::string separators) const;

  int _charge;
  int _spin;  // 2S+1
  std::string _memory;
  std::string _options;

  std::string _executable;
  std::string _input_file_name;
  std::string _log_file_name;
  std::string _mo_file_name;

  std::string _run_dir;

  std::string _basisset_name;
  std::string _auxbasisset_name;
  std::string _ecp_name;

  std::string _shell_file_name;
  std::string _scratch_dir;
  bool _is_optimization = false;

  std::string _cleanup = "";

  bool _get_charges = false;

  bool _write_guess = false;
  bool _write_charges = false;
  bool _write_basis_set = false;
  bool _write_auxbasis_set = false;
  bool _write_pseudopotentials = false;

  Logger* _pLog;

  std::vector<std::unique_ptr<StaticSite> > _externalsites;
  double _dpl_spacing = 0.1;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QM_PACKAGE_H
