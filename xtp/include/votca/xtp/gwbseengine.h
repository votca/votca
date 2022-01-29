/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#ifndef VOTCA_XTP_GWBSEENGINE_H
#define VOTCA_XTP_GWBSEENGINE_H

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "logger.h"

namespace votca {
namespace xtp {
class QMPackage;
class Orbitals;

/**
 * \brief Electronic Excitations via Density-Functional Theory
 *
 * Evaluates electronic ground state in molecular systems based on
 * density functional theory with Gaussian Orbitals.
 *
 */

class GWBSEEngine {
 public:
  std::string Identify() { return "gwbse_engine"; }

  void Initialize(tools::Property& options, std::string archive_filename);
  void ExcitationEnergies(Orbitals& orbitals);

  void setLog(Logger* pLog) { pLog_ = pLog; }

  void setQMPackage(QMPackage* qmpackage) { qmpackage_ = qmpackage; }

  std::string GetDFTLog() const { return dftlog_file_; };

  void setLoggerFile(std::string logger_file) { logger_file_ = logger_file; };

  const tools::Property& ReportSummary() const { return summary_; };

 private:
  QMPackage* qmpackage_;

  Logger* pLog_;

  // task options
  bool do_guess_ = false;
  bool do_dft_input_ = false;
  bool do_dft_run_ = false;
  bool do_dft_parse_ = false;
  bool do_gwbse_ = false;
  bool do_localize_ = false;
  bool do_dft_in_dft_ = false;

  // DFT log and MO file names
  std::string MO_file_;      // file containing the MOs from qmpackage...
  std::string dftlog_file_;  // file containing the Energies etc... from
                             // qmpackage...
  std::string logger_file_;
  std::string archive_file_;
  std::string guess_archiveA_;
  std::string guess_archiveB_;

  // Options for GWBSE module
  tools::Property gwbse_options_;
  tools::Property localize_options_;
  tools::Property dft_in_dft_options_;
  tools::Property summary_;

  void WriteLoggerToFile(Logger* pLog);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GWBSEENGINE_H
