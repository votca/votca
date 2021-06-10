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
#ifndef VOTCA_XTP_SPECTRUM_H
#define VOTCA_XTP_SPECTRUM_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {
class Orbitals;

class Spectrum final : public QMTool {
 public:
  Spectrum() = default;

  ~Spectrum() = default;

  std::string Identify() { return "spectrum"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string orbfile_;
  std::string output_file_ = "spectrum.dat";

  Logger log_;

  void CheckContent(const Orbitals& orbitals_);

  double evtonm(double eV);
  double evtoinvcm(double eV);
  double nmtoinvcm(double nm);
  double invcmtonm(double invcm);
  double nmtoev(double nm);

  double lower_ = 0.0;  // in eV
  double upper_;        // in eV
  Index n_pt_ = 100;

  Index minexc_ = 0;
  Index maxexc_ = std::numeric_limits<Index>::max();

  double fwhm_;  // in eV
  double shiftby_ = 0.0;

  std::string spectrum_type_ = "energy";
  // lineshape functions
  double Gaussian(double x, double center, double fwhm);
  double Lorentzian(double x, double center, double fwhm);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SPECTRUM_H
