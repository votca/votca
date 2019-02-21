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

#ifndef _VOTCA_XTP_SPECTRUM_H
#define _VOTCA_XTP_SPECTRUM_H

#include <stdio.h>

#include <votca/xtp/logger.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/qmtool.h>

namespace votca {
namespace xtp {
class Orbitals;

class Spectrum : public QMTool {
 public:
  Spectrum(){};

  ~Spectrum(){};

  std::string Identify() { return "spectrum"; }

  void Initialize(tools::Property& options);
  bool Evaluate();

 private:
  std::string _orbfile;
  std::string _output_file;

  Logger _log;

  void CheckContent(const Orbitals& _orbitals);

  double evtonm(double eV);
  double evtoinvcm(double eV);
  double nmtoinvcm(double nm);
  double invcmtonm(double invcm);
  double nmtoev(double nm);

  double _lower;  // in eV
  double _upper;  // in eV
  int _n_pt;

  int _minexc;     // in eV
  int exc_lambda;  // in eV

  double _fwhm;  // in eV
  double _shiftby;

  std::string _spectrum_type;
  // lineshape functions
  double Gaussian(double x, double center, double fwhm);
  double Lorentzian(double x, double center, double fwhm);
  double TruncatedLorentzian(double x, double center, double fwhm);
};

}  // namespace xtp
}  // namespace votca

#endif
