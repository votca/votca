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

#ifndef _VOTCA_XTP_DFTGWBSE_H
#define _VOTCA_XTP_DFTGWBSE_H

#include <stdio.h>

#include <votca/tools/constants.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/geometry_optimization.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/qmtool.h>
#include <votca/xtp/segment.h>

namespace votca {
namespace xtp {

class DftGwBse : public QMTool {
 public:
  DftGwBse(){};

  ~DftGwBse(){};

  std::string Identify() { return "dftgwbse"; }

  void Initialize(tools::Property *options);
  bool Evaluate();

 private:
  std::string _guess_file;
  bool _do_guess;

  std::string _mpsfile;
  bool _do_external;
  double _dipole_spacing;

  std::string _xyzfile;
  std::string _xml_output;  // .xml output
  std::string _package;
  std::string _archive_file;  // .orb file to parse to
  std::string _reporting;
  std::string _guess_orbA;
  std::string _guess_orbB;

  tools::Property _package_options;
  tools::Property _gwbseengine_options;
  tools::Property _geoopt_options;

  Logger _log;

  bool _do_optimize;
};

}  // namespace xtp
}  // namespace votca

#endif
