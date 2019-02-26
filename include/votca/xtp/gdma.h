/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __VOTCA_XTP_GDMA_H
#define __VOTCA_XTP_GDMA_H

#include <fstream>
#include <map>
#include <string>
#include <votca/ctp/logger.h>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {
/**
    \brief information about an element

    The Atom class stores atom id, name, type, mass, charge, residue number

*/
class GDMA {
 public:
  void WriteInputFile();
  void RunExternal();
  void ParseOutputFile();
  void Initialize(tools::Property& options);

  // functions to override Initialize
  void SetLimit(double rank) { _limit = rank; };
  void SetRunDir(std::string dir) { _runFolder = dir; }
  void SetChkFile(std::string file) { _chkFile = file; };
  void SetExecutable(std::string exec) { _executable = exec; };
  void SetDensity(std::string density) { _density = density; };
  void SetRadius(double radius) { _radius = radius; };
  void SetSwitch(double sw) { _switch = sw; };

  void setLog(ctp::Logger* pLog) { _log = pLog; };
  std::vector<std::vector<double> >& GetMultipoles() { return _multipoles; };

 private:
  std::vector<std::vector<double> > _multipoles;

  std::string _runFolder;
  std::string _chkFile;
  std::string _executable;

  std::string _density;
  int _limit;
  double _radius;
  double _switch;
  std::string _outFile;
  ctp::Logger* _log;
};
}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_GDMA_H */
