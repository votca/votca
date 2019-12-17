/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _dlpolytrajectoryreader_H
#define _dlpolytrajectoryreader_H

#include <fstream>
#include <iostream>
#include <string>
#include <votca/csg/trajectoryreader.h>
#include <votca/tools/unitconversion.h>

namespace votca {
namespace csg {

/**
    \brief class for reading dlpoly trajectory and configuration files

    This class encapsulates the dlpoly trajectory and configuration reading
   function and provides an interface to fill a topology class

*/

class DLPOLYTrajectoryReader : public TrajectoryReader {
 public:
  /// open original trajectory file
  bool Open(const std::string &file) override;
  /// read in the first frame
  bool FirstFrame(Topology &conf) override;
  /// read in the next frame
  bool NextFrame(Topology &conf) override;
  /// close original trajectory file
  void Close() override;

  /// set/get the original configuration or trajectory file name:
  /// <name>.dlpc/<name>.dlph (convention: ".dlpc"="CONFIG", ".dlph"="HISTORY")
  void setFname(std::string name) {
    _fname = name;
    return;
  }
  std::string getFname() { return _fname; }

  /// set/check the flag for the read-in file as configuration, i.e. not
  /// trajectory format
  void setIsConfig(bool isConf) {
    _isConfig = isConf;
    return;
  }
  bool getIsConfig() { return _isConfig; }

 private:
  std::ifstream _fl;
  std::string _fname;
  bool _first_frame;
  bool _isConfig;
};

}  // namespace csg
}  // namespace votca

#endif /* _dlpolytrajectoryreader_H */
