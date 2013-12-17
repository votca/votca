/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
#define	_dlpolytrajectoryreader_H

#include <string>
#include <iostream>
#include <fstream>
#include <votca/csg/trajectoryreader.h>

namespace votca { namespace csg {
using namespace votca::tools;

using namespace std;

/**
    \brief class for reading dlpoly trajectory files

    This class provides the TrajectoryReader interface and encapsulates the trajectory reading function of dl_poly

*/
class DLPOLYTrajectoryReader 
   : public TrajectoryReader
{
 public:
  // open original trajectory file
  bool Open(const string &file);
  // read in the first frame
  bool FirstFrame(Topology &top);
  // read in the next frame
  bool NextFrame(Topology &top);
  // close original trajectory file
  void Close();

  void   setFname(string name) { _fname = name; return; }
  string getFname()            { return _fname; }
  
  void setIsConfig(bool isConf) { _isConfig=isConf; return; }
  bool getIsConfig()            { return _isConfig; }
  
 private:
  ifstream _fl;
  string _fname;
  bool _first_frame;
  bool _isConfig;
};

}}

#endif	/* _dlpolytrajectoryreader_H */

