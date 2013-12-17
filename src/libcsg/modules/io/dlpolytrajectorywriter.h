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

#ifndef _DLPOLYTRAJECTORYWRITER_H
#define _DLPOLYTRAJECTORYWRITER_H

#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>

namespace votca { namespace csg {
using namespace votca::tools;

class DLPOLYTrajectoryWriter
   : public TrajectoryWriter
{
 public:
  // open transformed trajectory file
  void Open(string file, bool bAppend=false);
  // close transformed trajectory file
  void Close();
  // write a frame into transformed trajectory file
  void Write(Topology *conf);
  
  void   setFname(string name) { _fname = name; return; }
  string getFname()            { return _fname; }
  
  void setIsConfig(bool isConf) { _isConfig=isConf;  return; }
  bool getIsConfig()            { return _isConfig; }
  
 private:
  ofstream _fl;
  string _fname;
  bool _isConfig;
};

}}

#endif  /* _DLPOLYTRAJECTORYWRITER_H */
