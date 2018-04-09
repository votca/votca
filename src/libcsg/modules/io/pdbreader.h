/*
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_CSG_PDBREADER_H
#define __VOTCA_CSG_PDBREADER_H

#include <string>
#include <iostream>
#include <fstream>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/tools/elements.h>

namespace votca {
namespace csg {
using namespace votca::tools;

using namespace std;

/**
    brief class for reading pdb files

    This class provides the Trajectory and Topology reader interface
    for pdb files

*/
class PDBReader : public TopologyReader,
                  public TrajectoryReader,
                  public Elements {
  public:
    /// Constuctor
    PDBReader() {}
    /// Destructor
    ~PDBReader() {}
    /// open a topology file
    bool ReadTopology(string file, Topology &top);
    /// open a trajectory file
    bool Open(const string &file);
    /// read in the first frame
    bool FirstFrame(Topology &top);
    /// read in the next frame
    bool NextFrame(Topology &top);
    void Close();

  private:
    ifstream _fl;
    bool _topology;
};
}
}

#endif  // __VOTCA_CSG_PDBREADER_H
