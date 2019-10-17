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

#ifndef _VOTCA_CSG_GROREADER_H
#define _VOTCA_CSG_GROREADER_H

#include <fstream>
#include <iostream>
#include <string>
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>

namespace votca {
namespace csg {

/**
    \brief reader for gro files

    This class provides the TrajectoryReader + Topology reader interface
    for gro files

*/
class GROReader : public TrajectoryReader, public TopologyReader {
 public:
  GROReader() {}
  ~GROReader() override {}

  /// open a topology file
  bool ReadTopology(std::string file, Topology &top) override;

  /// open a trejectory file
  bool Open(const std::string &file) override;
  /// read in the first frame
  bool FirstFrame(Topology &top) override;
  /// read in the next frame
  bool NextFrame(Topology &top) override;

  void Close() override;

 private:
  std::ifstream _fl;
  bool _topology;
};

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_GROREADER_H */
