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

#ifndef VOTCA_CSG_TRAJECTORYWRITER_H
#define VOTCA_CSG_TRAJECTORYWRITER_H

// Standard includes
#include <fstream>
#include <iostream>
#include <string>

// Local VOTCA includes
#include "fileformatfactory.h"
#include "topology.h"

namespace votca {
namespace csg {

class TrajectoryWriter {
 public:
  TrajectoryWriter() = default;
  virtual ~TrajectoryWriter() = default;

  virtual void Open(std::string, bool = false) {}
  virtual void Close() {};

  virtual void Write(Topology *) {}

  static void RegisterPlugins(void);
};

// important - singleton pattern, make sure factory is created before accessed
inline FileFormatFactory<TrajectoryWriter> &TrjWriterFactory() {
  static FileFormatFactory<TrajectoryWriter> TrjWriterFactory_;
  return TrjWriterFactory_;
}

}  // namespace csg
}  // namespace votca

#endif /*  VOTCA_CSG_TRAJECTORYWRITER_H */
