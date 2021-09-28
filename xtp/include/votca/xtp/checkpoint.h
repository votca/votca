/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#pragma once
#ifndef VOTCA_XTP_CHECKPOINT_H
#define VOTCA_XTP_CHECKPOINT_H

// Standard includes
#include <fstream>

// Third party includes
#include <H5Cpp.h>

// Local VOTCA includes
#include "checkpoint_utils.h"
#include "checkpointreader.h"
#include "checkpointtable.h"
#include "checkpointwriter.h"

namespace votca {
namespace xtp {

enum class CheckpointAccessLevel {
  READ = 0,    // only read no write access
  MODIFY = 1,  // if file exists, change it
  CREATE = 2   // create new file
};

std::ostream& operator<<(std::ostream& s, CheckpointAccessLevel l);

class CheckpointFile {
 public:
  CheckpointFile(std::string fN);
  CheckpointFile(std::string fN, CheckpointAccessLevel access);

  std::string getFileName();
  std::string getVersion();

  H5::H5File getHandle();

  CheckpointWriter getWriter();
  CheckpointWriter getWriter(const std::string path_);
  CheckpointReader getReader();
  CheckpointReader getReader(const std::string path_);

 private:
  std::string fileName_;
  H5::H5File fileHandle_;
  CptLoc rootLoc_;
  CheckpointAccessLevel accessLevel_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_H
