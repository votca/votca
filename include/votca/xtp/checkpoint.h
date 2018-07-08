/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_XTP_CHECKPOINT_H
#define _VOTCA_XTP_CHECKPOINT_H


#include <H5Cpp.h>
#include <votca/xtp/checkpoint_utils.h>


namespace votca {
namespace xtp {

class CheckpointFile {
 public:
  CheckpointFile(std::string fileName, bool overWrite);
  CheckpointFile(std::string fileName);

  std::string getFileName();
  std::string getVersion();

  H5::H5File getHandle();

 private:
  std::string _fileName;
  H5::H5File _fileHandle;

};


}  // namespace xtp
}  // namespace votca
#endif  // _VOTCA_XTP_CHECKPOINT_H
