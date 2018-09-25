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


#include <iostream>
#include <fstream>
#include <stdexcept>
#include <votca/xtp/votca_config.h>

#include <string>
#include <typeinfo>
#include <vector>
#include <type_traits>
#include <votca/tools/vec.h>
#include <votca/xtp/checkpoint_utils.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpoint.h>
namespace votca {
namespace xtp {

using namespace checkpoint_utils;

CheckpointFile::CheckpointFile(std::string fN): CheckpointFile(fN, true) {};

CheckpointFile::CheckpointFile(std::string fN, bool overWrite)
    : _fileName(fN){

  try {
      bool fileExists = false;

      if (!overWrite){
      // Check if file exists
          std::ifstream file(_fileName);
          fileExists = (bool)file;
      }

      H5::Exception::dontPrint();

      if (fileExists){
          _fileHandle = H5::H5File(_fileName, H5F_ACC_RDWR);
      }
      else {
          _fileHandle = H5::H5File(_fileName, H5F_ACC_TRUNC);
      }
      
  } catch (H5::Exception& error) {
    throw std::runtime_error(error.getDetailMsg());
  }
};

std::string CheckpointFile::getFileName() { return _fileName; };

H5::H5File CheckpointFile::getHandle() { return _fileHandle; };

}  // namespace xtp
}  // namespace votca
