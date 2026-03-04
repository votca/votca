/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Standard includes
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#include <filesystem>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/checkpointreader.h"
#include "votca/xtp/checkpointwriter.h"
#include "votca/xtp/votca_xtp_config.h"

namespace votca {
namespace xtp {

using namespace checkpoint_utils;
namespace fs = std::filesystem;

std::ostream& operator<<(std::ostream& s, CheckpointAccessLevel l) {

  switch (l) {
    case CheckpointAccessLevel::READ:
      s << "read";
      break;
    case CheckpointAccessLevel::MODIFY:
      s << "modify";
      break;
    case CheckpointAccessLevel::CREATE:
      s << "create";
      break;
  }

  return s;
}

bool FileExists(const std::string& fileName) { return fs::exists(fileName); }

CheckpointFile::CheckpointFile(std::string fN)
    : CheckpointFile(fN, CheckpointAccessLevel::MODIFY) {}

CheckpointFile::CheckpointFile(std::string fN, CheckpointAccessLevel access)
    : fileName_(fN), accessLevel_(access) {

  std::lock_guard<std::recursive_mutex> lock(checkpoint_utils::Hdf5Mutex());

  try {
    H5::Exception::dontPrint();
    hid_t fcpl_id = H5Pcreate(H5P_FILE_CREATE);
    H5::FileCreatPropList fcpList(fcpl_id);
    switch (accessLevel_) {
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
      case CheckpointAccessLevel::READ:
        fileHandle_ = H5::H5File(fileName_, H5F_ACC_RDONLY);
        break;
      case CheckpointAccessLevel::CREATE:
        fileHandle_ = H5::H5File(fileName_, H5F_ACC_TRUNC, fcpList);
        break;
      case CheckpointAccessLevel::MODIFY:
        if (!FileExists(fileName_)) {
          fileHandle_ = H5::H5File(fileName_, H5F_ACC_TRUNC, fcpList);
        } else {
          fileHandle_ = H5::H5File(fileName_, H5F_ACC_RDWR, fcpList);
        }
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    }

  } catch (H5::Exception&) {
    std::stringstream message;
    message << "Could not access file " << fileName_;
    message << " with permission to " << accessLevel_ << "." << std::endl;

    throw std::runtime_error(message.str());
  }
}

std::string CheckpointFile::getFileName() { return fileName_; }

H5::H5File CheckpointFile::getHandle() { return fileHandle_; }

CheckpointWriter CheckpointFile::getWriter(const std::string path_) {
  if (accessLevel_ == CheckpointAccessLevel::READ) {
    throw std::runtime_error("Checkpoint file opened as read only.");
  }
  std::lock_guard<std::recursive_mutex> lock(checkpoint_utils::Hdf5Mutex());

  try {
    return CheckpointWriter(fileHandle_.createGroup(path_), path_);
  } catch (H5::Exception&) {
    try {
      return CheckpointWriter(fileHandle_.openGroup(path_), path_);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not create or open " << fileName_ << ":" << path_
              << std::endl;

      throw std::runtime_error(message.str());
    }
  }
}

CheckpointWriter CheckpointFile::getWriter() { return getWriter("/"); }

CheckpointReader CheckpointFile::getReader(const std::string path_) {
    std::lock_guard<std::recursive_mutex> lock(checkpoint_utils::Hdf5Mutex());

  try {
    return CheckpointReader(fileHandle_.openGroup(path_), path_);
  } catch (H5::Exception&) {
    std::stringstream message;
    message << "Could not open " << fileName_ << ":" << path_ << std::endl;

    throw std::runtime_error(message.str());
  }
}

CheckpointReader CheckpointFile::getReader() { return getReader("/"); }

}  // namespace xtp
}  // namespace votca
