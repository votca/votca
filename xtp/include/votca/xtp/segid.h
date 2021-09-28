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

#pragma once
#ifndef VOTCA_XTP_SEGID_H
#define VOTCA_XTP_SEGID_H

// Standard includes
#include <string>

// VOTCA includes
#include <votca/tools/tokenizer.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "qmstate.h"

/**
 * \brief Small wrapper for a segment id and the corresponding QMState or
 * filename
 *
 *
 *
 */

namespace votca {
namespace xtp {

class SegId {
 public:
  SegId(std::string input) {
    tools::Tokenizer tok(input, ":");
    std::vector<std::string> results = tok.ToVector();
    if (results.size() != 2) {
      throw std::runtime_error("Malformed string '" + input + "' for segment");
    }
    id_ = std::stoi(results[0]);
    TestStringForQMState(results[1]);
  }

  SegId(Index id, std::string geometry) : id_(id) {
    TestStringForQMState(geometry);
  }

  Index Id() const { return id_; }
  bool hasFile() const { return hasfilename_; }
  std::string FileName() const { return filename_; }
  QMState getQMState() const { return state_; }

  bool operator==(const SegId& other) {
    if (this->id_ == other.Id()) {
      return true;
    }
    return false;
  }


 private:
  void TestStringForQMState(const std::string& result) {
    std::string extension = tools::filesystem::GetFileExtension(result);
    if (extension == "pdb" || extension == "xyz" || extension == "mps") {
      hasfilename_ = true;
      filename_ = result;
    } else {
      try {
        state_ = QMState(result);
        hasfilename_ = false;
      } catch (std::runtime_error&) {
        throw std::runtime_error("'" + result +
                                 "' is neither a QMState nor a filename. Did "
                                 "you maybe forget the fileending");
      }
    }
  }
  bool hasfilename_ = false;
  Index id_;
  std::string filename_ = "";
  QMState state_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGID_H
