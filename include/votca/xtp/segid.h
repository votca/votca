/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <string>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/qmstate.h>
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
    _id = std::stoi(results[0]);
    TestStringForQMState(results[1]);
  }

  SegId(Index id, std::string geometry) : _id(id) {
    TestStringForQMState(geometry);
  }

  Index Id() const { return _id; }
  bool hasFile() const { return _hasfilename; }
  std::string FileName() const { return _filename; }
  QMState getQMState() const { return _state; }

 private:
  void TestStringForQMState(const std::string& result) {
    std::string extension = tools::filesystem::GetFileExtension(result);
    if (extension == "pdb" || extension == "xyz" || extension == "mps") {
      _hasfilename = true;
      _filename = result;
    } else {
      try {
        _state = QMState(result);
        _hasfilename = false;
      } catch (std::runtime_error&) {
        throw std::runtime_error("'" + result +
                                 "' is neither a QMState nor a filename. Did "
                                 "you maybe forget the fileending");
      }
    }
  }
  bool _hasfilename = false;
  Index _id;
  std::string _filename = "";
  QMState _state;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGID_H
