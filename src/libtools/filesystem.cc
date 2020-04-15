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

#include "../../include/votca/tools/filesystem.h"
#include <boost/algorithm/string.hpp>
#include <fstream>

namespace votca {
namespace tools {

namespace filesystem {

std::string GetFileExtension(const std::string& filename) {
  size_t i = filename.rfind('.', filename.length());
  if (i != std::string::npos) {
    return (filename.substr(i + 1, filename.length() - i));
  }
  return ("");
}

std::string GetFileBase(const std::string& filename) {
  size_t i = filename.rfind('.', filename.length());
  if (i != std::string::npos) {
    return (filename.substr(0, i));
  }
  return filename;
}

bool FileExists(const std::string& filename) {
  std::ifstream infile(filename);
  return infile.good();
}

}  // namespace filesystem
}  // namespace tools
}  // namespace votca
