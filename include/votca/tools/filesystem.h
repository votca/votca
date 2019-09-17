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

#ifndef VOTCA_TOOLS_FILESYSTEM_H
#define VOTCA_TOOLS_FILESYSTEM_H

#include <string>

namespace votca {
namespace tools {

namespace filesystem {

// return the file ending like .jpg .gro etc.., returns an empty string
// otherwise
std::string GetFileExtension(const std::string& filename);

// return the filename without the file extension
std::string GetFileBase(const std::string& filename);

// returns true if file exists otherwise false
bool FileExists(const std::string& filename);

}  // namespace filesystem
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_FILESYSTEM_H
