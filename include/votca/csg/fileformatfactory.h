/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_FILEFORMATFACTORY_H
#define VOTCA_CSG_FILEFORMATFACTORY_H

// Standard includes
#include <string>

// VOTCA includes
#include <votca/tools/filesystem.h>
#include <votca/tools/objectfactory.h>

namespace votca {
namespace csg {

template <typename T>
class FileFormatFactory : public tools::ObjectFactory<std::string, T> {
 public:
  FileFormatFactory() = default;

  std::unique_ptr<T> Create(const std::string &file);
};

template <typename T>
std::unique_ptr<T> FileFormatFactory<T>::Create(const std::string &file) {
  std::string filetype = tools::filesystem::GetFileExtension(file);
  try {
    return tools::ObjectFactory<std::string, T>::Create(filetype);
  } catch (std::exception &) {
    throw std::runtime_error("Error '" + filetype +
                             "' file format of file "
                             "'" +
                             file + "' cannot be read or written");
  }
  return nullptr;
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_FILEFORMATFACTORY_H
