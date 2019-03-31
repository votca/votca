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

#ifndef _VOTCA_CSG_FILEFORMATFACTORY_H
#define _VOTCA_CSG_FILEFORMATFACTORY_H

#include <string>
#include <votca/tools/filesystem.h>
#include <votca/tools/objectfactory.h>

namespace votca {
namespace csg {

namespace TOOLS = votca::tools;

template <typename T>
class FileFormatFactory : public TOOLS::ObjectFactory<std::string, T> {
 public:
  FileFormatFactory() {}

  T *Create(const std::string &file);
};

template <typename T>
T *FileFormatFactory<T>::Create(const std::string &file) {
  std::string filetype = tools::filesystem::GetFileExtension(file);
  try {
    return TOOLS::ObjectFactory<std::string, T>::Create(filetype);
  } catch (std::exception &error) {
    throw std::runtime_error("Error '" + filetype +
                             "' file format of file "
                             "'" +
                             file + "' cannot be read or written");
  }
  return NULL;
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_FILEFORMATFACTORY_H */
