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
/// For earlier commit history see ctp commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_OPTIONS_H
#define VOTCA_XTP_OPTIONS_H

#include <votca/tools/property.h>

namespace votca {
namespace xtp {
class Options {
 public:
  Options(const std::string& defaultspath) : _defaultspath(defaultspath){};

  void ResolveLinks(tools::Property& prop) const;

 private:
  std::string _defaultspath = "";
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_OPTIONS_H
