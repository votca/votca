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
#ifndef VOTCA_XTP_INDEXPARSER_H
#define VOTCA_XTP_INDEXPARSER_H

#include <string>
#include <vector>
#include <votca/tools/types.h>

/**
 * \brief Parser to read strings containing indexes in the format "1 2:7 8"
 * and returning a sorted expanded std::vector<Index> with only unique entries
 * It can also do the opposite
 */

namespace votca {
namespace xtp {

class IndexParser {
 public:
  std::vector<Index> CreateIndexVector(const std::string& Ids) const;

  std::string CreateIndexString(const std::vector<Index>& indeces) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ATOMCONTAINER_H
