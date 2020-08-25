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

// Standard inlcudes
#include <numeric>
#include <set>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// VOTCA includes
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/IndexParser.h"

namespace votca {
namespace xtp {

std::vector<Index> IndexParser::CreateIndexVector(
    const std::string& Ids) const {
  std::vector<Index> result;
  tools::Tokenizer tok(Ids, " ,\n\t");
  std::vector<std::string> results;
  tok.ToVector(results);
  const std::string delimiter = ":";
  for (std::string s : results) {
    if (s.find(delimiter) != std::string::npos) {
      std::string copy_s = s;
      try {
        Index start =
            boost::lexical_cast<Index>(s.substr(0, s.find(delimiter)));
        Index stop = boost::lexical_cast<Index>(
            s.erase(0, s.find(delimiter) + delimiter.length()));
        for (Index i = start; i <= stop; i++) {
          result.push_back(i);
        }
      } catch (boost::bad_lexical_cast&) {
        throw std::runtime_error("Could not convert " + copy_s +
                                 " to range of integers.");
      }
    } else {
      try {
        result.push_back(boost::lexical_cast<Index>(s));
      } catch (boost::bad_lexical_cast&) {
        throw std::runtime_error("Could not convert " + s + " to integer.");
      }
    }
  }
  // Eliminates duplicates and sorts the vector
  std::set<Index> s(result.begin(), result.end());
  result.assign(s.begin(), s.end());
  return result;
}

std::string IndexParser::CreateIndexString(
    const std::vector<Index>& indeces) const {
  std::set<Index> s(indeces.begin(), indeces.end());
  std::vector<Index> sorted_unique(s.begin(), s.end());
  std::string result = "";
  if (sorted_unique.empty()) {
    return result;
  }

  std::vector<Index> difference(sorted_unique.size() + 1);
  difference.back() = 0;
  std::adjacent_difference(sorted_unique.begin(), sorted_unique.end(),
                           difference.begin());

  bool range_started = false;
  Index startindex;
  for (Index i = 0; i < Index(sorted_unique.size()); i++) {
    if (difference[i + 1] == 1) {
      if (range_started) {
        continue;
      } else {
        range_started = true;
        startindex = sorted_unique[i];
      }
    } else {
      if (range_started) {
        result += std::to_string(startindex) + ":" +
                  (std::to_string(sorted_unique[i])) + " ";
        range_started = false;
      } else {
        result += (std::to_string(sorted_unique[i])) + " ";
      }
    }
  }
  boost::trim(result);
  return result;
}

}  // namespace xtp
}  // namespace votca
