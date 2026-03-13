/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_EXTENDEDHUCKEL_H
#define VOTCA_XTP_EXTENDEDHUCKEL_H

// Local VOTCA includes

/**
 * \brief Parameters for Extended Hueckel calculation
 */

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>

namespace votca {
namespace xtp {

class ExtendedHuckelParameters {
 public:
  ExtendedHuckelParameters();

  // Strict lookup, throws if exact parameter is missing
  double Get(const std::string& element, int l) const;

  // Lookup with fallback for missing higher-l shells
  double GetWithFallback(const std::string& element, int l,
                         int* used_l = nullptr) const;

  bool Has(const std::string& element, int l) const;
  bool HasElement(const std::string& element) const;

 private:
  std::map<std::pair<std::string, int>, double> eps_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_EXTENDEDHUECKEL_H
