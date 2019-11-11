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

#ifndef VOTCA_TOOLS_CROSSCORRELATE_H
#define VOTCA_TOOLS_CROSSCORRELATE_H

#include "datacollection.h"
#include <iostream>
#include <vector>
#include <votca/tools/eigen.h>
#include <votca/tools/types.h>

namespace votca {
namespace tools {

/**
    \brief class to calculate cross correlations and autocorrelations
*/
class CrossCorrelate {
 public:
  CrossCorrelate() = default;
  ~CrossCorrelate() = default;

  void AutoCorrelate(DataCollection<double>::selection& data);

  std::vector<double>& getData() { return _corrfunc; }
  const std::vector<double>& getData() const { return _corrfunc; }

 private:
  std::vector<double> _corrfunc;
};

inline std::ostream& operator<<(std::ostream& out, const CrossCorrelate& c) {
  for (Index i = 0; i < Index(c.getData().size()); i++) {
    out << i << " " << c.getData()[i] << std::endl;
  }
  return out;
}

}  // namespace tools
}  // namespace votca

#endif /* VOTCA_TOOLS_CROSSCORRELATE_H */
