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

namespace votca {
namespace tools {

/**
    \brief class to calculate cross correlkations and autocorrelations

    This class is relatively outdated and only used in csg_boltzmann!

    \todo implementation
*/
class CrossCorrelate {
 public:
  /// constructor
  CrossCorrelate(){};
  /// destructor
  ~CrossCorrelate(){};

  /**
      calculate the cross correlation
   */
  // void CrossCorrelate(DataCollection<double>::selection *data1,
  //    DataCollection<double>::selection *data2, bool average = false);

  /**
      calculate the auto correlation

   */
  void AutoCorrelate(DataCollection<double>::selection* data,
                     bool average = false);

  // Calculates only the Fourier trafo
  void FFTOnly(std::vector<double>& ivec);

  // Calculates only the Discrete Cosine trafo
  void DCTOnly(std::vector<double>& ivec);

  // Calculates Fourier trafo and then auto correlation
  void AutoFourier(std::vector<double>& ivec);

  // Calculates Discrete Cosine trafo and then auto correlation
  void AutoCosine(std::vector<double>& ivec);

  // Calculates auto correlation via two Fourier trafos
  void AutoCorr(std::vector<double>& ivec);

  std::vector<double>& getData() { return _corrfunc; }

 private:
  std::vector<double> _corrfunc;
};

inline std::ostream& operator<<(std::ostream& out, CrossCorrelate& c) {
  std::vector<double>& data = c.getData();
  for (size_t i = 0; i < data.size(); i++) {
    out << i << " " << c.getData()[i] << std::endl;
  }
  return out;
}

}  // namespace tools
}  // namespace votca

#endif /* VOTCA_TOOLS_CROSSCORRELATE_H */
