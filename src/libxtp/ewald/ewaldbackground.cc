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

#include "ewaldbackground.h"
#include "kspace.h"
#include "rspace.h"
#include "votca/xtp/tictoc.h"
#include <fstream>
#include <iostream>
#include <vector>

namespace votca {
namespace xtp {


void EwaldBackground::ImportBackgroundFromHdf5(std::string filename){
  CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("polar_background");
  for (Index i =0 ; i < r.getNumDataSets(); ++i){
    CheckpointReader rr = r.openChild("background_" + std::to_string(i));
    background.push_back(EwdSegment(rr));
  }

  // for(const auto& site : background[0]){
  //   std::cout << site << std::endl;
  // }
}
}
}  // namespace votca