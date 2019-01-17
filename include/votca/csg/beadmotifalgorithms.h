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

#ifndef VOTCA_CSG_BEADMOTIFALGORITHMS_H
#define VOTCA_CSG_BEADMOTIFALGORITHMS_H

#include "beadstructure.h"

namespace votca {
namespace csg {

template <class T>
T breakIntoMotifs(BeadStructure& beadstructure) {
  T bead_motifs;
  std::vector<BeadStructure> structures = breakIntoStructures(beadstructure);
  for (BeadStructure& structure : structures) {
    BeadMotif bead_motif;
    bead_motif.BeadStructure::operator=(structure);
    bead_motifs.push_back(bead_motif);
  }
  return bead_motifs;
}
/**
 * \brief This function will take a beadmotif and break it into its elemental
 *motifs
 **/
// std::unordered_map<int,BeadMotif> breakIntoSimpleMotifs(BeadMotif beadmotif);

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMOTIFALGORITHMS_H
