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

#include "beadmotif.h"
#include "beadmotifconnector.h"
#include "beadstructure.h"
#include "beadstructurealgorithms.h"

namespace votca {
namespace csg {

/**
 * \brief breaks a beadstructure into individual motifs.
 *
 * Will essentially take structures that are indenpendent of each other and
 * break them into separate beadmotifs.
 *
 * @param[in] - reference to beadstructure
 * @return - a container of beadmotifs
 **/
template <class T>
T breakIntoMotifs(BeadStructure<BaseBead> &beadstructure) {
  T bead_motifs;
  std::vector<BeadStructure<BaseBead>> structures =
      breakIntoStructures(beadstructure);
  for (BeadStructure<BaseBead> &structure : structures) {
    bead_motifs.push_back(BeadMotif(structure));
  }
  return bead_motifs;
}

/**
 * \brief This function will take a beadmotif and break it into its simple
 *motifs.
 *
 * A simple motif is one of four types:
 * single_bead
 * line
 * loop
 * fused_ring
 *
 * So given a beadmotif like this:
 *
 *       H1               H1
 *       |
 *       C1          =>   C1
 *      /  \
 *     H2  H3           H2   H3
 *
 * This structure which is originally of type 'single_structure' will be broken
 * up into 4 singles of type 'single_bead'
 *
 * Something like this:
 *
 *    C1 - C2                     C1 - C2
 *    |    |                =>     |    |
 *    C3 - C4 - C5 - H1           C3 - C4     C5 - H1
 *
 * This structure is of type 'single_structure' will be broken up into two
 * separate structures 1 of type 'loop' and the other of type line.
 **/

std::pair<std::unordered_map<int, BeadMotif>, BeadMotifConnector>
    breakIntoSimpleMotifs(BeadMotif beadmotif);

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMOTIFALGORITHMS_H
