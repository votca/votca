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

#include <votca/csg/beadmotif.h>
#include <votca/csg/beadstructurealgorithms.h>

using namespace std;

namespace votca {
namespace csg {

vector<BeadMotif> breakInToMotifs(BeadStructure beadstructure) {
  vector<BeadMotif> bead_motifs;

  vector<BeadStructure> structures = breakInToStructures(beadstruccture);
  for (const BeadStructure& structure : structures) {
    BeadMotif bead_motif;
    bead_motif.BeadStructure::operator= structure;
    bead_motifs.push_back(bead_motif);
  }
  return bead_motifs;
}

/**
 * \brief This function will take a beadmotif and break it into its elemental
 *motifs
 *
 * It will return the elemental motif and the edges connecting each motif with
 *other motifs as well as the edges describing which vertices are connecting the
 *motifs
 **/
pair<map<int, BeadMotif>, vector<pair<Edge, Edge>>> breakInToElementalMotifs(
    BeadMotif beadmotif) {}

}  // namespace csg
}  // namespace votca
