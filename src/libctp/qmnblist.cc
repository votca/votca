/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#include <votca/ctp/qmnblist.h>

namespace votca { namespace ctp {

QMPair *QMNBList::Add(Segment* seg1, Segment* seg2) {

    if (this->FindPair(seg1, seg2) != NULL) {
        throw std::runtime_error("Critical bug: pair already exists");
    }

    int id = this->size();

    QMPair *pair = new QMPair(id, seg1, seg2);

    this->AddPair(pair);

    return pair;
    
}




void QMNBList::PrintInfo(FILE *out) {

    QMNBList::iterator nit;

    for (nit = this->begin();
         nit != this->end();
         nit++) {

        QMPair *pair = *nit;

        int ghost;
        if (pair->HasGhost()) { ghost = 1; }
        else { ghost = 0; }

        fprintf(out, "PairID %5d  | Seg1 %4d Seg2 %4d dR %2.4f Ghost? %1d | lOuter %1.4 J %2.4f r12 %2.4f r21 %2.4f \n",
                pair->getId(),
                pair->first->getId(),
                pair->second->getId(),
                pair->Dist(),
                ghost,
                0.0, // pair->getLambdaO(),
                0.0, // pair->calcJeff2(),
                0.0, // pair->getRate12(),
                0.0 ); // pair->getRate21() );
    }
}

}}
