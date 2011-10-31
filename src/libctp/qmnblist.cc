/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/ctp/qmnblist.h>
#include <votca/ctp/qmbead.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include <votca/ctp/qmtopology.h>

namespace votca { namespace ctp {

void QMNBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    //Cleanup();
    
    _father = dynamic_cast<QMTopology*> (list1.getTopology());
    
    NBListGrid nb;
    nb.setCutoff(_cutoff);
    nb.SetMatchFunction(this, &QMNBList::Match);
    nb.Generate(list1, list2, do_exclusions);

}

bool QMNBList::Match(Bead *ab1, Bead *ab2, const vec &r, const double notused)
{
    QMBead *b1=dynamic_cast<QMBead*>(ab1);
    QMBead *b2=dynamic_cast<QMBead*>(ab2);

    if(b1->getMolecule() == b2->getMolecule()) return false;
    if(b1->GetCrgUnit() == NULL || b2->GetCrgUnit() == NULL) return false;

    QMCrgUnit *crg1, *crg2;
    crg1 = b1->GetCrgUnit();
    crg2 = b2->GetCrgUnit();

    if(crg1->getId() > crg2->getId())
        swap(crg1, crg2);
        
    if(!FindPair(crg1, crg2))
        AddPair(new QMPair(crg1, crg2, _father));

    return false;
}

}}