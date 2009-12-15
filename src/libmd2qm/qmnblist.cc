#include "qmnblist.h"
#include "qmbead.h"
#include <votca/csg/nblist.h>

void QMNBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    Cleanup();
    NBList nb;
    nb.setCutoff(_cutoff);
    nb.Generate(list1, list2);
    nb.Generate(list1, list2, do_exclusions);

    for(NBList::iterator iter=nb.begin(); iter!=nb.end();++iter) {
        QMBead *b1=dynamic_cast<QMBead*>((*iter)->first);
        QMBead *b2=dynamic_cast<QMBead*>((*iter)->second);

        if(b1->getMolecule() == b2->getMolecule()) continue;

        if(!FindPair(b1->GetCrgUnit(), b2->GetCrgUnit()))
            AddPair(new QMPair(b1->GetCrgUnit(), b2->GetCrgUnit(), (*iter)->r()));
    }
}

