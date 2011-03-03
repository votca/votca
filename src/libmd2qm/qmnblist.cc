#include "qmnblist.h"
#include "qmbead.h"
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>
#include "qmtopology.h"

void QMNBList::Generate(BeadList &list1, BeadList &list2, bool do_exclusions)
{
    Cleanup();
    
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

