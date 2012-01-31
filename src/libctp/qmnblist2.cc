#include <votca/ctp/qmnblist2.h>

namespace votca { namespace ctp {

QMPair2 *QMNBList2::Add(Segment* seg1, Segment* seg2) {

    if (this->FindPair(seg1, seg2) != NULL) {
        ; // throw std::runtime_error("Critical bug: pair already exists");
    }

    Topology *top = seg1->getTopology();
    int id = this->size();

    QMPair2 *pair = new QMPair2(id, seg1, seg2);

    this->AddPair(pair);
    
}




void QMNBList2::PrintInfo(ostream &out) {

    QMNBList2::iterator nit;

    for (nit = this->begin();
         nit != this->end();
         nit++) {

        QMPair2 *pair = *nit;

        int ghost;
        if (pair->HasGhost()) { ghost = 1; }
        else { ghost = 0; }

        std::printf("PairID %5d Seg1ID %4d Seg2ID %4d dR %2.4f Ghost? %1d \n",
                pair->getId(),
                pair->first->getId(),
                pair->second->getId(),
                pair->Dist(),
                ghost);

    }

}

}}