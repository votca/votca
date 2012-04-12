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
