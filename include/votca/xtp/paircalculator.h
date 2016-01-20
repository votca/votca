#ifndef _PAIRCALCULATOR2_H
#define	_PAIRCALCULATOR2_H

#include "qmcalculator.h"

namespace votca { namespace xtp {

class PairCalculator2 : public QMCalculator
{
public:

    PairCalculator2() {};
    virtual ~PairCalculator2() {};

    bool EvaluateFrame(Topology *top);
    virtual void EvaluatePair(Topology *top, QMPair *pair) { };
};

bool PairCalculator2::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;

    QMNBList &nblist = top->NBList();

    QMNBList::iterator pit;
    for (pit = nblist.begin(); pit != nblist.end(); pit++) {

        EvaluatePair(top, *pit);

        if ( (*pit)->getId() == -1 ) {
            
            string pairId = boost::lexical_cast<string>((*pit)->getId());
            string pdbname = "Pair" + pairId + ".pdb";
            (*pit)->WritePDB(pdbname);
        }

    }
    
    return 1;
}

}}

#endif	/* _PAIRCALCULATOR2_H */

 
