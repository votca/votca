#ifndef _PAIRCALCULATOR2_H
#define	_PAIRCALCULATOR2_H

#include "qmcalculator2.h"

namespace votca { namespace ctp {

class PairCalculator2 : public QMCalculator2
{
public:

    PairCalculator2() {};
    virtual ~PairCalculator2() {};

    bool EvaluateFrame(Topology *top);
    virtual void EvaluatePair(Topology *top, QMPair2 *pair) { };
};

bool PairCalculator2::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;

    QMNBList2 &nblist = top->NBList();

    QMNBList2::iterator pit;
    for (pit = nblist.begin(); pit != nblist.end(); pit++) {
        EvaluatePair(top, *pit);
        // return 1; // OVERRIDE
    }
    
    return 1;
}

}}

#endif	/* _PAIRCALCULATOR2_H */

 
