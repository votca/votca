#ifndef SANDBOX2_H
#define SANDBOX2_H

#include <votca/ctp/parallelpaircalc.h>

namespace votca { namespace ctp {

class Sandbox2 : public ParallelPairCalculator
{

public:
    
    Sandbox2() { };
   ~Sandbox2() { };

    string  Identify() { return "Sandbox"; }
    void    Initialize(Topology *top, Property *options);
    void    EvalPair(Topology *top, QMPair2 *qmpair, int slot);

};


void Sandbox2::Initialize(Topology *top, Property *options) {

    _nThreads = 1;

    cout << endl << "... ... Initialize with " << _nThreads << " threads ";    

}


void Sandbox2::EvalPair(Topology *top, QMPair2 *qmpair, int slot) {
  
    this->LockCout();
    cout << "\r... ... Overloading pair " << qmpair->getId() << ". " << flush;
    this->UnlockCout();

    int ij;
    for (int i = 0; i < 2000; i++) {
        for (int j = 0; j < 2000; j++) {
            ij = i+j;
        }
    }

}













}} /* exit namespace votca::ctp */

#endif /* SANDBOX2_H */
