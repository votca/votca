#ifndef PARALLELSITECALC_H
#define PARALLELSITECALC_H


#include <votca/ctp/qmcalculator2.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace ctp {

class ParallelSiteCalculator : public QMCalculator2
{

public:

    ParallelSiteCalculator() : _nThreads(2) { };
   ~ParallelSiteCalculator() { };

    string Identify() { return "Parallel Site Calculator"; }

    bool EvaluateFrame(Topology *top);


private:

    int _nThreads;


};

bool ParallelSiteCalculator::EvaluateFrame(Topology *top) {


    cout << endl << "... ... Start evaluating...";


    


}
    




}}




#endif
