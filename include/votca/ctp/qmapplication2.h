#ifndef _QMAPPLICATION2_H
#define	_QMAPPLICATION2_H

#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/application.h>
#include <votca/ctp/calculatorfactory2.h>

#include <votca/ctp/topology.h>
#include "statesaversqlite2.h"
#include "qmcalculator2.h"


using namespace std;
namespace votca { namespace ctp {

class QMApplication2 : public Application
{
public:
    QMApplication2();
   ~QMApplication2() { };

   void Initialize();
   bool EvaluateOptions();
   void Run(void);

   void ShowHelpText(std::ostream &out);





   virtual void BeginEvaluate(int nThreads);
   virtual bool EvaluateFrame();
   virtual void EndEvaluate();

   void AddCalculator(QMCalculator2 *calculator);

protected:

    CTP::Topology           _top;
    Property                _options;
    list< QMCalculator2* >  _calculators;

    void ReadData() {};
    void WriteData() {};
    void LoadOptions() {};

};

}}









#endif /* _QMAPPLICATION2_H */














