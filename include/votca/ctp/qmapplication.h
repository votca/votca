#ifndef _QMApplication_H
#define	_QMApplication_H

#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/application.h>
#include <votca/ctp/calculatorfactory.h>

#include <votca/ctp/topology.h>
#include "statesaversqlite.h"
#include "qmcalculator.h"


using namespace std;
namespace votca { namespace ctp {

class QMApplication : public Application
{
public:
    QMApplication();
   ~QMApplication() { };

   void Initialize();
   bool EvaluateOptions();
   void Run(void);

   void ShowHelpText(std::ostream &out);





   virtual void BeginEvaluate(int nThreads);
   virtual bool EvaluateFrame();
   virtual void EndEvaluate();

   void AddCalculator(QMCalculator *calculator);

protected:

    CTP::Topology           _top;
    Property                _options;
    list< QMCalculator* >   _calculators;

    void ReadData() {};
    void WriteData() {};
    void LoadOptions() {};

};

}}









#endif /* _QMApplication_H */














