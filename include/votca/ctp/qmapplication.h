/* 
 * File:   qmapplication.h
 * Author: vehoff
 *
 * Created on January 29, 2010, 10:38 AM
 */

#ifndef _QMAPPLICATION_H
#define	_QMAPPLICATION_H

#include <votca/tools/application.h>
#include "statesaversqlite.h"
#include "qmcalculator.h"
#include "qmtopology.h"

namespace votca { namespace ctp {

class QMApplication : public Application
{
public:
    QMApplication();
    ~QMApplication();

    void Initialize();
    bool EvaluateOptions();

    void Run(void);

    void ShowHelpText(std::ostream &out);

    // print neighbor list to file in human-readable format
    //this function is obsolate, please use the writexml calculator
    // void PrintNbs(string filename);
    /// get string of nearest neighbor names (NNnames)
    vector <string>& get_nnnames() {return _nnnames;}

    /// return true if evaluation should be continued, abort only if something important is missing
    virtual void BeginEvaluate();
    /// called for each frame, return true if evaluation should be continued
    virtual bool EvaluateFrame();
    /// stop evaluation & do final analysis if necessary
    virtual void EndEvaluate();

    /// void add a calculator for later use (compare: cg_engine -> AddObserver)
    void AddCalculator(QMCalculator *calculator);

protected:
    /// QM topology containing all relevant system information
    QMTopology _qmtop;
    /// Property object to parse xml files elegantly
    Property _options;
    /// List of strings that the concatenation of the two molnames must match to be analyzed
    vector <string> _nnnames;
    /// List of CTP observers for easy access of calculators
    list<QMCalculator *> _calculators;

    /// load system information from statesaver
    void ReadData();
    /// write information to statesaver
    void WriteData();
    /// loads the options in from the options file
    void LoadOptions();

};

}}

#endif	/* _QMAPPLICATION_H */

