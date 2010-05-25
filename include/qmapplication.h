/* 
 * File:   qmapplication.h
 * Author: vehoff
 *
 * Created on January 29, 2010, 10:38 AM
 */

#ifndef _QMAPPLICATION_H
#define	_QMAPPLICATION_H

#include <boost/program_options.hpp>
#include "qmtopology.h"
#include <votca/tools/property.h>
#include "statesaver.h"
#include "qmcalculator.h"

class QMApplication
{
public:
    QMApplication();
    virtual ~QMApplication();

    /// executes the program
    void Run(int argc, char **argv);
    /// print neighbor list to file in human-readable format
    void PrintNbs(string filename);
    /// get string of nearest neighbor names (NNnames)
    vector <string>& get_nnnames() {return _nnnames;}

    /// parse program options from command line
    boost::program_options::options_description_easy_init
        AddProgramOptions() { return _op_desc_specific.add_options(); }
    /// get available program options & descriptions
    boost::program_options::variables_map &OptionsMap() { return _op_vm; }
    boost::program_options::options_description &OptionsDesc() { return _op_desc; }

    /// function implementations in child classes
    virtual void HelpText();
    /// define and add program specific parameters if necessary
    virtual void AddSpecificOptions() {}
    /// initialize variables of child class etc
    virtual void Initialize() {}
    /// check whether required input is present and correct
    virtual void CheckInput() {}
    /// return true if evaluation should be continued, abort only if something important is missing
    virtual bool BeginEvaluate();
    /// called for each frame, return true if evaluation should be continued
    virtual bool EvaluateFrame(int nr, int nframes);
    /// stop evaluation & do final analysis if necessary
    virtual void EndEvaluate();
    /// void add a calculator for later use (compare: cg_engine -> AddObserver)
    void AddCalculator(QMCalculator *calculator);

protected:
    /// Variable map containing all program options
    boost::program_options::variables_map _op_vm;
    /// Program options required by child classes
    boost::program_options::options_description _op_desc_specific;
    /// program options required by all applications
    boost::program_options::options_description _op_desc;
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

private:
    /// get input parameters from file, location may be specified in command line
    void ParseCommandLine(int argc, char **argv);
};

#endif	/* _QMAPPLICATION_H */

