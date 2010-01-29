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

class QMApplication
{
public:
    QMApplication();
    ~QMApplication();

    /// executes the program
    void Run(int argc, char **argv);

    /// parse program options from command line
    boost::program_options::options_description_easy_init
        AddProgramOptions() { return _op_desc_specific.add_options(); }
    /// get available program options & descriptions
    boost::program_options::variables_map &OptionsMap() { return _op_vm; }
    boost::program_options::options_description &OptionsDesc() { return _op_desc; }

    /// function implementations in child classes
    virtual void HelpText();
    /// define and add program specific parameters if necessary
    virtual void Initialize() {}
    /// check whether required input is present and correct
    virtual void CheckInput() {}
    /// return true if evaluation should be started or continued
    virtual bool BeginEvaluate() {return true;}
    /// called for each frame, return true if evaluation should be continued
    virtual bool EvaluateFrame() { return true; }
    /// stop evaluation & do final analysis if necessary
    virtual void EndEvaluate() {}

protected:
    /// QM topology containing all relevant system information
    QMTopology *_qmtop;    
    Property _options;

    /// load system information from statesaver
    void ReadData();

    /// write information to statesaver
    void WriteData();

private:
    /// get input parameters from file, location may be specified in command line
    void ParseCommandLine(int argc, char **argv);

    boost::program_options::options_description _op_desc;
    boost::program_options::options_description _op_desc_specific;
    boost::program_options::variables_map _op_vm;
};

#endif	/* _QMAPPLICATION_H */

