/* 
 * File:   ctp_run.cc
 * Author: schrader
 *
 * Created on July 12, 2010, 3:11 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include <calculatorfactory.h>
/*
 *
 */
class QMAppRun : public QMApplication
{
public:
    void HelpText() {}

    void AddSpecificOptions(){
        _op_desc_specific.add_options()
        ("exec", boost::program_options::value<string>(), "execution list");
    }

    //TODO: Support for XML-File based options
    void CheckInput() {
        if (_op_vm.count("exec")) {
            Tokenizer tok(_op_vm["exec"].as<string>(), " ,\n\t");
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n)
                AddCalculator(Calculators().Create((*n).c_str()));
        }
    }
};

int main(int argc, char** argv) {
    QMAppRun qmapprun;
    qmapprun.Run(argc, argv);
    return (EXIT_SUCCESS);
}