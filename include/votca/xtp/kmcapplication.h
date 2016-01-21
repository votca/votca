/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_KMC_APPLICATION_H
#define	__VOTCA_KMC_APPLICATION_H

#include <votca/tools/application.h>
#include <votca/xtp/version.h>
#include "kmccalculator.h"
#include <votca/xtp/graphsql.h>

namespace votca { namespace xtp {

class KMCApplication : public Application
{
public:
    /// constructor
    KMCApplication();
    /// destructor
   ~KMCApplication();
    /// adds program options
    void Initialize(void);
    /// print options from the XML file
    bool EvaluateOptions();
    /// print help 
    void ShowHelpText(std::ostream &out);
    /// add a calculator to the list
    void AddCalculator(KMCCalculator* calculator);
    /// run all calculators
    void Run(void);
     /// return true if evaluation should be continued, abort only if something important is missing
    virtual void BeginEvaluate();
    /// called for each frame, return true if evaluation should be continued
    // TO DO - having filename here is a hack. Shall be changed to a pointer to a GRAPH object
    virtual bool EvaluateFrame();
    /// stop evaluation & do final analysis if necessary
    virtual void EndEvaluate();
   
protected:
    /// List of calculators
    list<KMCCalculator *> _calculators;
    /// program options from the xml file
    Property _options;
    /// sql database file
    string _filename;
    string _outputfile;
    
    int _nThreads;

    /// load system information from statesaver
    void ReadData();
    /// write information to statesaver
    void WriteData();
    /// loads the options in from the options file
    void LoadOptions();

private:
    /// application reads-in a Graph object from an sql file and provides it to all calculators
    GraphSQL<NodeSQL,LinkSQL>* graph;
};

}}

#endif	/* __VOTCA_KMC_APPLICATION_H */
