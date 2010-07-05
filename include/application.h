/*
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_APPLICATION_H
#define	__VOTCA_APPLICATION_H

#include <boost/program_options.hpp>
#include "property.h"

namespace votca { namespace tools {

class Application
{
public:
    Application();
    virtual ~Application();

    /// executes the program
    int Run(int argc, char **argv);

    /// Initialize application data
    /// 
    /// The initialize function is called by run before parsing the command line.
    /// All necesassary command line arguments can be added here
    virtual void Initialize() {};
    
    /// parse program options from command line
    boost::program_options::options_description_easy_init
        AddProgramOptions() { return _op_desc.add_options(); }
    /// get available program options & descriptions
    boost::program_options::variables_map &OptionsMap() { return _op_vm; }
    boost::program_options::options_description &OptionsDesc() { return _op_desc; }

    /// function implementations in child classes
    virtual string HelpText() { return ""; }

    /// initialize variables of child class etc
    ///virtual void Initialize() {}
    /// check whether required input is present and correct
    ///virtual void CheckInput() {}


protected:
    /// Variable map containing all program options
    boost::program_options::variables_map _op_vm;

    /// program options required by all applications
    boost::program_options::options_description _op_desc;
   
private:
    /// get input parameters from file, location may be specified in command line
    void ParseCommandLine(int argc, char **argv);
};

}}

#endif	/* __VOTCA_APPLICATION_H */

