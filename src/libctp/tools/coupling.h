/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_CTP_COUPLINGH_H
#define _VOTCA_CTP_COUPLINGH_H

#include <votca/ctp/logger.h>

namespace votca { namespace ctp {
    using namespace std;
    
class Coupling : public QMTool
{
public:

    Coupling() { };
   ~Coupling() { };

    string Identify() { return "coupling"; }

    void   Initialize(Property *options);
    bool   Evaluate();

 

private:
    string      _input_file;
    string      _output_file;
    Logger      _log;

};

void Coupling::Initialize(Property* options) 
{
    
    // options reading
    _input_file  = options->get("options.pdb2map.file").as<string> ();
    _output_file = options->get("options.pdb2map.outfile").as<string> ();
    
    
}



bool Coupling::Evaluate() {
    
    LOG( logINFO, _log ) << "Loading" << _input_file << flush;
    
    cout << _log;
    
    return true;
}


}}


#endif