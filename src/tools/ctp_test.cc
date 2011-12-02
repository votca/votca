/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <stdlib.h>
#include <string>
#include <votca/tools/application.h>
#include <votca/ctp/atom.h>
#include <votca/ctp/molecule.h>

using namespace std;
using namespace votca::tools;
using namespace votca::ctp;

class CTPTest : public Application 
{
    public:
        string ProgramName() { return "ctp_test"; }

        void HelpText(std::ostream &out) {
             out << "Runs tests of CTP objects" << endl;
        }

        void Initialize(void) {
                AddProgramOptions()
                ("options,o", boost::program_options::value<string>(), "  program options");
        }

        bool EvaluateOptions() { return true; }
            
        void Run()
        {
            cout << "Start of the test" << endl;

            cout << " Molecule object" << endl;         
            Molecule molecule(1, "DCV2T");
            cout << "  name: " << molecule.getName() << endl;          
            
            molecule.Init( "coordinates.dat" );            
            cout << "  number of atoms: " << molecule.NumberOfAtoms() << endl;
            molecule.WritePDB( cout );
            
        }
        
        

};


int main(int argc, char** argv) {
    CTPTest ctptest;
    return ctptest.Exec(argc, argv);
}
