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

#include <stdlib.h>
#include <votca/xtp/kmcapplication.h>
#include <string>

using namespace votca::xtp;

//using namespace votca::xtp;

class KMCRun : public KMCApplication {
    public:
        string ProgramName() { return "kmc_run"; }

        void HelpText(std::ostream &out) {
             out << "Runs specified calculators" << endl;
        }

};


int main(int argc, char** argv) {
    
    KMCRun kmcrun;
    return kmcrun.Exec(argc, argv);
}
