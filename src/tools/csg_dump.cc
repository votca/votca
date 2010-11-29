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

#include <stdlib.h>
#include <csgapplication.h>

using namespace std;
using namespace votca::csg;

class CsgDumpApp
    : public CsgApplication
{
    string ProgramName() { return "csg_dump"; }
    void HelpText(ostream &out) { out << "Print atoms that are read from topology file to help"
        " debugging atom naming."; }

    bool EvaluateTopology(Topology *top, Topology *top_ref);
};

int main(int argc, char** argv)
{
    CsgDumpApp app;

    return app.Exec(argc, argv);
}

bool CsgDumpApp::EvaluateTopology(Topology *top, Topology *top_ref)
{
    MoleculeContainer::iterator mol;
    for (mol = top->Molecules().begin(); mol != top->Molecules().end(); ++mol) {
        cout << "molecule: " << (*mol)->getId() + 1 << " " << (*mol)->getName()
                << " beads: " << (*mol)->BeadCount() << endl;
        for (int i = 0; i < (*mol)->BeadCount(); ++i) {
            cout << (*mol)->getBeadId(i) << " " <<
                    (*mol)->getBeadName(i) << " " << (*mol)->getBead(i)->getType()->getName() << endl;
        }
    }
    return true;
}

