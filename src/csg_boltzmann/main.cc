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


// TODO: This code need lots of cleaning up! please do not look at anything in here!
//

#include <math.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <votca/tools/tokenizer.h>
#include <votca/tools/rangeparser.h>
#include "analysistool.h"
#include "bondedstatistics.h"
#include "tabulatedpotential.h"
#include "stdanalysis.h"
#include <csgapplication.h>

using namespace std;

class CsgBoltzmann
    : public CsgApplication
{
public:

    string ProgramName() { return "csg_boltzmann"; }
    void HelpText(ostream &out) {
        out << "Performs tasks that are needed for simple boltzmann\n"
            "inversion in an interactive environment.";
    }
    bool DoTrajectory() { return true; }
    bool DoMapping() { return true; }

    void Initialize();
    void Run();

    void InteractiveMode();
    bool EvaluateTopology(Topology *top, Topology *top_ref);

protected:

    ExclusionList *CreateExclusionList(Molecule &atomistic, Molecule &cg);
    BondedStatistics _bs;

};
void CsgBoltzmann::Initialize()
{
    CsgApplication::Initialize();
    AddProgramOptions("Special options")
        ("excl", boost::program_options::value<string>(), "write exclusion list to file");

    AddObserver(&_bs);
}

bool CsgBoltzmann::EvaluateTopology(Topology *top, Topology *top_ref)
{
    if (OptionsMap().count("excl")) {
        ExclusionList *ex;
        if (top_ref->MoleculeCount() > 1)
            cout << "WARNING: cannot create exclusion list for topology with"
                "multiple molecules, using only first molecule\n";

        cout << "Writing exclusion list for atomistic molecule "
                << top_ref->MoleculeByIndex(0)->getName()
                << " in coarse grained representation "
                << top_ref->MoleculeByIndex(0)->getName() << endl;

        ex = CreateExclusionList(*top_ref->MoleculeByIndex(0), *top->MoleculeByIndex(0));
        ofstream fl;
        fl.open(OptionsMap()["excl"].as<string > ().c_str());
        fl << "# atomistic: " << top_ref->MoleculeByIndex(0)->getName()
                << " cg: " << top_ref->MoleculeByIndex(0)->getName()
                << " cgmap: " << OptionsMap()["cg"].as<string > () << endl;
        fl << *ex;
        fl.close();
        delete ex;
        return false;
    }
    return true;
}

ExclusionList *CsgBoltzmann::CreateExclusionList(Molecule &atomistic, Molecule &cg)
{
    list<int> exclude;
   
    ExclusionList *ex = new ExclusionList();
    ex->ExcludeAll(atomistic.BeadCount());

    // reintroduce bead internal nonbonded interaction
    for(int i=0; i<cg.BeadCount(); ++i) {
        exclude.clear();
        
        vector<int> &v = cg.getBead(i)->ParentBeads();
        exclude.insert(exclude.begin(), v.begin(), v.end());
        ex->Remove(exclude);
    }

    Topology *top_cg = cg.getParent();
    InteractionContainer::iterator iter;
    // reintroduce nonbonded interactions for bonded beads
    for(iter = top_cg->BondedInteractions().begin();
            iter!=top_cg->BondedInteractions().end(); ++iter) {
        Interaction *ic = *iter;
        exclude.clear();
        for(int i=0; i<ic->BeadCount(); i++) {
            vector<int> &v = top_cg->getBead(ic->getBeadId(i))->ParentBeads();
            exclude.insert(exclude.end(), v.begin(), v.end());
        }
        ex->Remove(exclude);
    }
    return ex;
}

void CsgBoltzmann::Run()
{
    CsgApplication::Run();
    if (OptionsMap().count("excl"))
        return;
    InteractiveMode();
}

void CsgBoltzmann::InteractiveMode()
{    
    std::map<std::string, AnalysisTool *> cmds;
    TabulatedPotential tab;
    StdAnalysis std;
    tab.Register(cmds);
    std.Register(cmds);
    
    string help_text = 
        "Interactive mode, expecting commands:\n"
        "help: show this help\n"
        "q: quit\n"
        "list: list all available bonds\n"
    	"vals <file> <selection>: write values to file\n"
    	"hist <file> <selection>: create histogram\n"
    	"tab <file> <selection>: create tabulated potential\n"
    	"autocor <file> <selection>: calculate autocorrelation, only one row allowed in selection!\n" 
    	"cor <file> <selection>: calculate correlations, first row is correlated with all other rows";

    cout << help_text << endl;
    
    while(1) {
        string line;
        cout << "> ";
        getline(cin, line);

        boost::trim(line);
        vector<string> args;
        Tokenizer tok(line, " \t");
        tok.ToVector(args);

        if(args.size() == 0) continue;

        string cmd = args.front();
        args.erase(args.begin());

    
            if(cmd == "q") break;            

            std::map<string, AnalysisTool *>::iterator tool;
            if(cmd == "help") {
                if(args.size() == 0) {
                    cout << help_text << endl;
                    continue;
                }
                cmd = args.front();
                args.erase(args.begin());
                tool = cmds.find(cmd);
                if(tool == cmds.end()) {
                    cout << "error, no help item found" << endl;
                    continue;
                }
                tool->second->Help(cmd, args);
                cout << endl;
                continue;
            }

            tool = cmds.find(cmd);
            if(tool == cmds.end()) {
                cout << "error, command not found" << endl;
                continue;
            }
            
            tool->second->Command(_bs, cmd, args);
    }
}

int main(int argc, char **argv)
{
    CsgBoltzmann app;
    app.Exec(argc, argv);
}
