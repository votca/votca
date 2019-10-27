/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

// TODO: This code need lots of cleaning up! please do not look at anything in
// here!
//

#include "analysistool.h"
#include "bondedstatistics.h"
#include "stdanalysis.h"
#include "tabulatedpotential.h"
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <votca/csg/csgapplication.h>
#include <votca/tools/rangeparser.h>
#include <votca/tools/tokenizer.h>

using namespace std;
using namespace votca::csg;

class CsgBoltzmann : public CsgApplication {
 public:
  string ProgramName() override { return "csg_boltzmann"; }
  void HelpText(ostream &out) override {
    out << "Performs tasks that are needed for simple boltzmann\n"
           "inversion in an interactive environment.";
  }
  bool DoTrajectory() override { return true; }
  bool DoMapping() override { return true; }

  void Initialize() override;
  bool EvaluateOptions() override;
  void Run() override;

  void InteractiveMode();
  bool EvaluateTopology(Topology *top, Topology *top_ref) override;

 protected:
  ExclusionList *CreateExclusionList(Molecule &atomistic, Molecule &cg);
  BondedStatistics _bs;
};
void CsgBoltzmann::Initialize() {
  CsgApplication::Initialize();
  AddProgramOptions("Special options")(
      "excl", boost::program_options::value<string>(),
      "write atomistic exclusion list to file");

  AddObserver(&_bs);
}

bool CsgBoltzmann::EvaluateOptions() {
  CsgApplication::EvaluateOptions();
  if (OptionsMap().count("excl")) {
    CheckRequired("cg", "excl options needs a mapping file");
  }
  return true;
}

bool CsgBoltzmann::EvaluateTopology(Topology *top, Topology *top_ref) {
  if (OptionsMap().count("excl")) {
    ExclusionList *ex;
    if (top_ref->MoleculeCount() > 1) {
      cout << "WARNING: cannot create exclusion list for topology with"
              "multiple molecules, using only first molecule\n";
    }

    cout << "Writing exclusion list for atomistic molecule "
         << top_ref->MoleculeByIndex(0)->getName()
         << " in coarse grained representation "
         << top_ref->MoleculeByIndex(0)->getName() << endl;

    ex = CreateExclusionList(*top_ref->MoleculeByIndex(0),
                             *top->MoleculeByIndex(0));
    std::ofstream fl;
    fl.open(OptionsMap()["excl"].as<string>().c_str());
    fl << "# atomistic: " << top_ref->MoleculeByIndex(0)->getName()
       << " cg: " << top_ref->MoleculeByIndex(0)->getName()
       << " cgmap: " << OptionsMap()["cg"].as<string>() << endl;
    fl << *ex;
    fl.close();
    delete ex;
    return false;
  }
  return true;
}

ExclusionList *CsgBoltzmann::CreateExclusionList(Molecule &atomistic,
                                                 Molecule &cg) {
  ExclusionList *ex = new ExclusionList();
  // exclude all with all
  {
    list<Bead *> excl_list;
    for (int i = 0; i < atomistic.BeadCount(); ++i) {
      excl_list.push_back(atomistic.getBead(i));
    }
    ex->ExcludeList(excl_list);
  }

  // remove exclusions from inside a mapped bead
  Topology *at_top = atomistic.getParent();
  for (int i = 0; i < cg.BeadCount(); ++i) {
    const vector<long> &parent_beads = cg.getBead(i)->ParentBeads();
    list<Bead *> excl_list;

    for (const long &parent_bead_id : parent_beads) {
      excl_list.push_back(at_top->getBead(parent_bead_id));
    }
    ex->Remove(excl_list);
  }

  // remove exclusion which come from atomistic topology and hence bonds and
  // angles
  Topology *cg_top = cg.getParent();
  for (long i = 0; i < cg.BeadCount() - 1; ++i) {
    for (long j = i + 1; j < cg.BeadCount(); ++j) {
      if (cg_top->getExclusions().IsExcluded(cg.getBead(i), cg.getBead(j))) {
        const vector<long> &parent_beads_w = cg.getBead(i)->ParentBeads();
        const vector<long> &parent_beads_v = cg.getBead(j)->ParentBeads();

        for (const long parent_bead_id_w : parent_beads_w) {
          for (const long parent_bead_id_v : parent_beads_v) {
            ex->RemoveExclusion(at_top->getBead(parent_bead_id_w),
                                at_top->getBead(parent_bead_id_v));
          }
        }
      }
    }
  }
  return ex;
}

void CsgBoltzmann::Run() {
  CsgApplication::Run();
  if (OptionsMap().count("excl")) {
    return;
  }
  InteractiveMode();
}

void CsgBoltzmann::InteractiveMode() {
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
      "autocor <file> <selection>: calculate autocorrelation, only one row "
      "allowed in selection!\n"
      "cor <file> <selection>: calculate correlations, first row is correlated "
      "with all other rows";

  cout << help_text << endl;

  while (1) {
    string line;
    cout << "> ";
    getline(cin, line);

    boost::trim(line);

    votca::tools::Tokenizer tok(line, " \t");
    vector<string> args = tok.ToVector();

    if (args.size() == 0) {
      continue;
    }

    string cmd = args.front();
    args.erase(args.begin());

    if (cmd == "q") {
      break;
    }

    std::map<string, AnalysisTool *>::iterator tool;
    if (cmd == "help") {
      if (args.size() == 0) {
        cout << help_text << endl;
        continue;
      }
      cmd = args.front();
      args.erase(args.begin());
      tool = cmds.find(cmd);
      if (tool == cmds.end()) {
        cout << "error, no help item found" << endl;
        continue;
      }
      tool->second->Help(cmd, args);
      cout << endl;
      continue;
    }

    tool = cmds.find(cmd);
    if (tool == cmds.end()) {
      cout << "error, command not found" << endl;
      continue;
    }

    tool->second->Command(_bs, cmd, args);
  }
}

int main(int argc, char **argv) {
  CsgBoltzmann app;
  app.Exec(argc, argv);
}
