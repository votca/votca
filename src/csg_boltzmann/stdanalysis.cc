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

#include "stdanalysis.h"
#include "analysistool.h"
#include "bondedstatistics.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <votca/tools/correlate.h>
#include <votca/tools/crosscorrelate.h>

namespace votca {
namespace csg {

void StdAnalysis::Register(map<string, AnalysisTool *> &lib) {
  lib["list"] = this;
  lib["vals"] = this;
  lib["cor"] = this;
  lib["autocor"] = this;
}

void StdAnalysis::Command(BondedStatistics &bs, string cmd,
                          vector<string> &args) {
  if (cmd == "vals") WriteValues(bs, args);
  if (cmd == "cor") WriteCorrelations(bs, args);
  if (cmd == "autocor") WriteAutocorrelation(bs, args);
  if (cmd == "list") {
    DataCollection<double>::selection *sel = bs.BondedValues().select("*");
    DataCollection<double>::selection::iterator i;
    cout << "Available bonded interactions:" << endl;
    for (i = sel->begin(); i != sel->end(); ++i)
      cout << (*i)->getName()
           << " ";  // << "[" << (*i).second->size() << "]" << " ";
    cout << endl;
    delete sel;
  }
}

void StdAnalysis::Help(string cmd, vector<string> &args) {
  if (cmd == "vals") {
    cout << "vals <file> <selection>\n"
         << "write values to file. The first row is the frame number, then one "
         << "row for each interaction specified. The output can be used to "
            "generate "
         << "2D correlation plots.\n\n"
         << "example: vals angle *angle*\n";
  }
  if (cmd == "cor") {
    cout << "cor <file> <selection>\n"
         << "Calculate linear correlation coefficient of the first item in "
            "selection with all the other items\n"
         << "WARNING: for evaluating correlations in the system, it is not "
            "sufficient to calculate the "
         << "linear correlation coefficient, 2D histograms with data from the "
            "vals command should be used instead!\n";
  }
  if (cmd == "autocor") {
    cout << "autocor <file> <interaction>\n"
         << "calculate autocorrelation function of first item in selection. "
            "The output is periodic since FFTW3 is used to "
            "calcualte correlations.\n";
  }
  if (cmd == "list") {
    cout << "list\nlists all available interactions\n";
  }
}

void StdAnalysis::WriteValues(BondedStatistics &bs, vector<string> &args) {
  ofstream out;

  DataCollection<double>::selection *sel = NULL;

  for (size_t i = 1; i < args.size(); i++)
    sel = bs.BondedValues().select(args[i], sel);

  out.open(args[0].c_str());
  out << *sel << endl;
  out.close();
  cout << "written " << sel->size() << " data rows to " << args[0] << endl;
  delete sel;
}

void StdAnalysis::WriteAutocorrelation(BondedStatistics &bs,
                                       vector<string> &args) {
  ofstream out;
  DataCollection<double>::selection *sel = NULL;

  for (size_t i = 1; i < args.size(); i++)
    sel = bs.BondedValues().select(args[i], sel);

  CrossCorrelate c;
  c.AutoCorrelate(sel, false);
  out.open(args[0].c_str());
  out << c << endl;
  out.close();
  cout << "calculated autocorrelation for " << sel->size()
       << " data rows, written to " << args[0] << endl;
  delete sel;
}

void StdAnalysis::WriteCorrelations(BondedStatistics &bs,
                                    vector<string> &args) {
  ofstream out;
  DataCollection<double>::selection *sel = NULL;

  for (size_t i = 1; i < args.size(); i++)
    sel = bs.BondedValues().select(args[i], sel);

  Correlate c;
  c.CalcCorrelations(sel);
  out.open(args[0].c_str());
  out << c << endl;
  out.close();
  cout << "calculated correlations for " << sel->size() << " rows, written to "
       << args[0] << endl;
  delete sel;
}

}  // namespace csg
}  // namespace votca
