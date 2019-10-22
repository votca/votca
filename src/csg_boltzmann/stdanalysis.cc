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

void StdAnalysis::Register(std::map<std::string, AnalysisTool *> &lib) {
  lib["list"] = this;
  lib["vals"] = this;
  lib["cor"] = this;
  lib["autocor"] = this;
}

void StdAnalysis::Command(BondedStatistics &bs, const std::string &cmd,
                          std::vector<std::string> &args) {
  if (cmd == "vals") {
    WriteValues(bs, args);
  }
  if (cmd == "cor") {
    WriteCorrelations(bs, args);
  }
  if (cmd == "autocor") {
    WriteAutocorrelation(bs, args);
  }
  if (cmd == "list") {
    votca::tools::DataCollection<double>::selection *sel =
        bs.BondedValues().select("*");
    std::cout << "Available bonded interactions:" << std::endl;
    for (auto &array : *sel) {
      std::cout << array->getName() << " " << std::endl;
    }
    delete sel;
  }
}

void StdAnalysis::Help(const std::string &cmd, std::vector<std::string> &args) {
  if (cmd == "vals") {
    std::cout
        << "vals <file> <selection>\n"
        << "write values to file. The first row is the frame number, then one "
        << "row for each interaction specified. The output can be used to "
           "generate "
        << "2D correlation plots.\n\n"
        << "example: vals angle *angle*\n";
  }
  if (cmd == "cor") {
    std::cout
        << "cor <file> <selection>\n"
        << "Calculate linear correlation coefficient of the first item in "
           "selection with all the other items\n"
        << "WARNING: for evaluating correlations in the system, it is not "
           "sufficient to calculate the "
        << "linear correlation coefficient, 2D histograms with data from the "
           "vals command should be used instead!\n";
  }
  if (cmd == "autocor") {
    std::cout
        << "autocor <file> <interaction>\n"
        << "calculate autocorrelation function of first item in selection. "
           "The output is periodic since FFTW3 is used to "
           "calcualte correlations.\n";
  }
  if (cmd == "list") {
    std::cout << "list\nlists all available interactions\n";
  }
}

void StdAnalysis::WriteValues(BondedStatistics &bs,
                              std::vector<std::string> &args) {
  std::ofstream out;

  votca::tools::DataCollection<double>::selection *sel = nullptr;

  for (size_t i = 1; i < args.size(); i++) {
    sel = bs.BondedValues().select(args[i], sel);
  }

  out.open(args[0]);
  out << *sel << std::endl;
  out.close();
  std::cout << "written " << sel->size() << " data rows to " << args[0]
            << std::endl;
  delete sel;
}

void StdAnalysis::WriteAutocorrelation(BondedStatistics &bs,
                                       std::vector<std::string> &args) {
  std::ofstream out;
  votca::tools::DataCollection<double>::selection *sel = nullptr;

  for (size_t i = 1; i < args.size(); i++) {
    sel = bs.BondedValues().select(args[i], sel);
  }

  votca::tools::CrossCorrelate c;
  c.AutoCorrelate(sel);
  out.open(args[0]);
  out << c << std::endl;
  out.close();
  std::cout << "calculated autocorrelation for " << sel->size()
            << " data rows, written to " << args[0] << std::endl;
  delete sel;
}

void StdAnalysis::WriteCorrelations(BondedStatistics &bs,
                                    std::vector<std::string> &args) {
  std::ofstream out;
  votca::tools::DataCollection<double>::selection *sel = nullptr;

  for (size_t i = 1; i < args.size(); i++) {
    sel = bs.BondedValues().select(args[i], sel);
  }

  votca::tools::Correlate c;
  c.CalcCorrelations(sel);
  out.open(args[0]);
  out << c << std::endl;
  out.close();
  std::cout << "calculated correlations for " << sel->size()
            << " rows, written to " << args[0] << std::endl;
  delete sel;
}

}  // namespace csg
}  // namespace votca
