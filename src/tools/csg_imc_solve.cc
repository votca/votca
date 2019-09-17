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

#include "csg_imc_solve.h"
#include <fstream>
#include <votca/tools/table.h>

int main(int argc, char **argv) {
  CG_IMC_solve app;
  return app.Exec(argc, argv);
}

void CG_IMC_solve::Initialize(void) {
  namespace propt = boost::program_options;

  AddProgramOptions()("regularization,r",
                      propt::value<double>()->default_value(0.0),
                      "regularization factor");

  AddProgramOptions()("imcfile,i", boost::program_options::value<string>(),
                      "imc statefile");

  AddProgramOptions()("gmcfile,g", boost::program_options::value<string>(),
                      "gmc statefile");

  AddProgramOptions()("outputfile,o", boost::program_options::value<string>(),
                      "outputfile");
}

bool CG_IMC_solve::EvaluateOptions() {
  CheckRequired("imcfile", "Missing imcfile");
  CheckRequired("gmcfile", "Missing gmcfile");
  CheckRequired("outputfile", "Missing outputfile");
}

Eigen::VectorXd LoadFile(std::string filename) {
  std::ifstream intt;
  intt.open(filename);
  if (!intt) throw runtime_error(string("error, cannot open file ") + filename);

  std::string line;
  std::vector<double> result;
  while (getline(intt, line)) {
    if (line[0] == '#') {
      continue;
    }
    Tokenizer tok(line, " \t");
    vector<string> tokens = ToVector();
    result.append(std::stod(tokens[0]));
  }
  intt.close();

  return Eigen::Map<Eigen::VectorXd>(result.data(), result.size());
}

void CG_IMC_solve::Run() {
  string imcfile = _op_vm["imcfile"].as<string>();
  string gmcfile = _op_vm["gmcfile"].as<string>();
  string outputfile = _op_vm["outputfile"].as<string>();

  double reg = _op_vm["regularization"].as<double>();

  Eigen::VectorXd A = LoadFile(imcfile);
  Eigen::VectorXd B = LoadFile(gmcfile);
}
