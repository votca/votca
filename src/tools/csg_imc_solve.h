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

#ifndef _VOTCA_CSG_IMC_SOLVE_H
#define _VOTCA_CSG_IMC_SOLVE_H

#include <votca/csg/version.h>
#include <votca/tools/application.h>
#include <votca/tools/property.h>

/**
    \brief Solves linear system for IMCS
 *
 **/

class CG_IMC_solve : public votca::tools::Application {
 public:
  std::string ProgramName() { return "csg_imc_solve"; }
  void HelpText(std::ostream &out) {
    out << "Solves the linear system for IMCs";
  }

  void ShowHelpText(std::ostream &out) {
    std::string name = ProgramName();
    if (VersionString() != "") name = name + ", version " + VersionString();

    votca::csg::HelpTextHeader(name);
    HelpText(out);

    out << "\n\n" << VisibleOptions() << std::endl;
  }

  bool EvaluateOptions();
  void Initialize();
  void Run();

 private:
  Eigen::MatrixXd LoadMatrix(const std::string &filename) const;
};

#endif /* _VOTCA_CSG_IMC_SOLVE_H */
