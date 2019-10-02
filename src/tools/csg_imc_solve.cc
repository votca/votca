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
#include <votca/csg/imcio.h>
#include <votca/tools/table.h>

int main(int argc, char** argv) {
  CG_IMC_solve app;
  return app.Exec(argc, argv);
}

void CG_IMC_solve::Initialize(void) {
  namespace propt = boost::program_options;

  AddProgramOptions()("regularization,r",
                      propt::value<double>()->default_value(0.0),
                      "regularization factor");

  AddProgramOptions()("imcfile,i", boost::program_options::value<std::string>(),
                      "imc statefile");

  AddProgramOptions()("gmcfile,g", boost::program_options::value<std::string>(),
                      "gmc statefile");

  AddProgramOptions()("idxfile,n", boost::program_options::value<std::string>(),
                      "[OPTIONAL] idx statefile");

  AddProgramOptions()("outputfile,o",
                      boost::program_options::value<std::string>(),
                      "outputfile");
}

bool CG_IMC_solve::EvaluateOptions() {
  CheckRequired("imcfile", "Missing imcfile");
  CheckRequired("gmcfile", "Missing gmcfile");
  CheckRequired("outputfile", "Missing outputfile");
  return true;
}

void CG_IMC_solve::Run() {
  std::string imcfile = _op_vm["imcfile"].as<std::string>();
  std::string gmcfile = _op_vm["gmcfile"].as<std::string>();
  std::string outputfile = _op_vm["outputfile"].as<std::string>();

  double reg = _op_vm["regularization"].as<double>();

  Eigen::MatrixXd A = votca::csg::imcio_read_matrix(gmcfile);
  votca::tools::Table B;
  B.Load(imcfile);

  votca::tools::Table x;
  x.resize(B.size());
  x.x() = B.x();

  // Calculating https://en.wikipedia.org/wiki/Tikhonov_regularization

  Eigen::MatrixXd ATA = A.transpose() * A;

  // We could add and invert but for stability reasons we will diagonalize the
  // matrix which is symmetric by construction

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ATA);

  Eigen::VectorXd inv_diag = Eigen::VectorXd::Zero(es.eigenvalues().size());
  double etol = 1e-12;

  int numerics_unstable = 0;
  // Constructing pseudoinverse if not stable
  // https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
  for (int i = 0; i < es.eigenvalues().size(); i++) {
    if (std::abs(es.eigenvalues()[i] + reg) < etol) {
      numerics_unstable++;
    } else {
      inv_diag[i] = 1 / (es.eigenvalues()[i] + reg);
    }
  }

  if (numerics_unstable > 0) {
    std::cout
        << "Regularisation parameter was to small, a pseudo inverse was "
           "constructed instead.\n Use a larger regularisation parameter R."
           " Smallest (eigenvalue+R)="
        << (es.eigenvalues().array() + reg).abs().minCoeff() << " Found "
        << numerics_unstable << " eigenvalues of " << es.eigenvalues().size()
        << " below " << etol << std::endl;
  }

  Eigen::MatrixXd inverse =
      es.eigenvectors() * inv_diag.asDiagonal() * es.eigenvectors().transpose();

  x.y() = -inverse * A.transpose() * B.y();

  x.Save(outputfile);

  if (OptionsMap().count("idxfile")) {
    std::string idxfile = _op_vm["idxfile"].as<std::string>();
    std::vector<std::pair<std::string, votca::tools::RangeParser> > ranges =
        votca::csg::imcio_read_index(idxfile);
    for (std::pair<std::string, votca::tools::RangeParser>& range : ranges) {
      votca::tools::Table tbl;
      for (int r : range.second) {
        tbl.push_back(x.x(r - 1), x.y(r - 1), 'i');
      }
      tbl.Save(range.first + ".dpot.imc");
    }
  }
}
