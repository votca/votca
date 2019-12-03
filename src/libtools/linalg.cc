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
#include <votca/tools/linalg.h>

#include <iostream>
#include <sstream>
namespace votca {
namespace tools {

Eigen::VectorXd linalg_constrained_qrsolve(const Eigen::MatrixXd &A,
                                           const Eigen::VectorXd &b,
                                           const Eigen::MatrixXd &constr) {
  // check matrix for zero column

  bool nonzero_found = false;
  for (Index j = 0; j < A.cols(); j++) {
    nonzero_found = A.col(j).isApproxToConstant(0.0, 1e-9);
    if (nonzero_found) {
      throw std::runtime_error("constrained_qrsolve_zero_column_in_matrix");
    }
  }

  const long int NoVariables = A.cols();
  const Index NoConstrains =
      constr.rows();  // number of constraints is number of rows of constr
  const long int deg_of_freedom = NoVariables - NoConstrains;

  Eigen::HouseholderQR<Eigen::MatrixXd> QR(constr.transpose());

  // Calculate A * Q and store the result in A
  auto A_new = A * QR.householderQ();
  // A_new = [A1 A2], so A2 is just a block of A
  // [A1 A2] has N rows. A1 has ysize columns
  // A2 has 2*ngrid-ysize columns
  Eigen::MatrixXd A2 = A_new.rightCols(deg_of_freedom);
  // now perform QR-decomposition of A2 to solve the least-squares problem A2 *
  // z = b A2 has N rows and (2*ngrid-ysize) columns ->
  Eigen::HouseholderQR<Eigen::MatrixXd> QR2(A2);
  Eigen::VectorXd z = QR2.solve(b);

  // Next two steps assemble vector from y (which is zero-vector) and z
  Eigen::VectorXd result = Eigen::VectorXd::Zero(NoVariables);
  result.tail(deg_of_freedom) = z;
  // To get the final answer this vector should be multiplied by matrix Q
  return QR.householderQ() * result;
}

EigenSystem linalg_eigenvalues(Eigen::MatrixXd &A, Index nmax) {

  EigenSystem result;
#ifdef MKL_FOUND

  Index n = A.rows();
  std::vector<MKL_INT> ifail(n);
  MKL_INT lda = MKL_INT(n);
  // make sure that containers for eigenvalues and eigenvectors are of correct
  // size
  result.eigenvalues().resize(nmax);
  result.eigenvectors().resize(n, nmax);

  MKL_INT il = 1;
  MKL_INT iu = MKL_INT(nmax);
  double abstol = 0.0;  // use default
  double vl = 0.0;
  double vu = 0.0;
  // make a pointer to the EIGEN matrix so that LAPACK understands it
  double *pA = A.data();
  double *pV = result.eigenvectors().data();
  double *pE = result.eigenvalues().data();

  MKL_INT info;
  MKL_INT m;

  // call LAPACK via C interface
  info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, 'V', 'I', 'U', lda, pA, lda, vl, vu,
                        il, iu, abstol, &m, pE, pV, lda, ifail.data());
  if (info == 0) {
    result.info() = Eigen::Success;
  } else if (info < 0) {
    result.info() = Eigen::NumericalIssue;
  } else {
    result.info() = Eigen::NoConvergence;
  }

#else
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  result.eigenvectors() = es.eigenvectors().leftCols(nmax);
  result.eigenvalues() = es.eigenvalues().head(nmax);
  result.info() = es.info();
#endif
  return result;
}

}  // namespace tools
}  // namespace votca
