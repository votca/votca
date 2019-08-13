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

void linalg_constrained_qrsolve(Eigen::VectorXd &x, Eigen::MatrixXd &A,
                                const Eigen::VectorXd &b,
                                const Eigen::MatrixXd &constr) {
  // check matrix for zero column

  bool nonzero_found = false;
  for (int j = 0; j < A.cols(); j++) {
    nonzero_found = A.col(j).isApproxToConstant(0.0, 1e-9);
    if (nonzero_found) {
      throw std::runtime_error("constrained_qrsolve_zero_column_in_matrix");
    }
  }

  const int NoVariables = x.size();
  const int NoConstrains =
      constr.rows();  // number of constraints is number of rows of constr

  Eigen::HouseholderQR<Eigen::MatrixXd> QR(constr.transpose());
  Eigen::MatrixXd Q = QR.householderQ();

  // Calculate A * Q and store the result in A
  A = A * Q;
  // A = [A1 A2], so A2 is just a block of A
  // [A1 A2] has N rows. A1 has ysize columns
  // A2 has 2*ngrid-ysize columns
  Eigen::MatrixXd A2 =
      A.block(0, NoConstrains, A.rows(), NoVariables - NoConstrains);
  // now perform QR-decomposition of A2 to solve the least-squares problem A2 *
  // z = b A2 has N rows and (2*ngrid-ysize) columns ->
  Eigen::HouseholderQR<Eigen::MatrixXd> QR2(A2);
  Eigen::VectorXd z = QR2.solve(b);

  // Next two steps assemble vector from y (which is zero-vector) and z
  Eigen::VectorXd result = Eigen::VectorXd::Zero(NoVariables);
  for (int i = NoConstrains; i < NoVariables; i++) {
    result[i] = z(i - NoConstrains);
  }
  // To get the final answer this vector should be multiplied by matrix Q
  x = Q * result;
  return;
}

EigenSystem linalg_eigenvalues(Eigen::MatrixXd &A, int nmax) {

  EigenSystem result;
#if defined(MKL)
  double wkopt;
  double *work;
  double abstol, vl, vu;

  MKL_INT lda;
  MKL_INT info;
  MKL_INT lwork;
  MKL_INT il, iu, m, ldz;

  int n = A.rows();
  MKL_INT ifail[n];
  lda = n;
  ldz = nmax;

  // make sure that containers for eigenvalues and eigenvectors are of correct
  // size
  result.eigenvalues().resize(nmax);
  result.eigenvectors().resize(n, nmax);

  lwork = -1;
  il = 1;
  iu = nmax;
  abstol = 0.0;  // use default
  vl = 0.0;
  vu = 0.0;
  // make a pointer to the EIGEN matrix so that LAPACK understands it
  double *pA = A.data();
  double *pV = result.eigenvectors().data();
  double *pE = result.eigenvalues().data();

  // call LAPACK via C interface
  info = LAPACKE_dsyevx(LAPACK_COL_MAJOR, 'V', 'I', 'U', n, pA, lda, vl, vu, il,
                        iu, abstol, &m, pE, pV, n, ifail);
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
