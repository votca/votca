
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE davidson_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/bseoperator_btda.h>
#include <votca/xtp/davidsonsolver.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>

using namespace votca::xtp;
using namespace std;

Eigen::MatrixXd symm_matrix(int N, double eps) {
  Eigen::MatrixXd matrix;
  matrix = eps * Eigen::MatrixXd::Random(N, N);
  Eigen::MatrixXd tmat = matrix.transpose();
  matrix = matrix + tmat;
  return matrix;
}

Eigen::MatrixXd init_matrix(int N, double eps) {
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = i; j < N; j++) {
      if (i == j) {
        matrix(i, i) = std::sqrt(static_cast<double>(1 + i));
      } else {
        matrix(i, j) = eps / std::pow(static_cast<double>(j - i), 2);
        matrix(j, i) = eps / std::pow(static_cast<double>(j - i), 2);
      }
    }
  }
  return matrix;
}

BOOST_AUTO_TEST_SUITE(davidson_test)

BOOST_AUTO_TEST_CASE(davidson_full_matrix) {

  int size = 100;
  int neigen = 10;
  double eps = 0.01;
  Eigen::MatrixXd A = init_matrix(size, eps);
  Logger log;
  DavidsonSolver DS(log);
  DS.set_ortho("QR");
  DS.solve(A, neigen);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
}

BOOST_AUTO_TEST_CASE(davidson_full_matrix_large) {

  int size = 400;
  int neigen = 10;
  double eps = 0.01;
  Eigen::MatrixXd A = init_matrix(size, eps);
  Logger log;
  DavidsonSolver DS(log);
  DS.solve(A, neigen);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
}

BOOST_AUTO_TEST_CASE(davidson_full_matrix_fail) {

  int size = 100;
  int neigen = 10;
  double eps = 0.01;
  Eigen::MatrixXd A = init_matrix(size, eps);

  Logger log;
  DavidsonSolver DS(log);
  DS.set_iter_max(1);
  DS.solve(A, neigen);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 0);
}

class TestOperator : public MatrixFreeOperator {
 public:
  TestOperator() = default;
  ;
  Eigen::RowVectorXd OperatorRow(int index) const override;

 private:
};

//  get a col of the operator
Eigen::RowVectorXd TestOperator::OperatorRow(int index) const {
  int lsize = this->size();
  Eigen::RowVectorXd row_out = Eigen::RowVectorXd::Zero(lsize);
  for (int j = 0; j < lsize; j++) {
    if (j == index) {
      row_out(j) = std::sqrt(static_cast<double>(index + 1));
    } else {
      row_out(j) = 0.01 / std::pow(static_cast<double>(j - index), 2);
    }
  }
  return row_out;
}

BOOST_AUTO_TEST_CASE(davidson_matrix_free) {

  int size = 100;
  int neigen = 10;

  // Create Operator
  TestOperator Aop;
  Aop.set_size(size);

  Logger log;
  DavidsonSolver DS(log);
  DS.set_tolerance("normal");
  DS.set_size_update("safe");
  DS.set_ortho("QR");
  DS.solve(Aop, neigen);

  Eigen::MatrixXd A = Aop.get_full_matrix();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(davidson_matrix_free_large) {

  int size = 400;
  int neigen = 10;

  // Create Operator
  TestOperator Aop;
  Aop.set_size(size);

  Logger log;
  log.setReportLevel(logDEBUG);
  DavidsonSolver DS(log);
  DS.set_tolerance("normal");
  DS.set_size_update("safe");
  DS.solve(Aop, neigen);
  std::cout << log << std::endl;
  Eigen::MatrixXd A = Aop.get_full_matrix();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }
}

class BlockOperator : public MatrixFreeOperator {
 public:
  BlockOperator() = default;
  ;
  Eigen::MatrixXd OperatorBlock(int row, int col) const override;

  bool useRow() const override { return false; }
  bool useBlock() const override { return true; }
  int getBlocksize() const override { return size() / 10; }

 private:
};

//  get a block of the operator
Eigen::MatrixXd BlockOperator::OperatorBlock(int row, int col) const {
  int blocksize = getBlocksize();
  Eigen::MatrixXd block = Eigen::MatrixXd::Zero(blocksize, blocksize);
  int blocdisttodiagonal = std::abs(row - col) * blocksize;
  for (int i_col = 0; i_col < blocksize; i_col++) {
    for (int i_row = 0; i_row < blocksize; i_row++) {
      block(i_row, i_col) =
          0.01 / std::pow(static_cast<double>(std::abs(i_row - i_col) +
                                              blocdisttodiagonal),
                          2);
    }
  }
  if (blocdisttodiagonal == 0) {
    for (int i = 0; i < blocksize; i++) {
      block(i, i) = std::sqrt(static_cast<double>(row * blocksize + i + 1));
    }
  }

  return block;
}

BOOST_AUTO_TEST_CASE(davidson_matrix_free_block) {

  int size = 100;
  int neigen = 10;

  // Create Operator
  BlockOperator Aop;
  Aop.set_size(size);

  Logger log;
  DavidsonSolver DS(log);
  DS.set_tolerance("normal");
  DS.set_ortho("QR");
  DS.set_size_update("safe");

  Eigen::MatrixXd A = Aop.get_full_matrix();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  DS.solve(Aop, neigen);
  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(davidson_matrix_free_block_large) {

  int size = 400;
  int neigen = 10;

  // Create Operator
  BlockOperator Aop;
  Aop.set_size(size);

  Logger log;
  DavidsonSolver DS(log);
  DS.set_tolerance("normal");
  DS.set_size_update("safe");

  Eigen::MatrixXd A = Aop.get_full_matrix();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  DS.solve(Aop, neigen);
  auto lambda = DS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
  if (!check_eigenvalues) {
    std::cout << "ref" << std::endl;
    std::cout << es.eigenvalues().head(neigen).transpose() << std::endl;
    std::cout << "result" << std::endl;
    std::cout << DS.eigenvalues().transpose() << std::endl;
  }
}

Eigen::ArrayXi index_eval(Eigen::VectorXd ev, int neigen) {

  int nev = ev.rows();
  int npos = nev / 2;

  Eigen::ArrayXi idx = Eigen::ArrayXi::Zero(npos);
  int nstored = 0;

  // get only positives
  for (int i = 0; i < nev; i++) {
    if (ev(i) > 0) {
      idx(nstored) = i;
      nstored++;
    }
  }

  // sort the epos eigenvalues
  std::sort(idx.data(), idx.data() + idx.size(),
            [&](int i1, int i2) { return ev[i1] < ev[i2]; });
  return idx.head(neigen);
}

Eigen::MatrixXd extract_eigenvectors(const Eigen::MatrixXd &V,
                                     const Eigen::ArrayXi &idx) {
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(V.rows(), idx.size());
  for (int i = 0; i < idx.size(); i++) {
    W.col(i) = V.col(idx(i));
  }
  return W;
}

class HermitianBlockOperator : public MatrixFreeOperator {
 public:
  HermitianBlockOperator() = default;
  ;

  void attach_matrix(const Eigen::MatrixXd &mat);
  Eigen::RowVectorXd OperatorRow(int index) const override;
  void set_diag(int diag);
  Eigen::VectorXd diag_el;

 private:
  Eigen::MatrixXd _mat;
};

void HermitianBlockOperator::attach_matrix(const Eigen::MatrixXd &mat) {
  _mat = mat;
}

//  get a col of the operator
Eigen::RowVectorXd HermitianBlockOperator::OperatorRow(int index) const {
  return _mat.row(index);
}

BOOST_AUTO_TEST_CASE(davidson_hamiltonian_matrix_free) {

  int size = 60;
  int neigen = 5;
  Logger log;

  // Create Operator
  HermitianBlockOperator Rop;
  Rop.set_size(size);
  Eigen::MatrixXd rmat = init_matrix(size, 0.01);
  Rop.attach_matrix(rmat);

  HermitianBlockOperator Cop;
  Cop.set_size(size);
  Eigen::MatrixXd cmat = symm_matrix(size, 0.01);
  Cop.attach_matrix(cmat);

  // create Hamiltonian operator
  HamiltonianOperator<HermitianBlockOperator, HermitianBlockOperator> Hop(Rop,
                                                                          Cop);

  DavidsonSolver DS(log);
  DS.set_tolerance("normal");
  DS.set_size_update("max");
  DS.set_ortho("QR");
  DS.set_matrix_type("HAM");
  DS.solve(Hop, neigen);

  auto lambda = DS.eigenvalues().real();
  std::sort(lambda.data(), lambda.data() + lambda.size());
  Eigen::MatrixXd H = Hop.get_full_matrix();

  Eigen::EigenSolver<Eigen::MatrixXd> es(H);
  Eigen::ArrayXi idx = index_eval(es.eigenvalues().real(), neigen);
  Eigen::VectorXd lambda_ref = idx.unaryExpr(es.eigenvalues().real());

  bool check_eigenvalues = lambda.isApprox(lambda_ref.head(neigen), 1E-6);
  if (!check_eigenvalues) {
    cout << "Davidson not converged after " << DS.num_iterations()
         << " iterations" << endl;
    cout << "Reference eigenvalues" << endl;
    cout << lambda_ref.head(neigen) << endl;
    cout << "Davidson eigenvalues" << endl;
    cout << lambda << endl;
    cout << "Residue norms" << endl;
    cout << DS.residues() << endl;
  }
  BOOST_CHECK_EQUAL(check_eigenvalues, 1);

  Eigen::MatrixXd evect_dav = DS.eigenvectors().real();
  Eigen::MatrixXd evect = es.eigenvectors().real();
  Eigen::MatrixXd evect_ref = extract_eigenvectors(evect, idx);

  bool check_eigenvectors =
      evect_ref.cwiseAbs2().isApprox(evect_dav.cwiseAbs2(), 0.001);
  BOOST_CHECK_EQUAL(check_eigenvectors, 1);
}

BOOST_AUTO_TEST_SUITE_END()
