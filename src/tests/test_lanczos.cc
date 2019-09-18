
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE davidson_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/lanczossolver.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(lanczos_test)

class BlockOperator: public MatrixFreeOperator {
 public:
  BlockOperator() {};
  Eigen::RowVectorXd row(int index) const;
  void set_diag(int diag);
  Eigen::VectorXd diag_el;

 private:
  int _diag;
};

// constructors
void BlockOperator::set_diag(int diag) {
  int lsize = this->size();
  diag_el = Eigen::VectorXd::Zero(lsize);
  if (diag == 1) {
    for (int i = 0; i < lsize; i++) {
      diag_el(i) = static_cast<double>(1. + (std::rand() % 1000) / 10.);
    }
  }
}

//  get a col of the operator
Eigen::RowVectorXd BlockOperator::row(int index) const {
  int lsize = this->size();
  Eigen::RowVectorXd row_out = Eigen::RowVectorXd::Zero(lsize);
  for (int j = 0; j < lsize; j++) {
    if (j == index) {
      row_out(j) = diag_el(j);
    } else {
      row_out(j) = 0.01 / std::pow(static_cast<double>(j - index), 2);
    }
  }
  return row_out;
}

class HamiltonianOperator : public MatrixFreeOperator {
  /* construct an hamiltonian matrix from sub operator
  H = [ R   C 
       -C  -R ]
  we assume that R and C are symmetric
  */

 public:

  HamiltonianOperator(const BlockOperator &R, const BlockOperator &C);
  Eigen::RowVectorXd row(int index) const;
  void set_diag();
  Eigen::VectorXd diag_el;

 private:
  BlockOperator _R;
  BlockOperator _C;
};

// constructors
HamiltonianOperator::HamiltonianOperator(const BlockOperator &R, const BlockOperator &C)
{
  BlockOperator _R = R;
  BlockOperator _C = C;
  HamiltonianOperator::set_size( 2*_R.size() );
}

//  get a col of the operator
Eigen::RowVectorXd HamiltonianOperator::row(int index) const {
  int lsize = this->size();
  Eigen::RowVectorXd row_out = Eigen::RowVectorXd::Zero(lsize);

  Eigen::RowVectorXd r = this->_R.row(index);
  Eigen::RowVectorXd c = this->_C.row(index);

  if (index < lsize/2)
  {
    for (int j = 0; j < lsize/2; j++) {
        row_out(j) = r(j);
    }
    for (int j = lsize/2; j < lsize; j++){
        row_out(j) = c(j);
    }
  } else {
    for (int j = 0; j < lsize/2; j++) {
        row_out(j) = -c(j);
    }
    for (int j = lsize/2; j < lsize; j++) {
        row_out(j) = -r(j);
    }
  }
}

BOOST_AUTO_TEST_CASE(lanczos_matrix_free) {

  int size = 100;
  int neigen = 10;

  // Create Operator
  int diag = 1;
  BlockOperator Rop;
  Rop.set_size(size);
  Rop.set_diag(diag);

  diag=0;
  BlockOperator Cop;
  Cop.set_size(size);
  Cop.set_diag(0);

  // create Hamiltonian matrix
  HamiltonianOperator Hop(Rop,Cop);

  Logger log;
  LanczosSolver LS(log);
  LS.solve(Hop, neigen);

  Eigen::MatrixXd H = Hop.get_full_matrix();
  Eigen::EigenSolver<Eigen::MatrixXd> es(H);

  auto lambda = LS.eigenvalues();
  auto lambda_ref = es.eigenvalues().head(neigen);
  bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
}

BOOST_AUTO_TEST_SUITE_END()
