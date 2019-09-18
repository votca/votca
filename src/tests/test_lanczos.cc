
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE lanczos_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/lanczossolver.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/hamiltonianoperator.h>

using namespace votca::xtp;
using namespace std;



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


BOOST_AUTO_TEST_SUITE(lanczos_test)


BOOST_AUTO_TEST_CASE(lanczos_matrix_free) {

  int size = 10;
  int neigen = 1;

  Eigen::VectorXd b(size);
  b.setRandom();

  // Create Operator
  int diag = 1;
  BlockOperator Rop;
  Rop.set_size(size);
  Rop.set_diag(1);
  Eigen::MatrixXd R = Rop.get_full_matrix();
  //std::cout << "R matrix" << R << std::endl;

  diag=0;
  BlockOperator Cop;
  Cop.set_size(size);
  Cop.set_diag(0);
  Eigen::MatrixXd C = Cop.get_full_matrix();
  //std::cout << "C matrix" << C << std::endl;

  // create Hamiltonian matrix
  HamiltonianOperator<BlockOperator> Hop(Rop,Cop);
  std::cout << "Hop done" << std::endl;
  Eigen::MatrixXd H = Hop.get_full_matrix();

  Logger log;
  LanczosSolver LS(log);

  LS.solve(Hop, neigen);
  std::cout << "Solve done" << std::endl;

  // Eigen::MatrixXd H = Hop.get_full_matrix();
  // Eigen::EigenSolver<Eigen::MatrixXd> es(H);
  // std::cout << "check done" << std::endl;


  // auto lambda = LS.eigenvalues().real();
  // auto lambda_ref = es.eigenvalues().head(neigen);
  // bool check_eigenvalues = lambda.isApprox(lambda_ref, 1E-6);

  // BOOST_CHECK_EQUAL(check_eigenvalues, 1);
} 

BOOST_AUTO_TEST_SUITE_END()
