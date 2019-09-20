
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE lanczos_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/lanczossolver.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/bseoperator_btda.h>

using namespace votca::xtp;
using namespace std;

Eigen::VectorXd sort_ev(Eigen::VectorXd ev)
{

  int nev = ev.rows();
  int npos = nev/2;

  Eigen::VectorXd ev_pos = Eigen::VectorXd::Zero(npos);
  int nstored = 0;

  // get only positives
  for (int i=0; i<nev; i++) {
    if (ev(i) > 0) {
      ev_pos(nstored) = ev(i);
      nstored ++;
    }
  }
  
  // sort the epos eigenvalues
  std::sort(ev_pos.data(), ev_pos.data() + ev_pos.size());
  return ev_pos;
}

class BlockOperator: public MatrixFreeOperator {
 public:
  BlockOperator(double alpha) : _alpha(alpha) {};
  Eigen::RowVectorXd row(int index) const;
  void set_diag(int diag);
  Eigen::VectorXd diag_el;

 private:
  int _diag;
  int _alpha;
};

// constructors
void BlockOperator::set_diag(int diag) {
  int lsize = this->size();
  diag_el = Eigen::VectorXd::Zero(lsize);
  if (diag == 1) {
    for (int i = 0; i < lsize; i++) {
      diag_el(i) = static_cast<double>(1. + (std::rand() % 1000) / 10. );
    }
  }

  if (diag == 0) {
    diag_el = Eigen::VectorXd::Random(lsize);
    // for (int i = 0; i < lsize; i++) {
    //   diag_el(i) = static_cast<double>( 1. + (std::rand() % 1000) / 10. ) * 0.1;
    // }
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
      row_out(j) = _alpha / std::pow(static_cast<double>(j - index), 2);
    }
  }
  return row_out;
}


BOOST_AUTO_TEST_SUITE(lanczos_test)

BOOST_AUTO_TEST_CASE(lanczos_matrix_free) {

  int size = 60;
  int neigen = 5;
  Logger log;

  // Create Operator
  BlockOperator Rop(0.01);
  Rop.set_size(size);
  Rop.set_diag(1);

  BlockOperator Cop(0.005);
  Cop.set_size(size);
  Cop.set_diag(0);

  // create Hamiltonian operator
  HamiltonianOperator<BlockOperator,BlockOperator> Hop(Rop,Cop);
  
  // Lanczos solver
  LanczosSolver LS(log);
  LS.solve(Hop, neigen);
  auto lambda = LS.eigenvalues().real();
  std::sort(lambda.data(), lambda.data() + lambda.size());  

  Eigen::MatrixXd H = Hop.get_full_matrix();
  Eigen::EigenSolver<Eigen::MatrixXd> es(H);
  auto lambda_ref = sort_ev(es.eigenvalues().real());
  
  //std::cout << "\nlanczos : \n" << lambda << std::endl;
  //std::cout << "\neigen : \n" << lambda_ref.head(neigen) << std::endl;

  bool check_eigenvalues = lambda.isApprox(lambda_ref.head(neigen), 1E-6);

  if (!check_eigenvalues) {
    cout << "Reference eigenvalues" << endl;
    cout << lambda_ref << endl;
    cout << "Lanczos eigenvalues" << endl;
    cout << lambda << endl;
  }

  BOOST_CHECK_EQUAL(check_eigenvalues, 1);
} 

BOOST_AUTO_TEST_SUITE_END()
