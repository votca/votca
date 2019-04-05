
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE davidson_test


#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/davidsonsolver.h>
#include <votca/xtp/eigen.h>

using namespace votca::xtp;
using namespace std;

Eigen::MatrixXd init_matrix(int N, double eps)
{
    Eigen::MatrixXd matrix;
    matrix =  eps * Eigen::MatrixXd::Random(N,N);
    Eigen::MatrixXd tmat = matrix.transpose();
    matrix = matrix + tmat; 

    for (int i = 0; i<N; i++) {
        matrix(i,i) =  static_cast<double> (1. + (std::rand() %1000 ) / 10.);
    }
    return matrix;
}


BOOST_AUTO_TEST_SUITE(davidson_test)

BOOST_AUTO_TEST_CASE(davidson_full_matrix) {

    int size = 100;
    int neigen = 10;
    double eps = 0.01;
    Eigen::MatrixXd A = init_matrix(size,eps);

    votca::ctp::Logger log;
    DavidsonSolver DS(log);
    DS.solve(A,neigen);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    auto lambda = DS.eigenvalues();
    auto lambda_ref = es.eigenvalues().head(neigen);
    bool check_eigenvalues = lambda.isApprox(lambda_ref,1E-6);

    BOOST_CHECK_EQUAL(check_eigenvalues,1);

}

BOOST_AUTO_TEST_CASE(davidson_full_matrix_fail) {

    int size = 100;
    int neigen = 10;
    double eps = 0.01;
    Eigen::MatrixXd A = init_matrix(size,eps);

    votca::ctp::Logger log;
    DavidsonSolver DS(log);
    DS.set_iter_max(1);
    DS.solve(A,neigen);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    auto lambda = DS.eigenvalues();
    auto lambda_ref = es.eigenvalues().head(neigen);
    bool check_eigenvalues = lambda.isApprox(lambda_ref,1E-6);
    
    BOOST_CHECK_EQUAL(check_eigenvalues,0);

}


class TestOperator : public MatrixFreeOperator
{
    public : 
        TestOperator() {};
        Eigen::VectorXd col(int index) const;
        void set_diag();
        Eigen::VectorXd diag_el;

    private:
};

// constructors
void TestOperator::set_diag()
{
    int lsize = this->size();
    diag_el = Eigen::VectorXd::Zero(lsize);
    for (int i=0; i<lsize;i++){
        diag_el(i) = static_cast<double> (1. + (std::rand() %1000 ) / 10.);
    }
    
}

//  get a col of the operator
Eigen::VectorXd TestOperator::col(int index) const
{
    int lsize = this->size();
    Eigen::VectorXd col_out = Eigen::VectorXd::Zero(lsize,1);    
    for (int j=0; j < lsize; j++)
    {
        if (j==index) {
            col_out(j) = diag_el(j); 
        }
        else{
            col_out(j) = 0.01 / std::pow( static_cast<double>(j-index),2) ;
        }
    }
    return col_out;
}


BOOST_AUTO_TEST_CASE(davidson_matrix_free) {

    int size = 100;
    int neigen = 10;
    
    // Create Operator
    TestOperator Aop;
    Aop.set_size(size);
    Aop.set_diag();

    votca::ctp::Logger log;
    DavidsonSolver DS(log);
    DS.set_tolerance("normal");
    DS.set_size_update("safe");
    DS.solve(Aop,neigen);

    Eigen::MatrixXd A = Aop.get_full_matrix();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    auto lambda = DS.eigenvalues();
    auto lambda_ref = es.eigenvalues().head(neigen);
    bool check_eigenvalues = lambda.isApprox(lambda_ref,1E-6);

    BOOST_CHECK_EQUAL(check_eigenvalues,1);

}

BOOST_AUTO_TEST_SUITE_END()