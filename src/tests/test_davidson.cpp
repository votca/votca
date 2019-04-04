
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE davidson_test


#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/matrixfreeoperator.h>
#include <votca/xtp/davidsonsolver.h>
#include <votca/xt/eigen.h>

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

class TestOperator : public MatrixFreeOperator
{
    public : 
    TestOperator(int n) {_size = n;}
    Eigen::VectorXd col(int index) const;
};

//  get a col of the operator
Eigen::VectorXd TestOperator::col(int index) const
{
    Eigen::VectorXd col_out = Eigen::VectorXd::Zero(_size,1);    
    for (int j=0; j < _size; j++)
    {
        if (j==index) {
            col_out(j) = static_cast<double> (j+1); 
        }
        else{
            col_out(j) = 0.01 / std::pow( static_cast<double>(j-index),2) ;
        }
    }
    return col_out;
}


BOOST_AUTO_TEST_SUITE(davidson_test)

BOOST_AUTO_TEST_CASE(davidson_full_matrix) {

    int size = 100;
    int neigen = 10;
    double eps = 0.01;
    Eigen::MatrixXd A = init_matrix(size,eps);

    DavidsonSolver DS;
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

    DavidsonSolver DS;
    DS.set_iter_max(1);
    DS.solve(A,neigen);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    auto lambda = DS.eigenvalues();
    auto lambda_ref = es.eigenvalues().head(neigen);
    bool check_eigenvalues = lambda.isApprox(lambda_ref,1E-6);
    
    BOOST_CHECK_EQUAL(check_eigenvalues,0);

}



BOOST_AUTO_TEST_CASE(davidson_matrix_free) {

    int size = 100;
    int neigen = 10;
    
    // Create Operator
    TestOperator Aop(size);
    DavidsonSolver DS;
    DS.solve(Aop,neigen);

    Eigen::MatrixXd A = Aop.get_full_mat();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);

    auto lambda = DS.eigenvalues();
    auto lambda_ref = es.eigenvalues().head(neigen);
    bool check_eigenvalues = lambda.isApprox(lambda_ref,1E-6);
    
    BOOST_CHECK_EQUAL(check_eigenvalues,1);

}

BOOST_AUTO_TEST_SUITE_END()