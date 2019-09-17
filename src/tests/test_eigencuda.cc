
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE eigencuda_test

#include <boost/test/unit_test.hpp>
#include <votca/xtp/eigen.h>
#include <votca/xtp/eigencuda.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(eigecuda_test)

BOOST_AUTO_TEST_CASE(right_matrix_multiplication) {
  // Call the class to handle GPU resources
  EigenCuda<double> EC;

  // Call matrix multiplication GPU
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(3, 2);

  // Define matrices
  A << 1., 2., 3., 4.;
  B << 5., 6., 7., 8., 9., 10.;
  C << 9., 10., 11., 12., 13., 14.;
  D << 13., 14., 15., 16., 17., 18.;
  X << 23., 34., 31., 46., 39., 58.;
  Y << 39., 58., 47., 70., 55., 82.;
  Z << 55., 82., 63., 94., 71., 106.;

  std::vector<Eigen::MatrixXd> tensor{B, C, D};
  std::vector<Eigen::MatrixXd> rs = EC.right_matrix_tensor(A, tensor);

  // Expected results
  BOOST_TEST(X.isApprox(rs[0]));
  BOOST_TEST(Y.isApprox(rs[1]));
  BOOST_TEST(Z.isApprox(rs[2]));
}

BOOST_AUTO_TEST_CASE(triple_tensor_product) {
  // Call the class to handle GPU resources
  EigenCuda<double> EC;

  // Define matrices
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 3);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, 3);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2, 2);
  Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(2, 2);

  A << 1., 2., 3., 4., 5., 6.;
  B << 5., 6., 7., 8., 9., 10.;
  C << 9., 10., 11., 12., 13., 14., 15., 16., 17.;
  D << 13., 14., 15., 16., 17., 18., 19., 20., 21.;
  X << 1788., 2040., 4281., 4884.;
  Y << 2292., 2616., 5541., 6324.;

  std::vector<Eigen::MatrixXd> tensor{C, D};
  std::vector<Eigen::MatrixXd> rs = EC.triple_tensor_product(A, B, tensor);

  // Expected results
  BOOST_TEST(X.isApprox(rs[0]));
  BOOST_TEST(Y.isApprox(rs[1]));
}

BOOST_AUTO_TEST_SUITE_END()
