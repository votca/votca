
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE eigencuda_test

#include <boost/test/unit_test.hpp>
#include <string>
#include <votca/xtp/eigen.h>
#include <votca/xtp/eigencuda.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(eigecuda_test)

BOOST_AUTO_TEST_CASE(right_matrix_multiplication) {
  // Call the class to handle GPU resources
  EigenCuda EC;

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
  EC.right_matrix_tensor_mult(tensor, A);

  // Expected results
  BOOST_TEST(X.isApprox(tensor[0]));
  BOOST_TEST(Y.isApprox(tensor[1]));
  BOOST_TEST(Z.isApprox(tensor[2]));
}

BOOST_AUTO_TEST_CASE(wrong_shape_cublas) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(2, 2);
  Eigen::MatrixXd B = Eigen::MatrixXd::Random(5, 5);

  EigenCuda EC;
  std::vector<Eigen::MatrixXd> tensor{B};

  BOOST_REQUIRE_THROW(EC.right_matrix_tensor_mult(tensor, A),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(triple_matrix_multiplication) {
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 2);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2, 3);
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3, 2);
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(2, 2);
  A << 1., 2., 3., 4.;
  B << 5., 6., 7., 8., 9., 10.;
  C << 9., 10., 11., 12., 13., 14.;
  X << 804., 876., 1810., 1972.;

  EigenCuda cuda_handle;
  const cudaStream_t stream = cuda_handle.get_stream();
  CudaMatrix matrixA{2, 2};
  CudaMatrix matrixC{3, 2};
  matrixA.copy_matrix_to_gpu(A, stream);
  matrixC.copy_matrix_to_gpu(C, stream);

  Eigen::MatrixXd result = cuda_handle.triple_matrix_mult(matrixA, B, matrixC);
  BOOST_TEST(X.isApprox(result));
}

BOOST_AUTO_TEST_SUITE_END()
