
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE dipoledipoleinteraction_test

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <votca/xtp/dipoledipoleinteraction.h>
#include <votca/xtp/eigen.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(dipoledipoleinteraction_test)

BOOST_AUTO_TEST_CASE(dipoledipoleinteraction_test) {

  PolarSite zero(0, "H", Eigen::Vector3d::UnitX());
  PolarSegment seg0("zero", 0);
  seg0.push_back(zero);
  PolarSite one(1, "C", Eigen::Vector3d::UnitZ());
  PolarSegment seg1("one", 1);
  seg1.push_back(one);

  std::vector<PolarSegment> segs;
  segs.push_back(seg0);
  segs.push_back(seg1);
  double tholedamp = 0.39;
  eeInteractor interactor(tholedamp);

  DipoleDipoleInteraction dipdip(interactor, segs);

  // building reference
  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(6, 6);
  PolarSegment seg_ref("ref", 100);
  seg_ref.push_back(zero);
  seg_ref.push_back(one);
  for (Index i = 1; i < seg_ref.size(); i++) {
    for (Index j = 0; j < i; j++) {
      ref.block<3, 3>(3 * i, 3 * j) =
          interactor.FillTholeInteraction(seg_ref[i], seg_ref[j]);
      ref.block<3, 3>(3 * j, 3 * i) =
          interactor.FillTholeInteraction(seg_ref[j], seg_ref[i]);
    }
  }

  for (Index i = 0; i < seg_ref.size(); i++) {
    ref.block<3, 3>(3 * i, 3 * i) = seg_ref[i].getPInv();
  }
  // building matrix via (i,j) operator
  Eigen::MatrixXd elementwise = Eigen::MatrixXd::Zero(6, 6);
  for (Index i = 0; i < 6; i++) {
    for (Index j = 0; j < 6; j++) {
      double value = dipdip(i, j);  // do not remove gcc debug needs this
      elementwise(i, j) = value;
    }
  }
  bool op_check = elementwise.isApprox(ref, 1e-6);
  BOOST_CHECK_EQUAL(op_check, 1);
  if (!op_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "operator" << std::endl;
    std::cout << elementwise << std::endl;
  }

  // building matrix via matrix vector product
  Eigen::MatrixXd gemv = Eigen::MatrixXd::Zero(6, 6);
  Eigen::MatrixXd ident = Eigen::MatrixXd::Identity(6, 6);
  for (Index i = 0; i < 6; i++) {
    gemv.col(i) = dipdip * ident.col(i);
  }
  bool gemv_check = gemv.isApprox(ref, 1e-6);
  BOOST_CHECK_EQUAL(gemv_check, 1);
  if (!gemv_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "gemv" << std::endl;
    std::cout << gemv << std::endl;
  }
  // building matrix via iterator product
  Eigen::MatrixXd iterator = Eigen::MatrixXd::Zero(6, 6);
  for (Index k = 0; k < dipdip.outerSize(); ++k) {
    for (DipoleDipoleInteraction::InnerIterator it(dipdip, k); it; ++it) {
      iterator(it.row(), k) = it.value();
    }
  }
  bool iterator_check = iterator.isApprox(ref, 1e-6);
  BOOST_CHECK_EQUAL(iterator_check, 1);
  if (!iterator_check) {
    std::cout << "ref" << std::endl;
    std::cout << ref << std::endl;
    std::cout << "iterator" << std::endl;
    std::cout << gemv << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
