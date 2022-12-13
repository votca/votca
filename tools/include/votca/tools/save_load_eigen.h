#ifndef BOOST_SERIALIZATION_IO_SAVELOAD_EIGEN_H
#define BOOST_SERIALIZATION_IO_SAVELOAD_EIGEN_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>

// Serialization methods for fixed or dynamic-size Eigen::Matrix type
namespace boost {
namespace serialization {

//  friend class boost::serialization::access;

template <class Archive, typename _Scalar, int _Rows, int _Cols, int _Options,
          int _MaxRows, int _MaxCols>
inline void serialize(
    Archive& ar,
    Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix,
    const unsigned int /* aVersion */) {
  Eigen::Index rows = matrix.rows();
  Eigen::Index cols = matrix.cols();
  ar&(rows);
  ar&(cols);
  if (rows != matrix.rows() || cols != matrix.cols()) matrix.resize(rows, cols);
  if (matrix.size() != 0)
    ar& boost::serialization::make_array(matrix.data(), rows * cols);
}
}  // namespace serialization
}  // namespace boost

#endif  // BOOST_SERIALIZATION_IO_SAVELOAD_EIGEN_H */