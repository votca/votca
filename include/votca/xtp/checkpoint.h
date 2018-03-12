/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <Eigen/Eigen>
#include <H5Cpp.h>
#include <string>
#include <typeinfo>
#include <vector>

namespace votca {
namespace xtp {

// Define the types needed

namespace hdf5_utils {

H5::IntType int_type(H5::PredType::NATIVE_INT);
H5::StrType str_type(0, H5T_VARIABLE);
H5::FloatType double_type(H5::PredType::NATIVE_DOUBLE);

H5::DataSpace str_scalar(H5::DataSpace(H5S_SCALAR));

inline H5::DataSpace StrScalar() { return H5::DataSpace(H5S_SCALAR); }

// Declare some HDF5 data type inference stuff:
// Adapted from
// https://github.com/garrison/eigen3-hdf5/blob/2c782414251e75a2de9b0441c349f5f18fe929a2/eigen3-hdf5.hpp#L18

template <typename T>
struct InferDataType;

template <>
struct InferDataType<float> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_FLOAT; }
};

template <>
struct InferDataType<double> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_DOUBLE; }
};

template <>
struct InferDataType<std::string> {
  static const H5::DataType* get(void) {
    static const H5::StrType strtype(0, H5T_VARIABLE);
    return &strtype;
  }
};

template <typename T>
void WriteData(H5::Group loc, const Eigen::MatrixBase<T>& matrix,
               const std::string name) {

  hsize_t dims[2] = {matrix.rows(),
                     matrix.cols()};  // eigen vectors are n,1 matrices

  const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();

  H5::DataSet dataset;
  H5::DataSpace dp(2, dims);

  dataset = loc.createDataSet(name.c_str(), *dataType, dp);

  dataset.write(matrix.derived().data(), *dataType);
}

}  // namespace hdf5_utils

class CheckpointFile {
 public:
  CheckpointFile(std::string fileName);
  ~CheckpointFile();

  std::string getName();

 private:
  std::string _fileName;
  H5::H5File _fileHandle;
  H5::Attribute _version;

  H5::Group _child1;
};
}  // namespace xtp
}  // namespace votca
#endif  // CHECKPOINT_H
