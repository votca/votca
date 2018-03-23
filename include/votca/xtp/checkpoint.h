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

#define CHECKPOINT_WRITE_SCAL(loc, var) WriteScalarAttribute(loc, var, #var)

namespace votca {
namespace xtp {

// Define the types needed

namespace hdf5_utils {

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
struct InferDataType<int> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_INT; }
};

template <>
struct InferDataType<unsigned int> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_UINT; }
};

template <>
struct InferDataType<std::string> {
  static const H5::DataType* get(void) {
    static const H5::StrType strtype(0, H5T_VARIABLE);
    return &strtype;
  }
};

template <typename T>
void WriteScalarAttribute(const H5::H5Location& loc, const T& value,
                          const std::string& name) {

  hsize_t dims[1] = {1};
  H5::DataSpace dp(1, dims);
  const H5::DataType* dataType = InferDataType<T>::get();

  H5::Attribute attr = loc.createAttribute(name.c_str(), *dataType, dp);
  attr.write(*dataType, &value);
}

void WriteScalarAttribute(const H5::H5Location& obj, const std::string& value,
                          const std::string& name);

class ScalarAttrWriter {
 public:
  ScalarAttrWriter(H5::Group& loc) : _loc(loc){};

  template <typename T>
  void operator()(const T& value, const std::string& name) {
    WriteScalarAttribute(_loc, value, name);
  };

 private:
  H5::Group _loc;
};

template <typename T>
void WriteData(const H5::Group& loc, const Eigen::MatrixBase<T>& matrix,
               const std::string& name) {

  hsize_t dims[2] = {matrix.rows(),
                     matrix.cols()};  // eigen vectors are n,1 matrices

  const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();

  H5::DataSet dataset;
  H5::DataSpace dp(2, dims);

  dataset = loc.createDataSet(name.c_str(), *dataType, dp);

  dataset.write(matrix.derived().data(), *dataType);
}

template <typename T>
void WriteData(const H5::Group& loc, const std::vector<T> v,
               const std::string& name) {
  hsize_t dims[2] = {(hsize_t)v.size(), 1};

  const H5::DataType* dataType = InferDataType<T>::get();
  H5::DataSet dataset;
  H5::DataSpace dp(2, dims);

  dataset = loc.createDataSet(name.c_str(), *dataType, dp);
  dataset.write(&(v[0]), *dataType);
}

template <typename T1, typename T2>
void WriteData(const H5::Group& loc, const std::map<T1, std::vector<T2>> map,
               const std::string& name) {

  size_t c = 0;
  std::string r;
  std::cout << name << std::endl;
  // Iterate over the map and write map as a number of vectors with T1 as index
  for (auto const& x : map) {
    r = std::to_string(c);
    H5::Group tempGr = loc.createGroup("/" + name);
    hdf5_utils::ScalarAttrWriter w(tempGr);
    w(x.first, "index" + r);
    WriteData(tempGr, x.second, "val" + r);
    ++c;
  }
}
}  // namespace hdf5_utils

class CheckpointFile {
 public:
  CheckpointFile(std::string fileName);

  std::string getFileName();
  std::string getVersion();

  H5::H5File getHandle();

 private:
  std::string _fileName;
  H5::H5File _fileHandle;
  H5::Group _child1;
  std::string _version;
};
}  // namespace xtp
}  // namespace votca
#endif  // CHECKPOINT_H
