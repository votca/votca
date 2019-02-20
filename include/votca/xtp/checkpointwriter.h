/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_XTP_CHECKPOINT_WRITER_H
#define VOTCA_XTP_CHECKPOINT_WRITER_H

#include <H5Cpp.h>
#include <map>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>
#include <votca/tools/vec.h>
#include <votca/xtp/eigen.h>

#include <votca/xtp/checkpoint_utils.h>

namespace votca {
namespace xtp {

using namespace checkpoint_utils;

class CheckpointWriter {
 public:
  CheckpointWriter(const CptLoc& loc) : CheckpointWriter(loc, "/"){};

  CheckpointWriter(const CptLoc& loc, const std::string& path)
      : _loc(loc), _path(path){};

  // see the following links for details
  // https://stackoverflow.com/a/8671617/1186564
  template <typename T>
  typename std::enable_if<!std::is_fundamental<T>::value>::type operator()(
      const T& data, const std::string& name) const {
    try {
      WriteData(_loc, data, name);
    } catch (H5::Exception& error) {
      std::stringstream message;
      message << "Could not write " << name << " to " << _loc.getFileName()
              << ":" << _path;

      throw std::runtime_error(message.str());
    }
  }

  // Use this overload if T is a fundamental type
  // int, double, unsigned int, etc, but not bool
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value &&
                          !std::is_same<T, bool>::value>::type
      operator()(const T& v, const std::string& name) const {
    try {
      WriteScalar(_loc, v, name);
    } catch (H5::Exception& error) {
      std::stringstream message;
      message << "Could not write " << name << " to " << _loc.getFileName()
              << ":" << _path << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  void operator()(const bool& v, const std::string& name) const {
    int temp = static_cast<int>(v);
    try {
      WriteScalar(_loc, temp, name);
    } catch (H5::Exception& error) {
      std::stringstream message;
      message << "Could not write " << name << " to " << _loc.getFileName()
              << ":" << _path << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  void operator()(const std::string& v, const std::string& name) const {
    try {
      WriteScalar(_loc, v, name);
    } catch (H5::Exception& error) {
      std::stringstream message;
      message << "Could not write " << name << " to " << _loc.getFileName()
              << ":" << _path << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  CheckpointWriter openChild(const std::string& childName) const {
    try {
      return CheckpointWriter(_loc.openGroup(childName),
                              _path + "/" + childName);
    } catch (H5::Exception& e) {
      try {
        return CheckpointWriter(_loc.createGroup(childName),
                                _path + "/" + childName);
      } catch (H5::Exception& e) {
        std::stringstream message;
        message << "Could not open or create" << _loc.getFileName() << ":/"
                << _path << "/" << childName << std::endl;

        throw std::runtime_error(message.str());
      }
    }
  }

 private:
  const CptLoc _loc;
  const std::string _path;
  template <typename T>
  void WriteScalar(const CptLoc& loc, const T& value,
                   const std::string& name) const {

    hsize_t dims[1] = {1};
    H5::DataSpace dp(1, dims);
    const H5::DataType* dataType = InferDataType<T>::get();
    H5::Attribute attr;
    try {
      attr = loc.createAttribute(name, *dataType, dp);
    } catch (H5::AttributeIException& error) {
      attr = loc.openAttribute(name);
    }
    attr.write(*dataType, &value);
  }

  void WriteScalar(const CptLoc& loc, const std::string& value,
                   const std::string& name) const {

    hsize_t dims[1] = {1};
    H5::DataSpace dp(1, dims);
    const H5::DataType* strType = InferDataType<std::string>::get();

    H5::Attribute attr;

    try {
      attr = loc.createAttribute(name, *strType, dp);
    } catch (H5::AttributeIException& error) {
      attr = loc.openAttribute(name);
    }
    attr.write(*strType, &value);
  }

  template <typename T>
  void WriteData(const CptLoc& loc, const Eigen::MatrixBase<T>& matrix,
                 const std::string& name) const {

    hsize_t matRows = hsize_t(matrix.rows());
    hsize_t matCols = hsize_t(matrix.cols());

    hsize_t dims[2] = {matRows, matCols};  // eigen vectors are n,1 matrices

    if (dims[1] == 0) dims[1] = 1;

    H5::DataSpace dp(2, dims);
    const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();
    H5::DataSet dataset;
    try {
      dataset = loc.createDataSet(name.c_str(), *dataType, dp);
    } catch (H5::GroupIException& error) {
      dataset = loc.openDataSet(name.c_str());
    }

    hsize_t matColSize = matrix.derived().outerStride();

    hsize_t fileRows = matCols;

    hsize_t fStride[2] = {1, fileRows};
    hsize_t fCount[2] = {1, 1};
    hsize_t fBlock[2] = {1, fileRows};

    hsize_t mStride[2] = {matColSize, 1};
    hsize_t mCount[2] = {1, 1};
    hsize_t mBlock[2] = {matCols, 1};

    hsize_t mDim[2] = {matCols, matColSize};
    H5::DataSpace mspace(2, mDim);

    for (hsize_t i = 0; i < matRows; i++) {
      hsize_t fStart[2] = {i, 0};
      hsize_t mStart[2] = {0, i};
      dp.selectHyperslab(H5S_SELECT_SET, fCount, fStart, fStride, fBlock);
      mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart, mStride, mBlock);
      dataset.write(matrix.derived().data(), *dataType, mspace, dp);
    }
  }

  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value>::type WriteData(
      const CptLoc& loc, const std::vector<T> v,
      const std::string& name) const {
    hsize_t dims[2] = {(hsize_t)v.size(), 1};

    const H5::DataType* dataType = InferDataType<T>::get();
    H5::DataSet dataset;
    H5::DataSpace dp(2, dims);
    try {
      dataset = loc.createDataSet(name.c_str(), *dataType, dp);
    } catch (H5::GroupIException& error) {
      dataset = loc.openDataSet(name.c_str());
    }
    dataset.write(&(v[0]), *dataType);
  }

  void WriteData(const CptLoc& loc, const std::vector<Eigen::Vector3d>& v,
                 const std::string& name) const {

    size_t c = 0;
    std::string r;
    CptLoc parent;
    try {
      parent = loc.createGroup(name);
    } catch (H5::GroupIException& error) {
      parent = loc.openGroup(name);
    }
    for (auto const& x : v) {
      r = std::to_string(c);
      WriteData(parent, x, "ind" + r);
      ++c;
    }
  }

  template <typename T1, typename T2>
  void WriteData(const CptLoc& loc, const std::map<T1, std::vector<T2>> map,
                 const std::string& name) const {

    size_t c = 0;
    std::string r;
    // Iterate over the map and write map as a number of vectors with T1 as
    // index
    for (auto const& x : map) {
      r = std::to_string(c);
      CptLoc tempGr;
      try {
        tempGr = loc.createGroup(name);
      } catch (H5::GroupIException& error) {
        tempGr = loc.openGroup(name);
      }
      WriteData(tempGr, x.second, "index" + r);
      ++c;
    }
  }
};
}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_WRITER_H
