/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#pragma once
#ifndef VOTCA_XTP_CHECKPOINTWRITER_H
#define VOTCA_XTP_CHECKPOINTWRITER_H

// Standard includes
#include <map>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif

// Third party includes
#include <H5Cpp.h>

// VOTCA includes
#include <votca/tools/linalg.h>

// Local VOTCA includes
#include "checkpoint_utils.h"
#include "checkpointtable.h"
#include "eigen.h"

namespace votca {
namespace xtp {

using namespace checkpoint_utils;

class CheckpointWriter {
 public:
  CheckpointWriter(const CptLoc& loc) : CheckpointWriter(loc, "/"){};

  CheckpointWriter(const CptLoc& loc, const std::string& path)
      : loc_(loc), path_(path){};

  // see the following links for details
  // https://stackoverflow.com/a/8671617/1186564
  template <typename T>
  typename std::enable_if<!std::is_fundamental<T>::value>::type operator()(
      const T& data, const std::string& name) const {
    try {
      WriteData(loc_, data, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name << " to " << loc_.getFileName()
              << ":" << path_;

      throw std::runtime_error(message.str());
    }
  }

  // Use this overload if T is a fundamental type
  // int, double, unsigned, etc, but not bool
  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value &&
                          !std::is_same<T, bool>::value>::type
      operator()(const T& v, const std::string& name) const {
    try {
      WriteScalar(loc_, v, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name << " to " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  void operator()(const bool& v, const std::string& name) const {
    Index temp = static_cast<Index>(v);
    try {
      WriteScalar(loc_, temp, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name << " to " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  void operator()(const std::string& v, const std::string& name) const {
    try {
      WriteScalar(loc_, v, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name << " to " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  CheckpointWriter openChild(const std::string& childName) const {
    try {
      return CheckpointWriter(loc_.openGroup(childName),
                              path_ + "/" + childName);
    } catch (H5::Exception&) {
      try {
        return CheckpointWriter(loc_.createGroup(childName),
                                path_ + "/" + childName);
      } catch (H5::Exception&) {
        std::stringstream message;
        message << "Could not open or create" << loc_.getFileName() << ":/"
                << path_ << "/" << childName << std::endl;

        throw std::runtime_error(message.str());
      }
    }
  }

  template <typename T>
  CptTable openTable(const std::string& name, std::size_t nRows,
                     bool compact = false) {
    CptTable table;
    try {
      table = CptTable(name, sizeof(typename T::data), loc_);
      T::SetupCptTable(table);
    } catch (H5::Exception&) {
      try {
        table = CptTable(name, sizeof(typename T::data), nRows);
        T::SetupCptTable(table);
        table.initialize(loc_, compact);
      } catch (H5::Exception&) {
        std::stringstream message;
        message << "Could not open table " << name << " in "
                << loc_.getFileName() << ":" << path_ << std::endl;
        throw std::runtime_error(message.str());
      }
    }

    return table;
  }

 private:
  const CptLoc loc_;
  const std::string path_;
  template <typename T>
  void WriteScalar(const CptLoc& loc, const T& value,
                   const std::string& name) const {

    hsize_t dims[1] = {1};
    H5::DataSpace dp(1, dims);
    const H5::DataType* dataType = InferDataType<T>::get();
    H5::Attribute attr;
    try {
      attr = loc.createAttribute(name, *dataType, dp);
    } catch (H5::AttributeIException&) {
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
    } catch (H5::AttributeIException&) {
      attr = loc.openAttribute(name);
    }

    const char* c_str_copy = value.c_str();
    attr.write(*strType, &c_str_copy);
  }

  template <typename T>
  void WriteData(const CptLoc& loc, const Eigen::MatrixBase<T>& matrix,
                 const std::string& name) const {

    hsize_t matRows = hsize_t(matrix.rows());
    hsize_t matCols = hsize_t(matrix.cols());

    hsize_t dims[2] = {matRows, matCols};  // eigen vectors are n,1 matrices

    if (dims[1] == 0) {
      dims[1] = 1;
    }

    H5::DataSpace dp(2, dims);
    const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();
    H5::DataSet dataset;
    try {
      dataset = loc.createDataSet(name.c_str(), *dataType, dp);
    } catch (H5::GroupIException&) {
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
    } catch (H5::GroupIException&) {
      dataset = loc.openDataSet(name.c_str());
    }
    dataset.write(v.data(), *dataType);
  }

  void WriteData(const CptLoc& loc, const std::vector<std::string>& v,
                 const std::string& name) const {

    hsize_t dims[1] = {(hsize_t)v.size()};

    std::vector<const char*> c_str_copy;
    c_str_copy.reserve(v.size());
    for (const std::string& s : v) {
      c_str_copy.push_back(s.c_str());
    }
    const H5::DataType* dataType = InferDataType<std::string>::get();
    H5::DataSet dataset;
    H5::DataSpace dp(1, dims);
    try {
      dataset = loc.createDataSet(name.c_str(), *dataType, dp);
    } catch (H5::GroupIException&) {
      dataset = loc.openDataSet(name.c_str());
    }
    dataset.write(c_str_copy.data(), *dataType);
  }

  void WriteData(const CptLoc& loc, const std::vector<Eigen::Vector3d>& v,
                 const std::string& name) const {

    size_t c = 0;
    std::string r;
    CptLoc parent;
    try {
      parent = loc.createGroup(name);
    } catch (H5::GroupIException&) {
      parent = loc.openGroup(name);
    }
    for (auto const& x : v) {
      r = std::to_string(c);
      WriteData(parent, x, "ind" + r);
      ++c;
    }
  }

  void WriteData(const CptLoc& loc, const tools::EigenSystem& sys,
                 const std::string& name) const {

    CptLoc parent;
    try {
      parent = loc.createGroup(name);
    } catch (H5::GroupIException&) {
      parent = loc.openGroup(name);
    }

    WriteData(parent, sys.eigenvalues(), "eigenvalues");
    WriteData(parent, sys.eigenvectors(), "eigenvectors");
    WriteData(parent, sys.eigenvectors2(), "eigenvectors2");
    WriteScalar(parent, Index(sys.info()), "info");
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
      } catch (H5::GroupIException&) {
        tempGr = loc.openGroup(name);
      }
      WriteData(tempGr, x.second, "index" + r);
      ++c;
    }
  }
};
}  // namespace xtp
}  // namespace votca

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif  // VOTCA_XTP_CHECKPOINTWRITER_H
