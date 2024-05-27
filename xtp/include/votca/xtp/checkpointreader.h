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
#ifndef VOTCA_XTP_CHECKPOINTREADER_H
#define VOTCA_XTP_CHECKPOINTREADER_H

// Standard includes
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

class CheckpointReader {
 public:
  CheckpointReader(const CptLoc& loc) : CheckpointReader(loc, "/"){};

  CheckpointReader(const CptLoc& loc, const std::string path)
      : loc_(loc), path_(path){};

  template <typename T>
  typename std::enable_if<!std::is_fundamental<T>::value>::type operator()(
      T& var, const std::string& name) const {
    try {
      ReadData(loc_, var, name);
    } catch (H5::Exception&) {
      std::stringstream message;

      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value &&
                          !std::is_same<T, bool>::value>::type
      operator()(T& var, const std::string& name) const {
    try {
      ReadScalar(loc_, var, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << "/" << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  void operator()(bool& v, const std::string& name) const {
    Index temp = Index(v);
    try {
      ReadScalar(loc_, temp, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
    v = static_cast<bool>(temp);
  }

  void operator()(std::string& var, const std::string& name) const {
    try {
      ReadScalar(loc_, var, name);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  CheckpointReader openChild(const std::string& childName) const {
    try {
      return CheckpointReader(loc_.openGroup(childName),
                              path_ + "/" + childName);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not open " << loc_.getFileName() << ":/" << path_ << "/"
              << childName << std::endl;

      throw std::runtime_error(message.str());
    }
  }

  std::vector<std::string> getChildGroupNames() const {
    std::vector<std::string> result;
    char pStr[128];
    for (hsize_t i = 0; i < loc_.getNumObjs(); i++) {
      memset(pStr, 0, 128);
      loc_.getObjnameByIdx(i, pStr, 128);
      try {
        loc_.openGroup(pStr);
        result.push_back(pStr);
      } catch (H5::Exception&) {
        ;
      }
    }
    return result;
  }

  Index getNumDataSets() const { return loc_.getNumObjs(); }

  CptLoc getLoc() { return loc_; }

  template <typename T>
  CptTable openTable(const std::string& name) {
    try {
      CptTable table = CptTable(name, sizeof(typename T::data), loc_);
      T::SetupCptTable(table);
      return table;
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not open table " << name << " in " << loc_.getFileName()
              << ":" << path_ << std::endl;
      throw std::runtime_error(message.str());
    }
  }

 private:
  const CptLoc loc_;
  const std::string path_;
  void ReadScalar(const CptLoc& loc, std::string& var,
                  const std::string& name) const {

    H5::Attribute attr = loc.openAttribute(name);
    H5::StrType stype = attr.getStrType();

    attr.read(stype, var);
  }

  template <typename T>
  void ReadScalar(const CptLoc& loc, T& value, const std::string& name) const {

    H5::Attribute attr = loc.openAttribute(name);
    const H5::DataType* dataType = InferDataType<T>::get();

    attr.read(*dataType, &value);
  }

  template <typename T>
  void ReadData(const CptLoc& loc, Eigen::MatrixBase<T>& matrix,
                const std::string& name) const {

    const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();

    H5::DataSet dataset = loc.openDataSet(name);

    H5::DataSpace dp = dataset.getSpace();

    hsize_t dims[2];
    dp.getSimpleExtentDims(dims, nullptr);  // ndims is always 2 for us

    hsize_t matRows = dims[0];
    hsize_t matCols = dims[1];

    matrix.derived().resize(matRows, matCols);
    if (matrix.size() == 0) {
      return;
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
      dataset.read(matrix.derived().data(), *dataType, mspace, dp);
    }
  }

  template <typename T>
  typename std::enable_if<std::is_fundamental<T>::value>::type ReadData(
      const CptLoc& loc, std::vector<T>& v, const std::string& name) const {

    H5::DataSet dataset = loc.openDataSet(name);
    H5::DataSpace dp = dataset.getSpace();

    const H5::DataType* dataType = InferDataType<T>::get();

    hsize_t dims[2];
    dp.getSimpleExtentDims(dims, nullptr);

    v.resize(dims[0]);
    if (v.empty()) {
      return;
    }
    try {
      dataset.read(&(v[0]), *dataType);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << std::endl;
      throw std::runtime_error(message.str());
    }
  }

  void ReadData(const CptLoc& loc, std::vector<std::string>& v,
                const std::string& name) const {

    H5::DataSet dataset = loc.openDataSet(name);
    H5::DataSpace dp = dataset.getSpace();

    const H5::DataType* dataType = InferDataType<std::string>::get();

    hsize_t dims[2];
    dp.getSimpleExtentDims(dims, nullptr);

    std::vector<char*> temp(dims[0]);
    if (temp.empty()) {
      return;
    }
    try {
      dataset.read(temp.data(), *dataType);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name << " from " << loc_.getFileName()
              << ":" << path_ << std::endl;
      throw std::runtime_error(message.str());
    }
    v.reserve(dims[0]);
    for (char* s : temp) {
      v.push_back(std::string(s));
      free(s);
    }
  }

  void ReadData(const CptLoc& loc, tools::EigenSystem& sys,
                const std::string& name) const {

    CptLoc parent = loc.openGroup(name);
    ReadData(parent, sys.eigenvalues(), "eigenvalues");
    ReadData(parent, sys.eigenvectors(), "eigenvectors");
    ReadData(parent, sys.eigenvectors2(), "eigenvectors2");
    Index info;
    ReadScalar(parent, info, "info");
    sys.info() = static_cast<Eigen::ComputationInfo>(info);
  }

  void ReadData(const CptLoc& loc, std::vector<Eigen::Vector3d>& v,
                const std::string& name) const {

    CptLoc parent = loc.openGroup(name);
    size_t count = parent.getNumObjs();

    v.resize(count);

    size_t c = 0;
    for (auto& vec : v) {
      ReadData(parent, vec, "ind" + std::to_string(c));
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

#endif  // VOTCA_XTP_CHECKPOINTREADER_H
