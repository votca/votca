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
#ifndef VOTCA_XTP_CHECKPOINTTABLE_H
#define VOTCA_XTP_CHECKPOINTTABLE_H

// Standard includes
#include <cstddef>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif

// Third party includes
#include <H5Cpp.h>

// Local VOTCA includes
#include "checkpoint_utils.h"

#define CPT_MEM_FROM_STRUCT(m, s) HOFFSET(s, m)
#define CPT_MEMBER(m, s) HOFFSET(s, m)

namespace votca {
namespace xtp {
using namespace checkpoint_utils;

class CptTable {
 public:
  CptTable() = default;
  CptTable(const std::string& name, const std::size_t& rowSize,
           const std::size_t& nRows)
      : name_(name),
        rowStructure_(rowSize),
        nRows_(nRows),
        props_(H5::DSetCreatPropList(H5::DSetCreatPropList::DEFAULT)){};

  CptTable(const std::string& name, const std::size_t& rowSize,
           const CptLoc& loc)
      : name_(name), loc_(loc), inited_(true), rowStructure_(rowSize) {
    dataset_ = loc_.openDataSet(name_);
    dp_ = dataset_.getSpace();
    hsize_t dims[2];
    dp_.getSimpleExtentDims(dims, nullptr);
    nRows_ = dims[0];
  }

  template <typename U>
  void addCol(const std::string& name, const size_t& offset);

  void initialize(const CptLoc& loc, bool compact) {
    // create the dataspace...
    if (inited_) {
      std::stringstream message;

      message << "Checkpoint tables cannot be reinitialized. " << name_
              << " has either already been initialized or already exists."
              << std::endl;

      throw std::runtime_error(message.str());
    }
    dims_[0] = nRows_;
    dims_[1] = 1;

    dp_ = H5::DataSpace(2, dims_);
    loc_ = loc;

    if (compact) {
      props_.setLayout(H5D_layout_t::H5D_COMPACT);
    }

    try {
      dataset_ = loc_.createDataSet(name_.c_str(), rowStructure_, dp_, props_);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name_ << " from " << loc_.getFileName();
      throw std::runtime_error(message.str());
    }

    inited_ = true;
  }

  void write(void* buffer, const std::size_t& startIdx,
             const std::size_t& endIdx) {

    if (!inited_) {
      std::stringstream message;
      message << "Checkpoint table uninitialized." << std::endl;
      throw std::runtime_error(message.str());
    }

    hsize_t s = (hsize_t)(startIdx);
    hsize_t e = (hsize_t)(endIdx);
    hsize_t l = e - s;

    hsize_t fStart[2] = {s, 0};
    hsize_t fCount[2] = {l, 1};

    hsize_t mStart[2] = {s, 0};
    hsize_t mCount[2] = {l, 1};

    hsize_t mDim[2] = {l, 1};

    H5::DataSpace mspace(2, mDim);

    dp_.selectHyperslab(H5S_SELECT_SET, fCount, fStart);
    mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart);
    try {
      dataset_.write(buffer, rowStructure_, mspace, dp_);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not write " << name_ << " from " << loc_.getFileName();
      throw std::runtime_error(message.str());
    }
  }

  void writeToRow(void* buffer, const std::size_t idx) {
    write(buffer, idx, idx + 1);
  }

  template <typename T>
  void write(std::vector<T>& dataVec) {
    write(dataVec.data(), 0, dataVec.size());
  }

  void read(void* buffer, const std::size_t& startIdx,
            const std::size_t& endIdx) {

    if (!inited_) {
      std::stringstream message;
      message << "Checkpoint table uninitialized." << std::endl;
      throw std::runtime_error(message.str());
    }

    hsize_t s = (hsize_t)(startIdx);
    hsize_t e = (hsize_t)(endIdx);
    hsize_t l = e - s;

    hsize_t fStart[2] = {s, 0};
    hsize_t fCount[2] = {l, 1};

    hsize_t mStart[2] = {s, 0};
    hsize_t mCount[2] = {l, 1};

    hsize_t mDim[2] = {l, 1};

    H5::DataSpace mspace(2, mDim);

    dp_.selectHyperslab(H5S_SELECT_SET, fCount, fStart);
    mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart);
    try {
      dataset_.read(buffer, rowStructure_, mspace, dp_);
    } catch (H5::Exception&) {
      std::stringstream message;
      message << "Could not read " << name_ << " from " << loc_.getFileName();
      throw std::runtime_error(message.str());
    }
  }

  void readFromRow(void* buffer, const std::size_t& idx) {
    read(buffer, idx, idx + 1);
  }

  template <typename T>
  void read(std::vector<T>& dataVec) {
    read(dataVec.data(), 0, dataVec.size());
  }

  std::size_t numRows() { return nRows_; }

  static const std::size_t MaxStringSize = 512;

 private:
  std::string name_;
  CptLoc loc_;
  bool inited_ = false;
  H5::CompType rowStructure_;
  std::size_t nRows_;
  hsize_t dims_[2];
  H5::DataSpace dp_;
  H5::DataSet dataset_;
  H5::DSetCreatPropList props_;
};

template <typename U>
inline void CptTable::addCol(const std::string& name, const size_t& offset) {
  static_assert(
      std::is_fundamental<U>::value,
      "Columns can only be added for fundamental types and 'const char*'");
  rowStructure_.insertMember(name, offset, *InferDataType<U>::get());
}

template <>
inline void CptTable::addCol<std::string>(const std::string& name,
                                          const size_t& offset) {
  rowStructure_.insertMember(name, offset, *InferDataType<std::string>::get());
}

}  // namespace xtp
}  // namespace votca

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif  // VOTCA_XTP_CHECKPOINTTABLE_H
