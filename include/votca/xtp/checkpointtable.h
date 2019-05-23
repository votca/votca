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
#pragma once
#ifndef VOTCA_XTP_CHECKPOINT_TABLE_H
#define VOTCA_XTP_CHECKPOINT_TABLE_H

#include <H5Cpp.h>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <votca/xtp/checkpoint_utils.h>

#define CPT_MEM_FROM_STRUCT(m, s) HOFFSET(s, m)
#define CPT_MEMBER(m, s) HOFFSET(s, m)

namespace votca {
namespace xtp {
using namespace checkpoint_utils;

class CptTable {
 public:
  CptTable(){};
  CptTable(const std::string& name, const std::size_t& rowSize,
           const std::size_t& nRows)
      : _name(name),
        _rowStructure(rowSize),
        _nRows(nRows),
        _props(H5::DSetCreatPropList(H5::DSetCreatPropList::DEFAULT)){};

  CptTable(const std::string& name, const std::size_t& rowSize,
           const CptLoc& loc)
      : _name(name), _loc(loc), _inited(true), _rowStructure(rowSize) {

    _dataset = _loc.openDataSet(_name);
    _dp = _dataset.getSpace();
    hsize_t dims[2];
    _dp.getSimpleExtentDims(dims, NULL);
    _nRows = dims[0];
  }

  template <typename U>
  typename std::enable_if<std::is_fundamental<U>::value>::type addCol(
      const U& item, const std::string& name, const size_t& offset) {
    _rowStructure.insertMember(name, offset, *InferDataType<U>::get());
  }

  void addCol(const std::string& item, const std::string& name,
              const size_t& offset) {

    _rowStructure.insertMember(name, offset,
                               *InferDataType<std::string>::get());
  }

  void addCol(const char* item, const std::string& name, const size_t& offset) {

    H5::DataType fixedWidth(H5T_STRING, MaxStringSize);

    _rowStructure.insertMember(name, offset, fixedWidth);
  }

  void initialize(const CptLoc& loc, bool compact) {
    // create the dataspace...
    if (_inited) {
      std::stringstream message;

      message << "Checkpoint tables cannot be reinitialized. " << _name
              << " has either already been initialized or already exists."
              << std::endl;

      throw std::runtime_error(message.str());
    }
    _dims[0] = _nRows;
    _dims[1] = 1;

    _dp = H5::DataSpace(2, _dims);
    _loc = loc;

    if (compact) {
      _props.setLayout(H5D_layout_t::H5D_COMPACT);
    }

    _dataset = _loc.createDataSet(_name.c_str(), _rowStructure, _dp, _props);
    _inited = true;
  }

  void write(void* buffer, const std::size_t& startIdx,
             const std::size_t& endIdx) {

    if (!_inited) {
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

    _dp.selectHyperslab(H5S_SELECT_SET, fCount, fStart);
    mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart);
    _dataset.write(buffer, _rowStructure, mspace, _dp);
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

    if (!_inited) {
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

    _dp.selectHyperslab(H5S_SELECT_SET, fCount, fStart);
    mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart);
    _dataset.read(buffer, _rowStructure, mspace, _dp);
  }

  void readFromRow(void* buffer, const std::size_t& idx) {
    read(buffer, idx, idx + 1);
  }

  template <typename T>
  void read(std::vector<T>& dataVec) {
    read(dataVec.data(), 0, dataVec.size());
  }

  std::size_t numRows() { return _nRows; }

  static const std::size_t MaxStringSize = 512;

 private:
  std::string _name;
  CptLoc _loc;
  bool _inited = false;
  H5::CompType _rowStructure;
  std::size_t _nRows;
  hsize_t _dims[2];
  H5::DataSpace _dp;
  H5::DataSet _dataset;
  H5::DSetCreatPropList _props;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_TABLE_H
