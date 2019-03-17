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

#ifndef VOTCA_XTP_CHECKPOINT_UTILS_H
#define VOTCA_XTP_CHECKPOINT_UTILS_H

#include <H5Cpp.h>
#include <string>
#include <cstddef>

namespace votca {
namespace xtp {

typedef H5::Group CptLoc;

namespace checkpoint_utils {

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

H5::DataSpace str_scalar(H5::DataSpace(H5S_SCALAR));

struct CptTable{
CptTable(const std::string& name, std::size_t rowSize, std::size_t nRows)
:_name(name), _rowStructure(rowSize), _nRows(){};
    ~CptTable(){};

    template<typename U>
    //typename std::enable_if <std::is_fundamental<U>::value>::type
    void addCol(const U& item, const std::string& name, const size_t& offset){
        _rowStructure.insertMember(name, offset, *InferDataType<U>::get());
    }

    void initialize(CptLoc loc){
        // create the dataspace...
        _dims[0] = _nRows;
        _dims[1] = 1;
        _dp = H5::DataSpace(2, _dims);

        _dataset = loc.createDataSet(_name.c_str(), _rowStructure,
                                    _dp);
    }

    void writeToRow(void* buffer, std::size_t idx){

        hsize_t fStart[2] = {idx};
        hsize_t fStride[2] = {1};
        hsize_t fCount[2] = {1};
        hsize_t fBlock[2] = {1};

        hsize_t mStart[2] = {0};
        hsize_t mStride[2] = {1};
        hsize_t mCount[2] = {1};
        hsize_t mBlock[2] = {1};

        hsize_t mDim[2] = {1,1};

        H5::DataSpace mspace(2, mDim);

        /* _dp.selectHyperslab(H5S_SELECT_SET, fCount, fStart, fStride, fBlock); */
        /* mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart, mStride, mBlock); */
        /* _dataset.write(buffer, _rowStructure, mspace, _dp); */

    }

    H5::CompType _rowStructure;
    std::size_t _nRows;
    hsize_t _dims[2];
    H5::DataSpace _dp;
    H5::DataSet _dataset;
    std::string _name;
};

} // checkpoint_utils
} // xtp
} // votca
#endif //VOTCA_XTP_CHECKPOINT_UTILS_H
