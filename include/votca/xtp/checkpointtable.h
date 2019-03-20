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
#ifndef VOTCA_XTP_CHECKPOINT_TABLE_H
#define VOTCA_XTP_CHECKPOINT_TABLE_H

#include <H5Cpp.h>
#include <string>
#include <cstddef>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <votca/xtp/checkpoint_utils.h>



namespace votca {
namespace xtp {
    using namespace checkpoint_utils;

struct CptTable{
CptTable(const std::string& name, std::size_t rowSize, std::size_t nRows)
:_name(name), _rowStructure(rowSize), _nRows(nRows), _maxStringSize(256){};
    ~CptTable(){};

    template<typename U>
    typename std::enable_if <std::is_fundamental<U>::value>::type
    addCol(const U& item, const std::string& name, const size_t& offset){
        _rowStructure.insertMember(name, offset, *InferDataType<U>::get());
    }

    void addCol(const std::string& item, const std::string& name,
                const size_t& offset){

        if (item.size() > _maxStringSize){
            std::stringstream message;
            message << "String too long this table." << std::endl
                    << "Maximum allowed string size is" << _maxStringSize
                    << std::endl;

            throw std::runtime_error(message.str());
        }

        /* std::cout << "String is " << item.size() << " characters long." */
        /*           << std::endl; */

        _rowStructure.insertMember(name, offset, *InferDataType<std::string>::get());
    }

    void initialize(const CptLoc& loc){
        // create the dataspace...
        _dims[0] = _nRows;
        _dims[1] = 1;
        _dp = H5::DataSpace(2, _dims);

        _dataset = loc.createDataSet(_name.c_str(), _rowStructure,
                                    _dp);
    }

    void writeToRow(void* buffer, std::size_t idx){

        hsize_t fStart[2] = {idx, 0};
        hsize_t fCount[2] = {1, 1};

        hsize_t mStart[2] = {0, 0};
        hsize_t mCount[2] = {1, 1};

        hsize_t mDim[2] = {1,1};

        H5::DataSpace mspace(2, mDim);

        _dp.selectHyperslab(H5S_SELECT_SET, fCount, fStart);
        mspace.selectHyperslab(H5S_SELECT_SET, mCount, mStart);
        _dataset.write(buffer, _rowStructure, mspace, _dp);

    }

    H5::CompType _rowStructure;
    std::size_t _nRows;
    hsize_t _dims[2];
    H5::DataSpace _dp;
    H5::DataSet _dataset;
    std::string _name;
    long int _maxStringSize;
};

} // xtp
} // votca
#endif // VOTCA_XTP_CHECKPOINT_TABLE_H
