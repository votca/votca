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

#ifndef _VOTCA_XTP_CHECKPOINT_WRITER_H
#define _VOTCA_XTP_CHECKPOINT_WRITER_H

#include <votca/xtp/eigen.h>
#include <map>
#include <H5Cpp.h>
#include <string>
#include <typeinfo>
#include <vector>
#include <type_traits>
#include <votca/tools/vec.h>

#include <votca/xtp/checkpoint_utils.h>

namespace votca { namespace xtp {

using namespace checkpoint_utils;

class CheckpointWriter {
public:
CheckpointWriter(const CptLoc& loc) : _loc(loc){};

    // see the following links for details
    // https://stackoverflow.com/a/8671617/1186564
    template <typename T>
        typename std::enable_if<!std::is_fundamental<T>::value>::type
        operator()(const T& data, const std::string& name){
        WriteData(_loc, data, name);
    }

    // Use this overload iff T is a fundamental type
    // int, double, unsigned int, etc.
    template<typename T>
        typename std::enable_if<std::is_fundamental<T>::value>::type
        operator()(const T& v, const std::string& name){
        WriteScalar(_loc, v, name);
    }

    void operator()(const std::string& v, const std::string& name){
        WriteScalar(_loc, v, name);
    }

private:
    CptLoc _loc;
    template <typename T>
        void WriteScalar(const CptLoc& loc, const T& value,
                         const std::string& name) {

        hsize_t dims[1] = {1};
        H5::DataSpace dp(1, dims);
        const H5::DataType* dataType = InferDataType<T>::get();

        H5::Attribute attr = loc.createAttribute(name, *dataType, dp);
        attr.write(*dataType, &value);
    }

    void WriteScalar(const CptLoc& loc, const std::string& value,
                     const std::string& name) {

        hsize_t dims[1] = {1};
        H5::DataSpace dp(1, dims);
        const H5::DataType* strType = InferDataType<std::string>::get();

        H5::Attribute attr = loc.createAttribute(name, *strType, StrScalar());
        attr.write(*strType, &value);
    }

    template <typename T>
        void WriteData(const CptLoc& loc, const Eigen::MatrixBase<T>& matrix,
                       const std::string& name) {

        hsize_t dims[2] = {size_t(matrix.rows()),
                           size_t(matrix.cols())};  // eigen vectors are n,1 matrices

        if (dims[1] == 0) dims[1] = 1;

        const H5::DataType* dataType = InferDataType<typename T::Scalar>::get();

        H5::DataSet dataset;
        H5::DataSpace dp(2, dims);

        dataset = loc.createDataSet(name.c_str(), *dataType, dp);

        dataset.write(matrix.derived().data(), *dataType);
    }

    template <typename T>
        typename std::enable_if<std::is_fundamental<T>::value>::type
        WriteData(const CptLoc& loc, const std::vector<T> v,
                  const std::string& name) {
        hsize_t dims[2] = {(hsize_t)v.size(), 1};

        const H5::DataType* dataType = InferDataType<T>::get();
        H5::DataSet dataset;
        H5::DataSpace dp(2, dims);

        dataset = loc.createDataSet(name.c_str(), *dataType, dp);
        dataset.write(&(v[0]), *dataType);
    }

    void WriteData(const CptLoc& loc, const votca::tools::vec& v,
                   const std::string& name) {

        // store tools::vec as n vector of three elements
        std::vector<double> data = {v.getX(), v.getY(), v.getZ()};
        WriteData(loc, data, name);
    }

    void WriteData(const CptLoc& loc, const std::vector<votca::tools::vec>& v,
                   const std::string& name) {

        size_t c = 0;
        std::string r;
        CptLoc parent = loc.createGroup(name);
        for (auto const& x: v){
            r = std::to_string(c);
            WriteData(parent, x, "ind"+r);
            ++c;
        }
    }

    template <typename T1, typename T2>
        void WriteData(const CptLoc& loc, const std::map<T1, std::vector<T2>> map,
                       const std::string& name) {

        size_t c = 0;
        std::string r;
        // Iterate over the map and write map as a number of vectors with T1 as index
        for (auto const& x : map) {
            r = std::to_string(c);
            CptLoc tempGr = loc.createGroup("/" + name);
            WriteData(tempGr, x.second, "index" + r);
            ++c;
        }
    }

};
}  // namespace xtp
}  // namespace votca
#endif  // _VOTCA_XTP_CHECKPOINT_WRITER_H
