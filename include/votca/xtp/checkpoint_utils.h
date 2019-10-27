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
#ifndef VOTCA_XTP_CHECKPOINT_UTILS_H
#define VOTCA_XTP_CHECKPOINT_UTILS_H

#include <H5Cpp.h>
#include <cstddef>
#include <cstring>
#include <string>

namespace votca {
namespace xtp {

using CptLoc = H5::Group;

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
struct InferDataType<long> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_LONG; }
};

template <>
struct InferDataType<unsigned> {
  static const H5::DataType* get(void) { return &H5::PredType::NATIVE_UINT; }
};

template <>
struct InferDataType<std::string> {
  static const H5::DataType* get(void) {

#if (defined(__GNUC__) && defined(__clang__))
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#elif (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif
    static const H5::StrType strtype(H5T_C_S1, H5T_VARIABLE);
#if (defined(__GNUC__) && defined(__clang__))
#pragma clang diagnostic pop
#elif (defined(__GNUC__) && !defined(__INTEL_COMPILER))
#pragma GCC diagnostic pop
#endif

    return &strtype;
  }
};

H5::DataSpace str_scalar(H5::DataSpace(H5S_SCALAR));

}  // namespace checkpoint_utils
}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_CHECKPOINT_UTILS_H
